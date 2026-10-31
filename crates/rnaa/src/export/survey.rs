//! Frozen survey consistency checks and biological-sample metadata summary.
use super::Entry;
use anyhow::{Context, Result, bail};
use rnaa_core::{
    Database, ManifestArtifact,
    generations::{GenerationHead, OutputGroup},
};
use serde::Serialize;
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    path::{Path, PathBuf},
};

#[derive(Serialize)]
pub(super) struct SurveyProof {
    matrix_sha256: String,
    matrix_bytes: u64,
    correlation_input_generation: Option<String>,
    shared_dependencies: BTreeMap<String, String>,
    #[serde(skip)]
    pub(super) design_contract: Value,
}

pub(super) fn compatibility(
    db: &Database,
    de: &GenerationHead,
    corr: &GenerationHead,
    files: &[(GenerationHead, ManifestArtifact, bool)],
) -> Result<SurveyProof> {
    let de_request = db
        .scientific_request(&de.request_sha256)?
        .context("analysis has no verified scientific request")?;
    let mut corr_request = db
        .scientific_request(&corr.request_sha256)?
        .context("correlation has no verified scientific request")?;
    let mut shared = BTreeMap::new();
    for key in ["active_cohort", "eligibility", "reference"] {
        let left = de_request
            .dependencies
            .get(key)
            .with_context(|| format!("analysis request lacks {key} identity"))?;
        let right = corr_request
            .dependencies
            .get(key)
            .with_context(|| format!("correlation request lacks {key} identity"))?;
        if left != right {
            bail!("correlation and analysis {key} identities differ");
        }
        shared.insert(key.into(), left.clone());
    }
    drop(de_request);
    let vst: Vec<_> = files
        .iter()
        .filter(|(head, file, provenance)| {
            head.output_group == OutputGroup::De
                && !provenance
                && file.kind == "VST"
                && Path::new(&file.path)
                    .file_name()
                    .is_some_and(|name| name == "vst.tsv")
        })
        .collect();
    if vst.len() != 1 {
        bail!("survey requires exactly one selected VST matrix");
    }
    let vst = &vst[0].1;
    let matrix: Vec<_> = corr_request
        .inputs
        .iter()
        .filter(|input| input.role == "MATRIX")
        .collect();
    if matrix.len() != 1 || matrix[0].sha256 != vst.checksum || matrix[0].bytes != vst.bytes {
        bail!("correlation input does not match the selected analysis VST matrix");
    }
    drop(matrix);
    let design = std::mem::take(
        corr_request
            .parameters
            .get_mut("design_contract")
            .context("correlation design contract is missing")?,
    );
    if !design.is_null() && !design.is_object() {
        bail!("correlation design contract is invalid");
    }
    drop(corr_request);
    let lineage = db.generation_lineage(corr)?;
    if !lineage.recorded {
        bail!("correlation generation has unresolved input lineage");
    }
    let producers: Vec<_> = lineage
        .dependencies
        .iter()
        .filter(|edge| {
            edge.input.kind == "MATRIX"
                && edge.input.checksum_type == "sha256"
                && edge.input.checksum == vst.checksum
                && edge.input.bytes == vst.bytes
        })
        .collect();
    let external = lineage
        .external_inputs
        .iter()
        .filter(|input| {
            input.kind == "MATRIX"
                && input.checksum_type == "sha256"
                && input.checksum == vst.checksum
                && input.bytes == vst.bytes
        })
        .count();
    if producers.len() + external != 1 {
        bail!("correlation matrix identity lacks one verified lineage witness");
    }
    Ok(SurveyProof {
        matrix_sha256: vst.checksum.clone(),
        matrix_bytes: vst.bytes,
        correlation_input_generation: producers
            .first()
            .map(|edge| edge.producer_generation_id.clone()),
        shared_dependencies: shared,
        design_contract: design,
    })
}

fn table(root: &Path, entries: &[Entry], name: &str) -> Result<PathBuf> {
    let matches: Vec<_> = entries
        .iter()
        .filter(|entry| entry.output == name)
        .collect();
    if matches.len() != 1 {
        bail!("survey requires exactly one {name} table");
    }
    Ok(root.join(&matches[0].output))
}

fn reader(path: &Path) -> Result<csv::Reader<File>> {
    Ok(csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?)
}

pub(super) fn summarize(
    root: &Path,
    entries: &[Entry],
    stratify: Option<&str>,
    design_contract: &Value,
    check: &impl Fn() -> Result<()>,
) -> Result<Value> {
    let metadata_path = table(root, entries, "de/sample_metadata.tsv")?;
    let mut metadata = reader(&metadata_path)?;
    let headers = metadata.headers()?.clone();
    let sample_index = headers
        .iter()
        .position(|name| name == "sample_id")
        .context("frozen sample metadata has no sample_id")?;
    let field = stratify.map(|name| {
        if name.trim().is_empty() || name != name.trim() { bail!("stratification column must be nonempty and trimmed"); }
        let indices: Vec<_> = headers.iter().enumerate().filter(|(_, column)| *column == name).map(|(index, _)| index).collect();
        if indices.len() != 1 { bail!("stratification column {name:?} is missing or ambiguous in frozen sample metadata"); }
        Ok(indices[0])
    }).transpose()?;
    let mut samples = Vec::new();
    let mut unique = BTreeSet::new();
    let mut levels = BTreeMap::<String, u64>::new();
    let mut missing = 0_u64;
    let mut design_rows = Vec::new();
    for record in metadata.records() {
        if samples.len().is_multiple_of(1024) {
            check()?;
        }
        let record = record?;
        let sample = record
            .get(sample_index)
            .context("metadata row lacks sample ID")?;
        if sample.trim().is_empty() || !unique.insert(sample.to_string()) {
            bail!("frozen sample IDs must be nonempty and unique");
        }
        samples.push(sample.to_string());
        if !design_contract.is_null() {
            design_rows.push(record.clone());
        }
        if let Some(index) = field {
            let value = record
                .get(index)
                .context("metadata row lacks stratification value")?;
            if value.is_empty() {
                missing += 1;
            } else {
                *levels.entry(value.to_string()).or_default() += 1;
            }
        }
    }
    if samples.is_empty() {
        bail!("survey sample metadata is empty");
    }
    if !design_contract.is_null() {
        let expected_samples = design_contract["samples"]
            .as_array()
            .context("correlation design samples missing")?;
        if !expected_samples
            .iter()
            .map(Value::as_str)
            .eq(samples.iter().map(|sample| Some(sample.as_str())))
        {
            bail!("correlation design sample IDs differ from frozen sample metadata");
        }
        for (name, column) in design_contract["columns"]
            .as_object()
            .context("correlation design columns missing")?
        {
            let index = headers
                .iter()
                .position(|header| header == name)
                .with_context(|| {
                    format!("correlation covariate {name:?} is missing from frozen sample metadata")
                })?;
            let values = column["values"]
                .as_array()
                .context("correlation design values missing")?;
            if values.len() != design_rows.len() {
                bail!("correlation design values have the wrong length");
            }
            for (row, expected) in design_rows.iter().zip(values) {
                let actual = row
                    .get(index)
                    .context("frozen sample metadata row is incomplete")?;
                let same = match column["type"].as_str() {
                    Some("factor") => expected.as_str() == Some(actual),
                    Some("numeric") => actual
                        .parse::<f64>()
                        .ok()
                        .zip(expected.as_f64())
                        .is_some_and(|(a, b)| a.is_finite() && a == b),
                    _ => bail!("unsupported correlation covariate type"),
                };
                if !same {
                    bail!("correlation covariate {name:?} differs from frozen sample metadata");
                }
            }
        }
    }
    let mut matrices = [
        reader(&table(root, entries, "de/gene_counts.tsv")?)?,
        reader(&table(root, entries, "de/gene_norm_counts.tsv")?)?,
        reader(&table(root, entries, "de/vst.tsv")?)?,
    ];
    for matrix in &mut matrices {
        let header = matrix.headers()?;
        if header.get(0) != Some("gene_id")
            || !header.iter().skip(1).eq(samples.iter().map(String::as_str))
        {
            bail!("survey matrix sample IDs differ from frozen sample metadata");
        }
    }
    let mut genes = BTreeSet::new();
    let mut iterators = matrices
        .iter_mut()
        .map(|matrix| matrix.byte_records())
        .collect::<Vec<_>>();
    loop {
        if genes.len().is_multiple_of(1024) {
            check()?;
        }
        let next = iterators
            .iter_mut()
            .map(|records| records.next().transpose())
            .collect::<std::result::Result<Vec<_>, _>>()?;
        if next.iter().all(Option::is_none) {
            break;
        }
        let ids = next
            .iter()
            .map(|record| record.as_ref().and_then(|row| row.get(0)))
            .collect::<Vec<_>>();
        let first = ids[0].context("survey matrices have different gene counts")?;
        if first.is_empty()
            || ids.iter().any(|id| *id != Some(first))
            || !genes.insert(first.to_vec())
        {
            bail!("survey matrices contain inconsistent, duplicate or empty gene IDs");
        }
    }
    if genes.is_empty() {
        bail!("survey matrices contain no genes");
    }
    let mut adjusted = reader(&table(root, entries, "corr/adjusted_matrix.tsv")?)?;
    let header = adjusted.headers()?;
    if header.get(0) != Some("gene_id")
        || !header.iter().skip(1).eq(samples.iter().map(String::as_str))
    {
        bail!("correlation matrix sample IDs differ from frozen sample metadata");
    }
    let mut correlation_genes = BTreeSet::new();
    for row in adjusted.byte_records() {
        if correlation_genes.len().is_multiple_of(1024) {
            check()?;
        }
        let row = row?;
        let id = row.get(0).context("correlation matrix row lacks gene ID")?;
        if !genes.contains(id) || !correlation_genes.insert(id.to_vec()) {
            bail!("correlation genes are duplicated or absent from the analysis matrix");
        }
    }
    if correlation_genes.is_empty() {
        bail!("correlation matrix contains no genes");
    }
    for entry in entries
        .iter()
        .filter(|entry| entry.source.kind == "DE_TABLE")
    {
        let mut result = reader(&root.join(&entry.output))?;
        let gene_column = result
            .headers()?
            .iter()
            .position(|name| name == "gene_id")
            .context("DE result lacks gene_id")?;
        let mut seen = BTreeSet::new();
        for row in result.byte_records() {
            if seen.len().is_multiple_of(1024) {
                check()?;
            }
            let row = row?;
            let id = row
                .get(gene_column)
                .context("DE result row lacks gene ID")?;
            if !genes.contains(id) || !seen.insert(id.to_vec()) {
                bail!("DE result genes are duplicated or absent from the analysis matrix");
            }
        }
        if seen.is_empty() {
            bail!("DE result contains no genes");
        }
    }
    for entry in entries
        .iter()
        .filter(|entry| entry.source.kind == "CORR_EDGES")
    {
        let mut edges = reader(&root.join(&entry.output))?;
        let columns = edges.headers()?;
        let first = columns
            .iter()
            .position(|name| matches!(name, "gene_a" | "gene_a_id"))
            .context("correlation edges lack gene_a")?;
        let second = columns
            .iter()
            .position(|name| matches!(name, "gene_b" | "gene_b_id"))
            .context("correlation edges lack gene_b")?;
        for (index, row) in edges.byte_records().enumerate() {
            if index.is_multiple_of(1024) {
                check()?;
            }
            let row = row?;
            if !row
                .get(first)
                .is_some_and(|id| correlation_genes.contains(id))
                || !row
                    .get(second)
                    .is_some_and(|id| correlation_genes.contains(id))
            {
                bail!("correlation edges reference genes absent from the selected matrix");
            }
        }
    }
    for entry in entries
        .iter()
        .filter(|entry| entry.source.kind == "CORR_DENSE")
    {
        let mut dense = reader(&root.join(&entry.output))?;
        let header = dense.byte_headers()?.clone();
        if header.get(0) != Some(b"gene_id".as_slice()) {
            bail!("dense correlation lacks a gene_id header");
        }
        let columns: Vec<Vec<u8>> = header.iter().skip(1).map(<[u8]>::to_vec).collect();
        let ids: BTreeSet<&[u8]> = columns.iter().map(Vec::as_slice).collect();
        if ids.len() != columns.len()
            || ids.len() != correlation_genes.len()
            || ids.iter().any(|id| !correlation_genes.contains(*id))
        {
            bail!("dense correlation columns differ from selected genes");
        }
        let mut rows = 0_usize;
        for row in dense.byte_records() {
            if rows.is_multiple_of(1024) {
                check()?;
            }
            let row = row?;
            if row.get(0) != columns.get(rows).map(Vec::as_slice) {
                bail!("dense correlation row IDs differ from its columns");
            }
            rows += 1;
        }
        if rows != columns.len() {
            bail!("dense correlation has the wrong number of rows");
        }
    }
    for entry in entries.iter().filter(|entry| {
        let name = Path::new(&entry.output)
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("");
        matches!(
            name,
            "modules_top_genes.tsv"
                | "modules_top_genes_annotated.tsv"
                | "module_membership.tsv"
                | "module_t_values.tsv"
                | "module_t_values_annotated.tsv"
        )
    }) {
        let mut table = reader(&root.join(&entry.output))?;
        let gene = table
            .headers()?
            .iter()
            .position(|name| name == "gene")
            .context("module table lacks gene column")?;
        for (index, row) in table.byte_records().enumerate() {
            if index.is_multiple_of(1024) {
                check()?;
            }
            let row = row?;
            if !row
                .get(gene)
                .is_some_and(|id| correlation_genes.contains(id))
            {
                bail!("module table references genes absent from the selected matrix");
            }
        }
    }
    for entry in entries
        .iter()
        .filter(|entry| entry.output == "corr/module_eigengenes.tsv")
    {
        let mut scores = reader(&root.join(&entry.output))?;
        let header = scores.headers()?;
        if header.get(0) != Some("module_id")
            || !header.iter().skip(1).eq(samples.iter().map(String::as_str))
        {
            bail!("module score sample IDs differ from frozen sample metadata");
        }
        let mut modules = BTreeSet::new();
        for row in scores.byte_records() {
            if modules.len().is_multiple_of(1024) {
                check()?;
            }
            let row = row?;
            let id = row.get(0).context("module score row lacks module ID")?;
            if id.is_empty() || !modules.insert(id.to_vec()) {
                bail!("module score IDs must be nonempty and unique");
            }
        }
    }
    check()?;
    Ok(
        json!({"schema_version":1,"sample_count":samples.len(),"gene_count":genes.len(),
        "correlation_gene_count":correlation_genes.len(),"metadata_columns":headers.iter().collect::<Vec<_>>(),
        "differential_tables":entries.iter().filter(|entry| entry.source.kind == "DE_TABLE").count(),
        "correlation_edge_tables":entries.iter().filter(|entry| entry.source.kind == "CORR_EDGES").count(),
        "correlation_dense_tables":entries.iter().filter(|entry| entry.source.kind == "CORR_DENSE").count(),
        "module_tables":entries.iter().filter(|entry| entry.source.kind == "CORR_MODULES").count(),
        "enrichment_tables":entries.iter().filter(|entry| entry.source.kind == "ENRICHMENT_TABLE").count(),
        "stratification":stratify.map(|column| json!({"column":column,"levels":levels,"missing":missing})),
        "scope":"frozen biological-sample metadata from the selected analysis generation"}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    fn fixture(root: &Path) -> Result<Vec<Entry>> {
        let mut entries = Vec::new();
        for (name, text) in [
            (
                "de/sample_metadata.tsv",
                "sample_id\tcondition\tage\nS1\tA\t21.0\nS2\tB\t32\n",
            ),
            (
                "de/gene_counts.tsv",
                "gene_id\tS1\tS2\nG1\t1\t2\nG2\t3\t4\n",
            ),
            (
                "de/gene_norm_counts.tsv",
                "gene_id\tS1\tS2\nG1\t1\t2\nG2\t3\t4\n",
            ),
            ("de/vst.tsv", "gene_id\tS1\tS2\nG1\t1\t2\nG2\t3\t4\n"),
            ("corr/adjusted_matrix.tsv", "gene_id\tS1\tS2\nG2\t3\t4\n"),
            ("de/de_test.tsv", "gene_id\tpvalue\nG1\t0.5\nG2\t0.2\n"),
            ("corr/edges_topk.tsv", "gene_a\tgene_b\tr\nG2\tG2\t1\n"),
            ("corr/corr_dense.tsv", "gene_id\tG2\nG2\t1\n"),
            ("corr/module_membership.tsv", "module_id\tgene\nM1\tG2\n"),
            (
                "corr/module_eigengenes.tsv",
                "module_id\tS1\tS2\nM1\t0.5\t1\n",
            ),
        ] {
            let path = root.join(name);
            fs::create_dir_all(path.parent().context("fixture parent")?)?;
            fs::write(&path, text)?;
            entries.push(Entry {
                generation_id: "fixture".into(),
                source: ManifestArtifact::snapshot(
                    if name.starts_with("de/de_") {
                        "DE_TABLE"
                    } else if name.contains("edges_") {
                        "CORR_EDGES"
                    } else if name.contains("dense") {
                        "CORR_DENSE"
                    } else if name.contains("module_") {
                        "CORR_MODULES"
                    } else {
                        "TEST"
                    },
                    &path,
                )?,
                output: name.into(),
                output_sha256: String::new(),
                output_bytes: 0,
                rows: None,
                columns: None,
            });
        }
        Ok(entries)
    }

    #[test]
    fn summary_checks_design_values_and_gene_identity() -> Result<()> {
        let temp = tempfile::tempdir()?;
        let entries = fixture(temp.path())?;
        let design = json!({"samples":["S1","S2"],"columns":{
            "condition":{"type":"factor","values":["A","B"]},
            "age":{"type":"numeric","values":[21.0,32.0]}
        }});
        let summary = summarize(
            temp.path(),
            &entries,
            Some("condition"),
            &design,
            &|| Ok(()),
        )?;
        assert_eq!(summary["stratification"]["levels"]["A"], 1);
        assert_eq!(summary["sample_count"], 2);
        assert_eq!(summary["gene_count"], 2);
        assert_eq!(summary["correlation_gene_count"], 1);
        assert_eq!(summary["correlation_dense_tables"], 1);
        assert_eq!(summary["module_tables"], 2);
        fs::write(
            temp.path().join("corr/module_membership.tsv"),
            "module_id\tgene\nM1\tG1\n",
        )?;
        assert!(
            format!(
                "{:#}",
                summarize(temp.path(), &entries, None, &design, &|| Ok(())).unwrap_err()
            )
            .contains("module table references genes")
        );
        fs::write(
            temp.path().join("corr/module_membership.tsv"),
            "module_id\tgene\nM1\tG2\n",
        )?;
        fs::write(
            temp.path().join("corr/module_eigengenes.tsv"),
            "module_id\tS2\tS1\nM1\t0.5\t1\n",
        )?;
        assert!(
            format!(
                "{:#}",
                summarize(temp.path(), &entries, None, &design, &|| Ok(())).unwrap_err()
            )
            .contains("module score sample IDs")
        );
        fs::write(
            temp.path().join("corr/module_eigengenes.tsv"),
            "module_id\tS1\tS2\nM1\t0.5\t1\n",
        )?;
        fs::write(
            temp.path().join("corr/corr_dense.tsv"),
            "gene_id\tG2\nG1\t1\n",
        )?;
        assert!(
            format!(
                "{:#}",
                summarize(temp.path(), &entries, None, &design, &|| Ok(())).unwrap_err()
            )
            .contains("dense correlation row IDs")
        );
        fs::write(
            temp.path().join("corr/corr_dense.tsv"),
            "gene_id\tG2\nG2\t1\n",
        )?;
        fs::write(
            temp.path().join("corr/edges_topk.tsv"),
            "gene_a\tgene_b\tr\nG1\tG2\t1\n",
        )?;
        assert!(
            format!(
                "{:#}",
                summarize(temp.path(), &entries, None, &design, &|| Ok(())).unwrap_err()
            )
            .contains("edges reference genes")
        );
        fs::write(
            temp.path().join("corr/edges_topk.tsv"),
            "gene_a\tgene_b\tr\nG2\tG2\t1\n",
        )?;
        fs::write(
            temp.path().join("de/sample_metadata.tsv"),
            "sample_id\tcondition\tage\nS1\tA\t22\nS2\tB\t32\n",
        )?;
        assert!(
            format!(
                "{:#}",
                summarize(temp.path(), &entries, None, &design, &|| Ok(())).unwrap_err()
            )
            .contains("covariate")
        );
        fs::write(
            temp.path().join("de/sample_metadata.tsv"),
            "sample_id\tcondition\tage\nS1\tA\t21\nS2\tB\t32\n",
        )?;
        fs::write(
            temp.path().join("de/vst.tsv"),
            "gene_id\tS1\tS2\nG1\t1\t2\nG3\t3\t4\n",
        )?;
        assert!(
            format!(
                "{:#}",
                summarize(temp.path(), &entries, None, &design, &|| Ok(())).unwrap_err()
            )
            .contains("gene IDs")
        );
        Ok(())
    }
}
