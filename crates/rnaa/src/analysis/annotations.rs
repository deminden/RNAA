//! Gene annotation sidecars for analysis outputs.
use anyhow::Context;
use anyhow::Result;
use std::path::{Path, PathBuf};

fn load_gene_annotation_map(
    path: &Path,
) -> Result<std::collections::HashMap<String, (String, String)>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to read gene annotation {}", path.display()))?;
    let headers = reader
        .headers()
        .with_context(|| format!("failed to parse gene annotation header {}", path.display()))?
        .iter()
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let gene_id_idx = headers
        .iter()
        .position(|col| col == "gene_id")
        .ok_or_else(|| anyhow::anyhow!("gene annotation missing gene_id column"))?;
    let gene_name_idx = headers
        .iter()
        .position(|col| col == "gene_name")
        .ok_or_else(|| anyhow::anyhow!("gene annotation missing gene_name column"))?;
    let biotype_idx = headers
        .iter()
        .position(|col| col == "gene_biotype")
        .ok_or_else(|| anyhow::anyhow!("gene annotation missing gene_biotype column"))?;
    let mut map = std::collections::HashMap::new();
    for row in reader.records() {
        let row = row?;
        let gene_id = normalize_stable_gene_id(row.get(gene_id_idx).unwrap_or_default());
        if gene_id.is_empty() {
            continue;
        }
        let gene_name = row.get(gene_name_idx).unwrap_or_default().to_string();
        let biotype = row.get(biotype_idx).unwrap_or_default().to_string();
        map.entry(gene_id).or_insert((gene_name, biotype));
    }
    Ok(map)
}

fn normalize_stable_gene_id(value: &str) -> String {
    rnaa_core::util::normalize_feature_id(value)
}

pub(crate) fn annotate_corr_outputs(
    outputs: &[PathBuf],
    annotation_path: &Path,
) -> Result<Vec<PathBuf>> {
    let annotations = load_gene_annotation_map(annotation_path)?;
    let mut annotated_outputs = Vec::new();
    for output in outputs {
        let Some(file_name) = output.file_name().and_then(|v| v.to_str()) else {
            continue;
        };
        if file_name == "edges_topk.tsv" || file_name == "edges_threshold.tsv" {
            let annotated_path = output.with_file_name(file_name.replace(".tsv", "_annotated.tsv"));
            annotate_edges_table(output, &annotated_path, &annotations)?;
            annotated_outputs.push(annotated_path);
        } else if file_name == "modules_top_genes.tsv" {
            let annotated_path = output.with_file_name("modules_top_genes_annotated.tsv");
            annotate_single_gene_table(
                output,
                &annotated_path,
                "gene",
                &annotations,
                Some(("gene_symbol", "gene_biotype")),
            )?;
            annotated_outputs.push(annotated_path);
        } else if file_name == "module_t_values.tsv" {
            let annotated_path = output.with_file_name("module_t_values_annotated.tsv");
            annotate_single_gene_table(
                output,
                &annotated_path,
                "gene",
                &annotations,
                Some(("gene_symbol", "gene_biotype")),
            )?;
            annotated_outputs.push(annotated_path);
        }
    }
    Ok(annotated_outputs)
}

fn annotate_edges_table(
    input: &Path,
    output: &Path,
    annotations: &std::collections::HashMap<String, (String, String)>,
) -> Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input)
        .with_context(|| format!("failed to read {}", input.display()))?;
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    writer.write_record([
        "gene_a_id",
        "gene_a_symbol",
        "gene_a_biotype",
        "gene_b_id",
        "gene_b_symbol",
        "gene_b_biotype",
        "r",
    ])?;
    for row in reader.records() {
        let row = row?;
        let gene_a = row.get(0).unwrap_or_default();
        let gene_b = row.get(1).unwrap_or_default();
        let r = row.get(2).unwrap_or_default();
        let (gene_a_symbol, gene_a_biotype) = annotations
            .get(gene_a)
            .cloned()
            .unwrap_or_else(|| (String::new(), String::new()));
        let (gene_b_symbol, gene_b_biotype) = annotations
            .get(gene_b)
            .cloned()
            .unwrap_or_else(|| (String::new(), String::new()));
        writer.write_record([
            gene_a,
            gene_a_symbol.as_str(),
            gene_a_biotype.as_str(),
            gene_b,
            gene_b_symbol.as_str(),
            gene_b_biotype.as_str(),
            r,
        ])?;
    }
    writer.flush()?;
    Ok(())
}

fn annotate_single_gene_table(
    input: &Path,
    output: &Path,
    gene_column: &str,
    annotations: &std::collections::HashMap<String, (String, String)>,
    extra_columns: Option<(&str, &str)>,
) -> Result<()> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(input)
        .with_context(|| format!("failed to read {}", input.display()))?;
    let headers = reader.headers()?.clone();
    let gene_idx = headers
        .iter()
        .position(|col| col == gene_column)
        .ok_or_else(|| anyhow::anyhow!("{} missing {} column", input.display(), gene_column))?;
    let mut out_headers = headers.iter().map(ToString::to_string).collect::<Vec<_>>();
    if let Some((symbol_col, biotype_col)) = extra_columns {
        out_headers.insert(gene_idx + 1, symbol_col.to_string());
        out_headers.insert(gene_idx + 2, biotype_col.to_string());
    }
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_path(output)
        .with_context(|| format!("failed to create {}", output.display()))?;
    writer.write_record(out_headers)?;
    for row in reader.records() {
        let row = row?;
        let gene_id = row.get(gene_idx).unwrap_or_default();
        let (symbol, biotype) = annotations
            .get(gene_id)
            .cloned()
            .unwrap_or_else(|| (String::new(), String::new()));
        let mut out = row.iter().map(ToString::to_string).collect::<Vec<_>>();
        if extra_columns.is_some() {
            out.insert(gene_idx + 1, symbol);
            out.insert(gene_idx + 2, biotype);
        }
        writer.write_record(out)?;
    }
    writer.flush()?;
    Ok(())
}
