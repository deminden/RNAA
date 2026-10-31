use super::*;

#[derive(Debug, Clone)]
pub(super) struct OverrideRow {
    run_accession: Option<String>,
    sample_accession: Option<String>,
    values: BTreeMap<String, String>,
    merge_strategy: MetadataMergeStrategy,
}

pub(super) fn load_override_tables(
    overrides: &[MetadataOverrideRecord],
) -> Result<Vec<OverrideRow>> {
    let mut rows = Vec::new();
    for override_record in overrides {
        let path = Path::new(&override_record.stored_path);
        rows.extend(load_override_table(path, override_record.merge_strategy)?);
    }
    Ok(rows)
}

fn load_override_table(
    path: &Path,
    merge_strategy: MetadataMergeStrategy,
) -> Result<Vec<OverrideRow>> {
    let content =
        fs::read_to_string(path).with_context(|| format!("failed to read {}", path.display()))?;
    let delimiter = if path.extension().and_then(|ext| ext.to_str()) == Some("csv") {
        b','
    } else if content
        .lines()
        .next()
        .is_some_and(|line| line.contains('\t'))
    {
        b'\t'
    } else {
        b','
    };

    let mut reader = ReaderBuilder::new()
        .delimiter(delimiter)
        .from_reader(content.as_bytes());
    let headers = reader.headers()?.clone();
    let run_idx = headers.iter().position(|name| name == "run_accession");
    let sample_idx = headers.iter().position(|name| name == "sample_accession");
    if run_idx.is_none() && sample_idx.is_none() {
        bail!(
            "metadata override {} must contain run_accession or sample_accession",
            path.display()
        );
    }

    let mut overrides = Vec::new();
    for record in reader.records() {
        let record = record?;
        let mut values = BTreeMap::new();
        for (index, field) in record.iter().enumerate() {
            let header = headers.get(index).unwrap_or_default();
            if field.trim().is_empty() || header == "run_accession" || header == "sample_accession"
            {
                continue;
            }
            values.insert(header.to_string(), field.to_string());
        }
        overrides.push(OverrideRow {
            run_accession: run_idx
                .and_then(|index| record.get(index))
                .map(ToString::to_string),
            sample_accession: sample_idx
                .and_then(|index| record.get(index))
                .map(ToString::to_string),
            values,
            merge_strategy,
        });
    }
    Ok(overrides)
}

pub(super) fn apply_overrides(
    run: &mut RunRecord,
    overrides: &[OverrideRow],
    _default_strategy: MetadataMergeStrategy,
) {
    for override_row in overrides.iter().filter(|override_row| {
        override_row.run_accession.as_deref() == Some(run.run_accession.as_str())
            || override_row.sample_accession.as_deref() == run.sample_accession.as_deref()
    }) {
        merge_metadata(
            &mut run.metadata,
            &override_row.values,
            override_row.merge_strategy,
        );
    }
}

pub fn merge_metadata(
    target: &mut BTreeMap<String, String>,
    incoming: &BTreeMap<String, String>,
    strategy: MetadataMergeStrategy,
) {
    for (key, value) in incoming {
        match strategy {
            MetadataMergeStrategy::Override | MetadataMergeStrategy::PreferLocal => {
                target.insert(key.clone(), value.clone());
            }
            MetadataMergeStrategy::PreferRemote => {
                target.entry(key.clone()).or_insert_with(|| value.clone());
            }
            MetadataMergeStrategy::UnionWithPrefix => {
                if let Some(existing) = target.get(key).cloned() {
                    if existing != *value {
                        target.insert(format!("remote_{key}"), existing);
                        target.insert(format!("local_{key}"), value.clone());
                    }
                } else {
                    target.insert(key.clone(), value.clone());
                }
            }
        }
    }
}
