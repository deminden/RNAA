use super::*;

pub(super) fn parse_row_map(project_id: &str, row: Map<String, Value>) -> Result<RunRecord> {
    let run_accession = required_string(&row, "run_accession")?;
    let primary_study = opt_string(&row, "study_accession");
    let secondary_study = opt_string(&row, "secondary_study_accession");
    let study_accession = secondary_study
        .clone()
        .or(primary_study.clone())
        .ok_or_else(|| anyhow!("ENA row missing study accession for {}", run_accession))?;
    let source_study_accession = match (&primary_study, &secondary_study) {
        (Some(primary), Some(_)) => Some(primary.clone()),
        _ => None,
    };
    let library_layout = opt_string(&row, "library_layout")
        .as_deref()
        .map(LibraryLayout::from_str)
        .transpose()?
        .unwrap_or_default();
    let fastq_files = expand_remote_files(
        &opt_string(&row, "fastq_ftp").unwrap_or_default(),
        &opt_string(&row, "fastq_md5").unwrap_or_default(),
        &opt_string(&row, "fastq_bytes").unwrap_or_default(),
        library_layout,
        false,
    );
    let sra_files = expand_remote_files(
        &opt_string(&row, "submitted_ftp").unwrap_or_default(),
        &opt_string(&row, "submitted_md5").unwrap_or_default(),
        &opt_string(&row, "submitted_bytes").unwrap_or_default(),
        library_layout,
        true,
    );
    let remote_files = if fastq_files.is_empty() {
        sra_files
    } else {
        fastq_files.into_iter().chain(sra_files).collect()
    };

    let metadata = row
        .iter()
        .filter_map(|(key, value)| {
            if matches!(
                key.as_str(),
                "run_accession"
                    | "study_accession"
                    | "secondary_study_accession"
                    | "sample_accession"
                    | "experiment_accession"
                    | "library_layout"
                    | "instrument_platform"
                    | "instrument_model"
                    | "fastq_ftp"
                    | "fastq_md5"
                    | "fastq_bytes"
                    | "submitted_ftp"
                    | "submitted_md5"
                    | "submitted_bytes"
            ) {
                return None;
            }
            let stringified = json_value_to_string(value);
            (!stringified.is_empty()).then(|| (format!("ena_{key}"), stringified))
        })
        .chain(
            primary_study
                .iter()
                .map(|value| ("ena_primary_study_accession".to_string(), value.to_string())),
        )
        .chain(secondary_study.iter().map(|value| {
            (
                "ena_secondary_study_accession".to_string(),
                value.to_string(),
            )
        }))
        .collect::<BTreeMap<_, _>>();

    Ok(RunRecord {
        project_id: project_id.to_string(),
        study_accession,
        source_study_accession,
        run_accession,
        sample_accession: opt_string(&row, "sample_accession")
            .or_else(|| opt_string(&row, "secondary_sample_accession")),
        experiment_accession: opt_string(&row, "experiment_accession"),
        library_layout,
        instrument_platform: opt_string(&row, "instrument_platform"),
        instrument_model: opt_string(&row, "instrument_model"),
        remote_files,
        metadata,
        state: rnaa_core::state::RunState::Resolved,
        last_error: None,
        updated_at: now_rfc3339(),
    })
}

fn expand_remote_files(
    urls_text: &str,
    md5_text: &str,
    bytes_text: &str,
    layout: LibraryLayout,
    submitted: bool,
) -> Vec<RemoteFile> {
    let urls = split_field(urls_text);
    let md5s = split_field(md5_text);
    let bytes = split_field(bytes_text);

    urls.into_iter()
        .enumerate()
        .filter_map(|(index, url)| {
            if url.is_empty() {
                return None;
            }
            let normalized = normalize_url(&url);
            let basename = normalized
                .rsplit('/')
                .next()
                .map(sanitize_basename)
                .unwrap_or_else(|| format!("file_{index}"));
            let kind = if submitted && basename.ends_with(".sra") {
                RemoteFileKind::Sra
            } else if layout == LibraryLayout::Paired {
                if basename.contains("_1.") || basename.contains("R1") {
                    RemoteFileKind::FastqR1
                } else if basename.contains("_2.") || basename.contains("R2") {
                    RemoteFileKind::FastqR2
                } else if submitted {
                    RemoteFileKind::Submitted
                } else {
                    RemoteFileKind::Fastq
                }
            } else if submitted && basename.ends_with(".sra") {
                RemoteFileKind::Sra
            } else {
                RemoteFileKind::Fastq
            };

            Some(RemoteFile {
                url: normalized,
                md5: md5s.get(index).cloned().filter(|value| !value.is_empty()),
                bytes: bytes.get(index).and_then(|value| value.parse::<u64>().ok()),
                kind,
                basename,
            })
        })
        .collect()
}

fn split_field(value: &str) -> Vec<String> {
    value
        .split(';')
        .map(str::trim)
        .filter(|field| !field.is_empty())
        .map(ToString::to_string)
        .collect()
}

fn required_string(row: &Map<String, Value>, key: &str) -> Result<String> {
    opt_string(row, key).ok_or_else(|| anyhow!("ENA row missing required field {key}"))
}

fn opt_string(row: &Map<String, Value>, key: &str) -> Option<String> {
    row.get(key)
        .map(json_value_to_string)
        .filter(|value| !value.is_empty())
}

fn json_value_to_string(value: &Value) -> String {
    match value {
        Value::Null => String::new(),
        Value::String(text) => text.trim().to_string(),
        Value::Number(number) => number.to_string(),
        Value::Bool(boolean) => boolean.to_string(),
        other => other.to_string(),
    }
}
