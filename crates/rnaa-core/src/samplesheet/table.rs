//! Deterministic run-keyed three-way merging of editable sample metadata.
use super::publication::ConflictPolicy;
use anyhow::{Context, Result, bail};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone)]
pub(super) struct Table {
    headers: Vec<String>,
    rows: BTreeMap<String, BTreeMap<String, String>>,
}
impl Table {
    pub(super) fn overlay(&self, newer: &Self) -> Self {
        let mut combined = self.clone();
        combined.headers.extend(
            newer
                .headers
                .iter()
                .filter(|name| !self.headers.contains(name))
                .cloned(),
        );
        combined.rows.extend(newer.rows.clone());
        combined
    }
    pub(super) fn project(&self, active: &Self) -> Result<Vec<u8>> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());
        writer.write_record(&self.headers)?;
        for id in active.rows.keys() {
            let row = self
                .rows
                .get(id)
                .context("active metadata row missing from merged history")?;
            writer.write_record(
                self.headers
                    .iter()
                    .map(|name| row.get(name).map(String::as_str).unwrap_or_default()),
            )?;
        }
        writer.flush()?;
        Ok(writer.into_inner()?)
    }
    pub(super) fn parse(bytes: &[u8]) -> Result<Self> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(bytes);
        let headers = reader
            .headers()?
            .iter()
            .map(str::to_string)
            .collect::<Vec<_>>();
        if headers
            .iter()
            .any(|name| name.trim().is_empty() || name.trim() != name)
            || headers.iter().collect::<BTreeSet<_>>().len() != headers.len()
        {
            bail!("sample metadata columns must be unique, trimmed and nonempty");
        }
        let key = headers
            .iter()
            .position(|name| name == "run_accession")
            .context("sample metadata missing run_accession")?;
        let mut rows = BTreeMap::new();
        for row in reader.records() {
            let row = row?;
            let id = row[key].trim();
            if id.is_empty() || rows.contains_key(id) {
                bail!("sample metadata has duplicate or empty run ID '{id}'");
            }
            rows.insert(
                id.into(),
                headers
                    .iter()
                    .cloned()
                    .zip(row.iter().map(str::to_string))
                    .collect(),
            );
        }
        Ok(Self { headers, rows })
    }
    pub(super) fn merge(
        base: &Self,
        local: &Self,
        remote: &Self,
        policy: ConflictPolicy,
    ) -> Result<Vec<u8>> {
        let mut headers = remote.headers.clone();
        headers.extend(
            local
                .headers
                .iter()
                .filter(|name| !remote.headers.contains(name))
                .cloned(),
        );
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(Vec::new());
        writer.write_record(&headers)?;
        let mut conflicts = Vec::new();
        for (run, fresh) in &remote.rows {
            let empty = BTreeMap::new();
            let previous = base.rows.get(run).unwrap_or(&empty);
            let edited = local.rows.get(run).unwrap_or(&empty);
            let mut row = Vec::new();
            for name in &headers {
                let old = previous.get(name);
                let ours = edited.get(name);
                let theirs = fresh.get(name);
                let chosen = if operational(name) || ours.is_none() {
                    theirs
                } else if protected(name) {
                    if ours != old && ours != theirs {
                        bail!(
                            "sample metadata identity '{name}' for {run} was edited; use a registered override and resolve again"
                        );
                    }
                    theirs
                } else if ours == old || ours == theirs {
                    theirs
                } else if theirs == old {
                    ours
                } else {
                    match policy {
                        ConflictPolicy::Fail => {
                            conflicts.push(format!("{run}/{name}"));
                            ours
                        }
                        ConflictPolicy::KeepLocal => ours,
                        ConflictPolicy::UseSource => theirs,
                    }
                };
                row.push(chosen.cloned().unwrap_or_default());
            }
            writer.write_record(&row)?;
        }
        if !conflicts.is_empty() {
            bail!(
                "metadata refresh conflicts in {}; samplesheet and catalog are unchanged. Review retained versions; retry resolve with --metadata-conflicts keep-local/use-source, or metadata refresh with --conflicts keep-local/use-source",
                conflicts.join(", ")
            );
        }
        writer.flush()?;
        Ok(writer.into_inner()?)
    }
}
fn operational(name: &str) -> bool {
    matches!(name, "read1_path" | "read2_path")
}
fn protected(name: &str) -> bool {
    matches!(
        name,
        "project_id"
            | "run_accession"
            | "study_accession"
            | "sample_accession"
            | "experiment_accession"
            | "library_layout"
            | "instrument_platform"
            | "instrument_model"
            | "biological_sample_id"
            | "identity_source"
    )
}
