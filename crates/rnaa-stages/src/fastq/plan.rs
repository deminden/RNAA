use super::{FastqLimits, FastqSummary, PreparedFastqs, layout, reader};
use anyhow::{Context, Result, bail};
use rnaa_core::{LibraryLayout, ManifestArtifact, model::VerifiedFile};
use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

#[derive(Debug)]
pub struct FastqPlan {
    layout: layout::Layout,
    limits: FastqLimits,
    inputs: Vec<ManifestArtifact>,
    summary: FastqSummary,
}
impl FastqPlan {
    /// Validate and count the full decoded group without creating temporary files.
    pub fn inspect(
        layout: LibraryLayout,
        files: &[VerifiedFile],
        limits: FastqLimits,
        check: &impl Fn() -> Result<()>,
    ) -> Result<Self> {
        limits.validate()?;
        check()?;
        let mut group = layout::inspect(layout, files)?;
        let inputs = files
            .iter()
            .map(|file| ManifestArtifact::snapshot(file.artifact_kind.as_str(), &file.path))
            .collect::<Result<Vec<_>>>()?;
        let summary = scan(&group, limits, None, None, check)?;
        group.staging = summary.needs_staging;
        let value = Self {
            layout: group,
            limits,
            inputs,
            summary,
        };
        value.verify_sources(check)?;
        Ok(value)
    }
    pub fn summary(&self) -> &FastqSummary {
        &self.summary
    }
    pub fn inputs(&self) -> &[ManifestArtifact] {
        &self.inputs
    }
    pub fn verify_sources(&self, check: &impl Fn() -> Result<()>) -> Result<()> {
        for input in &self.inputs {
            check()?;
            if !input.matches_file()? {
                bail!(
                    "FASTQ source changed after decoded preflight: {}",
                    input.path
                );
            }
        }
        Ok(())
    }
    /// Materialize within the preflight byte promise. Resource schedulers can
    /// reserve summary.staging_bytes() at this parent filesystem before calling.
    pub fn prepare_in(
        &self,
        parent: &Path,
        check: &impl Fn() -> Result<()>,
    ) -> Result<PreparedFastqs> {
        self.verify_sources(check)?;
        let staging = if self.layout.staging {
            Some(
                tempfile::Builder::new()
                    .prefix("rnaa-fastq-")
                    .tempdir_in(parent)?,
            )
        } else {
            None
        };
        let first = staging
            .as_ref()
            .map(|directory| directory.path().join("read1.fastq"))
            .unwrap_or_else(|| self.layout.lanes[0].first.clone());
        let second = self.layout.lanes[0].second.as_ref().map(|source| {
            staging
                .as_ref()
                .map(|directory| directory.path().join("read2.fastq"))
                .unwrap_or_else(|| source.clone())
        });
        if staging.is_some() {
            let mut limits = self.limits;
            limits.max_decoded_bytes = Some(self.summary.decoded_bytes);
            limits.max_staging_bytes = Some(self.summary.staging_bytes()?);
            let mut a = BufWriter::new(File::create(&first)?);
            let mut b = second
                .as_ref()
                .map(|path| File::create(path).map(BufWriter::new))
                .transpose()?;
            let observed = scan(
                &self.layout,
                limits,
                Some(&mut a),
                b.as_mut().map(|value| value as &mut dyn Write),
                check,
            )?;
            a.flush()?;
            if let Some(b) = &mut b {
                b.flush()?;
            }
            if observed != self.summary {
                bail!("FASTQ decoded shape changed after preflight");
            }
            self.verify_sources(check)?;
        }
        Ok(PreparedFastqs {
            read1: first,
            read2: second,
            read1_kind: self.layout.first_kind,
            summary: self.summary.clone(),
            _staging: staging,
        })
    }
}

fn scan(
    layout: &layout::Layout,
    limits: FastqLimits,
    mut first: Option<&mut dyn Write>,
    mut second: Option<&mut dyn Write>,
    check: &impl Fn() -> Result<()>,
) -> Result<FastqSummary> {
    let mut control = reader::ScanControl::new(check);
    let mut summary = FastqSummary {
        lanes: layout.lanes.len(),
        paired: layout.lanes[0].second.is_some(),
        needs_staging: layout.staging,
        min_record_bytes: u64::MAX,
        ..Default::default()
    };
    for lane in &layout.lanes {
        let mut a = reader::open(&lane.first)?;
        let mut b = lane
            .second
            .as_ref()
            .map(|path| reader::open(path))
            .transpose()?;
        let mut count = 0_u64;
        loop {
            let left = reader::record(&mut *a, limits, &mut summary, &mut control)
                .with_context(|| format!("invalid FASTQ {}", lane.first.display()))?;
            let right = b
                .as_mut()
                .map(|reader| reader::record(&mut **reader, limits, &mut summary, &mut control))
                .transpose()?
                .flatten();
            match (left, right) {
                (None, None) => break,
                (Some(a), b) if !summary.paired || b.is_some() => {
                    if let Some(b) = &b
                        && reader::read_id(&a.lines[0]) != reader::read_id(&b.lines[0])
                    {
                        bail!("paired FASTQ read identifiers differ");
                    }
                    summary.fragments = summary
                        .fragments
                        .checked_add(1)
                        .context("FASTQ fragment count overflows")?;
                    summary.read_records = summary
                        .read_records
                        .checked_add(if summary.paired { 2 } else { 1 })
                        .context("FASTQ record count overflows")?;
                    summary.r1_normalized_bytes = summary
                        .r1_normalized_bytes
                        .checked_add(a.normalized_bytes)
                        .context("normalized FASTQ byte count overflows")?;
                    if let Some(b) = &b {
                        summary.r2_normalized_bytes = summary
                            .r2_normalized_bytes
                            .checked_add(b.normalized_bytes)
                            .context("normalized FASTQ byte count overflows")?;
                    }
                    let normalized = summary.normalized_bytes()?;
                    if summary.needs_staging
                        && limits
                            .max_staging_bytes
                            .is_some_and(|limit| normalized > limit)
                    {
                        bail!(
                            "decoded FASTQ staging exceeds its byte limit {}",
                            limits.max_staging_bytes.unwrap()
                        );
                    }
                    if let Some(out) = &mut first {
                        for line in &a.lines {
                            out.write_all(line)?;
                        }
                    }
                    if let (Some(out), Some(b)) = (&mut second, b) {
                        for line in &b.lines {
                            out.write_all(line)?;
                        }
                    }
                    count += 1;
                }
                _ => bail!("paired FASTQ lane contains unequal mate counts"),
            }
        }
        if count == 0 {
            bail!("FASTQ lane {} contains no reads", lane.first.display());
        }
    }
    Ok(summary)
}
