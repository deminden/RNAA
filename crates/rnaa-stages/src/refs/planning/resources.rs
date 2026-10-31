//! Structural phase/output estimates, bounded staging and native shape checks.
use super::*;

impl Plan {
    pub fn stage_fasta(
        &self,
        source: &Path,
        target: &Path,
        check: &impl Fn() -> Result<()>,
    ) -> Result<()> {
        use std::io::Write;
        let mut reader = open(source)?.take(self.fasta.decoded_bytes);
        let mut output = File::create(target)?;
        let mut buffer = [0_u8; 65536];
        let (mut written, mut next) = (0_u64, 0_u64);
        loop {
            if written >= next {
                check()?;
                next = written.saturating_add(8 * MIB);
            }
            let n = reader.read(&mut buffer)?;
            if n == 0 {
                break;
            }
            output.write_all(&buffer[..n])?;
            written += n as u64;
        }
        if written != self.fasta.decoded_bytes || reader.into_inner().read(&mut [0_u8; 1])? != 0 {
            bail!("decoded FASTA changed beyond its staging promise");
        }
        drop(output);
        if rnaa_core::util::compute_sha256(target)? != self.fasta.decoded_sha256 {
            bail!("decoded FASTA changed after planning");
        }
        check()
    }
    pub fn validate_outputs(
        &self,
        bundle: &rnaa_core::ReferenceBundle,
        report: &kallistors::index::IndexBuildReport,
    ) -> Result<()> {
        if report.transcripts as u64 != self.fasta.transcripts
            || report.kmers as u64 > self.fasta.unique_kmers_upper
            || report.unitigs as u64 > self.fasta.unique_kmers_upper
        {
            bail!("native reference build exceeded its planned graph shape");
        }
        if std::fs::metadata(&bundle.kallisto_index_path)?.len() > self.index_bytes {
            bail!("native reference index exceeded its serialized output promise");
        }
        let table_bytes = std::fs::metadata(&bundle.tx2gene_path)?
            .len()
            .checked_add(
                std::fs::metadata(
                    bundle
                        .gene_annotation_path
                        .as_ref()
                        .context("missing annotation output")?,
                )?
                .len(),
            )
            .context("reference table bytes overflow")?;
        if table_bytes > self.table_bytes {
            bail!("reference tables exceeded their output promise");
        }
        Ok(())
    }
    pub fn inspect(
        cdna: &ManifestArtifact,
        gtf: &ManifestArtifact,
        config: &RefsConfig,
        context: &ExecutionResources,
        workspace: &Path,
        check: &impl Fn() -> Result<()>,
    ) -> Result<Self> {
        let limits = Limits::try_from(config)?;
        let scan_memory = (limits.line as u64)
            .checked_mul(16)
            .and_then(|n| n.checked_add(64 * MIB))
            .context("reference scan memory overflows")?;
        let mut unique = UniqueWords::new(context, scan_memory, workspace, check)?;
        let mut fasta = scan_fasta(
            Path::new(&cdna.path),
            limits,
            config.index_k,
            &mut unique,
            check,
        )?;
        fasta.clean_unique_kmers = unique.finish()?;
        fasta.unique_kmers_upper = fasta
            .clean_unique_kmers
            .saturating_add(fasta.uncertain_kmers)
            .min(fasta.kmer_positions);
        if fasta.unique_kmers_upper > u32::MAX as u64 {
            bail!("reference graph upper bound exceeds the engine's node width");
        }
        let gtf_limits = Limits {
            decoded: limits
                .decoded
                .map(|n| n.saturating_sub(fasta.decoded_bytes)),
            ..limits
        };
        let guard = context.reserve(
            "gtf_plan",
            ResourceRequest {
                cpu_threads: 1,
                memory_bytes: scan_memory,
                ..Default::default()
            },
            check,
        )?;
        let gtf_plan = scan_gtf(Path::new(&gtf.path), gtf_limits, check)?;
        if !cdna.matches_file()? || !gtf.matches_file()? {
            bail!("reference source changed during decoded planning");
        }
        drop(guard);
        let u = u128::from(fasta.unique_kmers_upper);
        let p = u128::from(fasta.kmer_positions);
        let t = u128::from(fasta.transcripts);
        let names = u128::from(fasta.name_bytes);
        let rows = u128::from(gtf_plan.rows);
        let strings = u128::from(gtf_plan.string_bytes);
        let threads = config
            .index_threads
            .min(usize::try_from(context.limits().cpu_threads).unwrap_or(usize::MAX))
            .max(1);
        let inventory = serde_json::to_vec(&(cdna, gtf))?.len() as u128;
        let metadata = 16 * u128::from(MIB) + inventory * 32;
        let native = 64 * u128::from(MIB)
            + u * (512 + config.index_k as u128 * 2)
            + p * 64
            + u128::from(fasta.decoded_bytes) * 4
            + t * 1024
            + names * 8
            + (threads as u128) * 4 * u128::from(MIB);
        let annotation = 64 * u128::from(MIB)
            + rows * 512
            + strings * 8
            + t * 512
            + names * 8
            + gtf_plan.max_line as u128 * 16;
        Ok(Self {
            model: "kallistors-reference-build-v1",
            staging_bytes: if crate::refs::source::is_gzip(Path::new(&cdna.path))? {
                fasta.decoded_bytes
            } else {
                0
            },
            memory_bytes: u64::try_from(
                native.max(annotation).max(u128::from(scan_memory)) + metadata * 4,
            )?,
            copy_bytes: cdna
                .bytes
                .checked_add(gtf.bytes)
                .context("reference copy bytes overflow")?,
            index_bytes: u64::try_from(
                64 * u128::from(MIB)
                    + u * (config.index_k as u128 + 256)
                    + p * 64
                    + t * 256
                    + names * 8,
            )?,
            table_bytes: u64::try_from(65536 + strings * 4 + rows * 128)?,
            metadata_bytes: u64::try_from(metadata)?,
            threads,
            fasta,
            gtf: gtf_plan,
        })
    }
    pub fn request(
        &self,
        generations: &Path,
        catalog: &Path,
        workspace: &Path,
    ) -> Result<ResourceRequest> {
        let future = self
            .copy_bytes
            .checked_add(self.index_bytes)
            .and_then(|n| n.checked_add(self.table_bytes))
            .and_then(|n| n.checked_add(self.metadata_bytes))
            .context("reference disk promises overflow")?;
        let mut disks = vec![
            DiskRequest {
                path: generations.into(),
                bytes: future,
            },
            DiskRequest {
                path: catalog.into(),
                bytes: self.metadata_bytes,
            },
        ];
        if self.staging_bytes > 0 {
            disks.push(DiskRequest {
                path: workspace.into(),
                bytes: self.staging_bytes,
            });
        }
        Ok(ResourceRequest {
            cpu_threads: self.threads as u64,
            memory_bytes: self.memory_bytes,
            disks,
        })
    }
}
