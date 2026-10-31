//! Export sealed analysis generations as independent, provenance-bound tables.
mod admission;
mod cohort;
mod design;
mod publication;
mod recovery;
mod request_guard;
mod survey;
use crate::cli::{ExportArgs, SurveyArgs};
use anyhow::{Context, Result, bail};
use rnaa_core::{
    Database, ManifestArtifact, ProjectConfig, ProjectPaths,
    generations::{GenerationHead, OutputGroup},
    jobs::{JobLease, JobOutcome},
    resources::{DiskRequest, ExecutionResources, ReservationOwner, ResourceRequest},
};
use serde::Serialize;
use std::{
    fs,
    path::{Path, PathBuf},
};

const MIB: u64 = 1024 * 1024;

#[derive(Clone, Copy)]
enum Selection {
    All,
    Counts,
    Normalized,
    Vst,
    De,
    Correlation,
}
impl Selection {
    fn parse(value: Option<&str>) -> Result<Self> {
        Ok(match value.unwrap_or("all") {
            "all" => Self::All,
            "counts" => Self::Counts,
            "normalized" => Self::Normalized,
            "vst" => Self::Vst,
            "de" => Self::De,
            "correlation" => Self::Correlation,
            other => bail!(
                "unsupported export selection {other}; choose all, counts, normalized, vst, de or correlation"
            ),
        })
    }
    fn includes(self, kind: &str) -> bool {
        match self {
            Self::All => matches!(
                kind,
                "COUNTS"
                    | "NORM_COUNTS"
                    | "VST"
                    | "DE_TABLE"
                    | "SIZE_FACTORS"
                    | "SAMPLESHEET"
                    | "CORR_EDGES"
                    | "CORR_DENSE"
                    | "CORR_MODULES"
                    | "ENRICHMENT_TABLE"
                    | "ADJUSTED_MATRIX"
            ),
            Self::Counts => kind == "COUNTS",
            Self::Normalized => kind == "NORM_COUNTS",
            Self::Vst => kind == "VST",
            Self::De => kind == "DE_TABLE",
            Self::Correlation => matches!(
                kind,
                "CORR_EDGES"
                    | "CORR_DENSE"
                    | "CORR_MODULES"
                    | "ENRICHMENT_TABLE"
                    | "ADJUSTED_MATRIX"
            ),
        }
    }
    fn groups(self) -> &'static [OutputGroup] {
        match self {
            Self::Correlation => &[OutputGroup::Corr],
            Self::All => &[OutputGroup::De, OutputGroup::Corr],
            _ => &[OutputGroup::De],
        }
    }
}

#[derive(Serialize)]
struct Entry {
    generation_id: String,
    source: ManifestArtifact,
    output: String,
    output_sha256: String,
    output_bytes: u64,
    rows: Option<u64>,
    columns: Option<u64>,
}

#[derive(Clone, Copy)]
enum Mode<'a> {
    Export,
    Survey { stratify: Option<&'a str> },
}

struct Request<'a> {
    stage: &'static str,
    selection: Selection,
    format: &'static str,
    delimiter: u8,
    mode: Mode<'a>,
    require_current_cohort: bool,
    require_current_admission: bool,
    require_current_design: bool,
    require_current: bool,
}

pub(crate) fn run(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    args: ExportArgs,
) -> Result<()> {
    let selection = Selection::parse(args.what.as_deref())?;
    let (format, delimiter) = match args.format.as_deref().unwrap_or("tsv") {
        "tsv" => ("tsv", b'\t'),
        "csv" => ("csv", b','),
        other => bail!("unsupported export format {other}; choose tsv or csv"),
    };
    run_request(
        paths,
        db,
        config,
        Request {
            stage: "export",
            selection,
            format,
            delimiter,
            mode: Mode::Export,
            require_current_cohort: args.require_current_cohort
                || args.require_current_admission
                || args.require_current_design
                || args.require_current,
            require_current_admission: args.require_current_admission
                || args.require_current_design
                || args.require_current,
            require_current_design: args.require_current_design || args.require_current,
            require_current: args.require_current,
        },
    )
}

pub(crate) fn survey(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    args: SurveyArgs,
) -> Result<()> {
    run_request(
        paths,
        db,
        config,
        Request {
            stage: "survey",
            selection: Selection::All,
            format: "tsv",
            delimiter: b'\t',
            mode: Mode::Survey {
                stratify: args.stratify.as_deref(),
            },
            require_current_cohort: args.require_current_cohort
                || args.require_current_admission
                || args.require_current_design
                || args.require_current,
            require_current_admission: args.require_current_admission
                || args.require_current_design
                || args.require_current,
            require_current_design: args.require_current_design || args.require_current,
            require_current: args.require_current,
        },
    )
}

fn run_request(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    request: Request<'_>,
) -> Result<()> {
    let mut lease = JobLease::try_acquire(db, request.stage, None, &[])?
        .context("another project analysis or export is active")?;
    let owned = lease.database();
    let result = recovery::reconcile(paths, &owned)
        .and_then(|()| execute(paths, &owned, config, &request, lease.attempt_id()));
    match result {
        Ok(path) => {
            lease.finish(JobOutcome::Succeeded, None)?;
            println!("{}\t{}", request.stage, path.display());
            println!(
                "{}_scope\t{}",
                request.stage,
                if request.require_current {
                    "current_scientific_requests_reconstructed_supported_inputs_only"
                } else if request.require_current_design {
                    "sealed_generations_current_design_checked_other_settings_not_revalidated"
                } else if request.require_current_admission {
                    "sealed_generations_current_admission_checked_other_settings_not_revalidated"
                } else if request.require_current_cohort {
                    "sealed_generations_current_cohort_checked_other_settings_not_revalidated"
                } else {
                    "sealed_generations_not_revalidated_against_current_settings"
                }
            );
            Ok(())
        }
        Err(error) => {
            lease.finish(JobOutcome::Failed, Some(&format!("{error:#}")))?;
            Err(error)
        }
    }
}

fn execute(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    request: &Request<'_>,
    attempt: &str,
) -> Result<PathBuf> {
    let check = || db.check_job();
    let resources = ExecutionResources::open(
        &config.resources,
        ReservationOwner {
            project_root: fs::canonicalize(&paths.root)?,
            attempt_id: attempt.into(),
            stage: request.stage.into(),
            run_accession: None,
        },
    )?;
    let inspect = resources.reserve(
        "inspect",
        ResourceRequest {
            cpu_threads: 1,
            memory_bytes: 64 * MIB,
            disks: vec![DiskRequest {
                path: paths.root.clone(),
                bytes: 4 * MIB,
            }],
        },
        &check,
    )?;
    let mut heads: Vec<GenerationHead> = Vec::new();
    let mut files = Vec::new();
    for group in request.selection.groups() {
        let Some(head) = db.generation_head(*group, db.project_id())? else {
            continue;
        };
        if !db.verify_generation_head(&head)? {
            bail!("selected {} generation failed verification", group.as_str());
        }
        let outputs = db.generation_output_artifacts(&head)?;
        let selected: Vec<_> = outputs
            .iter()
            .filter(|item| request.selection.includes(&item.kind) && item.path.ends_with(".tsv"))
            .cloned()
            .collect();
        if selected.is_empty() {
            continue;
        }
        for file in selected {
            files.push((head.clone(), file, false));
        }
        let provenance: Vec<_> = outputs
            .into_iter()
            .filter(|item| item.kind == "MANIFEST" && item.path.ends_with(".json"))
            .collect();
        if provenance.is_empty() {
            bail!(
                "selected {} generation has no registered provenance manifest",
                group.as_str()
            );
        }
        for file in provenance {
            files.push((head.clone(), file, true));
        }
        heads.push(head);
    }
    if heads.is_empty() {
        bail!("no sealed analysis outputs match the requested export");
    }
    let mut validation_heads = heads.clone();
    if request.require_current
        && config.corr.matrix.is_none()
        && heads
            .iter()
            .any(|head| head.output_group == OutputGroup::Corr)
        && !heads
            .iter()
            .any(|head| head.output_group == OutputGroup::De)
    {
        let de = db
            .generation_head(OutputGroup::De, db.project_id())?
            .context("current request verification requires a selected DE matrix")?;
        if !db.verify_generation_head(&de)? {
            bail!("selected DE generation failed verification");
        }
        validation_heads.push(de);
    }
    let cohort_guard = if request.require_current_cohort {
        Some(cohort::Guard::capture(
            paths,
            db,
            config,
            &validation_heads,
        )?)
    } else {
        None
    };
    let admission_guard = if request.require_current_admission {
        Some(admission::Guard::capture(db, config, &validation_heads)?)
    } else {
        None
    };
    let survey_proof = if matches!(request.mode, Mode::Survey { .. }) {
        let de = heads
            .iter()
            .find(|head| head.output_group == OutputGroup::De)
            .context("survey requires a selected analysis generation")?;
        let corr = heads
            .iter()
            .find(|head| head.output_group == OutputGroup::Corr)
            .context("survey requires a selected correlation generation")?;
        if !files
            .iter()
            .any(|(_, file, provenance)| !provenance && file.kind == "DE_TABLE")
        {
            bail!("survey requires differential-expression results; run deseq2 first");
        }
        Some(survey::compatibility(db, de, corr, &files)?)
    } else {
        None
    };
    // A single TSV record can span the entire source. Field-offset vectors can
    // exceed the input size substantially when a record has many tiny fields.
    // Disk output has a separate fourfold escaping allowance.
    let largest = files
        .iter()
        .map(|(_, file, _)| file.bytes)
        .max()
        .unwrap_or(0);
    let mut memory = largest
        .checked_mul(32)
        .and_then(|v| v.checked_add(64 * MIB))
        .and_then(|v| {
            (files.len() as u64)
                .checked_mul(64 * 1024)
                .and_then(|metadata| v.checked_add(metadata))
        })
        .context("export memory estimate overflow")?;
    // Each table has its own write cap. Reserve a separate, enforced budget for
    // the final manifest so it cannot consume space promised to those tables.
    let manifest_limit = (files.len() as u64)
        .checked_mul(64 * 1024)
        .and_then(|v| v.checked_add(4 * MIB))
        .context("export manifest estimate overflow")?;
    let mut disk = files.iter().try_fold(manifest_limit, |sum, (_, file, _)| {
        file.bytes
            .checked_mul(4)
            .and_then(|v| v.checked_add(64 * 1024))
            .and_then(|v| sum.checked_add(v))
            .context("export disk estimate overflow")
    })?;
    let mut scratch_disk = 4 * MIB;
    let mut external_matrix_bytes = None;
    if request.require_current_design {
        if validation_heads
            .iter()
            .any(|head| head.output_group == OutputGroup::Corr)
            && let Some(matrix) = &config.corr.matrix
        {
            let metadata = fs::metadata(matrix)
                .with_context(|| format!("cannot inspect external matrix {}", matrix.display()))?;
            if !metadata.is_file() {
                bail!(
                    "external matrix is not a regular file: {}",
                    matrix.display()
                );
            }
            let matrix_bytes = metadata.len();
            external_matrix_bytes = Some(matrix_bytes);
            memory = memory
                .checked_add(
                    matrix_bytes
                        .checked_mul(32)
                        .context("external matrix memory estimate overflow")?,
                )
                .context("external matrix memory estimate overflow")?;
        }
        let metadata_bytes = fs::metadata(paths.samplesheet_path())?.len();
        let runs = db.list_active_runs()?;
        let identities = runs.iter().try_fold(0_u64, |total, run| {
            let identity = rnaa_core::cohorts::BiologicalSample::for_run(run)?;
            total
                .checked_add(run.run_accession.len() as u64)
                .and_then(|v| v.checked_add(identity.id.len() as u64))
                .context("design identity estimate overflow")
        })?;
        let n = runs.len() as u64;
        let allowance = metadata_bytes
            .checked_add(identities)
            .context("design metadata estimate overflow")?;
        memory = memory
            .checked_add(
                allowance
                    .checked_mul(64)
                    .context("design memory estimate overflow")?,
            )
            .and_then(|v| {
                n.checked_mul(n)
                    .and_then(|cells| cells.checked_mul(128))
                    .and_then(|design| v.checked_add(design))
            })
            .context("design memory estimate overflow")?;
        scratch_disk = scratch_disk
            .checked_add(
                allowance
                    .checked_mul(16)
                    .context("design scratch estimate overflow")?,
            )
            .context("design scratch estimate overflow")?;
    }
    let summary_limit = if matches!(request.mode, Mode::Survey { .. }) {
        let metadata = files
            .iter()
            .find(|(_, file, provenance)| {
                !provenance
                    && file.kind == "SAMPLESHEET"
                    && Path::new(&file.path)
                        .file_name()
                        .is_some_and(|name| name == "sample_metadata.tsv")
            })
            .context("survey requires frozen sample metadata")?
            .1
            .bytes;
        memory = memory
            .checked_add(
                largest
                    .checked_mul(32)
                    .context("survey matrix estimate overflow")?,
            )
            .context("survey memory estimate overflow")?
            .checked_add(
                metadata
                    .checked_mul(64)
                    .context("survey metadata estimate overflow")?,
            )
            .context("survey memory estimate overflow")?;
        let limit = metadata
            .checked_mul(32)
            .and_then(|v| v.checked_add(64 * 1024))
            .context("survey summary estimate overflow")?;
        disk = disk
            .checked_add(limit)
            .context("survey disk estimate overflow")?;
        Some(limit)
    } else {
        None
    };
    drop(inspect);
    let _reservation = resources.reserve(
        "write",
        ResourceRequest {
            cpu_threads: 1,
            memory_bytes: memory,
            disks: vec![
                DiskRequest {
                    path: paths.root.clone(),
                    bytes: scratch_disk,
                },
                DiskRequest {
                    path: paths.exports_dir.clone(),
                    bytes: disk,
                },
            ],
        },
        &check,
    )?;
    let design_guard = if request.require_current_design {
        Some(design::Guard::capture(
            paths,
            db,
            config,
            &validation_heads,
            &files,
            admission_guard
                .as_ref()
                .context("design guard requires current admission")?
                .eligible(),
            external_matrix_bytes,
        )?)
    } else {
        None
    };
    let request_guard = if request.require_current {
        Some(request_guard::Guard::capture(
            paths,
            db,
            config,
            &validation_heads,
            admission_guard
                .as_ref()
                .context("request guard requires admission")?,
            design_guard
                .as_ref()
                .context("request guard requires design validation")?,
        )?)
    } else {
        None
    };
    if paths
        .exports_dir
        .symlink_metadata()
        .is_ok_and(|m| !m.is_dir() || m.file_type().is_symlink())
    {
        bail!("export directory must be a real directory");
    }
    fs::create_dir_all(&paths.exports_dir)?;
    let staged = tempfile::Builder::new()
        .prefix(&format!(".rnaa-{}-{attempt}-", request.stage))
        .tempdir_in(&paths.exports_dir)?;
    recovery::register(db, request.stage, staged.path())?;
    let mut entries = Vec::new();
    for (head, source, provenance) in files {
        check()?;
        let name = Path::new(&source.path)
            .file_name()
            .context("export source has no filename")?
            .to_str()
            .context("export filename is not UTF-8")?;
        let folder = if provenance {
            format!("provenance/{}", head.output_group.as_str())
        } else {
            head.output_group.as_str().into()
        };
        let name = if provenance {
            name.to_string()
        } else {
            format!(
                "{}.{}",
                name.strip_suffix(".tsv").context("expected TSV output")?,
                request.format
            )
        };
        let relative = Path::new(&folder).join(name);
        let target = staged.path().join(&relative);
        fs::create_dir_all(target.parent().context("export path has no parent")?)?;
        let limit = source
            .bytes
            .checked_mul(4)
            .and_then(|v| v.checked_add(64 * 1024))
            .context("export output limit overflow")?;
        let (rows, columns) = publication::write(
            &source,
            &target,
            request.delimiter,
            provenance,
            limit,
            &check,
        )?;
        let output = ManifestArtifact::snapshot("EXPORT", &target)?;
        entries.push(Entry {
            generation_id: head.generation_id.clone(),
            source,
            output: relative.to_string_lossy().into_owned(),
            output_sha256: output.checksum,
            output_bytes: output.bytes,
            rows,
            columns,
        });
    }
    if let Mode::Survey { stratify } = request.mode {
        let metadata_source = entries
            .iter()
            .find(|entry| {
                entry.source.kind == "SAMPLESHEET" && entry.output == "de/sample_metadata.tsv"
            })
            .context("survey metadata output missing")?
            .source
            .clone();
        let proof = survey_proof
            .as_ref()
            .context("missing survey compatibility proof")?;
        let summary = survey::summarize(
            staged.path(),
            &entries,
            stratify,
            &proof.design_contract,
            &check,
        )?;
        let output_path = staged.path().join("metadata_summary.json");
        publication::write_json(
            &output_path,
            &summary,
            summary_limit.context("missing survey output limit")?,
            &check,
        )?;
        let output = ManifestArtifact::snapshot("SURVEY_SUMMARY", &output_path)?;
        entries.push(Entry {
            generation_id: heads
                .iter()
                .find(|head| head.output_group == OutputGroup::De)
                .context("missing analysis head")?
                .generation_id
                .clone(),
            source: metadata_source,
            output: "metadata_summary.json".into(),
            output_sha256: output.checksum,
            output_bytes: output.bytes,
            rows: None,
            columns: None,
        });
    }
    for head in &validation_heads {
        check()?;
        if db
            .generation_head(head.output_group, &head.subject)?
            .as_ref()
            != Some(head)
            || !db.verify_generation_head(head)?
        {
            bail!("source generation changed during export; retry");
        }
    }
    if let Some(guard) = &cohort_guard {
        guard.verify(paths, db, config)?;
    }
    if let Some(guard) = &admission_guard {
        guard.verify_current(db, config)?;
    }
    if let Some(guard) = &design_guard {
        guard.verify_current(paths, db, config)?;
    }
    if let Some(guard) = &request_guard {
        guard.verify_current(
            paths,
            db,
            config,
            admission_guard
                .as_ref()
                .context("missing admission guard")?,
            design_guard.as_ref().context("missing design guard")?,
        )?;
    }
    let selection_scope = if request.require_current {
        "selected and required DE scientific requests reconstructed against current supported inputs; external files are not locked"
    } else {
        "sealed selected generations; current scientific settings are not revalidated"
    };
    let manifest = serde_json::json!({"schema_version":1,"project_id":db.project_id(),"attempt_id":attempt,"format":request.format,"kind":request.stage,
        "selection":selection_scope, "survey_compatibility":survey_proof,
        "current_cohort_validation":cohort_guard.as_ref().map(|guard| serde_json::json!({"membership_sha256":guard.digest(),"algorithm":"active_cohort_v2","full_freshness":false})),
        "current_admission_validation":admission_guard.as_ref().map(|guard| guard.proof()),
        "current_design_validation":design_guard.as_ref().map(|guard| guard.proof()),
        "current_request_validation":request_guard.as_ref().map(|guard| guard.proof()),
        "generations":heads,"files":entries});
    publication::write_json(
        &staged.path().join(if request.stage == "survey" {
            "survey-manifest.json"
        } else {
            "export-manifest.json"
        }),
        &manifest,
        manifest_limit,
        &check,
    )?;
    check()?;
    let id = staged
        .path()
        .file_name()
        .context("temporary export has no name")?
        .to_string_lossy()
        .trim_start_matches('.')
        .to_string();
    let destination = paths.exports_dir.join(id);
    publication::publish(db, staged.path(), &destination, |transaction| {
        if let Some(guard) = &cohort_guard {
            guard.verify(paths, transaction, config)?;
        }
        if let Some(guard) = &admission_guard {
            guard.verify_catalog(transaction, config)?;
        }
        if let Some(guard) = &design_guard {
            guard.verify_catalog(transaction, config)?;
        }
        if let Some(guard) = &request_guard {
            guard.verify_catalog(transaction, config)?;
        }
        transaction.mark_export_published()?;
        transaction.complete_owned_job(JobOutcome::Succeeded, None)
    })?;
    drop(staged);
    Ok(destination)
}
