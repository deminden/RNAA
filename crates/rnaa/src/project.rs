//! Project membership, metadata resolution and reference preparation commands.
use crate::bootstrap::init_logging;
use crate::cli::{
    AddArgs, InitArgs, RefsArgs, RefsCommand, RefsPrepareArgs, ResolveArgs, parse_merge_strategy,
    validate_where_predicate,
};
use anyhow::Context;
use anyhow::{Result, bail};
use rnaa_core::Database;
use rnaa_core::config::ProjectConfig;
use rnaa_core::model::InputType;
use rnaa_core::paths::ProjectPaths;
use rnaa_core::traits::ReferenceManager;
use rnaa_core::util::sanitize_basename;
use rnaa_stages::refs::EnsemblReferenceManager;
use rnaa_stages::resolver::EnaResolver;
use std::fs;
use std::path::Path;

pub(crate) fn cmd_init(root: &Path, args: InitArgs, log_file: Option<&Path>) -> Result<()> {
    let paths = ProjectPaths::new(root);
    if paths.db_path.exists() {
        bail!("RNAA project already exists at {}", root.display());
    }
    paths.ensure_layout()?;
    let _log_guard = init_logging(&paths, "init", log_file)?;

    let mut config = if let Some(path) = args.config {
        ProjectConfig::load(&path)?
    } else {
        ProjectConfig::default()
    };
    if config.project.name == "rnaa-project"
        && let Some(name) = root.file_name().and_then(|value| value.to_str())
    {
        config.project.name = name.to_string();
    }
    config.save(&paths.config_path)?;
    let db = Database::create(root, &config)?;

    println!("project_id\t{}", db.project_id());
    println!("root\t{}", root.display());
    println!("config\t{}", paths.config_path.display());
    println!("database\t{}", paths.db_path.display());
    Ok(())
}

pub(crate) fn cmd_add(
    paths: &ProjectPaths,
    db: &Database,
    _config: &ProjectConfig,
    args: AddArgs,
) -> Result<()> {
    let merge_strategy = parse_merge_strategy(&args.merge_strategy)?;
    let mut total_inputs = 0_usize;
    for input_id in &args.ids {
        let input_type = InputType::from_accession(input_id)?;
        db.add_input(input_id, input_type)?;
        total_inputs += 1;
    }
    for sample in &args.samples {
        db.add_input(sample, InputType::Sample)?;
        total_inputs += 1;
    }
    for sample in &args.include_samples {
        db.add_input(sample, InputType::SampleInclude)?;
        total_inputs += 1;
    }
    for sample in &args.exclude_samples {
        db.add_input(sample, InputType::SampleExclude)?;
        total_inputs += 1;
    }
    for predicate in &args.where_predicates {
        validate_where_predicate(predicate)?;
        db.add_input(predicate, InputType::FilterInclude)?;
        total_inputs += 1;
    }
    for predicate in &args.exclude_where_predicates {
        validate_where_predicate(predicate)?;
        db.add_input(predicate, InputType::FilterExclude)?;
        total_inputs += 1;
    }
    if total_inputs == 0 && args.metadata.is_empty() {
        bail!(
            "nothing to add; provide at least one --id/--sample/--include-sample/--exclude-sample/--where/--exclude-where or --metadata"
        );
    }

    for metadata_path in &args.metadata {
        let file_name = metadata_path
            .file_name()
            .and_then(|value| value.to_str())
            .ok_or_else(|| {
                anyhow::anyhow!("invalid metadata file name {}", metadata_path.display())
            })?;
        let observed = rnaa_core::ManifestArtifact::snapshot("METADATA_OVERRIDE", metadata_path)?;
        let staged = tempfile::Builder::new()
            .prefix("override-")
            .suffix(&format!("-{}", sanitize_basename(file_name)))
            .tempfile_in(&paths.overrides_dir)?;
        fs::copy(metadata_path, staged.path()).with_context(|| {
            format!(
                "failed to copy metadata override {} to {}",
                metadata_path.display(),
                staged.path().display()
            )
        })?;
        let copied = rnaa_core::ManifestArtifact::snapshot("METADATA_OVERRIDE", staged.path())?;
        if !observed.matches_file()?
            || observed.checksum != copied.checksum
            || observed.bytes != copied.bytes
        {
            bail!("metadata override changed while being registered; retry");
        }
        staged.as_file().sync_all()?;
        let (_file, dest) = staged.keep()?;
        #[cfg(unix)]
        fs::File::open(&paths.overrides_dir)?.sync_all()?;
        db.add_metadata_override(metadata_path, &dest, merge_strategy)?;
    }

    println!("added_inputs\t{}", total_inputs);
    println!("added_accessions\t{}", args.ids.len());
    println!("added_samples\t{}", args.samples.len());
    println!("added_include_samples\t{}", args.include_samples.len());
    println!("added_exclude_samples\t{}", args.exclude_samples.len());
    println!("added_where_predicates\t{}", args.where_predicates.len());
    println!(
        "added_exclude_where_predicates\t{}",
        args.exclude_where_predicates.len()
    );
    println!("metadata_overrides\t{}", args.metadata.len());
    Ok(())
}

pub(crate) fn cmd_resolve(
    paths: &ProjectPaths,
    db: &Database,
    config: &ProjectConfig,
    args: ResolveArgs,
) -> Result<()> {
    let resolver = EnaResolver::new(config.resolver.fields.clone())?.with_metadata_conflicts(
        crate::cli::metadata_conflict_policy(&args.metadata_conflicts)?,
    );
    db.set_project_config(config)?;
    if args.force {
        tracing::info!(
            "force flag set; resolver will refresh ENA metadata and upsert existing runs"
        );
    }
    let _resolved = resolver.resolve_project(db, paths, config)?;
    Ok(())
}

pub(crate) fn cmd_refs(
    paths: &ProjectPaths,
    db: &Database,
    config: ProjectConfig,
    args: RefsArgs,
) -> Result<()> {
    match args.command {
        RefsCommand::Prepare(prepare) => cmd_refs_prepare(paths, db, config, prepare),
    }
}

fn cmd_refs_prepare(
    paths: &ProjectPaths,
    db: &Database,
    mut config: ProjectConfig,
    args: RefsPrepareArgs,
) -> Result<()> {
    let mut changed = false;
    if let Some(organism) = args.organism {
        config.refs.organism = organism;
        changed = true;
    }
    if let Some(ensembl) = args.ensembl {
        config.refs.ensembl = ensembl;
        changed = true;
    }
    match (args.gtf, args.cdna) {
        (Some(gtf), Some(cdna)) => {
            config.refs.custom_gtf = gtf.display().to_string();
            config.refs.custom_cdna = cdna.display().to_string();
            changed = true;
        }
        (None, None) => {}
        _ => bail!("--gtf and --cdna must be provided together"),
    }

    if changed {
        config.save(&paths.config_path)?;
    }

    db.set_project_config(&config)?;
    let manager = EnsemblReferenceManager::new()?;
    let reference = manager.ensure_reference(paths, db, &config)?;
    println!("reference_id\t{}", reference.id);
    println!("cdna\t{}", reference.cdna_path.display());
    println!("gtf\t{}", reference.gtf_path.display());
    println!("index\t{}", reference.kallisto_index_path.display());
    println!("manifest\t{}", reference.manifest_path.display());
    Ok(())
}
