use super::*;
use rnaa_core::config::MetadataMergeStrategy;
use rnaa_core::model::InputType;

#[test]
fn metadata_conflicts_cannot_publish_new_cohort_and_explicit_resolution_preserves_edits()
-> Result<()> {
    let temp = tempfile::tempdir()?;
    let paths = ProjectPaths::new(temp.path());
    paths.ensure_layout()?;
    let config = ProjectConfig::default();
    config.save(&paths.config_path)?;
    let db = Database::create(temp.path(), &config)?;
    db.add_input("SRP114962", InputType::Study)?;
    let resolver = EnaResolver::new(config.resolver.fields.clone())?;
    let fetch = |condition: &str| -> Result<Vec<RunRecord>> {
        let mut runs = EnaResolver::parse_report_rows(
            include_str!("../../../../fixtures/ena_srp114962.json"),
            db.project_id(),
        )?;
        for run in &mut runs {
            run.metadata.insert("condition".into(), condition.into());
        }
        Ok(runs)
    };
    resolver.resolve_project_owned(&db, &paths, &config, |_| fetch("control"))?;
    let cohort = db.active_cohort()?.unwrap();
    let publication = db.metadata_publication()?.unwrap();
    let edited = fs::read_to_string(paths.samplesheet_path())?.replace("control", "edited");
    fs::write(paths.samplesheet_path(), &edited)?;
    let failure = resolver.resolve_project_owned(&db, &paths, &config, |_| fetch("source-update"));
    assert!(format!("{:#}", failure.unwrap_err()).contains("refresh conflicts"));
    assert_eq!(db.active_cohort()?.unwrap().cohort_id, cohort.cohort_id);
    assert_eq!(
        db.metadata_publication()?.unwrap().publication_id,
        publication.publication_id
    );
    assert!(
        db.list_active_runs()?
            .iter()
            .all(|run| run.metadata["condition"] == "control")
    );
    assert_eq!(fs::read_to_string(paths.samplesheet_path())?, edited);
    resolver
        .with_metadata_conflicts(ConflictPolicy::KeepLocal)
        .resolve_project_owned(&db, &paths, &config, |_| fetch("source-update"))?;
    assert_ne!(db.active_cohort()?.unwrap().cohort_id, cohort.cohort_id);
    assert!(
        db.list_active_runs()?
            .iter()
            .all(|run| run.metadata["condition"] == "source-update")
    );
    assert!(fs::read_to_string(paths.samplesheet_path())?.contains("edited"));
    publication.verify()?;
    db.metadata_publication()?.unwrap().verify()?;
    Ok(())
}

#[test]
fn project_resolution_persists_narrow_and_empty_membership_without_losing_runs() -> Result<()> {
    let temp = tempfile::tempdir()?;
    let paths = ProjectPaths::new(temp.path());
    paths.ensure_layout()?;
    let config = ProjectConfig::default();
    config.save(&paths.config_path)?;
    let db = Database::create(temp.path(), &config)?;
    db.add_input("SRP114962", InputType::Study)?;
    let resolver = EnaResolver::new(config.resolver.fields.clone())?;
    let fixture = include_str!("../../../../fixtures/ena_srp114962.json");
    let fetch = |_: &InputRecord| EnaResolver::parse_report_rows(fixture, db.project_id());
    let first = resolver.resolve_project_owned(&db, &paths, &config, fetch)?;
    assert_eq!(first.runs.len(), 2);
    let excluded = first.runs[0].sample_accession.as_deref().unwrap();
    db.add_input(excluded, InputType::SampleExclude)?;
    assert!(db.list_active_runs().is_err());
    let second = resolver.resolve_project_owned(&db, &paths, &config, fetch)?;
    assert_eq!(second.runs.len(), 1);
    assert_eq!(db.list_runs()?.len(), 2);
    assert_eq!(
        db.list_active_runs()?[0].run_accession,
        second.runs[0].run_accession
    );
    db.add_input(
        second.runs[0].sample_accession.as_deref().unwrap(),
        InputType::SampleExclude,
    )?;
    let empty = resolver.resolve_project_owned(&db, &paths, &config, fetch)?;
    assert!(empty.runs.is_empty());
    assert!(db.list_active_runs()?.is_empty());
    assert_eq!(db.list_runs()?.len(), 2);
    let mut sheet = ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(paths.samplesheet_path())?;
    assert_eq!(sheet.records().count(), 0);
    Ok(())
}

#[test]
fn changed_selectors_during_fetch_prevent_resolution_publication() -> Result<()> {
    let temp = tempfile::tempdir()?;
    let paths = ProjectPaths::new(temp.path());
    paths.ensure_layout()?;
    let config = ProjectConfig::default();
    config.save(&paths.config_path)?;
    let db = Database::create(temp.path(), &config)?;
    db.add_input("SRP114962", InputType::Study)?;
    let resolver = EnaResolver::new(config.resolver.fields.clone())?;
    let result = resolver.resolve_project_owned(&db, &paths, &config, |_: &InputRecord| {
        db.add_input("condition=control", InputType::FilterInclude)?;
        EnaResolver::parse_report_rows(
            include_str!("../../../../fixtures/ena_srp114962.json"),
            db.project_id(),
        )
    });
    assert!(
        result
            .unwrap_err()
            .to_string()
            .contains("selectors changed")
    );
    assert!(db.active_cohort()?.is_none());
    assert!(db.list_runs()?.is_empty());
    Ok(())
}

#[test]
fn resolver_lease_excludes_live_analysis_and_records_preflight_failure() -> Result<()> {
    use rnaa_core::jobs::{JobLease, JobOutcome};
    let temp = tempfile::tempdir()?;
    let config = ProjectConfig::default();
    let db = Database::create(temp.path(), &config)?;
    let paths = ProjectPaths::new(temp.path());
    let resolver = EnaResolver::new(config.resolver.fields.clone())?;
    let mut lease = JobLease::try_acquire(&db, "normalize", None, &[])?.unwrap();
    let blocked = resolver.resolve_project(&db, &paths, &config).unwrap_err();
    assert!(blocked.to_string().contains("another analysis or resolver"));
    lease.finish(JobOutcome::Succeeded, None)?;
    drop(lease);
    let preflight = resolver.resolve_project(&db, &paths, &config).unwrap_err();
    assert!(preflight.to_string().contains("no inputs found"));
    let attempts = db.list_job_attempts(10)?;
    assert!(
        attempts
            .iter()
            .any(|attempt| attempt.stage == "resolve" && attempt.state == "FAILED")
    );
    Ok(())
}

#[test]
fn parses_ena_json_fixture() {
    let fixture = include_str!("../../../../fixtures/ena_srp114962.json");
    let runs = EnaResolver::parse_report_rows(fixture, "project-1").unwrap();
    assert_eq!(runs.len(), 2);
    assert_eq!(runs[0].study_accession, "SRP114962");
    assert_eq!(
        runs[0].source_study_accession.as_deref(),
        Some("PRJNA396019")
    );
    assert!(
        runs[0]
            .remote_files
            .iter()
            .any(|file| file.kind == RemoteFileKind::FastqR1 || file.kind == RemoteFileKind::Fastq)
    );
}

#[test]
fn merges_metadata_prefer_local() {
    let mut remote = BTreeMap::from([
        ("condition".to_string(), "remote".to_string()),
        ("center".to_string(), "ebi".to_string()),
    ]);
    let incoming = BTreeMap::from([
        ("condition".to_string(), "local".to_string()),
        ("batch".to_string(), "b1".to_string()),
    ]);
    merge_metadata(&mut remote, &incoming, MetadataMergeStrategy::PreferLocal);
    assert_eq!(remote.get("condition").unwrap(), "local");
    assert_eq!(remote.get("batch").unwrap(), "b1");
}

#[test]
fn merges_metadata_union_with_prefix() {
    let mut remote = BTreeMap::from([("condition".to_string(), "remote".to_string())]);
    let incoming = BTreeMap::from([("condition".to_string(), "local".to_string())]);
    merge_metadata(
        &mut remote,
        &incoming,
        MetadataMergeStrategy::UnionWithPrefix,
    );
    assert_eq!(remote.get("remote_condition").unwrap(), "remote");
    assert_eq!(remote.get("local_condition").unwrap(), "local");
}

#[test]
fn filters_runs_by_sample_selectors() {
    let fixture = include_str!("../../../../fixtures/ena_srp114962.json");
    let runs = EnaResolver::parse_report_rows(fixture, "project-1").unwrap();
    let include_samples = BTreeSet::from_iter([String::from("SAMN09909193")]);
    let exclude_samples = BTreeSet::new();
    let include_filters = Vec::new();
    let exclude_filters = Vec::new();
    let included = filter_runs_by_selectors(
        runs.clone(),
        &include_samples,
        &exclude_samples,
        &include_filters,
        &exclude_filters,
    );
    assert_eq!(included.len(), 1);
    assert_eq!(
        included[0].sample_accession.as_deref(),
        Some("SAMN09909193")
    );

    let include_samples = BTreeSet::new();
    let exclude_samples = BTreeSet::from_iter([String::from("SAMN09909193")]);
    let excluded = filter_runs_by_selectors(
        runs,
        &include_samples,
        &exclude_samples,
        &include_filters,
        &exclude_filters,
    );
    assert_eq!(excluded.len(), 1);
    assert_eq!(
        excluded[0].sample_accession.as_deref(),
        Some("SAMN09909196")
    );
}

#[test]
fn splits_selector_inputs_from_query_inputs() {
    let inputs = vec![
        InputRecord {
            project_id: "p".to_string(),
            input_id: "SRP1".to_string(),
            input_type: InputType::Study,
            added_at: "now".to_string(),
        },
        InputRecord {
            project_id: "p".to_string(),
            input_id: "samn1".to_string(),
            input_type: InputType::SampleInclude,
            added_at: "now".to_string(),
        },
        InputRecord {
            project_id: "p".to_string(),
            input_id: "samn2".to_string(),
            input_type: InputType::SampleExclude,
            added_at: "now".to_string(),
        },
        InputRecord {
            project_id: "p".to_string(),
            input_id: "instrument_platform=illumina".to_string(),
            input_type: InputType::FilterInclude,
            added_at: "now".to_string(),
        },
        InputRecord {
            project_id: "p".to_string(),
            input_id: "sample_accession!=SAMN0".to_string(),
            input_type: InputType::FilterExclude,
            added_at: "now".to_string(),
        },
    ];
    let selectors = split_inputs(&inputs);
    assert_eq!(selectors.query_inputs.len(), 1);
    assert_eq!(selectors.query_inputs[0].input_id, "SRP1");
    assert!(selectors.include_samples.contains("SAMN1"));
    assert!(selectors.exclude_samples.contains("SAMN2"));
    assert_eq!(
        selectors.include_filters_raw,
        vec!["instrument_platform=illumina"]
    );
    assert_eq!(
        selectors.exclude_filters_raw,
        vec!["sample_accession!=SAMN0"]
    );
}

#[test]
fn parses_and_applies_where_predicates() {
    let fixture = include_str!("../../../../fixtures/ena_srp114962.json");
    let runs = EnaResolver::parse_report_rows(fixture, "project-1").unwrap();
    let include_filters =
        parse_filter_predicates(&[String::from("sample_accession=SAMN09909193")]).unwrap();
    let exclude_filters = Vec::new();
    let selected = filter_runs_by_selectors(
        runs.clone(),
        &BTreeSet::new(),
        &BTreeSet::new(),
        &include_filters,
        &exclude_filters,
    );
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].run_accession, "SRR7750979");

    let include_filters = Vec::new();
    let exclude_filters =
        parse_filter_predicates(&[String::from("sample_accession=SAMN09909193")]).unwrap();
    let selected = filter_runs_by_selectors(
        runs,
        &BTreeSet::new(),
        &BTreeSet::new(),
        &include_filters,
        &exclude_filters,
    );
    assert_eq!(selected.len(), 1);
    assert_eq!(selected[0].run_accession, "SRR7750980");
}
