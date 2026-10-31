use super::*;

#[derive(Debug, Default)]
pub(super) struct InputSelectors<'a> {
    pub(super) query_inputs: Vec<&'a InputRecord>,
    pub(super) include_samples: BTreeSet<String>,
    pub(super) exclude_samples: BTreeSet<String>,
    pub(super) include_filters_raw: Vec<String>,
    pub(super) exclude_filters_raw: Vec<String>,
}

pub(super) fn split_inputs(inputs: &[InputRecord]) -> InputSelectors<'_> {
    let mut query_inputs = Vec::new();
    let mut include_samples = BTreeSet::new();
    let mut exclude_samples = BTreeSet::new();
    let mut include_filters = Vec::new();
    let mut exclude_filters = Vec::new();
    for input in inputs {
        match input.input_type {
            InputType::SampleInclude => {
                include_samples.insert(normalize_sample_selector(&input.input_id));
            }
            InputType::SampleExclude => {
                exclude_samples.insert(normalize_sample_selector(&input.input_id));
            }
            InputType::FilterInclude => include_filters.push(input.input_id.clone()),
            InputType::FilterExclude => exclude_filters.push(input.input_id.clone()),
            _ => query_inputs.push(input),
        }
    }
    InputSelectors {
        query_inputs,
        include_samples,
        exclude_samples,
        include_filters_raw: include_filters,
        exclude_filters_raw: exclude_filters,
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PredicateOp {
    Eq,
    Neq,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct FieldPredicate {
    field: String,
    op: PredicateOp,
    value: String,
}

pub(super) fn parse_filter_predicates(raw: &[String]) -> Result<Vec<FieldPredicate>> {
    raw.iter()
        .map(|item| parse_filter_predicate(item))
        .collect()
}

fn parse_filter_predicate(spec: &str) -> Result<FieldPredicate> {
    let trimmed = spec.trim();
    if let Some((field, value)) = trimmed.split_once("!=") {
        return build_filter_predicate(field, PredicateOp::Neq, value, spec);
    }
    if let Some((field, value)) = trimmed.split_once('=') {
        return build_filter_predicate(field, PredicateOp::Eq, value, spec);
    }
    bail!(
        "invalid filter predicate '{}'; expected 'field=value' or 'field!=value'",
        spec
    )
}

fn build_filter_predicate(
    field: &str,
    op: PredicateOp,
    value: &str,
    original: &str,
) -> Result<FieldPredicate> {
    let field = field.trim().to_ascii_lowercase();
    let value = value.trim().to_ascii_lowercase();
    if field.is_empty() || value.is_empty() {
        bail!(
            "invalid filter predicate '{}'; field and value must be non-empty",
            original
        );
    }
    Ok(FieldPredicate { field, op, value })
}

pub(super) fn filter_runs_by_selectors(
    runs: Vec<RunRecord>,
    include_samples: &BTreeSet<String>,
    exclude_samples: &BTreeSet<String>,
    include_filters: &[FieldPredicate],
    exclude_filters: &[FieldPredicate],
) -> Vec<RunRecord> {
    runs.into_iter()
        .filter(|run| {
            let sample = run
                .sample_accession
                .as_deref()
                .map(normalize_sample_selector)
                .unwrap_or_default();
            (include_samples.is_empty() || include_samples.contains(&sample))
                && (exclude_samples.is_empty() || !exclude_samples.contains(&sample))
                && include_filters
                    .iter()
                    .all(|predicate| run_matches_predicate(run, predicate))
                && !exclude_filters
                    .iter()
                    .any(|predicate| run_matches_predicate(run, predicate))
        })
        .collect()
}

fn normalize_sample_selector(value: &str) -> String {
    value.trim().to_ascii_uppercase()
}

fn run_matches_predicate(run: &RunRecord, predicate: &FieldPredicate) -> bool {
    let Some(field_value) = run_field_value(run, &predicate.field) else {
        return false;
    };
    let normalized = field_value.trim().to_ascii_lowercase();
    match predicate.op {
        PredicateOp::Eq => normalized == predicate.value,
        PredicateOp::Neq => normalized != predicate.value,
    }
}

fn run_field_value<'a>(run: &'a RunRecord, field: &str) -> Option<&'a str> {
    match field {
        "run_accession" => Some(run.run_accession.as_str()),
        "sample_accession" => run.sample_accession.as_deref(),
        "study_accession" => Some(run.study_accession.as_str()),
        "source_study_accession" => run.source_study_accession.as_deref(),
        "experiment_accession" => run.experiment_accession.as_deref(),
        "library_layout" => Some(run.library_layout.as_str()),
        "instrument_platform" => run.instrument_platform.as_deref(),
        "instrument_model" => run.instrument_model.as_deref(),
        _ => run.metadata.get(field).map(String::as_str),
    }
}
