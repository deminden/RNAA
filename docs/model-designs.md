# Typed metadata and model designs

RNAA uses the same typed model builder for DE and correlation adjustment.
It does not infer a continuous variable from numeric-looking strings.
`condition` and `batch` are factors by default. Declare other covariates in
`rnaa.toml` before using them in a design:

```toml
[metadata]
minimum_replicates = 2
missing_values = ["NA", "not collected"]

[metadata.columns.condition]
type = "factor"
levels = ["control", "treated"]
reference = "control"

[metadata.columns.age]
type = "numeric"

[metadata.columns.donor]
type = "factor"

[metadata.columns.source_url]
type = "ignore"
```

Use the level labels in your own samplesheet. Unknown configuration sections/keys,
unsupported project schema versions (only version 1 is supported), and invalid
type/level/reference combinations are errors. Model column names use ASCII letters,
digits and underscores, starting with a letter or underscore. Identifier/ignored
columns can retain other source column names, using quoted TOML keys when needed.

## Types and missing values

- `factor` preserves string categories, including distinct labels `01` and `1`.
  Declared levels define the allowed values and their order. Unobserved nonreference
  levels are omitted from the actual model; the full declaration stays in provenance.
  An explicit reference must be observed after selection, grouping and QC.
- `numeric` requires finite floating-point values. Invalid numbers and infinities
  fail. Equivalent forms such as `21` and `21.0` agree across technical runs.
- `identifier` and `ignore` exclude a column from modeling. Biological/run IDs and
  read paths cannot be reassigned as numeric/factor covariates. Study or donor
  covariates can be explicitly declared when scientifically appropriate.

Values are trimmed before model interpretation. Empty/whitespace-only values always
fail for a referenced variable. Additional missing markers are explicit, exact and
case sensitive. `NA` is not silently made missing unless configured that way.
Missing values are neither imputed nor dropped automatically.

If factor levels are omitted, the observed levels are sorted and frozen in the
request. If its reference is omitted, RNAA uses the first declared level, or the
first sorted observed level. Declare a reference to prevent a selected subset from
choosing another baseline. A missing declared reference fails instead of silently
switching it. Model/contrast values must agree within a biological sample's runs.

## Formula and preflight

The native engine validates its supported formula syntax, including main effects,
interactions, intercept control and supported transforms. RNAA's variable discovery
distinguishes function names, quoted labels and named arguments from covariates.
Dot formulas expand only declared factor/numeric columns present in the metadata;
identity and operational columns do not become predictors through `.`.

An intercept-only `~ 1` model preserves the biological sample dimension. DE requires
at least one coefficient. Both DE and adjustment require a full-rank matrix and
positive residual degrees of freedom. Confounding and saturated designs fail
before a replacement output generation or expression fitting is started.

Categorical contrasts are resolved against that actual model, its observed levels
and reference. Every contrasted level must contain at least
`metadata.minimum_replicates` biological samples, default two. Technical runs do
not increase this count. Contrasts against numeric covariates, absent levels,
duplicate/unsafe output names and non-estimable effects fail before fitting.
Automatically generated contrast filenames are sanitized without changing the
scientific level labels. Supported DE formula offsets are frozen in the contract;
linear residualization rejects offsets explicitly.

Formula-level releveling is resolved against the resulting coefficients, including
transformed factor names. A full-indicator model such as `~ 0 + condition` has no
treatment baseline; its contrast subtracts the two group coefficients. A request
with multiple estimable views of the same source factor is rejected as ambiguous.
Categorical Wald results retain the engine's two-group zero rule: expression only
in a third group does not turn a zero-versus-zero comparison into an apparent effect.

## Saved design and publication

`model_design.json` is a registered, sealed output of DE and residualization. It
records the formula, schema, sample order, typed values, observed levels, references,
coefficient names, row-major matrix, offsets, rank, residual degrees of freedom and
resolved contrast vectors with source/actual factor names and per-sample group
labels. Its typed model values are part of the scientific request.
DE request revision 6 and correlation revision 5 use a [scientific metadata
projection](scientific-cache.md) while retaining full schema/source observations.

Publication checks scientific settings against both the current project catalog
and the saved configuration. Schema/design changes during a fit cannot publish an
old model under new settings. Catalog guards run inside the publication transaction;
file checks run outside it. Thread-count changes do not invalidate this observation
guard. Default correlation also refuses a DE matrix whose reconstructed current
model/membership or DE settings differ from its recorded scientific contract.
It also revalidates current count dependencies. DE/correlation publication observes
upstream quantification settings and source metadata so a policy edit during fitting
cannot silently attach old native counts to the new requested context.

Preflight failures leave existing result generations and valid run counts intact,
record a failed attempt and mark the requested project stage failed. Valid designs
can still fail later statistical fitting; structural checks do not establish power,
scientific parity or numerical success for every expression dataset.

Editable model values survive [controlled metadata refresh](metadata-refresh.md),
with retained source witnesses and explicit conflict policies. Identity/selection
corrections use registered overrides and resolution. Unused operational fields are
excluded from scientific cache identity; independent scientific oracle fixtures
remain open.
