# Run progress and cohort analysis

Each run records its own download, verification, preprocessing and quantification
progress. Normalization, differential expression and correlation operate on a
frozen cohort and record their attempts, errors and completion separately.
An analysis success, historical selection, publication failure or worker recovery
does not rewrite member run states, diagnostics or timestamps. An invalid count
input can still trigger the existing evidence-repair path; a progress label is
never sufficient proof that an artifact is usable.

Project-stage records give the current analysis summary. Durable attempts record
ownership and outcomes; sealed generation requests bind the specific cohort,
membership and scientific inputs. They are not independent per-run DE jobs.
Changing active membership invalidates the analysis summary while retaining the
excluded runs and their histories. An explicit compatible matrix can be correlated
with quantified cohort members without first inventing DE completion.

`rnaa status` prints `run_progress_model`, individual `cohort_stage` records,
`cohort_stage_error` diagnostics and the legacy migration receipt count. Ingest
counts, raw-byte totals and observed ingest throughput include active runs only.
Project-wide analysis is reported as zero or one completed stage. Empty cohorts
cannot report completed analysis. The live display calls analysis running only
when its project-stage record is `RUNNING`; failures and idle quantified cohorts
show waiting instead. A cohort's run count may exceed its biological sample count
when technical runs are pooled.

## Upgrading an existing project

Stop the project's analysis worker before opening it with the upgraded binary.
Migration `0021_run_progress` acquires the same operating-system writer lock used
by analysis. It waits up to two seconds for an active owner and otherwise returns
an idle-worker error without projecting labels. Once the migration is applied,
ordinary database opens can coexist with a live analysis owner.

The migration copies all six legacy `DE_*` and `CORR_*` run labels, original errors
and original update timestamps into `legacy_run_analysis_states`, together with a
projection timestamp. It then projects these runs to `QUANT_DONE` and clears their
analysis error from the run record. Original ingest states and diagnostics remain
untouched. The projection does not create artifact records, generations or
successful analysis stages. Normal artifact/admission checks still apply before
reuse. `CLEANED` remains a legacy retained-input progress label, not analysis proof.

The receipts are immutable and exposed by `Database::legacy_run_analysis_states`.
The SQL transaction includes both projection and migration version; an error
rolls back the entire change. Reopening or competing upgrade callers cannot
duplicate receipts. Database guards reject subsequent insertion or assignment of
legacy analysis labels, including by an older writer process. Use the upgraded
binary for subsequent commands; running old analysis code against the upgraded
database is unsupported.

## Recovery evidence

Migration tests cover all six labels, exact diagnostic preservation, rollback,
writer exclusion, concurrent upgrades and old-writer rejection. Native CLI tests
compare complete run records after successful analysis, historical reuse and
injected publication failures. They also test quantified-only custom matrices and
cohort narrowing. Real child workers are killed during normalization, DE and
correlation attempts, then recovered after lease expiry while an independent
ingest owner remains live. Recovery fails the interrupted analysis attempt and
preserves both the member records and the live ingest owner.

These local checks do not establish power-loss behavior, prompt cancellation
inside numerical engines or supported-host release qualification.
