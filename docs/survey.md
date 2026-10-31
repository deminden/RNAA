# Frozen survey bundles

```sh
rnaa survey
rnaa survey --stratify condition
rnaa survey --stratify condition --require-current-cohort
rnaa survey --stratify condition --require-current-admission
rnaa survey --stratify condition --require-current-design
rnaa survey --stratify condition --require-current
```

The command creates a new directory under `exports/` and prints its path. It
contains the selected sealed DE and correlation tables as TSV, their source
manifests, `metadata_summary.json`, and `survey-manifest.json` with generation
identities, source/output checksums, dimensions and compatibility evidence.
Sample metadata and source manifests may contain sensitive values or local paths;
inspect the bundle before sharing it.

Both a differential-expression result and a correlation result are required.
The command verifies the complete sealed inventories. Their scientific requests
must agree on active-cohort, eligibility and reference identities, and the
correlation input's recorded SHA-256 and byte count must match the selected VST
matrix with a verified lineage witness. The original correlation input may belong
to another generation when its matrix bytes are identical: this supports valid
historical reuse. The manifest records that input generation when available.

The survey checks sample IDs/order across the frozen metadata and matrices, gene
IDs across counts/normalized/VST tables, and membership of correlation, DE, edge,
dense and module gene IDs. Correlation design covariates are checked against the
frozen sample metadata when adjustment uses them. Enrichment tables and scientific
values retain their original sealed provenance; the survey does not independently
recompute them.

`--stratify` names a column in the **frozen biological-sample metadata**. The
summary counts samples by exact value, with empty values counted separately as
missing. It does not stratify runs or consult the current live samplesheet. A
missing/ambiguous column, inconsistent input, unavailable result or corrupt seal
fails without replacing a prior bundle.

Like `export`, a survey is a snapshot of selected frozen generations. It does not
establish that those generations match current project settings or metadata. It
can optionally require current cohort, admission, model-design or supported
scientific-request checks as described in [export](export.md). The survey
uses the same project job, resource pool, reader pin, per-file bounds, transactional
publication and cancellation behavior described in [export](export.md). Survey
memory estimates add conservative allowances for the input matrix and metadata;
large inputs can be rejected even when observed memory use would fit. Interrupted
tracked staging is preserved for inspection on the next export or survey;
automatic deletion and power-loss qualification are still outstanding.
Inspection of scientific request and lineage records is not fully memory-bounded.
