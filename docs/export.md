# Exporting analysis tables

```sh
rnaa export --what counts --format csv
rnaa export --what normalized --format tsv
rnaa export --what vst
rnaa export --what de
rnaa export --what correlation
rnaa export --what all --format csv
rnaa export --what all --require-current-cohort
rnaa export --what all --require-current-admission
rnaa export --what all --require-current-design
rnaa export --what all --require-current
```

The default selection is `all` and the default format is `tsv`. Every successful
command prints the path of a new directory under the project's `exports/`. The
`de/` and `corr/` subdirectories contain selected tables; `provenance/` retains
their registered source manifests. `export-manifest.json` records the exact source
generation and file checksums, output checksums, and table row/column counts.
Original provenance may contain local paths and analysis details; inspect it before
sharing an export.

An export is a **frozen snapshot of selected sealed generations**. It checks the
generation inventory before and after copying, verifies the bytes actually read
against each source checksum, and refuses a changed or corrupt source. It does
not reconstruct scientific requests or prove that the selected
results match the project's current cohort, metadata, reference or settings. The
DE and correlation heads in an `all` export may have different inputs. Rerun the
relevant analysis and inspect its provenance when current or mutually consistent
results are required. The command prints this scope and the manifest records it.

`--require-current-cohort` checks that every selected generation recorded the
same biological-sample/run membership as the currently resolved active cohort.
The check rejects unresolved or changed cohort selectors and rechecks the catalog
at publication. The manifest records the checked cohort digest. This option does
**not** validate current QC eligibility, quantification inputs, reference,
metadata values, scientific settings or compatibility between generations. It is
not a full freshness check. Without the option, frozen historical exports remain
available even if the current cohort has changed.

`--require-current-admission` also verifies the recorded QC eligibility and
quantified abundance identities for the active cohort, and the selected analysis
reference and its current source files. It implies the cohort check. It observes
existing QC decisions and refuses missing or changed evidence rather than
refreshing it. This may take time to hash current inputs, and a busy run prevents
the check. Scientific model options, metadata values and every other analysis
input are **not** fully revalidated; this is not a full freshness check. The
manifest records the eligibility and reference digests and their limited scope.

`--require-current-design` implies the admission check and reconstructs each
selected DE or correlation adjustment design from current eligible samples and
their used metadata. It checks DE design, transform and abundance-count options,
and correlation adjustment/model when applicable. It records which designs were
checked and a checksum of the current samplesheet. Unused metadata values can
change without invalidating a matching design. The recorded DE contrasts are
used for reconstruction; current contrast selections are not compared.
Correlation method, gene selection, pathway options and all other scientific
inputs or settings are also outside this check. It is not a full freshness check.
The check uses temporary metadata copies and memory/disk estimates, including an
allowance based on the size of a selected external matrix. These are not hard
limits; large or unusually complex designs may still exceed them. A killed
process can leave its temporary design metadata under `.rnaa_support`. The live
samplesheet is checked before publication but external edits are not locked.
For external matrices, sample-metadata preparation reads a size-bounded header
instead of loading every numeric value. Strict validation hashes the file and
reconstructs its model within the admitted file-size bound; an observed size change causes
the export to fail.

`--require-current` implies the preceding checks and reconstructs the supported
current scientific requests, including current contrast selections, correlation
method, matrix and gene selection, module settings and output options, engine
identities, and current native or explicitly declared external quantification
evidence. Correlation using the selected DE VST matrix verifies its producing DE
request, even when the export selects correlation alone. For an independently
supplied matrix, the selected file, current metadata and relevant options are
checked; correlation-only export does not require a DE result. Gene-list and
pathway files are checked against their recorded identities. The manifest identifies
the checked generations and request digests.

Module options from the latest correlation invocation are stored in project
configuration, including the resolved pathway-file path. Running correlation
without `--module-gmt` clears those options. Older module results that lack these
settings fail closed until correlation is rerun. The latest explicit external
matrix path is also stored in project configuration; running correlation without
`--matrix` selects the current DE matrix again. Older external-matrix results
without a saved selection may need correlation to be rerun. The external matrix
check attests its current bytes and metadata, not how its values were generated.
This mode requires recorded input
lineage for the selected analysis, reference and native quantification ancestry;
unresolved legacy lineage must be imported before this check can pass.
Other frozen-export modes remain
available for those results. The check observes current evidence at execution time;
external files and configuration are not locked against outside edits, so operators
must keep project inputs stable while it runs. External quantification declarations
attest supplied data and reference binding, not the provenance of upstream reads.

CSV conversion preserves field values and ordering, including quoted fields. It
does not change the canonical TSV artifacts. An explicit selection with no matching
sealed outputs fails. Unsupported formats/selections and capacity failures also
return a nonzero status. A project-owned export job and generation reader pin
coordinate with analysis and cleanup. The shared resource pool reserves estimated
CPU, memory and disk for inspection and writing; estimates are cooperative, not
operating-system limits. A very long single table record may delay cancellation,
and inspection inventory memory is not yet fully bounded by the estimate. The
write estimate conservatively budgets 32 times the largest source file plus
per-file metadata, so large inputs can be refused even when actual usage would fit.
The final manifest has a separate reserved size limit; exceptionally large
provenance metadata can cause the export to fail without publishing a bundle.
Manifest writing checks for cancellation during serialization.

Tables are staged in a private directory and published by a same-filesystem,
no-replacement rename after validation. Ordinary failures and cancellation remove
the temporary directory and preserve previous exports. A killed process can leave
staged data or a renamed bundle. On the next export or survey, after the previous
job lease has expired, RNAA checks its catalog record and the directory's filesystem
identity. It moves matching unfinished staging to `.rnaa-recovered-*` for
inspection and records renamed bundles as orphaned. Recovery never deletes these
directories or marks an interrupted job successful. Conflicting, replaced or
untracked paths are left untouched. Recovery events include the prior attempt ID
and the path relative to `exports/`.

Only directories registered by this version can be recovered automatically;
older leftovers and the small window before registration require manual
inspection. Power-loss durability remains unqualified. Publication and job
completion are fenced in one SQLite transaction, and a process-exit test covers
the rename/commit boundary. Inspect an orphaned bundle and its manifest before use.
