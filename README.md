# RNAA

Rust Nucleotide Alignment and Analytics.

RNAA is an RNA-seq survey pipeline for public cohorts. It resolves ENA accessions, downloads raw data with durable SQLite state, prepares transcriptome references, quantifies with kallistors, runs native Rust DE and normalization with rsdeseq2, and computes adjusted correlation outputs.

## Version

Current workspace version: `0.3.1` (unreleased).

The current development tree is being prepared for **0.3.1**. See the
[changelog](CHANGELOG.md), [release gates and upgrade guide](docs/releasing.md),
and [offline installation example](examples/README.md). Public release gates remain open. The current differential-expression engine is
`rsdeseq2 0.4.0`; see the changelog for compatibility and validation scope.

## Capabilities

- project init, add, resolve, status, doctor
- ENA resolution for run, study, project, and sample accessions
- sample include/exclude selectors and metadata include/exclude predicates
- restart-safe download worker with verification and manifests
- cached GENCODE preset references for human and mouse
  default presets use matching `transcripts.fa.gz` plus `annotation.gtf.gz`
- optional in-process `fasterp` preprocessing before quantification
- FASTQ retention after successful quantification: `both`, `original`, `trimmed`, or `none`
- threaded `fasterp` preprocessing with configurable per-run thread count
- normalized preprocess QC storage with conservative default gating before quant/downstream stages
- in-process kallistors reference indexing and quantification
- restart-safe quant reconciliation from existing artifacts and matching reference manifests
- normalization-only stage via `rnaa normalize`
- native Rust tximport-style abundance aggregation and DE through rsdeseq2
- adjusted correlation through Rust residualization plus `mincorr`
- TSV/CSV export of sealed analysis tables with generation provenance
- frozen survey bundles with cohort-compatibility checks and metadata summaries
- restart-safe DE/corr reconciliation from request manifests and existing outputs
- complete download, preprocessing, quantification, normalization/DE and correlation output generations with atomic database publication
- optional module-level enrichment from Spearman modules via Rust `rsfgsea`
- orchestration via `rnaa run` with download plus quant overlap
- quant run-level parallelism via `[quant].workers` or `rnaa quant --workers`

See [exporting analysis tables](docs/export.md) and [survey bundles](docs/survey.md)
for selections, provenance and the frozen-generation scope.

## Limitations

- shared blobs, unreachable project generations and count-backed FASTQ/SRA inputs support preview, trash, restore and expired purge
- `.sra` extraction fallback is not finished
- vendored `fasterp` remains patched for workspace quality gates

## Requirements

- Rust stable toolchain (Rust 1.97.1 or newer for the current kallistors/rsdeseq2 releases)
- See [Linux qualification and installation](docs/linux-support.md) for verified
  targets, native build prerequisites and the remaining platform gates.
- HTTP downloads, transcript indexing, quantification, preprocessing, differential expression, and correlation are provided by Rust libraries linked into RNAA; no R runtime or bioinformatics command-line tools are required.

## Quickstart

Build from the complete workspace and run the small offline installation check:

```bash
git clone https://github.com/deminden/RNAA
cd RNAA
cargo build --release --locked -p rnaa
python3 examples/offline_reference.py --binary target/release/rnaa
./target/release/rnaa --help
```

The offline check requires Python 3.11+; the RNAA executable itself does not
require Python. For a complete public-data walkthrough, including reference
checksums, sample metadata, quantification, differential expression, correlation,
export and survey, follow [the yeast workflow](examples/public_yeast.md). It
requires network access and several gigabytes of free disk space.

For your own cohort, prepare appropriate sample metadata and biological replicates
before fitting a model; do not assume a public dataset has columns such as `batch`
or `condition`. The [model design guide](docs/model-designs.md) explains the
requirements. `rnaa download` processes the available queue and exits;
`--forever` keeps a worker running and is intended for a separate long-lived
process, such as the [systemd example](docs/systemd/README.md).

FASTQ inputs are validated before processing. All single-end files are consumed;
paired multi-lane files must have matching lane basenames with R1/R2 or _1/_2
mate markers. Missing mates, unequal counts, mismatched read IDs, duplicate files,
and mixed raw/trimmed inputs fail explicitly. Multiple lanes use bounded-memory
temporary staging; temporary disk space must accommodate their decoded contents.

In `[quant]`, `fastq_max_line_mb` defaults to 1 MiB per decoded line and must be
positive; `fastq_max_decoded_mb` defaults to 0 (no aggregate cap). Preprocessing's
`preprocess_max_input_mb`, also available as `--preprocess-max-mb`, defaults to
1024 MiB across **all decoded mates and lanes**, replacing the old per-file
compressed-size check. Setting it to 0 disables that additional cap; the global
decoded cap still applies. Exact preflight bounds normalized temporary writes.
These limits govern new work; healthy cached results remain reusable after a
limit changes. See [FASTQ preparation](docs/fastq-preparation.md) for details.

References, downloads, preprocessing and quantification reserve CPU, estimated memory, staging and output
disk space from a common local pool before execution. Configure `[resources]` in
`rnaa.toml`: `pool_root` (empty uses `$HOME/.cache/rnaa/resources/v1`), `cpu_threads`
(0 detects available parallelism), `memory_mb` (0 uses 75% of detected capacity),
`disk_headroom_mb` (1024) and `wait_seconds` (300). All participating commands must
use the same local pool. Capacity failures cannot trigger preprocessing bypass;
healthy cache reuse does not open the pool. Fresh correlation adjustment and
network computation also make cooperative CPU, estimated-memory and disk
reservations, but other analysis and publication work is not fully admitted.
These estimates do not establish a pipeline-wide memory or disk ceiling. See
[resource reservations](docs/resource-reservations.md).

Download admission covers transfer buffers, remaining known-size writes or checked
unknown-size write windows, FASTQ validation and the extra generation publication
copy. A refused stream-window increase saves its verified prefix for retry. A
publication capacity failure preserves prior results and the verified transfer
cache. No partial generation becomes the active download.

Fresh quantification also checks a decoded index loading estimate against the
declared memory capacity before staging FASTQs or loading the graph. Compressed
bitmap expansion is counted without materializing the values. Cold loading and
the first run reserve capacity together; the resident graph retains its own memory
claim across worker waves. Idle graphs can be evicted to let waiting work progress.
Run estimates include packed reads, mapping caches and EC/EM state; conservative
ambiguity bounds and actual memory use still need production qualification. See
[quant index preflight](docs/quant-index-resources.md).

Quantification accepts `--strand unstranded|forward|reverse`,
`--single-fragment-length N`, and `--single-fragment-sd N`. These assumptions are
recorded in provenance and checked for cache reuse. Failed requested quantification
work returns a nonzero exit status.

## Runtime Behavior

`rnaa doctor --json` reports a versioned diagnostic object containing catalog,
configuration and resource-policy checks and the linked engine versions. An
uninitialized root is reported without creating a project. Invalid configuration,
an unreadable or corrupt catalog, and resource-detection failures produce a
nonzero exit status after the report is printed. When the configuration file is
missing, doctor validates the catalog copy used by normal project startup.
The report's `ok` means no diagnostic error was detected; it does not assert that
an uninitialized project is ready to run.
The catalog is opened read-only and checked with SQLite `quick_check`; this is not
a network, full integrity, scientific-validity or output-provenance audit. SQLite
may manage WAL shared-memory files while opening a live catalog.

`rnaa run` prints compact stage summaries such as:

```text
Project  <project_id>
Elapsed  00:02:14
Stage    quantification | download 6/6 | quant 2/6 | normalize 0/1 | corr 0/1
ETA      00:18:00 total | 00:00:00 download | 00:18:00 processing | limiting processing
Active   quantifying 2 run(s); 2/6 complete
Paths    /data/my-rnaa-project/metadata/samplesheet.tsv | /data/my-rnaa-project/logs
```

`rnaa status` uses the same project-level stage model, so `normalize` and `corr` are reported as `0/1` or `1/1` instead of fake per-run progress.

`rnaa status --json` emits a schema-versioned status object with active cohort
membership, per-run state and eligibility, stage outcomes, reference selection,
generation heads, progress, storage and recent errors. Durations are seconds;
unknown estimates are null. If the configured cohort cannot be verified, active
runs and cohort progress are null and the reason is reported. Like text status,
this command performs recovery/reconciliation and is an observational report,
not an atomic snapshot while other commands are running.

Normalization, DE and correlation keep their attempts and outcomes separate from
run ingest/quantification progress. Analysis retries, reuse and failures preserve
run diagnostics and completed counts. Progress totals include only the active
cohort; excluded runs remain in history. Status prints individual `cohort_stage`
records and errors. Existing projects receive an audited migration of legacy
analysis labels; see [analysis progress and upgrades](docs/analysis-progress.md).

Downloaded raw data, prepared references, quant outputs, and DE/corr outputs are persisted on disk and tracked in SQLite. Re-running commands reconciles existing artifacts where possible instead of redoing work blindly.

Downloads stage a complete run before publishing it. Interrupted HTTP transfers keep
durable prefix checkpoints under `.rnaa_transfers/`. Resume checks the URL, expected
size/MD5, local prefix SHA256, server validator and exact returned byte range. Without
an authoritative MD5 or a strong ETag, interrupted data is restarted. A full `200`
response replaces the partial file; a `416` response requires proof of completion.
Automatic HTTP decompression is disabled so recorded sizes and hashes describe the
downloaded bytes. Published manifests retain the HTTP source/version metadata.

Download time limits are configurable in `[download]`: `connect_timeout_seconds`
and `read_timeout_seconds` default to 30; `max_transfer_seconds` defaults to 86400.
All must be positive. Transfer deadlines and cancellation are checked between reads
and during retry backoff; a blocked read can take up to its configured timeout.
Retries are bounded by `retries`. Local SHA256 without a publisher checksum records
what was fetched; it does not independently attest the publisher's content.

When preprocessing is enabled, RNAA can prune FASTQ inputs after a successful quantification:

- `--fastq-retention both` keeps original and trimmed FASTQ files
- `--fastq-retention original` keeps only original FASTQ files
- `--fastq-retention trimmed` keeps only trimmed FASTQ files
- `--fastq-retention none` removes both original and trimmed FASTQ files and keeps only downstream quant outputs
- `--preprocess-threads N` sets per-run `fasterp` threads when preprocessing is enabled

Scheduling is budget-first:

- `--worker-budget N` sets the total CPU budget RNAA may spend on preprocessing and quantification
- `fasterp` is capped at 16 threads per run
- kallistors uses the remaining quant budget, deriving parallel run count automatically from the budget and per-run thread cap
- quantification workers in one phase share a lazy, verified read-only reference
  index; mapping and count state remain per run. See [index sharing](docs/quant-index-sharing.md).
- `--match-cache-mb N` on `quant` or `run` limits each run's exact-match lookup
  cache across its mapping threads (default 64 MiB). It does not change count identity.
- rsdeseq2/normalization is forced to a single thread
- correlation uses the full configured budget by default

Storage pressure can also gate downloads:

- set `storage.max_active_size = "100GB"` in `rnaa.toml`
- or pass `--storage-limit 100GB` to `rnaa download` or `rnaa run`
- when active project storage reaches the limit, RNAA pauses new downloads until quant cleanup frees space

Private transfer caches and retained generations count toward storage usage. Download
publication currently copies verified cache files into independent generation files,
so peak disk usage includes both copies until commit. Failed transfers may keep their
cache for retry. Downloads now reserve transfer and publication capacity. Automatic
cleanup of abandoned transfer/support data remains unfinished; the storage limit is an admission check,
not a hard quota. Tracked unpublished generation data has an explicit reversible
[cleanup workflow](docs/unpublished-cleanup.md).

`rnaa storage usage --details` reports observed logical bytes by category and,
on Unix, allocated blocks with hardlinked inodes counted once. Add `--include-shared`
to inspect shared storage as well. See [storage usage](docs/storage-usage.md) for
the concurrency flags and measurement limits.

Removed FASTQ files are permanently deleted after successful quantification and removed from the active artifact set. Their verified digests and quantification generation remain recorded so intentional retention can be distinguished from accidental deletion. Retention policy changes also apply to completed quantifications.

FASTQ retention runs as a separate owned attempt after quantification commits.
Cleanup failure preserves successful quantification and appears in attempt history;
rerunning quantification retries the policy. Cleanup verifies the committed output
generation before removing any inputs, and records removal evidence before unlinking.

Quant cache reuse checks the exact registered input set, content hashes, reference index, engine and scientific settings. CPU allocation changes preserve reuse. Legacy quant manifests without the typed request identity require re-quantification; existing output files remain available.

Preprocess QC is parsed into a normalized SQLite table and gated by default before quant/downstream analysis. The built-in defaults are conservative fail-fast thresholds intended to reject catastrophically bad files:

```toml
[quant.qc_gate]
min_pass_rate = 0.60
max_low_quality_rate = 0.30
max_too_many_n_rate = 0.10
max_too_short_rate = 0.40
```

You can override them in `rnaa.toml`, for example:

```toml
[quant.qc_gate]
min_pass_rate = 0.8
max_low_quality_rate = 0.2
max_too_short_rate = 0.2
```

Runs that fail the configured QC gate are recorded in `preprocess_qc` with a gate
reason and are excluded from quantification and downstream analysis. Eligibility
has a separate, append-only decision history; a new verdict preserves completed
run progress and previously published count files.

`--preprocess-bypass` allows an explicit fallback to raw FASTQ after an engine
error. It does not bypass QC rejection. When every read is rejected, RNAA preserves
the report, records a failed QC attempt and leaves the run excluded, including when
the QC rejection is reused from cache. Bypasses have their own report and cannot
reuse trimmed reads from an earlier preprocessing attempt.

If you change the gate after preprocessing has already been completed, rerun the gate without rerunning `fasterp`:

```bash
./target/release/rnaa qc reevaluate --root /data/my-rnaa-project
```

This recomputes eligibility from verified stored reports, including after original,
trimmed or all FASTQs have been intentionally removed. Changing a threshold or
strictness does not invalidate the count cache. New quantification requests use
revision 4; compatible revision 3 manifests retain their original bytes and retention
receipts. A stricter bypass policy still requires native preprocessing if the stored
report records an engine-error bypass.

Reports must match their registered hashes and preprocessing generation. Generated
reports also require unchanged raw sources or exact evidence of authorized cleanup.
Missing, changed or malformed evidence blocks admission and preserves the last valid
normalized QC record. Active run writers are skipped and counted as `qc_busy_runs`;
an active analysis writer must finish before manual reevaluation. A changed decision,
its history, normalized QC, event and downstream invalidation commit together.

Changed scientific quantification settings can now recover [pruned count inputs](docs/dependency-recovery.md):
RNAA prefers compatible retained trimmed reads, otherwise uses/restages a complete
original FASTQ group. Unavailable local sources block with a failing command while
prior outputs remain intact. Inspect exact producers with
`rnaa lineage --group quant --subject RUN1 --history`; add `--verify` to rehash
payloads and validate authorized pruning.

Normalization and DE admit only eligible runs. Analysis requests freeze the complete
quantified cohort, including excluded members, current count identities, policy fingerprints and report hashes;
publication rechecks them. Default correlation refuses a DE matrix from an outdated
QC cohort or policy and asks you to rerun normalize/deseq2. A custom matrix containing
an excluded project run is also rejected.

`rnaa status` includes `runs_skipped_qc` and per-run `qc_eligibility` decisions with
fingerprints and reasons.

Resolution records a durable active cohort. Resolving a narrower or empty selection
keeps previous runs, counts and analysis generations as history; those runs do not
reenter processing or analysis automatically. Changed selectors, metadata override
registrations or resolver settings require `rnaa resolve` again. Existing projects
with accession inputs need one resolution to establish membership. Projects with
explicitly supplied runs and no selectors retain their existing run selection.
Status distinguishes `runs_total` (active) from `runs_catalog_total` (history), reports
`cohort_ready`, and identifies the selected cohort and its biological sample mappings.

Biological identity comes from explicit `biological_sample_id` metadata, otherwise
`sample_accession`. Without either, the identity is labeled `RunFallback`; it does
not establish biological independence. Multiple admitted runs with the same
biological identity are pooled before gene summarization, scaling and rounding,
producing one model row per sample. Matrix columns use biological sample IDs, and
each DE generation includes `sample_metadata.tsv` with contributing run IDs.
Conflicting model metadata within one sample fails before fitting. See [the pooling
policy and its validation limits](docs/technical-runs.md) for the count/length treatment.

Analysis requests bind active membership, sample identities, canonical metadata and
every contributing run's counts/provenance. Publication rejects changes made during
computation. Default correlation checks the current DE model/membership against its
frozen producer contract, then uses current metadata for its own adjustment.
[Scientific cache identity](docs/scientific-cache.md) excludes unrelated metadata;
changes to DE's used covariates require a refreshed producer. Covariate types, factor levels,
references and missing-value rules are explicit; [typed model preflight](docs/model-designs.md)
records the actual design and rejects confounded, saturated or undersupported
contrasts before fitting. [Metadata refresh](docs/metadata-refresh.md) preserves
deliberate model edits, retains immutable source revisions and stops on conflicting
source updates. Use `rnaa metadata refresh` and `rnaa metadata show` to refresh and
inspect them; identity/selector corrections use registered overrides and resolution.

Correlation now uses bounded numerical tiles and output stripes, with bounded
top-k selection and configurable memory/output estimates checked before computing.
See [correlation limits and edge semantics](docs/correlation-resources.md) and
[independent scientific references](validation/scientific/README.md). The current
DESeq2 library remains in use; its replacement will be integrated at the end.

[Module enrichment](docs/module-enrichment.md) now scores full components with
standardized PCA, keeps full membership separately from top-gene reports, and
validates pathway overlap/options and resource estimates before graph computation.
Module graph, PCA, ranking and enrichment use the declared thread/memory limits.

For DE and correlation outputs, ENSG-keyed files remain the primary computational outputs. Annotated sidecar outputs are also written where applicable:

- `gene_counts_annotated.tsv`
- `gene_norm_counts_annotated.tsv`
- `vst_annotated.tsv`
- `edges_topk_annotated.tsv`
- `modules_top_genes_annotated.tsv`

## Configuration

See [`configs/rnaa.example.toml`](configs/rnaa.example.toml).

Shared artifact storage can be enabled by setting one of:

- `storage.shared_root` in project config
- `RNAA_SHARED_ROOT` environment variable

Configured relative shared-store paths are resolved from the project root. The
shared catalog records project/generation ownership and durable publication/reader
pins. Inspect it or preview orphan-blob collection:

```sh
rnaa storage status
rnaa storage gc
rnaa storage gc --apply
rnaa storage restore sha256:DIGEST
rnaa storage gc --apply --purge-expired
```

Collection preserves owned, busy, unreachable-project and legacy/unindexed data.
`rnaa storage import-legacy --project /data/project-a --project /data/project-b`
previews old shared ownership. Add `--apply` to import protection; releasing legacy
holds also requires `--complete-project-inventory` and verified link coverage.
See [shared ownership migration](docs/shared-legacy-import.md) for the inventory
contract and protected cases.
`--apply` moves eligible blobs to reversible trash; only `--purge-expired` deletes
trash past its recorded deadline. See [shared storage](docs/shared-storage.md) for
the ownership protocol, retained history, JSON output and remaining cleanup scope.

Project generations have a separate dependency-aware collector:

```sh
rnaa storage generations gc
rnaa storage generations gc --apply
rnaa storage generations keep GENERATION_ID --label publication
rnaa storage generations unkeep GENERATION_ID
rnaa storage generations restore GENERATION_ID
rnaa storage generations gc --apply --purge-expired
```

Current results, active artifact views, QC evidence, input bindings and pruning
proofs stay protected. Restore includes required retired dependencies. Workers and
pinned path snapshots block collection. Projects with unproven legacy ancestry
remain conservatively held. `rnaa storage generations import-lineage` previews
ancestry reconstruction from supported checksummed legacy seals/manifests; add
`--apply` to publish separate receipts without changing the original results.
See [legacy ancestry](docs/legacy-ancestry.md) and
[generation cleanup](docs/generation-cleanup.md).

Stopped jobs' tracked unpublished generation directories have their own preview,
trash, diagnostic restore and expired purge:

```sh
rnaa storage generations unpublished gc
rnaa storage generations unpublished gc --apply
rnaa storage generations unpublished restore GENERATION_ID
rnaa storage generations unpublished gc --apply --purge-expired
```

Restored files remain unpublished. Unindexed directories, referenced inputs and
unproven legacy ancestry remain held. See
[unpublished cleanup](docs/unpublished-cleanup.md).

Count-backed FASTQ/SRA inputs have reversible file cleanup:

```sh
rnaa storage files clean --mode fastq
rnaa storage files clean --mode fastq --apply
rnaa storage files status
rnaa storage files restore TRASH_ID
rnaa storage files gc --apply --purge-expired
rnaa quant --cleanup fastq --cleanup-on success --trash-days 7
```

Cleanup preserves sealed count results and their missing-input evidence. Restore
handles complete run groups; recovery can restore reads for a new fit before
trying a download. Automatic FASTQ trash takes precedence over permanent FASTQ
retention when enabled. See [input cleanup](docs/input-cleanup.md) for ownership,
scheduling, expiration and recovery contracts.

Project configuration currently uses the `[refs].ensembl` key for the preset release selector. In the current implementation this selects the matching GENCODE preset release, with `"latest"` mapping to the latest supported preset for the chosen organism.

## Worker ownership and recovery

Download, preprocessing, quantification, normalization/DE and correlation attempts
are recorded durably. Each claim has an owner, a monotonically increasing fencing
token, a 45-second lease and a heartbeat every five seconds. Opening `status` or a
second worker preserves live work. A stopped process becomes retryable once its
lease expires; legacy interrupted states without attempt records are also recovered.
Stop older RNAA processes before upgrading a project to this ownership schema.

Inspect the most recent 100 attempts, or request cancellation:

```sh
rnaa jobs list
rnaa jobs list --json
rnaa jobs cancel ATTEMPT_ID
```

Cancellation is cooperative: it rejects further worker database writes and successful
completion. An in-process scientific engine may finish its current computation
before the command observes cancellation. Writer locks still protect ingest outputs,
shared reference/configuration changes and other work awaiting full generation
conversion. They prevent a paused or expired process from overlapping a replacement
writer. Stop an unresponsive process before retrying its
work. Lock files under `.rnaa-locks/` are persistent coordination files; their
presence alone does not mean work is active, and normal recovery does not delete them.

Finite download and pipeline commands return a failure status when requested work
fails. Attempt records preserve failure details and cancellation/interruption state.

## Reference integrity

Reference preparation publishes the source FASTA, source GTF, annotation tables,
index and manifest as one verified generation. Cache reuse requires matching source
hashes, ID policy, index settings and engine identity. Custom references may be plain
or gzip files. Presets use one concrete GENCODE release with matching FASTA/GTF and
provider checksums. Existing reference folders without this provenance are rebuilt
once; previous files remain available.

`[refs]` supports `index_k = 31` and `index_threads = 1`. The k-mer length must be odd
and between 5 and 31. Changing k rebuilds the index; changing thread count preserves
cache identity. Reference preparation has its own durable job and can run inside an
analysis job. Failure or cancellation preserves the last committed reference bundle.

Fresh preparation now reserves source transfer/planning and complete build/publication
capacity from the shared pool. A bounded disk-backed scan counts distinct reference
k-mers before estimating native memory. It can need temporary disk and additional
planning time; conservative estimates still require production qualification.
`max_source_line_mb = 1` bounds decoded FASTA/GTF lines, and `max_decoded_mb = 0`
disables an optional aggregate decoded-source cap. These execution settings preserve
healthy reuse. See [reference resources](docs/resource-reservations.md#current-reference-integration).

Compressed FASTA and GTF include every gzip member. Reference request revision 2
records this complete decoding contract; older bundles rebuild once without removing
their files. Byte-identical indexes continue to use the same count-binding identity.

The project records one selected reference with separate index and analysis
fingerprints. Annotation changes refresh DE/correlation while counts remain reusable
when the index identity matches. A changed index or scientific quantification request
blocks incompatible counts from analysis. Original count files remain available.
If required reads were pruned, RNAA fails rather than declaring incompatible counts
current; recover the matching input reads before requantifying them. Automatic
reacquisition of pruned inputs remains planned.

Analysis publication and cache reuse recheck the selected reference and settings.
Default correlation refuses a DE matrix bound to a different reference. `status`
prints `selected_reference`, `reference_index_sha256`, `reference_ready` and a reason
when the prepared snapshot is stale. Status does not build references or resolve a
new upstream release; an active analysis keeps its writer ownership.

Advanced integrations supplying external transcript abundances must explicitly call
`rnaa_stages::quant::declare_external_quant(run, directory, reference, provenance)`.
The declaration validates kallisto-compatible TSV columns, complete transcript IDs,
lengths and finite nonnegative values, then binds the TSV and index hashes in
`rnaa_quant_import.json`. It records caller-supplied provenance, without claiming
RNAA produced the counts or reconstructing FASTQ history. The TSV is the computational
input for this mode; optional HDF5 sidecars are not its count source. Bare abundance
files are not admitted automatically. Revising supplied data requires a new explicit
declaration. This library API does not register runs, publish native quant outputs
or create a native quant cache hit.

## Output generations

Reference preparation, download, preprocessing, quantification, normalization/DE and correlation write complete output sets under
`.rnaa_generations/<generation-id>/payload/`. The database selects the current set
only after its required files, hashes and manifest verify and all files have been
flushed. Artifact registration, stage state, completion event and job success commit
with that selection. Engine failures and interruptions preserve the previous head.
QC rejection publishes its complete diagnostic generation with a failed attempt;
the run stays excluded. Prior output sets remain available by their permanent paths.
Normalization, DE and correlation can select an earlier verified generation when
the exact scientific request is repeated. Reuse preserves original files and
producer history while recording the new attempt. Inspect request variants with
`rnaa lineage --group de --all-requests --verify`; see the
[historical reuse contract](docs/scientific-cache.md).
Preprocessing and native counts also support historical selection after retention,
using immutable bindings when current equivalent reads have newer paths. Reuse
checks current source/processing content and QC admission; valid old deletion
receipts alone cannot make a stale source current. A compatible healthy older QC
report can repair current report evidence without redownloading reads.
If no healthy report survives, `rnaa quant` reconstructs it from a complete original
read group, downloading the selected source when necessary. Matching recovered
inputs reuse the original counts; changed inputs require new counts. QC rejection
still excludes the run. Without recoverable originals, the command reports the
missing dependency and preserves the count payload.

Downloads and reference preparation also select compatible verified history.
Download history with a publisher checksum can be selected offline. An older
unchecksummed URL snapshot requires a fresh matching strong HTTP validator;
otherwise RNAA fetches the source again. The currently selected verified download
remains usable offline. Reference reuse verifies source content, annotation policy,
index settings and engine identity before changing the project selection. Provider
metadata resolution still uses the provider's listing/checksum service.

`rnaa status` prints each current generation's group, subject, ID and payload location.
Downloaded files are inside `payload/raw/<run-accession>/`, with a manifest under
`payload/manifests/download/`.
Quant files are inside `payload/quant/<run-accession>/<engine>/`;
preprocessing files are inside `payload/preprocess/<run-accession>/fasterp/`.
DE files are inside `payload/de/<project-id>/`; correlation files are inside
`payload/corr/`, with project outputs in its `<project-id>/` subdirectory. The library
API `ProjectPaths::with_generation_heads(&db)` returns paths pinned to those sets.
Keep that `ProjectPaths` value alive while consuming its files; clones retain the
same maintenance pin. `GenerationHead::paths` and `Generation::paths` now return
`Result<ProjectPaths>` and acquire a pin before exposing paths.
Analysis inputs retain their metadata snapshots inside each generation. DE verifies
upstream quant generations when present; default correlation verifies its DE
generation. Quantification and DE verify download generations when present and
preprocessing generations when enabled.
Damaged outputs or missing seals are rejected. Legacy preprocessing files may remain
in `quant/<run-accession>/<engine>/fasterp/`; fresh work writes a generation.

Retention can remove a generated raw or trimmed FASTQ only with a durable record connecting
its producer and a committed quantification consumer of the exact bytes. The original
seal remains unchanged. Verification accepts the authorized absence while requiring
all remaining files to match; reports and manifests cannot be pruned this way.

The first normalization/DE or correlation after this schema upgrade recomputes legacy
analysis outputs into a generation; legacy files are kept. Run normalization/DE before
using the default correlation matrix after upgrading. Existing valid quant manifests,
including intentionally retained inputs, remain reusable; the next recomputation
creates a quant generation. Old and abandoned generations
remain on disk and count toward storage limits until explicitly collected. Shared
orphan collection and dependency-aware retirement of unreachable committed project
generations, supported sealed legacy ancestry import and observed logical/inode
storage accounting are implemented. Tracked unpublished generations also support
reversible diagnostic cleanup. Shared legacy ownership has an explicit verified
inventory/import workflow. Unresolved legacy data, unsealed result adoption,
unindexed/transfer/support cleanup and admission reservations remain open. FASTQ/SRA
file trash is implemented separately from generation collection. Reference
preparation publishes immutable generations.

## Documentation

- [Public-data tutorial](examples/public_yeast.md) and [offline installation check](examples/README.md)
- [Model designs](docs/model-designs.md), [technical runs](docs/technical-runs.md), and [metadata refresh](docs/metadata-refresh.md)
- [Export](docs/export.md), [survey bundles](docs/survey.md), and [analysis progress](docs/analysis-progress.md)
- [Resource reservations](docs/resource-reservations.md), [FASTQ preparation](docs/fastq-preparation.md), [quantification resources](docs/quant-index-resources.md), and [correlation resources](docs/correlation-resources.md)
- [Dependency recovery](docs/dependency-recovery.md), [generation cleanup](docs/generation-cleanup.md), [input cleanup](docs/input-cleanup.md), and [storage usage](docs/storage-usage.md)
- [Linux support](docs/linux-support.md), [architecture](docs/architecture.md), and [release preparation](docs/releasing.md)

## Repository Layout

- `crates/rnaa`: CLI binary
- `crates/rnaa-core`: schemas, config, DB, traits, state, manifests
- `crates/rnaa-stages`: resolver, downloader, refs, preprocess, quant, DESeq2, correlation
- `crates/rnaa-formats`: matrix readers and writers
- `crates/rnaa-stages/src/tximport.rs`: native kallisto abundance-to-gene aggregation
- `fixtures/`: offline test fixtures
- `crates/rnaa/tests`: CLI integration tests, shared fixtures and workflow scenarios
- `docs/systemd/`: example Linux services
- `validation/real-data-smoke/`: checksummed real-data run record for the native Rust workflow

## Contributing

Run before submitting:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

The workspace currently includes a vendored `fasterp` patch to keep the documented clippy gate passing under the workspace settings.

Direct dependency versions were checked against stable crates.io releases on 2026-11-04,
including the vendored preprocessing library's optional and development dependencies.
The audit is recorded in [validation/dependencies/2026-11-04-0.3.1.json](validation/dependencies/2026-11-04-0.3.1.json).
Run it again with Python 3.11+:

```sh
python3 scripts/check_dependencies.py
```

The audit checks both version requirements and resolved lockfile versions. Update
the root and vendor lockfiles with `cargo update` and
`cargo update --manifest-path vendor/fasterp/Cargo.toml`, then run the full checks.
Upstream requirements can retain older transitive versions: `crypto-common 0.1.7`,
for example, requires `generic-array =0.14.7`. These are not forced across incompatible
constraints. Dependabot is configured to propose weekly Cargo updates for both roots.

`mincorr` has no published stable release. Its exact Git revision is retained with
the Kendall tie-handling fixes absent from the older `v0.2.2-gpl` tag. The audit
flags a published stable GitHub release for review when one becomes available.
Engine versions are generated from Cargo.lock at build time. Scientific cache
requests include its digest, so dependency changes conservatively require fresh
quantification and analysis rather than reusing outputs from an older build.
See [vendor/fasterp/RNAA_PATCHES.md](vendor/fasterp/RNAA_PATCHES.md) for local changes
and vendor verification instructions.

## License

MIT.
