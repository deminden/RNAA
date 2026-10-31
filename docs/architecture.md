# Architecture

RNAA is organized as a contract-first Rust workspace.

## Crates

- `crates/rnaa`: CLI entrypoint and orchestration
- `crates/rnaa-core`: config, schemas, SQLite DB, manifests, state machine, traits
- `crates/rnaa-stages`: resolver, downloader, references, quantification, DESeq2, correlation
- `crates/rnaa-formats`: matrix readers and writers

## Workflow boundaries

The CLI entrypoint only opens the project, initializes logging and dispatches commands.
`crates/rnaa/src/cli.rs` owns arguments and conversion into typed settings. Ingest has
separate input/recovery, preprocessing, quantification, scheduling and retention
modules. Analysis separates normalization/DE, correlation, input admission, gene-set
selection and annotation sidecars. Shared resource policy, progress projections,
metadata projection and storage accounting have their own modules. Scheduler results
expose attempted run IDs and whether work occurred, keeping worker handles private.

Correlation separates reusable row preparation, bounded numerical tiles, output
selection/serialization and resource estimates. Top-k does not allocate a full
gene-by-gene matrix; dense and threshold writers hold one bounded row stripe.
Scientific reference generators and offline hash/numerical gates live under
`validation/scientific`; [resource limits](correlation-resources.md) are explicit
payload estimates, with global reservations still a separate scheduler concern.

The database module owns connection lifecycle, migrations, fenced reads/writes and
atomic contexts. Its `runs`, `artifacts`, `qc` and `analysis` submodules own their
respective queries and row decoding; every mutation still goes through the same
transaction/fencing boundary. Reference preparation separates source acquisition,
annotation derivation and immutable publication. Correlation separates adjustment,
network output and enrichment, retaining the existing public library imports.
Process-kill and ownership tests live beside their subsystems, outside production
implementation files; their test module paths remain stable.

## Design intent

- stage contracts live in Rust traits
- published artifacts live at permanent generation paths
- every stage writes or records enough provenance to be rerun safely
- library integrations run in-process; R is not a runtime dependency

## Native quantification and DE

- `kallistors` builds the transcript index and quantifies reads in-process.
- FASTQ and remote reference downloads use the native `reqwest` client; no curl or wget process is required.
- Quant stages keep kallisto-compatible `abundance.h5`, `abundance.tsv`, and `run_info.json` outputs. The HDF5 file is the canonical quant artifact; TSV is retained as the tximport-compatible interchange table.
- `rnaa-stages::tximport` maps transcript abundance rows through `tx2gene`, supports gene-level `no`, `scaledTPM`, and `lengthScaledTPM` count modes. `dtuScaledTPM` requires transcript-level output and is rejected. The importer normalizes GENCODE pipe-delimited headers and numeric versions while preserving `_PAR_Y`, rejects ambiguous transcript mappings and inconsistent transcript sets, and records mapping diagnostics. For `no` mode, size factors are estimated on length-corrected counts before completing the offset matrix.
- `rsdeseq2` performs size-factor estimation, NB model fitting, contrasts, normalization, and VST in Rust.
- Preprocessing, residualization, correlation, and module enrichment are also implemented through Rust crates. No R runtime or external bioinformatics CLI is required for the analysis path.
- Stage manifests record engine versions, input hashes, reference/index identity, options, and output hashes so old R/kallisto products are not reused as native outputs.

## Artifact integrity improvements

Import algorithm revision 2 fixes across-sample length scaling, zero-abundance
length imputation, offset estimation and ties-to-even rounding. Previous DE
requests are invalidated by the revision. DE reuse validates the engine and
algorithm revision plus SHA-256 digests and sizes of the samplesheet, mapping,
annotation, abundance tables and all declared output files.

JSON manifests and abundance HDF5 replacements are staged on fresh inodes and
atomically published. Shared blobs are copied into independent files; new
quantification cannot overwrite an existing blob through an aliased inode.
Reference preparation, download, preprocessing, quantification, normalization/DE and
correlation use complete generations as described below. Shared-blob ownership,
pins and reversible orphan collection are implemented. [Project generation retirement](generation-cleanup.md)
adds dependency/proof-aware GC, maintenance pins and reversible project generation
trash. [FASTQ/SRA input cleanup](input-cleanup.md) adds separate count-backed file
trash, restore and receipt-aware recovery. Migration 0019 adds
[verified legacy ancestry imports](legacy-ancestry.md) as separate immutable receipts;
unsupported or unsealed legacy history remains held. [Storage usage](storage-usage.md)
measures observed logical paths and deduplicated allocated inodes, with explicit
concurrency/error flags. Migration 0020 adds [tracked unpublished cleanup](unpublished-cleanup.md)
with immutable whole-container receipts and reversible diagnostic restore.
Shared schema 2 adds [legacy ownership migration](shared-legacy-import.md) from
explicit project inventory, read-only catalog inspection and immutable witnesses.
Unresolved legacy data, unsealed adoption, unindexed/transfer/support cleanup and
remaining producer admission remain open. Worker leases are described below.

Fresh CLI preprocessing and quantification now use a common host-local reservation
pool. Quantification separates bounded layout/input scanning, atomic cold index
plus first-run admission, resident graph ownership and per-run mapping/EM work.
Active `Arc` borrowers protect the graph; pressure-aware idle eviction and scheduler
release let independent projects and preprocessing make progress. Native packed
readers enforce preflight shape bounds, and joined mapper threads release private
thread-local caches before run guards end. Execution-only estimates and owner
descriptors are recorded separately from scientific identity. These conservative
models need production RSS qualification. Owned downloads also use separately
admitted resume/transfer/verification/FASTQ-scan/publication phases. HTTP responses
unwind before known-size readmission; unknown-size growth checkpoints and exits
on capacity refusal. Publication copies enforce their observed byte promise before
sealing and registration. References also admit source verification/transfer,
bounded decoded planning and complete native build/publication. A disk-backed
canonical-kmer counter supplies distinct-key observations without retaining all
keys in memory. Analysis and other final generation publication still need
admission. See [resource reservations](resource-reservations.md).

Quantification request revision 4 records typed semantic settings, reference index
identity and sorted input identities (including multiplicity). Execution resources
are excluded. The active registered input set must match exactly. Missing inputs
require a SQLite retention record matching the input snapshot and the specific
quantification manifest digest; a missing file alone never proves valid retention.
Original preprocessing provenance is preserved as input lineage when quantifying
retained trimmed reads. Schema migration statements and their version records now
commit in one immediate transaction, with rollback coverage for failed upgrades.

## Project stages

Normalization, DE and correlation are cohort-level operations. Their attempts,
project-stage records and frozen generation requests track analysis independently
of each run's ingest/quantification progress. Success, historical reuse, failure and
interrupted-worker recovery preserve member run records, including ingest errors
and timestamps. Cohort revisions persist explicit membership and canonical metadata.
Changing membership invalidates analysis progress without discarding retained runs.

Migration 0021 preserves old DE/correlation run labels, errors and timestamps in
immutable audit receipts before projecting those labels to `QUANT_DONE`. This is a
progress projection; existing input and generation validation still determines
whether results can be used. The migration holds the analysis writer lock and
rejects subsequent attempts to write analysis labels into runs. See
[analysis progress and upgrades](analysis-progress.md) for the operator contract.

## Durable local attempts

Migration 0007 adds job resources and attempt history. A claim checks run readiness,
increments the resource's fencing token and inserts the attempt in one immediate
SQLite transaction. Ingest claims also change the corresponding run's state;
analysis claims leave runs untouched. A partial unique index permits one live
attempt per resource. The resource covers all ingest stages for one run, or all
cohort analysis stages for one project. Reference construction has its own durable resource and writer lock so two commands
cannot build the same reference simultaneously.

Worker database handles carry their attempt identity. All database mutations now
run through a write transaction; worker mutations check owner, token, active state,
lease and cancellation inside that transaction. This includes failure handling,
artifact registration and events, so stale handles cannot overwrite a replacement
attempt. Download, preprocessing, quantification and analysis generation publication
group the complete registry update, stage state, completion event and attempt outcome
in one transaction. Reference selection publishes in the reference generation transaction.

Leases last 45 seconds and renew every five seconds. Recovery and reconciliation
use the same physical writer locks as claims. They do not infer process death from
a RUNNING state or a lock file's existence. An expired lease is recoverable only
when the previous process has released its writer lock. These locks protect remaining
direct output paths and shared reference/configuration changes while engines run
without cancellation callbacks. Converting all such writes and inputs to generations
is still needed before permitting replacement of an expired but live writer.
The implementation targets local workers sharing one project
filesystem; distributed scheduling and global resource admission are separate work.

The core library now includes a shared admission backend with atomic CPU/memory/
future-disk promises, bounded cancellation-aware waits and kernel-lock ownership.
Missing or replaced owner-file identities cannot release another worker's claim.
Its SQLite catalog supports dead-owner and orphan coordination-file recovery.
Fresh preprocessing now uses owned scan/engine reservations from a common local
pool, with explicit host policy and bounded output promises. Quantification also
uses scan/index/run claims. The operational milestone still requires phase estimates
and actual admission around analysis and publication outside downloads/references. Exact decoded
FASTQ preflight, capped normalized staging and the native preprocessing in-flight
window are delivered; see [FASTQ preparation](fastq-preparation.md). See
[resource reservations](resource-reservations.md) for the protocol
and the distinction between cooperative accounting and an operating-system ceiling.
The ingest driver cancels and joins remaining workers when it exits with an error,
preserving ownership and allocation lifetimes through shutdown.
Fresh quantification also performs a bounded decoded index-layout scan and local
loading-peak check before FASTQ staging. It carries the owned cancellation fence
through preflight/staging and between native phases. Atomic index/first-run
admission covers loading and mapping peaks; resident ownership and idle eviction
preserve the graph's allocation lifetime and progress. See [index preflight](quant-index-resources.md).

The downloader separates worker orchestration, private transfer caches, bounded
resource estimates, history selection and complete-generation publication. HTTP
resume and admission cases have separate test modules. FASTQ scanning shares one
request model across download validation, preprocessing and quantification.

Cancellation is a persisted request, checked at worker mutation/completion boundaries.
It can defer stopping a current native computation. Attempt shutdown records a
terminal outcome and restores interrupted running states transactionally. Killed
processes release physical locks automatically; their expired attempts are recovered
on subsequent status/worker checks. Legacy running states without attempt records
are handled under the same locks. Older binaries must be stopped before upgrading.

`rnaa jobs list --json` exposes recent attempts, owners, tokens, semantic request
digests, heartbeat/expiry times, terminal outcomes and cancellation status. Native
workflow tests verify request provenance for both fresh analysis and cache reuse.
Tests also cover competing threads/processes, coexistence with a downloader,
status during active work, a killed child, expired live writers, stale mutations,
heartbeats, cancellation and finite download failure exit status.

## Generation publication

Migration 0008 adds `artifact_generations` and `generation_heads`. Every fresh
reference preparation, download, preprocessing, quantification, normalization/DE and correlation attempt writes into its own UUID directory under
`.rnaa_generations/`. A generation moves through BUILDING, SEALED and COMMITTED;
unpublished work belonging to terminal attempts becomes ABANDONED during recovery.
No filesystem alias is the authority for publication. Payload paths are permanent,
so a reader holding a resolved generation keeps reading that set after a replacement.

The seal records the exact complete payload inventory with SHA-256 and size for
each file, plus stage, owner attempt and scientific request identity. Sealing checks
mandatory outputs, rejects symlinks and unsupported file types, and syncs files and
directories before writing the durable seal witness. Publication revalidates the
inventory and prepares all artifact digests outside the database transaction.
SQLite uses synchronous FULL mode. A short fenced transaction replaces the active
artifact registry, records the new generation head, marks stage/member state done,
appends completion events and records the owning attempt's outcome. The previous generation's
payload and registration receipt remain available; its records leave the live registry.
Errors and panics roll the transaction back, and escaped transaction handles close.

Generation cache reconciliation validates both the seal witness and complete inventory.
A damaged catalog head is an error, never a silent fallback to legacy paths. Successful
cache reuse updates state, registration, event and attempt atomically without changing
the generation. Default correlation verifies the current DE generation before reading
its matrix; corruption requires a DE/normalization repair. Legacy analyses require
one recomputation into generations after the upgrade; their original files are kept.
Valid legacy quant requests, including retention evidence, remain reusable without
conversion. A new computation creates a quant generation scoped to its run. Each
run resolves its own head without scanning every other run's catalog. A project
snapshot resolves all heads in one database query.

Quantification prepares shared HDF5 payloads and artifact hashes before sealing.
Its final transaction registers the shared blob, four required quant artifacts
(HDF5, TSV, run information and manifest), current head, run completion, event and
attempt success together. Existing generation files remain unchanged. DE verifies
quant and preprocessing generations before admitting runs when preprocessing is enabled.

Preprocessing request revision 2 captures raw input digests before lane staging and
verifies them after execution and before sealing. Read-only cache lookup never stages
lanes, creates directories or invokes an engine against a committed generation.
Cache hits preserve actual read counts. A fresh preprocessing generation contains
trimmed FASTQs, report and scientific manifest; its registry, QC row, run progress,
completion event and attempt outcome publish together. Empty trimmed outputs with
zero passing reads are a QC rejection, never an engine-error bypass. Rejection
publishes a complete diagnostic generation and fails the attempt atomically.
Explicitly permitted engine-error bypasses publish a report and manifest with no
trimmed files and retire prior trimmed records, so quantification consumes raw reads.

FASTQ retention acquires a separate run-resource lease after quantification commits,
including policy changes on reconciled results. Its request binds the policy and
quant manifest digest. Retention verifies the quant generation, matches every
selected input to provenance, records durable removal evidence before unlinking,
and removes the live artifact record afterward. Interrupted cleanup retries with
that evidence. Retention failure/recovery never demotes completed quant or analysis
progress. Migration 0009 adds `generation_pruned_files`: a receipt binds the producer,
committed consuming quant generation, quant manifest digest and exact input snapshot.
Quant seals now include frozen input snapshots. Retention writes both legacy evidence
and producer-aware pruning authorization in one transaction before unlinking.
Generation verification accepts absent trimmed files only with a matching registered
raw/trimmed FASTQ receipt and committed consumer seal; changed files, extra files, reports and
manifests never qualify. The original inventories/witnesses stay immutable. Historical
receipts remain valid after the consuming quant generation is superseded. Shared
blob ownership/pins and reversible orphan collection now use a separate store-wide
catalog and recoverable publication handoff; see [shared storage](shared-storage.md).
Project-generation retirement follows original/binding edges and pruning proofs;
file-specific cleanup uses its own reversible receipts and state machine, described
in [input cleanup](input-cleanup.md). Legacy ownership import remains open.

Download request revision 2 freezes the selected source descriptors and layout,
excluding execution allocation and retry policy. The worker validates the selected
basenames/layout and full FASTQ groups without composing decoded temporary files.
HTTP data and checkpoints remain in private request-scoped per-file cache directories;
complete payloads copy into independent inodes in an unpublished generation. SHA256
records, manifest, current head, VERIFIED progress, event and attempt success commit
together. Source descriptors are rechecked inside the final publication transaction.
Failed/cancelled replacements preserve the prior head and records. Completed heads
are checked for changed sources and corruption under analysis/run writer exclusion;
intentional raw pruning with a committed quant consumer remains valid.

The HTTP layer validates resume checkpoints (source descriptor, prefix length and
SHA256, strong ETag or expected MD5), discards uncommitted tails, sends Range/If-Range,
and validates Content-Range bounds, totals, response length and entity identity.
Checkpoints sync data before publishing prefix metadata at 8 MiB/one-second intervals
and ordinary interruption boundaries. A full response truncates the prior prefix;
416 succeeds only with consistent length and validator/checksum evidence. Requests
disable content decoding and require the byte representation used for verification.
Connect/read timeouts, a per-file transfer deadline, bounded retries and cancellable
backoff prevent silent infinite HTTP waits. Deadline/cancellation checks run between
reads, so a blocked read may take its configured timeout. Published manifests retain
the final effective URL, ETag, size and local byte digest. Without a publisher checksum,
this remains observed-source provenance rather than independent content attestation.

Successful publication removes its private transfer cache; interrupted caches and
abandoned generations remain for retry/diagnosis and count toward disk admission.
Tracked terminal unpublished generations can be moved into reversible diagnostic
trash; transfer caches retain their independent retry role. Copy publication
temporarily duplicates verified bytes. Atomic CPU/RAM/disk reservations and
unindexed/transfer/support cleanup remain unfinished work.

Reference preparation uses a distinct `reference:<reference-id>` resource, with its
own lease, heartbeat, cancellation and fencing. It can run while its caller owns an
analysis job; cancelled or expired parents cannot start reference work or complete
publication. Source/annotation/index payloads and their manifest publish together.
The request binds full source SHA-256 identities, normalized ID/annotation policy,
all index options and engine identity. Thread count is execution provenance. Reuse
pins a specific head and verifies its durable seal, complete inventory and request;
a loadable index or an existing annotation table alone never establishes reuse.

Custom files are hashed before copying and checked again before publication. Both
plain and gzip inputs are supported. Presets resolve one matching FASTA/GTF release,
use concrete release URLs and verify the provider's MD5SUMS before retaining source
SHA-256 snapshots. Native bounded HTTP transfers retain checkpoints on failure;
successful generations own independent copies. GTF derivation rejects malformed or
ambiguous normalized mappings; index-to-GTF mapping diagnostics remain in the
manifest. Legacy reference directories require one proven rebuild and remain intact.
Migration 0011 records the selected reference, its committed generation receipt and
stage-specific fingerprints. The index scope binds index bytes; the analysis scope
also binds transcript mapping and gene annotation. Paths and generation UUIDs are
provenance. Selection and downstream invalidation commit with reference publication;
an opaque compare-and-swap token prevents a slower producer from replacing a newer
choice, including a choice with scientifically equivalent data.
Descriptor schema, artifact roles and registrations must match the sealed inventory.
Config/source changes and corruption fail admission, preserving previous payloads.

Quant workers freeze selection before computation and verify it before publication;
the final transaction checks only catalog identity/settings. DE uses the same typed,
retention-aware quant reconciliation as quantification, under each run's writer lock.
It rejects mismatched index/engine/settings/input sets and atomically invalidates
progress/project stages with a diagnostic event. QC exclusions still preserve valid
count progress. DE requests bind the reference and per-run quant provenance files;
default correlation requires its upstream DE reference to match current selection.

External transcript abundances use an explicit `quant_import` declaration binding
TSV/index identity and caller provenance. Its parser validates complete target sets,
lengths and quantitative values. It is a supplied input contract, not a claim of
native quantification or FASTQ provenance. Mutable declarations and their TSVs are
hashed again before analysis publication. Unproven bare legacy files are rejected.

Status evaluates the prepared selection without fetching/building a new reference.
It reports readiness/reason and marks completed analysis stages pending when their
reference becomes stale, while respecting analysis writer ownership. Global shared
reference lookup, cross-project pins, complete artifact-DAG propagation/reacquisition
and generation garbage collection remain separate work.

Migration 0010 separates QC eligibility from scientific progress. Immutable decision
rows and current heads bind the enabled/strict/gate policy, verified report identity,
verdict and reason. Physical paths and producer UUIDs remain provenance and do not
change the semantic fingerprint. Report integrity and generated raw-source lineage
must verify; exact committed quant-consumption receipts authorize absent raw reads.
Manual reevaluation shares analysis/run writer exclusion and atomically updates the
decision/head, normalized QC when valid, event and downstream invalidation. It never
rewrites run progress or count artifacts. It skips busy run writers and reports their
count; an active analysis writer prevents manual reevaluation. Invalid evidence
records a blocked decision without replacing the normalized QC metrics.

Quant request revision 4 excludes gate policy and strictness from count identity.
Revision 3 reuse validates the original typed digest before projecting those two
fields out, preserving the original manifest and all retention receipts. Scientific
options and exact registered/retained inputs must still match. QC rejection does not
start another preprocessing attempt; genuine evidence/source repair may require one.
Preprocessing completion and admission are distinct: a rejected native result has
PREPROCESS_DONE progress, a rejected eligibility decision and a FAILED attempt.

DE and correlation requests bind a frozen eligibility-cohort digest, current abundance
identities and valid report snapshots for all quantified runs, including exclusions.
Before publication/reuse they recheck file hashes and policy while heartbeats can
renew. Inside the final transaction they recheck catalog membership, eligibility and
count generation IDs, keeping large file hashes outside the writer lock.
Normalization/DE filter eligibility before
fitting. Default correlation requires the upstream DE cohort and policy to match the
current snapshot; custom matrices cannot include known ineligible project runs.
Persisted user-selected cohorts and biological-sample aggregation remain unfinished.

DE uses a fresh filtered metadata snapshot for request checks instead of rewriting
the previous output directory. On computation it copies that snapshot into the new
payload. Correlation similarly freezes its metadata. Requests compare semantic input
identities while validating physical snapshot digests, allowing reuse across immutable
generation paths. Frozen biological sample/cohort identities and complete input
lineage still require the separate cohort/artifact-DAG work.

Child-process crash tests kill workers after payload writes, after sealing, inside
the publication transaction and immediately after commit. They verify one complete
current generation, consistent registry/job/stage state and preserved previous files.
Native workflow tests also reject partial replacement, damaged output inventories and
corrupted upstream matrices. The generation store includes all retained versions in
active disk accounting until explicit dependency-aware generation retirement.
Generation trash remains physical storage until its recorded window and purge complete.
