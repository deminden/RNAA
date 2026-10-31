# RNAA Development Plan

The implementation ledger below is a historical record. Release cleanup removed
superseded local logs and reports referenced by older entries; their recorded paths
are not promises that those files are present in the current tree. The remaining
implementation sequence and acceptance gates are retained below.

## Implementation ledger — 2026-10-29–30

The user authorized all improvements from the code review. Existing native migration
changes are retained. Completion requires evidence for every milestone below.

Scope-weighted progress estimate after batch thirty-nine (2026-10-30): approximately
**70% complete**, with a **65–75%** range. This includes partial delivery within the open milestones;
it is not a count of checked boxes, tests or completed implementation batches.
Independent reference gates, bounded correlation and full module PCA/resource
admission are delivered. Global resources, storage migration edge cases,
cross-project reuse, remaining CLI contracts and measured release gates account
for much of the unfinished work. Per-run/cohort progress separation also passes
its migration, native workflow and interrupted-worker acceptance checks.
Quant workers now share a verified index within one phase and use an explicit
per-run mapping-cache allowance. The shared reservation backend is verified;
pipeline-wide admission remains open until producers supply bounded phase requests.
Decoded-group preflight, exact staging bounds and the full native preprocessing
in-flight/result window are delivered. Fresh preprocessing now uses actual shared
scan/engine reservations. Owned CLI quantification also admits bounded scanning,
the cold index and first run atomically, independent resident ownership and later
run working state. Idle eviction preserves progress under memory contention.
Owned downloads now reserve cache verification, resume/transfer windows, FASTQ
validation and full extra publication copies. Reference acquisition, decoded
planning and full native build/publication are also admitted. Complete gzip FASTA
decoding has its own request revision. Analysis, publication outside downloads/
references and combined resource qualification remain open.

- [ ] Scientific correctness: import parity, shared ID policy and diagnostics,
      successful native DE/correlation fixtures, algorithm revision/invalidation.
      Implemented: importer, offsets, GENCODE IDs, mapping diagnostics and algorithm
      revision, successful native DE/correlation/enrichment workflow coverage.
      Implemented: pinned independent R import/offset, DE model/contrast/VST and
      five-method dense/edge reference fixtures with file/generator hashes and tolerances.
      Pending: final replacement DE library integration at the user's request and
      remaining scientific/release coverage outside these declared fixture cases.
- [ ] Input/configuration contracts: complete FASTQ groups, scientific quant options,
      validation, effective resource allocation.
      Implemented: complete decoded-group preflight, bounded lines, exact normalized
      staging promises and a bounded native preprocessing in-flight/reordering window.
      Implemented: common host policy and owned preprocessing scan/engine admission.
      Implemented: bounded decoded index layout preflight, local loading-peak
      rejection before materialization and ownership checks between quant phases.
      Implemented: shared quant scan/index/run admission, bounded packed input shape,
      independent resident ownership and pressure-aware idle eviction.
      Implemented: owned download phases, bounded unknown-length write windows,
      prefix-aware future promises and complete-copy publication admission.
      Implemented: complete decoded reference scans, bounded external distinct-key
      planning, source/annotation limits and full native build/publication admission.
      Pending: remaining producer admission and measured phase allocation qualification.
- [ ] Artifact validity: typed fingerprints, completed-run validation, dependency invalidation.
      Implemented: independent QC eligibility/history, policy-only reuse after FASTQ
      cleanup, frozen admission/count/report snapshots and stale-cohort rejection.
      Implemented: selected-reference fingerprints, native/imported count admission,
      reference-change publication guards and stale-reference status.
      Implemented: new sealed producer edges, file ownership indexes, supported
      pruned count-input recovery and current native count-policy admission/fencing.
      Implemented: indexed historical request lookup and verified atomic selection
      for normalization, DE and correlation, with current admission and reuse-attempt evidence.
      Implemented: full scientific request catalog, preprocessing/native count history,
      use-time input bindings and retention receipts, current source-content checks
      after pruning, and healthy historical QC report recovery.
      Implemented: historical download/reference selection with current source
      authority, immutable producer history and owned atomic reference choice updates.
      Implemented: QC report reconstruction without healthy history from original
      reads/supported sources, with independent count validity and current admission.
      Pending: generalized invalidation, legacy ancestry completion and cross-project reuse.
- [ ] Atomic publication: immutable generations/shared blobs, transactional registration,
      migrations and crash recovery drills.
      Implemented for fresh download, preprocessing, quantification, normalization/DE
      and correlation: sealed inventories, permanent generation paths, atomic
      head/registry/state/event/attempt commit, pinned readers and crash-boundary drills.
      Quant shared-blob registration commits with publication; retention has separate
      owned attempts and producer/consumer receipts for generated raw/trimmed FASTQs.
      Reference preparation also publishes a sealed generation under its own resource.
      Implemented: store-wide shared-blob ownership, pinned publication/readers,
      conservative orphan GC, recorded trash deadlines and restore/purge recovery.
      Implemented: project maintenance pins, dependency/binding/pruning-proof
      generation retirement, reversible trash/restore and receipt-backed purge.
      Implemented: count-backed FASTQ/SRA file trash, complete-group restore,
      persisted expiration, automatic timing and receipt-aware input recovery.
      Implemented: supported sealed legacy ancestry import with separate immutable
      receipts and observed project/shared logical and allocated-inode accounting.
      Implemented: tracked unpublished generation trash, diagnostic restore,
      recorded expiration and intent recovery for terminal attempts.
      Implemented: read-only legacy shared ownership import, explicit complete
      project census, immutable migration witnesses and current hardlink guards.
      Implemented: fresh download publication memory/disk admission before beginning
      a generation, capped copies/metadata and prior-head preservation on capacity failure.
      Implemented: reference generation creation follows a full copy/staging/index/
      annotation/publication claim, with original hashes and source-decoding revision 2.
      Pending: unresolved legacy evidence, unsealed result adoption,
      unindexed/transfer/support cleanup and accounting-backed admission reservations.
- [ ] Durable jobs: atomic claims, leases/owners, cancellation, separate run/cohort state,
      reliable outcomes and command exit codes.
      Implemented: persistent claims/attempt history, heartbeats, database fencing,
      writer exclusion, cooperative cancellation and finite ingest failure status.
      Reference jobs can run inside analysis ownership; publication remains fenced.
      Implemented: analysis attempts/publication/reuse/recovery preserve run ingest
      progress and diagnostics; audited legacy-state projection and database guards.
      Pending: prompt engine cancellation and broader release/crash qualification.
- [ ] Cohorts/design: persisted inclusion, biological sample grouping, typed metadata,
      rank/replication/contrast validation and frozen analysis snapshots.
      Implemented: durable active membership, independent retained run history,
      biological identity, technical-run pooling before gene import, per-sample
      model rows/outputs, frozen membership/canonical metadata and stale-publication
      fencing, with conflicting sample metadata and duplicate custom-matrix units rejected.
      Implemented: explicit factor/numeric/identifier/ignore roles, frozen shared
      model contracts, rank/residual-degree/replication/contrast preflight and
      scientific-settings publication guards.
      Implemented: controlled three-way metadata refresh, retained edit history across
      exclusion/reinclusion, atomic cohort/metadata heads and frozen source witnesses.
      Implemented: typed scientific metadata projection with separate source/config
      observations and current consumer covariates after producer-contract validation.
      Implemented: separate run/cohort progress, active-only counts/bytes/throughput,
      truthful running/waiting analysis display and quantified-only custom matrices.
      Pending: broader release/crash gates.
- [ ] Correlation: adjustment/formula correctness, reusable decomposition, bounded top-k,
      memory preflight, module score semantics and resource budgets.
      Implemented: reusable Pearson/Spearman/bicor preparation, bounded tiles,
      top-k heaps, streamed dense/threshold stripes, allocation-independent ties,
      memory/output payload estimates and independent edge references.
      Implemented: full-component standardized PCA, strict GMT/parameter/overlap
      admission, bounded mutual-neighbour module graph/connectivity, module
      thread/memory/output estimates and independent R enrichment references.
      Pending: global reservations and broader supported-host production-scale evidence.
- [ ] Operations: CPU/memory/storage reservations, index sharing, robust HTTP resumes,
      retry/stall limits and safe shared-store garbage collection.
      Implemented: one lazy immutable quant index per compatible worker group,
      bounded per-run exact-match cache and measured memory/semantic reuse checks.
      Implemented: shared admission backend with atomic CPU/memory/future-disk
      accounting, restrictive live-limit contracts, resizing and kernel-lock recovery.
      Implemented: exact decoded FASTQ/staging plans and bounded native preprocessing.
      Implemented: shared preprocessing CPU/memory/staging/output admission, typed
      non-bypassable admission errors, bounded source shape and scheduler shutdown joins.
      Implemented: decoded graph/bitmap/EC/transcript index planning and local
      peak checks without materializing compressed bitmap expansions.
      Implemented: atomic index/first-run admission, independently charged resident
      graphs and working state, joined mapper TLS teardown and constrained CLI progress.
      Implemented: download cache/resume/transfer/FASTQ-scan/publication reservations,
      capacity-aware HTTP windows, readmission without retained claims and source-byte caps.
      Implemented: admitted reference acquisition/planning/build/publication, temporary
      sorted runs with checked write windows, complete gzip staging and planning cancellation.
      Pending: analysis and other producer publication reservations and
      supported-host resource qualification.
- [ ] Maintainability/product: workflow extraction, strict matrix contracts, truthful storage
      documentation, useful doctor/status JSON, export and survey.
      Implemented: doctor/status JSON, frozen-generation TSV/CSV export, and
      compatible frozen survey bundles with metadata summaries and provenance.
      Pending: complete freshness validation, export/survey library workflows,
      and remaining CLI contracts.
- [ ] Release gates: scientific fixtures, corruption/invalidation tests, crash/concurrency
      tests, measured scale benchmarks and clean-install quickstart.
- [x] Current stable registry dependencies: update direct workspace and vendored
      requirements plus both lockfiles, check minimum Rust versions, verify engine
      provenance and the upgraded native workflow. One unpublished Git dependency
      remains an explicitly documented exception; older transitive versions required
      by upstream constraints are recorded, not silently forced incompatible.

Survey implementation (2026-10-31): sealed DE and correlation generations can be
combined into a frozen survey after checking their scientific identities, VST
lineage, sample and gene identities, and correlation design covariates. Summaries
use frozen biological-sample metadata; publication is resource-admitted, fenced,
and no-replacement. Current-settings freshness, library workflow extraction,
crash-boundary cleanup and broader platform qualification remain open. See
`docs/survey.md` and the corresponding validation evidence.

Export integrity follow-up (2026-10-31): each source stream is now hashed while
it is consumed, detecting temporary changes even when source contents are later
restored. A subprocess-exit test verifies that a rename followed by uncommitted
job completion leaves the database transaction uncommitted and the complete
directory available for inspection. Automatic reconciliation and power-loss
durability remain open.

Export reservation follow-up (2026-10-31): the final export/survey manifest now
has a distinct admitted disk allowance enforced during writing, and manifest and
summary serialization check cancellation while writing. Oversized manifests
fail before publication; inspection inventory memory is still not fully bounded.

Interrupted export follow-up (2026-10-31): new export/survey attempts record the
staging directory and filesystem identity. A later export or survey holding the
analysis lease preserves matching interrupted staging under a recovery name,
records already-renamed output as orphaned, and leaves conflicting or untracked
paths alone. Recovery is non-destructive and records an event. Legacy/unregistered
leftovers, automated disposal and power-loss qualification remain open.
The shared-store ownership compatibility boundary also recognizes migration
0022, which adds no shared-blob references. Failed ingest finalization retains
the original error when reporting a secondary job-state error.

Current-cohort guard (2026-10-31): export and survey can optionally require the
selected generations' recorded cohort digest to match the currently resolved
membership, with catalog revalidation at publication. This is a membership
check only; eligibility, metadata, reference and full scientific request
freshness remain open.

Current-admission guard (2026-10-31): export and survey can optionally compare
their recorded eligibility/abundance and reference identities with existing
current evidence, alongside cohort membership. Evidence is observed rather than
refreshed and the catalog is rechecked at publication. Metadata values, model
settings and full scientific request freshness remain open.

Current-design guard (2026-10-31): export and survey can additionally rebuild the
recorded DE and correlation adjustment designs against current used sample
metadata. Current contrast selections, correlation method, gene/pathway inputs,
other scientific settings and full request reconstruction remain open. Design
scratch and memory estimates are conservative but not fully bounded.

Current-request guard (2026-10-31): export and survey can reconstruct current
supported DE and correlation requests, validate native or declared external
quantification evidence, and compare current contrasts and correlation inputs.
Unpersisted module arguments and independently supplied matrices that do not match
the selected DE matrix fail closed. External files remain unlocked; support for
those unpersisted modes and stronger concurrency/resource guarantees remains open.

Verified first implementation batch (2026-10-29): all 45 workspace tests pass,
including seven import/offset regressions, local artifact corruption detection,
failed-publication and legacy-hardlink preservation, and end-to-end shared HDF5
snapshot preservation plus effective preprocessing thread allocation. Formatting,
strict workspace Clippy and `git diff --check` pass. DE requests now carry import
revision 2 and validate the scientific manifest's local input/output digests.
Shared blob publication uses separate inodes; full multi-file atomicity remains open.

Verified second implementation batch (2026-10-29): all 53 workspace tests pass;
formatting, strict Clippy and whitespace checks pass. Shared FASTQ preparation now
consumes every single-end lane and matches/validates paired lanes, counts and IDs.
Tests cover duplicates, missing mates, mismatched IDs, truncation, gzip members,
newlines at lane boundaries and native multi-lane processing. Multi-lane staging
uses bounded memory but requires decoded temporary disk space; streaming reader
composition and admission reservations remain operational work.

Quant scientific options now expose strand and single-end fragment mean/SD in
configuration and CLI; values are validated and included in reuse checks. Completed
quantification is revalidated and invalidated on changed available input hashes,
scientific settings or outputs, marking downstream stages pending. Preprocessing
has hash-verified source/output manifests with source mutation and output corruption
regressions. Failed requested quantification returns a nonzero command status,
verified with an incomplete paired input. Cross-project fingerprints, new-file-set
identity, retained/deleted input ownership, and QC-policy invalidation still need
completion; these incremental checks do not yet establish the full artifact DAG.

Verified third implementation batch (2026-10-29): all 57 workspace tests pass;
formatting, strict Clippy and whitespace checks pass. Quant request revision 3 uses
a typed semantic request and content fingerprint, excludes CPU allocations, and
compares the exact registered input set. Input snapshots are captured before work
and checked again before completion. Required quant outputs must match the manifest.
Historical preprocessing inputs remain recorded as lineage.

SQLite migration 0006 records intentional input removal against the exact quant
manifest generation and artifact snapshot. Retention records are written before
unlinking. Cache validation permits absent files only with matching generation
records; accidental deletion and records from another generation fail. Retention
policy changes also execute on completed/reconciled results. Migrations now run in
one immediate transaction, with an injected-failure rollback test.

Native integration coverage now checks newly registered lanes, reuse after raw/trimmed
retention, CPU-only reuse, changed raw reads, and tightening QC after completed
quantification. The latter reuses the stored preprocessing report and changes run
eligibility. Full eligibility fingerprints independent of quantification (including
QC policy changes after all FASTQs have been pruned) remain required.

Verified fourth implementation batch (2026-10-29): all 61 workspace tests pass;
formatting, strict workspace Clippy and whitespace checks pass. DE and correlation
requests now snapshot their inputs before computation and verify them again before
sealing their manifests. Required output coverage, every output digest, canonical
metadata, per-run abundance, annotation, gene-list content, pathway content, formula,
adjustment and enrichment options participate in reconciliation. Malformed request
descriptors and corrupted outputs cause recomputation. Selection, adjustment and
enrichment errors now enter correlation failure handling.

Correlation honors `adjust=none`. Residualization shares the DE formula parser,
supports interactions and intercept removal, performs one decomposition per design,
and rejects rank deficiency and insufficient residual degrees of freedom. Strict
matrix validation rejects duplicate/empty IDs, nonfinite values and mismatched
dimensions before writing. The 500-gene/eight-sample deterministic native fixture
runs DE, correlation and enrichment successfully and exercises cache reuse plus
mutation/corruption invalidation. It is not an independent statistical oracle.

Dependency verification (2026-10-29): the registry audit checks 49 packages across
58 direct requirements, including vendored optional/development dependencies, with
zero failures. `kallistors` is now 0.6.0, `rusqlite` 0.41.1 and `uuid` 1.29.0;
vendored `gzp` is 2.0.4 and `ureq` 3.5.0. Both lockfiles were refreshed. The 82
standalone preprocessing library tests pass, including local HTTP and parallel gzip
round trips. Multi-member gzip reading is fixed and verified through preprocessing,
preventing silent loss after the first member. Full Cargo metadata contains no
prerelease dependency versions and no declared dependency MSRV above Rust 1.97.1.

`mincorr` has no published stable release. The pinned revision with Kendall tie fixes
is retained rather than moving backward to the older license-qualified tag. The
audit records that exception and flags future published stable GitHub releases.
`crypto-common 0.1.7` requires `generic-array =0.14.7`, so that transitive version
cannot be updated independently. Engine versions now come from Cargo.lock; its
digest participates in scientific reuse keys to invalidate older dependency builds.
The timestamped audit is `validation/dependencies/2026-10-29.json`; the repeatable
check is `scripts/check_dependencies.py`. Dependabot is configured for both roots,
and CI verifies locked workspace builds plus vendor library tests.

Verified fifth implementation batch (2026-10-29): all 70 workspace tests pass;
formatting, strict workspace Clippy and whitespace checks pass. Migration 0007 adds
persistent job resources and attempt history. Readiness, a monotonic fencing token,
the attempt and a run's RUNNING state are claimed in one immediate transaction.
Leases last 45 seconds and heartbeat every five seconds. Fenced worker handles
check ownership, expiry and cancellation inside every database write transaction,
including failure handling, artifact registration and events. Failure details and
scientific request identities are preserved in attempt history, including analysis
cache reuse. `rnaa jobs list [--json]` and `rnaa jobs cancel ATTEMPT_ID` are available.

Download, preprocessing, quantification, normalization/DE and correlation now acquire
worker ownership. One run resource excludes concurrent ingest writers; one project
analysis resource excludes overlapping DE/correlation writers. Reference construction
has its own file lock. Recovery, quant validation and download reconciliation share
writer exclusion, so status and a second worker cannot reset or prune active work.
Claims occur before scheduler attempted-set insertion, avoiding permanently skipped
ready work after a competing claim. Finite download workers wait for concurrent
download claims and return errors for failed work; orchestration rejects ingest/quant
failures before reporting completion.

Executable evidence includes competing threads executing one claim, two actual CLI
download processes sharing one HTTP transfer, status while work is active, native
quantification coexisting with a separate downloader, heartbeat renewal, CLI
cancellation, and an actual child process killed before expiry/recovery/retry.
Expired live writers remain excluded, and stale database handles cannot alter
replacement run states, artifacts or events. The concurrency integration suite
passed ten consecutive runs. It exposed a startup metadata-backfill race; backfill
now acquires an immediate transaction before reading and rolls back earlier rows
if a later legacy record is invalid. Empty metadata needs no repeated backfill.

At this batch, physical writer exclusion remained necessary because engines still
wrote direct output paths. Cancellation rejects subsequent mutations and successful attempt
completion, but a native engine may finish its current computation before observing
it. Full generation publication is required before replacing an expired but live
writer; individual fenced mutations do not yet make the whole stage completion
atomic. The durable-jobs milestone therefore remains open. Old RNAA processes must
be stopped before upgrading to the ownership schema.

Verified sixth implementation batch (2026-10-29): all 80 workspace tests pass;
formatting, strict workspace Clippy and whitespace checks pass. Migration 0008 adds
generation catalogs and current output heads. Fresh normalization/DE and correlation
compute into unpublished UUID directories, seal the exact complete file inventory,
verify required files and artifact digests, sync payloads/directories and retain a
durable seal witness. SQLite now uses synchronous FULL mode. One fenced transaction
publishes the head, replaces live artifact records, completes stage/member state,
appends completion events and succeeds the attempt. A failed callback, panic or
publication error rolls back every database change. Transaction handles retained
past their callback close and reject later use.

Readers explicitly resolve and pin generation paths. Successful cache reuse keeps
the same generation and finalizes state/registration/event/job changes atomically.
Old payloads and their registration receipts remain available after replacement;
unpublished generations attached to terminal attempts become ABANDONED at recovery.
Missing seals, inventory corruption and incomplete catalog heads are rejected.
Default correlation verifies its upstream DE generation before reading the matrix.
Legacy analyses require one recomputation into generations after the upgrade;
original legacy files remain untouched. `rnaa status` exposes current payload paths.

DE request checks use temporary filtered metadata instead of rewriting a published
snapshot. Computation copies that snapshot into its generation; correlation also
freezes metadata. Cache identity uses semantic roles/names/hashes/sizes, while both
current inputs and prior physical snapshots still verify. Size factors now have an
explicit artifact kind and participate in publication. Disk admission accounting
includes the generation store exactly once, including retained and abandoned data.

Executable evidence includes successful normalization -> correlation and DE ->
correlation/enrichment workflows, reuse without rewriting prior generations, partial
DE/adjustment failure preserving previous results, damaged-upstream rejection and
repair, cancellation before publication and complete transaction rollback. Real
child processes are killed after payload writes, after sealing, inside the final
transaction and immediately after commit; each boundary exposes exactly one complete
generation with consistent registry/job/stage state and preserved previous files.
These checks do not establish independent statistical parity, power-loss simulation,
full artifact-DAG lineage, prompt engine cancellation or generation GC safety.

Verified seventh implementation batch (2026-10-29): all 85 workspace tests pass;
the complete mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. Quantification now produces a complete per-run generation
under `payload/quant/<run>/<engine>/`. Run subject and stage must match the owning
attempt; required HDF5, TSV, run information and manifest files must match the seal.
TSV abundance has an explicit artifact kind. Shared blob copying and digest checks
occur before publication; the shared-blob database record, artifact references,
current head, run completion, event and successful attempt commit together. A failed
registration/completion callback rolls the entire database change back. Retired
quant records leave the live registry while prior payloads and receipts remain.

Quant workers refresh their run's pinned head without scanning every other run's
catalog. Whole-project readers resolve all heads in one database query. Status
prints group, subject, generation ID and payload location. DE verifies current
quant generations when present. Incomplete/damaged seals invalidate quant reuse;
legacy typed quant manifests remain reusable, including after authorized input
retention, without requiring unavailable reads to be regenerated just for migration.
Preprocessing still writes its separate canonical directory, never beneath the
immutable quant payload. Its own generation conversion remains pending.

FASTQ retention now claims a separate fenced run-resource attempt after quant
publication, including policy changes on completed/reconciled outputs. The request
binds the policy and quant manifest digest. Removal verifies the committed quant
generation and input snapshot, writes durable retention evidence before unlinking,
then removes the live artifact record. Retention failure/cancellation/recovery never
changes completed quant/analysis progress; its own attempt preserves the outcome.
Empty cleanup has no new attempt. Invalidation of completed quant and downstream
stage states/event also now commits atomically.

New executable evidence kills quant writer child processes at payload, seal,
in-transaction and post-commit boundaries; each exposes a complete prior/new output
set, consistent run/job/event state and preserved pinned readers. Tests also reject
wrong-run registrations, roll back shared registration on completion failure, fence
retention against competing producers, preserve completed states after failed or
cancelled cleanup, and execute native failed/successful replacement followed by
retention and cache reuse. A missing seal after input retention prevents success.
The official registry/GitHub audit was rerun: 49 registry packages across 58 direct
requirements, zero failures, one documented unpublished Git dependency exception.

Verified eighth implementation batch (2026-10-29): all 97 workspace tests pass;
the complete mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. Preprocessing now publishes per-run generations under
`payload/preprocess/<run>/fasterp/`. One fenced transaction registers trimmed files,
report and scientific manifest, records QC, updates progress/events and terminalizes
the attempt. A legitimate QC rejection publishes its complete diagnostic generation
and FAILED outcome together; engine failure/interruption preserves the previous
head. Prior payloads and registration receipts remain available.

Preprocessing request revision 2 snapshots raw inputs before lane staging, verifies
them after native execution and checks the frozen snapshots again before sealing.
Read-only reconciliation never stages lanes, creates directories or invokes an
engine against published files. Cached results retain actual passed/failed counts.
Zero passing reads with empty trimmed outputs are a failed QC result, including
under an explicitly permitted engine-error bypass; they are no longer mislabeled
as a preprocessing exception and fed raw to quantification. Explicit engine-error
bypasses publish their own report/manifest with no trimmed files and retire prior
trimmed records. Changing back to strict policy rejects an oversized input without
rewriting either prior preprocessing or quant generations.

Migration 0009 records `generation_pruned_files`. New quant seals include frozen
input snapshots. Before unlinking a generated trimmed FASTQ, retention verifies
its registered producer receipt and a committed quant consumer of the exact input,
and atomically records both retained-input evidence and producer/consumer pruning
identity. Generation verification accepts that authorized absence while requiring
all remaining payloads and seal witnesses to match. It rejects changed or extra
files, wrong consumer/digest/kind, missing reports/manifests and legacy evidence
without a generated-producer receipt. Historical pruning evidence survives consumer
replacement. Inventories also reject duplicate paths and invalid digests.

Quantification and DE verify preprocessing generations when enabled. Manual QC
reevaluation verifies report registration hashes and preprocessing generations,
including authorized pruning. It shares analysis/run writer exclusion, reports
busy run writers, and atomically updates each run's QC, progress, events and
downstream invalidation. A passing reevaluation clears an earlier QC error.
Corrupt reports cannot overwrite prior QC verdicts; active writers are preserved.
This does not yet make eligibility independent of count fingerprints.

Executable evidence includes real child-process kills at preprocessing payload,
seal, publication-transaction and post-commit boundaries, plus retention intent,
unlink and registry-removal boundaries. Retry preserves complete prior/new heads,
consistent run/attempt/event state, pinned readers and authorized pruning. Tests
inject pruning-registration failure and verify both evidence tables roll back,
publish rejected QC atomically, reject changed source snapshots/duplicate inventories,
and exercise native cached rejection, engine-error bypass without stale trims,
strict failure, immutable replacement, retention/cache reuse and damaged-report
rejection. The restart fixture now supplies genuinely passing reads; the old
four-base/low-quality fixture had incorrectly depended on bypassing QC rejection.

Verified ninth implementation batch (2026-10-29): all 116 workspace tests pass;
the complete mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. Downloads now publish complete per-run generations under
`payload/raw/<run>/`, including a manifest under `payload/manifests/download/`.
Download request revision 2 binds source selection and layout, excludes execution
allocation/retry policy, and checks the latest source descriptors inside the final
publication transaction. Raw SHA256 records, manifest, head, VERIFIED progress,
completion event and attempt success commit together. Failed/cancelled replacements
preserve prior files and registry entries. Metadata changes or corruption invalidate
completed download heads under writer exclusion. Matching complete generations reuse
without network access. Invalid/incomplete FASTQ groups never become published runs.

Private request-scoped transfer caches remain outside generations. Each source has
its own directory, preventing payload basenames from colliding with another file's
partial/checkpoint/receipt files. Verified files copy into independent generation
inodes before sealing; success removes its private cache. Both generations and
transfer caches count toward disk usage. The temporary duplicate-byte peak and
remaining abandoned caches are documented; global CPU/RAM/disk reservations, cache
GC and hard-quota enforcement remain unfinished.

A shared HTTP layer now records durable prefix length/SHA256 checkpoints, binds
URL/expected size/MD5 and server identity, and discards uncommitted or corrupt tails.
Resume requires authoritative MD5 or a strong ETag. Range/If-Range requests validate
Content-Range start/end/total, Content-Length and entity identity. A full response
replaces the prefix; 416 requires matching complete length and validator/checksum
evidence. Truncated bodies resume, changed entities restart, and invalid encoded
representations cannot be silently decompressed into mismatched bytes. Manifests
retain effective URL, ETag, sizes and local digests. A local SHA256 without a publisher
checksum remains observed-source provenance, not independent content attestation.

Download settings expose positive connect/read timeouts (30 seconds each by default)
and a per-file transfer deadline (86400 seconds by default), with bounded retries
and cancellable backoff. Cancellation/deadline checks occur between reads, so a
blocked read may take its configured timeout. FASTQ validation checks ownership
periodically and consumes all lanes without decoded temporary staging. Reference
source transfers also use the bounded HTTP layer, but reference generation/ownership
and reference-specific request validation remain pending.

Generated raw FASTQs now use producer/consumer pruning receipts alongside generated
trimmed files. Quantification/preprocessing/DE verify download generations when
present. Native end-to-end coverage downloads actual FASTQ, preprocesses, quantifies,
removes both FASTQ sets, and reuses all verified heads with the HTTP server gone.
Additional evidence includes real download writer kills at payload, seal, transaction
and post-commit boundaries; a real HTTP child killed after an 8 MiB checkpoint and
resumed to the exact full payload; corrupt/unvalidated prefixes; partial/malformed
ranges; changed ETags/totals; ignored Range requests; truncated bodies; short ranges;
416; MD5 verification; encoded responses; stalls/slow trickles and cancellation.
CLI tests cover multi-file failure/retry preserving previous readers, competing
processes, cache reuse, same-size corruption, cancellation and source metadata changing
mid-transfer. A nonblocking accepted-socket race in the macOS test server was fixed;
the HTTP concurrency suite passed three consecutive runs before the full sequence.
The official dependency audit was refreshed after verification: 49 registry packages
across 58 direct requirements, zero failures, one documented Git-only exception;
both lockfiles still match the audit.

Verified tenth implementation batch (2026-10-29): all 123 workspace tests pass;
the full mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. Migration 0010 adds immutable eligibility decisions and
current per-run heads. Policy, verdict/reason and verified report identity have a
separate semantic fingerprint; physical report paths and producer generation IDs
remain provenance. Decision/head, normalized QC when valid, event and downstream
invalidation commit together. Repeating the same assessment adds no history row.
QC changes never demote or restore scientific progress or rewrite count artifacts.
A completed native preprocessing result has PREPROCESS_DONE progress even when its
QC decision rejects it and the producing attempt fails. Reusing a rejection starts
no new engine attempt. Engine failures still retain their own failure state.

Eligibility verifies registered report hashes, current preprocessing generations and
frozen raw-source lineage. Added/changed sources or unauthorized absence require
repair; exact committed quant-consumption evidence authorizes intentional cleanup.
Missing/malformed/changed reports block admission without overwriting the last valid
normalized metrics. Legacy report imports retain registered-digest verification;
they cannot reconstruct raw-source provenance that was never recorded.

Quant request revision 4 removes gate policy and strictness from count identity.
Revision 3 manifests first validate their original typed request digest, then project
only those policy fields out for compatibility. Original manifest bytes, generation
heads and retention receipts remain intact. Scientific options and exact input sets
still invalidate counts. Tightening/relaxing policy after each of all four retention
modes rejects/restores eligibility with unchanged counts and no new engine attempt,
including when every FASTQ is gone. Stricter bypass policy can still require native
preprocessing; admitting an engine-error bypass never admits zero-passing native QC.

Normalization/DE exclude ineligible runs before fitting, preserving completed counts
when QC evidence is damaged. Analysis requests freeze quantified membership, all
eligibility fingerprints (including excluded members), current abundance identities
and valid report hashes. Execution and publication/reuse transactions recheck policy,
membership, evidence and count identities; publication also rechecks request inputs.
A concurrent count replacement or newly quantified member cannot publish an analysis
of the previous cohort. Default correlation refuses a DE matrix from an old policy
or cohort; a custom matrix cannot include a known ineligible project run. Status
exposes each run's decision/fingerprint/reason and the skipped count.

Executable evidence adds decision/head/event rollback, corrupt-head rejection,
semantic identity under report relocation, strict-bypass/zero-read behavior, frozen
snapshot detection of changed counts/reports/policy/membership and busy writers,
legacy-request digest/unknown-field corruption rejection, policy-only reuse under
all retention modes, and native analysis exclusion/stale-matrix/relaxation workflows.
The official dependency audit was refreshed: 49 packages across 58 direct requirements,
zero failures, one documented unpublished Git exception; both lockfiles match it.

Verified eleventh implementation batch (2026-10-29): all 133 workspace tests pass;
the full mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. Reference preparation now has a distinct durable resource,
lease, heartbeat and cancellation/fencing boundary that works inside an owned
analysis. Source snapshots, GTF-derived tables, index and manifest publish as one
sealed reference generation. Head, registry, completion event and successful attempt
commit together. A failed build, cancelled parent or failed callback preserves the
prior head and pinned readers. Reuse verifies the exact pinned generation and request;
no attempt or manifest rewrite occurs on a normal cache hit.

Reference request revision 1 binds full source SHA-256 identities, normalized ID and
annotation policy, explicit index options and engine/dependency identity. Threads
remain execution provenance. `refs.index_k` (odd 5..31) and positive `index_threads`
are validated. Plain and gzip custom sources are supported and checked before/after
copying and before publication. GTF parsing rejects malformed rows and ambiguous
normalized transcript mappings; index-to-annotation mapping diagnostics are recorded.
Presets select a matched transcript/GTF pair, resolve latest to a concrete release,
use the correct concrete release directory and require the provider MD5SUMS. Source
transfers use the existing bounded/resumable HTTP path and retain failed checkpoints.
Old unproven reference directories require one rebuild and remain untouched.

New evidence includes two competing reference CLI processes publishing one complete
bundle, native reuse despite thread changes, option/source changes, repair of a
loadable index with the wrong settings, annotation corruption, malformed/conflicting
GTF rejection, source mutation during preparation, parent cancellation, independent
reference ownership and stale-handle fencing. Real child producers are killed after
payload writes, sealing, inside publication and after commit; every boundary keeps
one complete head, consistent registration/events/attempts and usable old readers.
A failed reference publication callback rolls all replacement changes back.

The user additionally requested smart refactoring of large files. The CLI entrypoint
is now 147 lines (previously 5,145), with separate command contracts, ingest/analysis
workflows, input admission, retention, resource policy, progress and storage modules.
The scheduler returns only attempted run IDs and whether work occurred instead of
exposing its internal worker state. Database lifecycle/transactions remain in the
391-line db module (previously 1,554); run metadata, artifacts, QC and analysis queries
now own their row decoding in explicit submodules without bypassing fenced access.
Reference service/source/annotation code is separated, and correlation separates
residualization, network output and enrichment while preserving public library
imports. Generation/job tests are separated from production implementation with
unchanged process-test paths. Top-level and quant command help match the previous
binary exactly. Existing native workflows and corruption/crash tests all still pass.

Large source/report hashes are checked before analysis publication while heartbeats
can renew. The final transaction rechecks catalog membership, eligibility and count
generation IDs without rehashing large payloads under SQLite's writer lock. Reference
requests also freeze before annotation/index computation, and source hashes stay
outside the final writer transaction. The registry audit was refreshed: 49 packages
across 58 direct requirements, zero failures, one documented Git-only exception.

Verified twelfth implementation batch (2026-10-29): all 141 workspace tests pass;
the full mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. Migration 0011 adds selected reference state. Selection
binds its exact committed reference generation and verifies descriptor schema,
artifact roles, registration receipts and sealed inventory. Separate index and
analysis fingerprints bind index bytes versus index/mapping/annotation bytes;
physical paths, producer UUIDs and execution threads remain provenance. Selection,
semantic-change invalidation and event commit atomically with fresh reference
publication. Cache selection uses the same catalog contract. An opaque choice token, separate from scientific fingerprints, prevents late
producers from overwriting newer choices, including scientifically equivalent data. Failed callbacks
roll selection, head, registry and invalidation back together.

Quantification pins the selected reference, records its binding in provenance,
checks source/settings/files before publication and checks catalog identity inside
the final transaction. Analysis freezes the full reference dependency. DE admission
now shares quantification's typed request, exact input/output and retention checks,
so changed references, engines or scientific settings cannot admit stale native
counts. Admission holds each run's writer lock, refreshes the queued record and
atomically invalidates incompatible progress/stages/events. QC exclusions remain
independent and preserve completed counts. Default correlation rejects an upstream
DE reference mismatch. Publication/reuse guards reject reference changes during
computation without hashing large source files inside the writer transaction.

Externally supplied abundances require an explicit `quant_import` declaration through
the library API, with caller provenance and exact TSV/index hashes. The constructor
validates unique columns, complete transcript IDs/lengths and finite nonnegative
values; it never invents native quantifier or FASTQ history. The TSV is its count
source; optional HDF5 attachments are not validated as equivalent count data. Bare
legacy files are rejected. The existing statistical fixtures now declare their
simulated/external inputs honestly, and changed input values require a renewed
explicit declaration. DE requests freeze the per-run provenance file as well as the
abundance bytes. A general import CLI and immutable external-input publication are
not supplied by this declaration API.

Status reports selected identity, index digest, readiness and failure reason without
building references or fetching latest provider data. Stale completed project stages
become pending under analysis exclusion; live analysis ownership is preserved.
Manual reference configuration changes are synchronized before `refs prepare`.

Executable evidence covers selection/invalidation rollback, invalid artifact roles,
unsupported/corrupt descriptors, stale selection races, source/settings mutation and
annotation-only scope compatibility. A reference job publishes a new choice while
an analysis owns its resource; the old analysis cannot publish a head or artifacts.
Native counts with all FASTQs pruned remain byte-for-byte reusable after annotation
changes, while index changes block analysis and preserve old count generations.
Native DE/correlation workflows reject stale reference matrices, refresh annotations
and reject incompatible supplied counts. External declarations reject changed data,
run IDs, reference indexes, invalid targets, nonfinite values and tampered digests.
The official dependency audit was refreshed with unchanged compatible lockfiles.

Verified thirteenth implementation batch (2026-10-29): all 154 workspace tests pass;
formatting, strict Clippy and whitespace checks pass. The official dependency audit
still reports 49 registry packages / 58 direct requirements current, with the same
documented Git-only dependency exception and unchanged compatible lockfiles.

Migration 0012 records immutable cohort descriptors and a selected head separately
from run/result history. Resolution publishes run metadata and the selected member
list in one immediate transaction with a compare-and-swap guard on prior selection
and selectors. Narrow and empty results remain narrow and empty; old runs, progress,
counts and analysis generations survive. The migration also preserves accession
roles separately, fixing silently ignored exclusions when an accession already
existed as a query or inclusion. Overlapping query responses are deduplicated, and
conflicting metadata is rejected. Registered override files are snapshotted and
rechecked during resolution; their hashes enter the resolver manifest.

Resolver attempts now share analysis writer exclusion and durable heartbeat/fencing
with fitting. Changed accession/selector intent, override registrations or resolver
settings require a new resolution. Download and preprocessing/quantification queues,
normalization/DE, correlation, samplesheet refresh, QC reevaluation and progress
reporting use active membership. Existing projects with accession inputs require
one resolution; selector-free projects with explicitly supplied runs retain their
previous selection. Inactive failed downloads are not retried, and inactive QC
records cannot invalidate a current analysis through ordinary status/worker refresh.
Status separates active run totals from the retained catalog and reports readiness,
selected cohort and biological identity mappings.

Biological identity uses explicit `biological_sample_id` metadata, then sample
accession, with an explicitly labeled run fallback when neither is available.
Multiple admitted sequencing runs with one biological identity fail before fitting;
they cannot silently become independent biological replicates. The statistical
runner independently checks identity columns in its supplied sheet. This is a
guard, not implemented technical-run aggregation or proof of independence for
fallback IDs. Analysis requests now bind a frozen active membership and stored
metadata projection, together with existing count/QC/reference inputs. Scientific
fingerprints exclude cohort revision UUIDs, but publication checks the concrete
head so equivalent later resolutions still fence superseded producers.

New executable evidence includes resolution from the local ENA fixture into full,
narrow and empty cohorts; changing selectors during a fetch; competing resolver/
analysis ownership; role preservation; resolution rollback; identity/design metadata
changes; and a changed selection rolling back an otherwise sealed analysis commit.
A native statistical workflow narrows eight supplied samples to six, refuses its
old matrix, rebuilds the selected result, preserves historical generations,
ignores inactive failures/QC history and rejects duplicate biological identities.
The larger native-analysis integration fixture now lives in a reusable test support
module; resolver tests also moved into their own module without changing test paths.

Verified fourteenth implementation batch (2026-10-29): all 162 workspace tests pass;
the full required formatting/strict-Clippy/test sequence and whitespace checks pass.
The refreshed official dependency audit still reports 49 registry packages / 58
direct requirements current, zero failures and the documented Git-only exception.

The provisional duplicate-run rejection is replaced with a typed `SamplePlan` and
one statistical unit per biological sample. Requests record identity sources,
every contributing run and the group metadata. Output matrices and size factors
use biological sample IDs; `sample_metadata.tsv` is a registered, sealed DE output
with JSON-encoded run membership. Model/contrast variables must agree across a
sample's contributions; other varying fields are blank in the sample projection
and remain available in the frozen run-level input. Single biological units cannot
be inflated into multiple observations by their sequencing-run count.

Import revision 3 implements the explicit `fractional_transcript_mass_v1` policy:
sum fractional transcript counts and count/effective-length mass, reconstruct
pooled lengths/TPM, then summarize to genes, correct/scale and round. Population
mean lengths are calculated across biological samples. Single-run samples retain
their original abundance estimates. Zero targets retain positive lengths without
inventing abundance. Invalid values, conflicting transcript sets/lengths,
normalized-ID collisions and reused contribution IDs fail. See
[the exact method and validation limits](docs/technical-runs.md). This is pooling
of estimated tables, not a claim of joint-FASTQ quantification or bootstrap parity.

DE now consumes the captured canonical metadata rows rather than rebuilding model
values silently from the database. Registered identities must agree with any
explicit biological IDs in that sheet. Both its source snapshot and the frozen
sample projection bind the request and are checked before publication. Default
correlation uses the sample metadata frozen with its DE matrix and refuses stale
metadata, QC/cohort or reference dependencies. Explicit paths to that selected
matrix retain the producer proof. Separately supplied matrices support known
biological IDs and legacy run IDs, reject ambiguous aliases and multiple columns
for one known sample, and enforce active/QC admission. Custom matrix admission and
metadata preparation now have their own module. Engine provenance also uses the
actual linked rsdeseq2 version instead of a hardcoded historical version.

New evidence covers closed-form fractional-count/unequal-length pooling, split-run
invariance in all three count modes, zero expression and malformed inputs; typed
group ordering, identity/fallback sources and contradictory model metadata; and
native fitting/correlation of eight supplied runs as seven biological samples.
The workflow checks all pooled gene counts, freezes exact contributors, reuses a
matching generation, rejects conflicting metadata without replacing prior output,
handles a separately supplied biological matrix, blocks duplicate run aliases and
accepts a proven matrix despite inactive history sharing its sample identity.
Metadata-edit tests confirm the edited condition reaches the actual model and
that the intervening stale-matrix refusal is a durable failed attempt.

Verified fifteenth implementation batch (2026-10-29): all 177 workspace tests pass;
the full mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. The refreshed official dependency audit reports 49 registry
packages across 58 direct requirements current, zero failures and one documented
unpublished Git dependency; the report matches both maintained lockfiles.

Metadata now has explicit factor/numeric/identifier/ignore roles, declared levels
and references, exact missing markers and a biological replication threshold.
Numeric-looking categories stay categorical; continuous values must be finite.
Protected identities/read paths cannot become predictors. Technical contributions
must agree after typed normalization. Missing or undeclared covariates, confounded
or saturated designs, insufficient biological replication, missing levels and
ambiguous/non-estimable contrasts fail before fitting or replacement-generation
creation. `~ 1` retains the sample dimension; `.` expands only declared predictors.

DE and correlation adjustment share one native model builder. The sealed
`model_design.json` records typed values, sample order, matrix, coefficients, offsets,
rank, residual degrees of freedom and resolved contrast weights. DE request revision
4 and correlation revision 3 invalidate older untyped requests. Contrast resolution
uses actual formula coding, including relevelled factors and models without an
intercept; full indicator models have no assumed treatment baseline. Categorical
Wald tests retain the two-group zero-count rule even when another group expresses
the gene. The metadata/formula parser, model/contrast preparation and settings
observations have focused modules rather than enlarging the DE runner.

Publication checks current scientific settings on disk and inside the catalog
transaction. Thread changes remain execution policy. Changed metadata types or
DE settings prevent stale default-correlation reuse. Failed preflight records a
failed attempt/project stage while preserving prior results and valid run counts.
All configuration sections reject unknown fields and unsupported project schema
versions; the shipped example is parsed by a regression test. That test also caught
and fixed its duplicate correlation table and misplaced DE transform setting.

Evidence includes shared typed DE/adjustment matrices, explicit reference ordering,
finite continuous values, dot/intercept/offset behavior, rank and replication
failures, catalog rollback under changed settings, thread-only reuse, native
relevelled/full-indicator effect equivalence, and a three-group zero-count fixture.
These establish internal contracts and invariants; independent scientific parity
remains a separate release gate. See [model design rules](docs/model-designs.md).

Verified sixteenth implementation batch (2026-10-29): all 189 workspace tests pass;
the complete mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. The official dependency audit was refreshed again: 49
registry packages / 58 direct requirements current, zero failures, the documented
unpublished Git exception, and exact hashes for both maintained lockfiles.

Metadata refresh now performs a run-keyed three-way merge of the previous generated
source, deliberate local edits and freshly generated values. Source-only changes
and local-only edits merge; differing concurrent cell changes fail before publication.
Explicit `keep-local`/`use-source` policies resolve conflicting cells while preserving
unrelated edits. `rnaa metadata refresh` updates from the stored catalog and
`rnaa metadata show` verifies/prints its retained publication. Resolver refreshes use
`--metadata-conflicts`; `--force` does not discard edits. Identities remain source-owned,
read paths refresh operationally, and selectors continue to use resolved source
metadata/registered overrides rather than silently consuming local model edits.

Migration 0013 stores immutable metadata revision descriptors and authoritative heads.
Each prepared revision retains generated, edited, merged and historical tables,
resolved run records, a column map, exact override source copies and a checked receipt.
Original source paths/hashes remain provenance. Cohort and metadata heads commit
together on resolution. Explicit refresh commits its head, event and successful
attempt together. Catalog checks run inside those transactions; large source hashes
remain outside. Candidates that fail validation or publication remain unselected.

Editable views install by atomic replacement and have a replayable installed marker.
Recovery handles stops before/after installation and merges later unrelated edits
against the exact pre-publication view. Conflicting later edits stay untouched.
Historical edit/source tables retain deliberately edited values across exclusion,
an empty active cohort and reinclusion. A missing working view can be restored.
Automatic download housekeeping defers refresh during live analysis/resolution;
explicit refresh respects the shared ownership boundary and records durable failure.

Override registration now uses unique, verified copies: same-named files registered
within one second cannot overwrite each other. DE/correlation retain a separate
canonical source copy inside their sealed generations, while mutable observations
still fence the running attempt. Later edits no longer mutate an earlier result's
source witness. DE request revision 5 and correlation revision 4 invalidate the prior
request form without claiming a statistical algorithm change.

New evidence covers native fitting after edited metadata survives refresh/download,
resolver conflicts preserving cohort/catalog/metadata heads, explicit resolution,
transaction rollback of metadata head/event/success, stale owners and candidates,
interrupted view installation with later edits, empty-cohort/reinclusion persistence,
same-name override registration, immutable original-source changes and corruption
rejection. See [metadata editing and recovery rules](docs/metadata-refresh.md).
These simulated transaction/installation boundaries do not establish power-loss or
complete process-kill release guarantees, and no automatic revision GC is claimed.

Verified seventeenth implementation batch (2026-10-29): all 196 workspace tests pass;
the mandatory formatting -> strict Clippy -> workspace test sequence and whitespace
checks pass. The refreshed official audit still reports 49 registry packages / 58
direct requirements current, zero failures and one documented unpublished Git
exception; its hashes match both maintained lockfiles.

Scientific metadata cache identity now projects the actual typed design and
run-to-biological-sample membership. Full metadata tables and unused schema
declarations stay in observations/provenance. Literal formulas, used typed values,
observed levels/references, actual X/offsets/coefficients, contrasts, samples,
count/reference/engine/transform inputs and QC admission remain scientific inputs.
Finite numeric spellings normalize; factor labels remain distinct. Dot expansion
still consumes declared predictors present in the sheet. Current schema/replication
preflight runs before reuse, even when a policy edit leaves the fitted X unchanged.

DE request revision 6 and correlation revision 5 separate those semantic values
from raw sheet hashes. Membership digest revision 2 ignores equivalent selector
expressions/source annotations after resolution. Full selectors, concrete cohort
head, catalog metadata, observed schema and raw source bytes still fence a live
attempt. Unrelated changes between attempts can reuse; changes during an attempt
conservatively require a retry. These request revisions change admission/cache
contracts without changing a statistical engine or the pooling method.

Default correlation reconstructs the current producer DE model/membership before
admitting its matrix. Its adjustment then uses current grouped metadata. With DE
`~ condition` and adjustment `~ age + condition`, changing age can reuse DE and refit
correlation using the new values; a producer-covariate change first requires DE.
An unchanged downstream matrix/request can reuse correlation even after a new DE
generation. This producer admission logic lives in a dedicated CLI module; typed
scientific design and source-observation handling also have focused modules.

Fresh sealed manifests/descriptors retain full configuration observations and raw
input witnesses outside request identity. Cache-hit attempts retain current raw
sheets under content-addressed `metadata/analysis_observations/` and attach those
witnesses/configuration observations to their atomic reconciliation event/outcome.
Original result generations and sidecars retain their original observation.
Existing observation files are verified without replacement; corruption stops reuse.
Hashing/copying remains outside the final SQLite writer transaction.

New evidence covers header/row ordering, read paths, unrelated catalog annotations,
unused numeric/ignored declarations and unused adjustment formulas preserving DE/
correlation heads; normalized numeric spellings; dot-formula sensitivity; effective
factor coding; replication preflight on equivalent projections; equivalent selectors
preserving membership identity while fencing old workers; current correlation-only
covariates changing fitted residuals with correct orthogonality; producer label
changes requiring DE; and retained raw witness immutability/corruption rejection.
The broader native/QC/reference/cohort/cancellation/crash tests remain green.
See [scientific cache contracts](docs/scientific-cache.md).

Scientific numerical release gates remain open. An exploratory exchange of RUN1
and RUN8 condition assignments in the synthetic native fixture triggered the
engine's dispersion nonconvergence path despite structural preflight passing.
This remains a numerical validation case for independent scientific fixtures;
preflight/reuse tests do not claim that every structurally valid model will fit.

Verified eighteenth implementation batch (2026-10-30): all 209 workspace tests pass;
the complete mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. The new daily official audit is
`validation/dependencies/2026-10-30.json`: 49 registry packages / 58 direct requirements
current, zero failures, the documented unpublished Git exception, and matching hashes
for both maintained lockfiles. The prior day's report is retained.

Migration 0014 indexes files owned by committed generations and stores exact
consumer/producer edges. New seals bind input roles, path/hash/size witnesses,
producer generation IDs and canonical producer-seal fingerprints. Files, edges,
head, registry and outcome publish atomically. Catalog changes before commit reject
publication; edge/file-index or producer-seal corruption rejects verification.
Managed inputs cannot silently downgrade to external sources, including aliases;
matrix paths are pinned canonically. Historical upstream heads remain readable.
Legacy file ownership is backfilled without rewriting old seals; missing recorded
ancestry is exposed as `recorded: false` and is not invented.

Fresh preprocessing, quantification, DE and correlation record their real producer
dependencies; quant also records its index. Recovery plans select compatible retained
trimmed reads, complete originals, full source restaging or an explicit blocked
outcome. Plans record operational target/source fingerprints and available input
observations. Actual count requests bind bytes after restaging. A pending recovery
can be reconstructed from run state, current settings and preserved committed
history after interruption. Restaging uses ordinary owned/resumable/atomic download
publication; failed replacements preserve prior counts, read producers and readers.
Disabling preprocessing recovers originals rather than treating trimmed data as raw.

`quant` restages previous stale-count dependencies while initial unrelated downloads
remain under download/run. Missing local or unsupported/incomplete remote sources
fail clearly. Blocked plans do not launch false replacement quant jobs. Competing
commands wait for one producer, and followers check outstanding previous-count
dependencies before claiming success. Explicit `run --no-download` reports required
recovery that conflicts with that option. Report admission remains independent of
count payload validity, including when QC evidence is corrupt and reads are pruned.

Count publication observes current layout, scientific options, selected source and
QC policy in the catalog and on disk, with execution allocations excluded. DE/corr
also observe upstream quantification policy and full remote-source metadata during
an attempt. Default correlation revalidates count dependencies after its reference
binding check, closing a stale-producer path when global count settings change.
Source selection is independent of progress state; regression coverage ensures
archive records cannot enter a FASTQ group when a run becomes complete.

`rnaa lineage` provides JSON inspection of selected generations, optional upstream
history and explicit payload rehashing. Failed payload verification retains the
report and returns failure; structural faults never print a falsely valid graph.
See [recovery/lineage contracts and limitations](docs/dependency-recovery.md).

Evidence includes native local-HTTP download/preprocess/quantification from all-pruned
and retained-trim states, durable pending recovery, preprocessing disablement,
failed restaging/retry, missing-source blocking, competing commands, exact and
transitive edges, legacy ownership migration, uncommitted/aliased producer refusal,
catalog rollback/corruption, lineage JSON failures and upstream setting fences.
Earlier native model/QC/reference/metadata/reuse/crash tests remain green.
Recovery-specific process-kill/power-loss/scale gates, generalized legacy lineage,
request reattachment, report-only repair after pruning and ownership/GC are still
open; successful synthetic/native scenarios do not establish those broader claims.

Verified nineteenth implementation batch (2026-10-30): all 218 workspace tests pass;
the full mandatory formatting -> strict Clippy -> workspace test sequence and
whitespace checks pass. The refreshed official dependency audit reports 49 registry
packages / 58 direct requirements current, zero failures and the same documented
unpublished Git exception; both lockfile hashes match the final tree. Migration
0015 indexes exact project/output-group/stage/subject/request lookup. A prepared,
verified reuse token binds the committed seal, complete original registrations and
observed selected head. Selection checks its owning attempt's scope/request,
candidate catalog/lineage and selected-head compare-and-swap, then restores the
registry and commits head/state/events/reuse-attempt success atomically. The producer
attempt, payload, seal and commit time remain unchanged. Fresh publication and reuse
share registry retirement, including legacy entries, without deleting files.

Normalization/DE/correlation now share current-request admission and historical
selection instead of checking only the selected result. Candidates require valid
payloads/lineage, descriptors/manifests, exact current scientific identity, required
outputs and current input/reference/QC/configuration checks. Corrupt historical
payloads are skipped; failed catalog/admission/ownership checks cannot claim reuse.
Current metadata/configuration observations are retained for each reuse attempt.
Read/reference/count historical selection and first-class typed request catalog
storage are explicitly still pending; this batch does not claim generalized DAG reuse.

`lineage --all-requests` inspects committed request variants for selected output
subjects, with a `selected` flag and optional ancestry/payload verification. The
tests exercise metadata/settings/covariate round trips, normalization versus
contrast contracts, exact restored registrations, unchanged old seals, pinned
intervening readers, corrupt-history recomputation/inspection, rollback after head
selection, changed catalogs/heads, cancellation, wrong scopes/requests and legacy
registry retirement. Existing metadata/preflight and native scientific workflows
remain required regression gates.

Verified twentieth implementation batch (2026-10-30): all 233 workspace tests pass,
with zero failures/ignored tests. The complete mandatory formatting -> strict
Clippy -> workspace test sequence and whitespace checks pass. Both history/recovery
concurrency scenarios also pass three consecutive focused runs after fixing stale
plan publication. The refreshed official dependency audit reports 49 registry
packages / 58 direct requirements current, zero failures and the documented
unpublished Git exception; both lockfile hashes match the verified tree. Migration
0016 stores complete typed scientific requests and policy indexes. Fresh reference,
download, preprocessing, native count and analysis publishers register requests in
the result publication transaction. Legacy native requests can be indexed only
after their full digest matches the committed head; old seals remain unchanged.

Preprocessing and native quantification now select verified historical results
under owned attempts before dependency reconstruction. Source snapshots bind the
current selected download, complete source identities and registry observations.
Local retained sources remain provable even when a trimmed-only count request no
longer records originals. Use-time bindings record exact input owners, request,
consumer and reuse attempt, and commit with head/registry/state/outcome. They do not
rewrite old manifests. Pruned FASTQs are not reintroduced as active registrations;
an engine execution worker still reconstructs missing trimmed reads for a new fit.

Bindings require exact scientific identity and full producer/pruning verification.
Trimmed-only request keys can retain current originals as provenance/retention
observations only through complete preprocessing ancestry. Dedicated binding
receipts authorize deletion of newer equivalent files without fabricating original
seal inputs. Integrity checks, current-source/processing comparisons after pruning,
QC admission and catalog/settings fences remain distinct. Fresh publication clears
the active binding while immutable historical receipts remain for retention proof.

Lineage inspection exposes original and use-time edges and separately verifies
bindings. Earlier healthy preprocessing reports can repair current QC evidence
without reads/network; no healthy evidence still requires reconstruction. This
does not close generalized report recovery or legacy ancestry. Request construction,
payload/scientific validation, physical input admission, source snapshots and
binding-retention evidence are separated into focused modules.
Recovery plans are rebuilt under run exclusion after historical lookup, so another
producer's completed read restoration cannot leave a stale download action to publish.
The local recovery fixture allows parallel native-work scheduling delays while
retaining exact transfer/producer-count assertions; dedicated HTTP tests still
verify finite stalls and transfer deadlines.

Regression coverage includes fully pruned history without HTTP, newer-path retention
and subsequent reuse, local/remote trimmed-only requests, preprocessing mode round
trips, corrupt-history failed recovery/retry, report repair, changed local/source
content and stale pruned heads, concurrent history commands, typed request atomicity,
lazy indexing without seal rewrites, exact ancestry and rejection of changed roles,
multiplicity/content, failed selection rollback and corrupt binding receipts.

Next implementation work after batch twenty: historical download/reference
selection (delivered in batch twenty-one below), generalized report reconstruction
and shared-store ownership/pins/GC. Cross-project sharing/pins,
immutable external inputs,
independent scientific oracle fixtures and measured resource/scale gates remain open.
Supported missing count inputs can now be restaged automatically. Unavailable
custom sources, archive extraction and generalized recovery remain explicit limits.

Verified twenty-first implementation batch (2026-10-30): all 248 workspace tests
pass, with zero failures/ignored tests, after the complete mandatory formatting ->
strict Clippy -> workspace test sequence. The official dependency audit again
reports 49 registry packages / 58 direct requirements current, zero failures and
the documented unpublished Git exception. Both lockfile hashes match the tree.

Download history selects earlier checksum-identified sources without origin I/O.
Older unchecksummed URL snapshots require fresh strong HTTP validation matching
the stored validator, effective URL, size and identity encoding; weak/missing or
changed validators, unsupported HEAD and failed validation fall back to a new GET.
A failed fetch preserves the current generation. The currently selected verified
snapshot retains its existing offline cache contract. Private unchecksummed file
caches left by failed multi-file attempts now require current origin validation
before reuse. Reuse proof, registry/head/state/event and successful attempt commit
together; producer seals and old files are unchanged. Pruned FASTQs cannot be
reattached as physical inputs.

Reference history validates full source content, annotation policy, scientific
index settings, engine identity and complete sealed inventories. Restoring index
options selects the original healthy generation. Restoring another previously
prepared source uses an owned atomic project-reference selection even when that
source's own resource head remains unchanged. Identical custom sources moved to
new paths reuse payloads while recording current observations. Unchanged access
creates no extra attempt. Source/configuration/parent-choice guards fence fresh
and historical publication; execution thread changes preserve reuse. Provider
reference resolution still consults provider listing/checksum metadata.

Coverage includes source/index round trips, corrupt-history replacement with
preserved evidence, offline publisher-checksum reuse, live strong-validator reuse,
changed/unversioned/weak/unsupported origins, failed validation/fetch followed by
verified retry, failed multi-file private-cache revalidation, source/configuration
change guards and cancellation. A real native download -> preprocessing -> counts
source round trip selects all three original generations without network access
when publisher checksums establish identity. No statistical engine/request revision
changes in this batch. An initial local HTTP integration failure was not reproduced
by the subsequent focused suite, complete mandatory sequence or three additional
consecutive focused-suite runs; failure output now includes the durable attempt
diagnostics and transfer log. The initial failure's cause remains unconfirmed.

Next implementation work after batch twenty-one: QC report reconstruction without
healthy history (delivered in batch twenty-two below), followed by shared-store
ownership/pins/GC. Cross-project reuse, complete
legacy lineage, per-run/cohort progress separation, independent scientific parity,
bounded correlation/resource allocation and product/release gates remain open.

Verified twenty-second implementation batch (2026-10-30): all 255 workspace tests
pass with zero failures/ignored tests after the complete mandatory formatting ->
strict Clippy -> workspace test sequence. Whitespace checks pass. The refreshed
official dependency audit reports 49 registry packages / 58 direct requirements
current, zero failures and the documented unpublished Git exception. Both lockfile
hashes match the verified tree.

QC evidence recovery now remains actionable when counts themselves reconcile.
When no healthy compatible historical report survives, `quant` plans report
reconstruction from complete original reads, downloading the supported selected
source when required. Trimmed reads do not stand in for originals. A blocked
report-only plan retains completed count progress and its missing-dependency
diagnostic. Operational recovery request revision 2 distinguishes `counts` from
`qc_report`; no statistical engine/scientific request revision changes.

Reconstruction publishes a fresh immutable preprocessing generation and reapplies
current QC admission. Matching recovered source/processing content selects the
original count payload without refitting, retaining its seal bytes and recording
new input bindings. Changed source content requires new counts. A failed fetch
preserves all current heads; rejected measured reads remain ineligible even with
engine-error bypass enabled. Existing damaged reports remain unchanged for
diagnosis. Fresh/cache preprocessing publication now checks complete current
original groups, frozen source/catalog/settings observations and the size policy,
using the same guards as historical selection. Thread allocation remains execution
policy. Shared source observations and reusable native test fixtures were extracted
into focused modules; the dependency recovery integration file shrank from over
1,000 lines to 821.

Seven new native integration tests cover missing/corrupt only reports after complete
pruning, retained local originals, retained trimmed-only refusal, failed transfer
and retry, changed unchecksummed origins, rejected reconstructed QC and competing
commands sharing one replacement download/report. They assert unchanged count
generations when scientific inputs match, preserved producer seals, renewed
admission and retention-compatible repeat access. Existing recovery, ownership,
retention and scientific workflow suites remain green.

Next implementation work after batch twenty-two: a store-wide ownership catalog,
durable publication/reader pins and shared orphan GC/trash/restore (delivered in
batch twenty-three below). Collection of project generations must also retain historical generation edges,
use-time bindings and pruning receipts, and protect in-flight producers/readers.
Cross-project request reuse, complete legacy ancestry, immutable external inputs,
per-run/cohort progress separation, independent scientific parity, bounded
correlation/resource allocation and the remaining product/release gates stay open.
The new recovery tests do not establish process-kill or power-loss guarantees for
every recovery boundary.

Verified twenty-third implementation batch (2026-10-30): all 269 workspace tests
pass with zero failures/ignored tests after the complete mandatory formatting ->
strict Clippy -> workspace test sequence; whitespace checks pass. The refreshed
official dependency audit reports 49 registry packages / 58 direct requirements
current, zero failures and the documented unpublished Git exception. Both lockfile
hashes match the verified tree. Tests use the already selected workspace rusqlite
release through an explicit CLI test dependency; the new lockfile edge changes
the conservative engine dependency fingerprint.

Shared HDF5 blob publication now uses a store-wide SQLite catalog and pinned
handoff. Durable PREPARED project/generation ownership is recorded before copying;
the exclusive blob lock spans verified copying, project publication and global
ownership finalization. The scientific project transaction remains unchanged in
scope. A crash between databases leaves a reconcilable reservation, never an
unprotected blob. Committed generation ownership survives active-head replacement.
Missing/unreadable/substituted project catalogs and changed committed ownership
proofs remain protected, with diagnostics. Normal managed owners use their exact
generation witness; legacy active references remain supported. Concurrent store
initialization/directory creation and publication are coordinated.

The library reader API combines durable pin rows with shared OS exclusion, so a
paused live reader/writer cannot lose its blob because a timer expired. Applied
collection can discard stale reader rows only after exclusive acquisition. New
shared copies use independent inodes. Existing unknown blobs receive legacy holds
when encountered, and unindexed files are never eligible for automatic deletion.

`storage status`, default-preview `storage gc`, `storage gc --apply`,
`storage restore` and explicit expired-trash purge print JSON. Collection reconciles
project ownership and records a transactional move/delete intent before filesystem
changes; interrupted trash/restore/purge operations resume under the same lock.
New trash deadlines are persisted and cannot be shortened by later GC options.
Corrupt/ambiguous payloads are retained for diagnosis. Catalog bytes, owner counts,
prepared owners, durable readers, held reasons and operation states remain explicit;
unindexed files and abandoned temporary copies are not counted as reclaimed.
Configured relative shared roots now follow the project root, independently of
the invoking directory; explicit store inspection can work without a project DB.

Ten core tests and four native CLI tests cover two-project ownership, eight
simultaneous publishers, live writer/reader exclusion, real killed writers on both
sides of project commit, killed readers, stale pin reconciliation, interrupted
move/restore intent, dry-run immutability, retained deadlines and purge, unavailable
projects, legacy holds, payload/symlink faults, failed native registration/retry,
historical ownership, concurrent native projects, corrupt ownership witnesses,
relative roots and JSON preview/trash/restore. Existing scientific/recovery/retention
and concurrency suites stay green. Focused modules separate catalog publication,
filesystem exclusion, project ownership proof and collection transitions.

Next implementation work after batch twenty-three: dependency-aware project-generation
retirement, project reader/publication pins and reversible generation cleanup
(delivered in batch twenty-four below). FASTQ/SRA-specific `TrashCleaner`, complete
legacy import, abandoned temporary cleanup and consistent global/project physical
accounting remain open. Cross-project scientific request reuse, resource budgets,
independent scientific parity, per-run/cohort progress separation, unfinished CLI
commands and measured release/scale/power-loss gates also remain open.

Verified twenty-fourth implementation batch (2026-10-30): all 280 workspace tests
pass with zero failures/ignored tests after the complete mandatory formatting ->
strict Clippy -> workspace test sequence. The earlier full maintenance regression
run and focused lifecycle/native suites also pass. Whitespace checks pass. The
refreshed official dependency audit reports 49 registry packages / 58 direct
requirements current, zero failures and the documented unpublished Git exception;
both lockfile hashes match the verified tree. Migration 0017 records project
maintenance observations, keep labels and generation retirement receipts.

Workers, reconciliation scopes, prepared reuse tokens, shared-blob publication and
generation path snapshots now hold shared project maintenance exclusion. Database
mutations acquire that boundary before their SQLite writer transaction. Cleanup
requires exclusive acquisition across graph validation and filesystem changes;
paused live readers/writers cannot be overridden by a timeout. Path snapshots and
clones retain a guard; head/generation path conversion now returns `Result` and
rejects retired data. Durable observations can be stale or absent under a containing
writer transaction; OS exclusion is authoritative. Lineage inspection pins its
complete read. Existing public source/path callers were updated to propagate errors.

Project cleanup roots include all selected heads, selected reference, active artifact
views/aliases, current QC, keep labels and retained-input quant witnesses. Traversal
includes original producer edges, every historical input binding and reverse
pruning-proof links to the exact count/binding consumer. Damaged unreachable results
keep their required dependencies; invalid catalog indexes stop collection. Unknown
legacy ancestry conservatively holds the project history rather than inventing
edges. Canonical ownership resolves active aliases and missing/dangling paths;
retired artifacts, QC evidence and references cannot be reactivated through aliases.

`storage generations status/gc/keep/unkeep/restore` provide project JSON operations
without a configured shared store. GC defaults to preview. Applied retirement records
inventory/registry/seal/observed-file receipts before moving complete generation
directories to project trash. Original metadata/seal bytes survive. Recorded deadlines
cannot be shortened by later options. Expired producers remain protected by other
live/restorable consumers. Restore admits the complete required dependency/proof
closure only after every directory is back and verified; it does not change heads.
Retired requests are excluded from reuse lookup until restoration finishes.

Durable intents resume interrupted moves, grouped restoration and partial recursive
purge without accepting unexpected remaining files. Purged metadata stays as an
inspectable tombstone. The shared collector can release matching generation
ownership only after a completed purge receipt and active-reference checks; project
trash retains its shared copy. Measured candidate payload sizes are distinguished
from unmeasured data. Collection is currently project-wide idle maintenance.

Eight core lifecycle tests and three native CLI tests cover transitive producers,
snapshot clones, prepared tokens, workers, explicit keep/unkeep, aliases, corruption,
legacy holds, deadlines, restore closure and unchanged original seals. Real killed
readers/collectors exercise intent-before-move, intent-after-move, interrupted restore
and partial purge. Native tests exercise restored-request selection, shared-owner
release after project purge, and binding/pruning-proof survival after another fit is
retired. Existing scientific, recovery, retention and shared-store suites remain green.
Graph planning, filesystem transitions, maintenance pins and ownership normalization
are separated into focused modules.

At the end of batch 24, next work was FASTQ/SRA-specific reversible file cleanup
and `TrashCleaner` (delivered in batch 25 below), followed by legacy ancestry/import
and abandoned data cleanup/accounting. File cleanup flags remain separate from
whole-generation GC.
Scientific parity, memory/resource admission, cross-project scientific request reuse,
per-run/cohort progress separation, prompt engine cancellation, unfinished product
commands and measured release/scale/power-loss gates remain open. This batch does
not establish full legacy or power-loss guarantees.

Metadata refresh/merge, immutable source witnesses and typed scientific metadata
cache projection are delivered; comprehensive process-kill/power-loss release gates
remain open. Run-level progress still
contains cohort completion states. Independent scientific oracle fixtures must
cover both the base importer and the new pooling policy before claiming parity.

Legacy ancestry/import, abandoned data cleanup and measured physical accounting
remain necessary before the entire storage roadmap closes. Cross-project request-to-artifact lookup, complete legacy DAG lineage,
per-run/cohort state separation, pinned independent scientific oracle fixtures,
bounded correlation/resource budgets and remaining product/release gates remain open.
Synthetic workflow and process-kill tests do not establish independent scientific
parity, power-loss guarantees or complete legacy generation garbage collection.

All commits and pushes require the full AGENTS.md verification sequence. Checked
items must point to executable tests or measured runtime evidence, not intent.

### Batch 25 — verified reversible FASTQ/SRA input cleanup (2026-10-30)

Migration 0018 adds immutable count-backed input trash descriptors and separate
physical witnesses. Manual `storage files status/clean/gc/restore` commands print
JSON and preview by default. Cleanup verifies exclusive project ownership and exact
native producer/consumer/input-binding proofs, holds complete selected run groups
when any member fails, records moves before removing registry entries, and restores
complete groups only after their physical bytes verify. Changed layout/source/
producer selections restore data without selecting stale inputs. External files,
hardlink aliases, extra files and conflicting bytes are preserved with diagnostics.

`TrashCleaner` now implements automatic `none|fastq|sra|all` cleanup and
`success|always|never` scheduling, including idle completed-count boundaries and
command outcomes. Reversible FASTQ trash takes precedence over permanent FASTQ
retention when enabled. Cleanup failure preserves count completion and the primary
processing error. Purge stays explicit, uses original deadlines, waits for the
complete operation to expire, and never expires newly trashed files in that call.

Trashed/purged receipts supply exact count-input/pruning evidence; restore retires
that operation's absence permission without editing scientific history. Project
GC protects restorable producers/consumers and retained absence-proof witnesses.
Recovery selects complete current trash groups before downloading, including local
inputs with no downloadable source. Ingest startup resumes pending transitions.
No SRA extraction or archive quantification is claimed.

Filesystem moves use exclusive destination publication, sync both locations and
recover both-name/copy-before-unlink states. Cross-device copy staging has a
reserved operation identity. Real killed processes exercise partial paired moves,
complete/partial copies, partial restore and unlink-before-purge-commit. Actual
multi-device and power-loss validation remains a release gate. Planning/catalog
transitions, ownership proofs, filesystem operations and restoration are separate
modules, keeping the new workflow reviewable.

Native tests also exercise historical count input bindings, offline restoration for
changed count settings, scheduling/retention precedence, source changes, collisions
and corruption. Broad recovery checks caught an unnecessary second reconciliation
after ordinary downloads; the extra pass now applies only to restored trash.
Recovery also refreshes QC admission under run exclusion before planning, rather
than assuming the earlier command-level observation is still current. Retry of a
completed/interrupted restore is idempotent. Selected FASTQ families must match the
complete input group consumed by the current count binding, including missing
registry members whose files still exist.

A downloader empty-queue/storage-threshold stall is fixed: only queued work waits
for capacity, so the final completed transfer can return. This does not establish
resource reservations or bounded storage admission. Idle `download --forever`
command shells no longer retain project-wide reader pins; individual download
attempts still pin their payloads. A native child-process test proves another
command can clean completed counts while the idle downloader stays alive.

Earlier full-suite attempts exposed intermittent HTTP fixture failures and apparent
extra recovery fetches. A local socket experiment confirmed that macOS accepted
sockets inherited the nonblocking listener flag: reads returned `WouldBlock` in
microseconds instead of waiting for request bytes. The shared recovery fixture
counted those accepts as requests, while the download-history fixture closed them
before a request arrived. Both now explicitly use blocking accepted streams;
recovery counts only complete HTTP requests. Delayed connection/header regressions
exercise both failure modes, and failure diagnostics retain request lines, events
and attempts. This fixes a confirmed fixture defect without relaxing production
transfer or concurrency assertions. Full measured concurrency/power-loss release
gates remain open.

The full mandatory sequence completed successfully on the final implementation:
`cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **298 passed, 0 failed, 0 ignored** across
26 test/doc-test suites. The 18 additions include seven core lifecycle/helper
checks, ten file-cleanup/HTTP fixture integration checks and one late-request
history-origin regression. The final execution log is
`/tmp/rnaa-file-cleanup-final-all-checks-sixth.log` (terminal exit 0); earlier
partial/failing development logs are not the verification result.

The fresh official-registry dependency audit checked 49 packages / 58 direct
requirements with zero failures and the documented pinned Git exception. The
report at `validation/dependencies/2026-10-30.json` matches both maintained lockfile
hashes. No dependency versions changed in this batch. `git diff --check` also
passed. See `docs/input-cleanup.md` for the complete operator contract.

Next work: complete legacy ancestry/shared ownership import, abandoned unpublished
cleanup and physical accounting; independent scientific parity, bounded
correlation/global resource budgets, cross-project scientific request reuse,
per-run/cohort progress separation, prompt engine cancellation, unfinished product
commands and measured release/scale/power-loss gates remain open. The overall goal
remains active; this batch does not close the storage or release roadmap.

### Batch 26 — sealed legacy ancestry and observed storage accounting (2026-10-30)

Migration 0019 adds additive lineage-import receipts and indexed producer edges.
`storage generations import-lineage` previews by default and applies only verified
ancestry from original seals, every registered manifest and supported physical
request witnesses. Catalog and filesystem witnesses bind the evidence without
rewriting original seal bytes, requests, payloads, active artifacts or selected heads.
Exact consumed inputs remain separate from broader provenance. Reference source
origins and quantification index witnesses protect their original producers without
expanding the read set used for count reuse. Missing, conflicting, cyclic or
unsupported evidence remains held. Imported ancestry participates in GC, restore,
input retention and offline recovery; unresolved history still conservatively
protects the project. Current source aliases add maintenance-only protection of
their present managed targets without altering immutable historical edges.

`storage usage [--details] [--include-shared]` reports observed logical bytes by
category and, on Unix, allocated blocks deduplicated by device/inode across included
roots. It observes symlinks without following them, skips repeated directory inode
aliases, reports hardlinks beyond the scan, and qualifies concurrent catalog/worker
activity and traversal errors. It creates no shared store and authorizes no cleanup.
Resource admission reservations and actual free-space/reclamation estimates remain
separate work. The implementations use focused legacy evidence/catalog/test modules
and separate accounting traversal/catalog modules instead of enlarging the existing
generation and storage routing files.

The final mandatory sequence passed on the updated lockfiles:
`cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **313 passed, 0 failed, 0 ignored** across
27 test/doc-test suites. The terminal exit was 0; the authoritative execution log is
`/tmp/rnaa-legacy-accounting-final-all-checks.log`. New coverage includes real killed
importers before/after catalog publication, rollback, corruption, pins, cycles,
reference origins, consumed/provenance separation, aliases, hardlink accounting and
native pipeline/analysis imports, retained count reuse and recovery with origins
unavailable. The native vendor library has **82 passing tests**, and its actual
`wasm32-unknown-unknown` library check passes with the optional WASM feature.

The fresh official-registry audit found and upgraded `js-sys` to **0.3.107** and
`wasm-bindgen` to **0.2.130**, synchronizing both vendor manifests, both lockfiles
and the related WASM dependency family. It then checked 49 packages / 58 direct
requirements with **zero failures**, retaining the documented pinned Git exception.
`validation/dependencies/2026-10-30.json` matches both current lockfile SHA-256
values. `git diff --check` also passes. Operator contracts are documented in
`docs/legacy-ancestry.md` and `docs/storage-usage.md`.

Remaining storage work includes complete legacy shared ownership migration,
unsealed result adoption, abandoned/unpublished cleanup and reservation-backed
admission. Independent scientific parity, bounded correlation/global resources,
cross-project scientific reuse, per-run/cohort progress separation, prompt engine
cancellation, remaining product commands and measured release/scale/power-loss
gates remain open. The overall goal remains active.

### Batch 27 — reversible unpublished generation cleanup and CI dependencies (2026-10-30)

Migration 0020 adds a separate cleanup catalog for tracked unpublished generation
containers. `storage generations unpublished gc` previews eligibility and logical
bytes; apply records immutable whole-container witnesses, marks terminal attempts'
generation rows ABANDONED and moves their data into reversible diagnostic trash.
`status`, `restore` and explicit expired purge complete the operator workflow.
Restore preserves bytes and empty directories without publishing a scientific
result, changing selected heads or making failed outputs reusable. Original
scientific inventory/request/registry metadata remains unchanged.

Live workers and readers exclude maintenance. Owning attempts must be terminal;
clock expiry alone does not authorize collection. Published catalog references,
committed ancestry and bindings, registered and retained input paths, and both
persisted/on-disk settings protect their containers, including aliases and current
trash paths. Unresolved legacy ancestry, unsealed legacy artifacts and unindexed
directories remain held. Invalid candidates, changed bytes/receipts, special files,
symlinks and ambiguous destinations never authorize deletion.

Trash deadlines are fixed at selection. A zero-day window still cannot purge newly
selected data during the same invocation. Durable move/restore/purge intents resume
after actual process termination; partial purge permits only missing inventory
entries while verifying every survivor. Current references are checked before each
maintenance pass. Restores complete before rebuilding protection for other pending
moves. The implementation keeps catalog, protection, traversal, recovery and tests
in separate focused modules. It indexes file ownership and partial-purge inventories
instead of quadratic membership scans, reads full catalog documents only for
unpublished rows and avoids rescanning the complete graph per candidate.

New core and native tests cover lifecycle/expiry, hardlink preservation, empty
directories, immutable seal bytes, rollback, new references after trash, aliases,
absolute paths containing URL-like filename components, legacy holds, live pins,
unindexed data, changed payload/catalog evidence, native count reuse and real killed
builders/collectors before/after trash and restore moves and during partial purge.
The final mandatory sequence passed on the latest implementation:
`cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **322 passed, 0 failed, 0 ignored** across
28 test/doc-test suites, terminal exit 0. The authoritative log is
`/tmp/rnaa-unpublished-final-all-checks-third.log`; earlier development runs do not
replace that final verification. `git diff --check` also passes. Nine new tests
extend the 313-test baseline from batch 26, including the real crash-worker helper.

The independent CI audit found `actions/checkout@v4` behind the current stable
release. The workflow now selects **checkout v7.0.1** and **rust-cache v2.9.3**, with
their release and action-manifest evidence in
`validation/dependencies/ci-actions-2026-10-30.json`. Both actions use Node 24;
the workflow uses hosted Ubuntu. The workflow and Dependabot YAML parse, selected
versions match the official latest non-prerelease releases and the recorded workflow
hash matches the edited file. Weekly Dependabot updates now include GitHub Actions.
The hosted workflow itself has not been executed locally. The renewed crate audit
checks 49 packages / 58 requirements with zero failures and matching lockfile
hashes; no crate or lockfile versions changed after batch 26.

See `docs/unpublished-cleanup.md` and `validation/dependencies/README.md`. Remaining
storage work includes complete legacy shared ownership migration, unsealed legacy
adoption, unindexed/transfer/support ownership and cleanup, admission reservations
and measured recovery/scale/power-loss release gates. Independent scientific parity,
bounded correlation/resources, cross-project scientific reuse, per-run/cohort state
separation, prompt engine cancellation and remaining product commands stay open.
The overall goal remains active.

### Batch 28 — verified legacy shared ownership import (2026-10-30)

Shared schema 2 adds explicit legacy project inventory, immutable import descriptors,
per-blob proof indexes/markers and a durable census of previously imported project
identities. `storage import-legacy --project ROOT ...` previews by default. Apply
imports verified ownership without changing scientific result files, manifests,
requests, selected heads or project schemas. Project catalogs are opened read-only,
including original shared-reference catalogs without generation tables. Library
and CLI output exposes migration receipt IDs for tracing the saved evidence.

Only an explicit `--complete-project-inventory` assertion, successful checks of
all known project identities, verified blob/view bytes and fully accounted Unix
hardlinks can release a legacy hold. Missing/invalid projects and unexplained links
retain protection. Canonically named unindexed blobs can be cataloged after actual
SHA-256/length verification; unexpected names, symlinks and temporary files stay
unindexed. A subset inventory cannot overwrite a previously released proof or
release ownerless blobs while another known project is omitted. Unsupported project
migration histories remain protected during both import and later ownership proof.
The shared reader explicitly supports migrations through 0020; future reference
contract changes require updating that boundary after compatibility inspection.

Project maintenance locks and deterministic blob exclusion cover evidence capture,
revalidation and store publication. Filesystem witnesses precede one atomic shared
catalog transaction. Failed transactions and killed pre-commit importers leave no
published ownership changes; retries never adopt orphan witness files by filename.
Known-project coverage is rechecked inside the store writer before releasing holds.
Concurrent store initialization now reads the schema version after acquiring that
writer, fixing a race caught by the existing simultaneous-publication regression.
Project creation/opening also acquires maintenance exclusion before connection
initialization, schema migration or metadata backfill. A real process test holds
exclusive legacy maintenance, proves an attempted normal open leaves the old
catalog unchanged, then releases the lock and verifies migration can finish.

GC verifies imported descriptor/catalog/filesystem evidence and active historical
claims, so a missing owner/index cannot authorize deletion. Current additional
hardlinks protect orphaned live and trashed blobs even after import. Pending trash/
purge intents recheck ownership and links before finishing. Completed purge clears
the live proof link while keeping immutable import history and the project census;
fresh later publication has its own ownership. Scientific adoption of unsealed
results remains separate from importing a shared ownership claim.

New coverage includes read-only original-schema projects, multiple-project claims,
explicit inventory/hold release, omitted known projects, preserved prior proofs,
unknown/new hardlinks, busy readers/writers, rollback, missing owner/index/witness
faults, unsupported future schemas and actual killed importers before/after commit.
Native CLI tests preserve original quantification seals and current count reuse,
import historical generation owners and release them only after project purge.
Explicit store import also bypasses unrelated old settings syntax, preserving the
invoking legacy project's configuration and catalog bytes.

The final mandatory sequence passed on this implementation:
`cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **336 passed, 0 failed, 0 ignored** across
29 test/doc-test suites, terminal exit 0. The authoritative final execution log is
`/tmp/rnaa-legacy-shared-final-all-checks-third.log`; earlier development and
pre-final inspection-output runs are not the final verification. Fourteen tests
extend the 322-test baseline, including two native CLI tests and crash/open helpers.
`git diff --check` also passes.

The fresh official-registry dependency audit checked 49 packages / 58 direct
requirements with zero failures and the documented pinned Git exception. Both
maintained lockfiles and the CI workflow still match their recorded evidence; this
batch adds no dependency or lockfile changes. `docs/shared-legacy-import.md` records
the full operator contract and platform/coverage limits.

The overall goal remains active. Unresolved legacy evidence and unavailable project
identities stay held; unsealed result adoption, unindexed/transfer/support cleanup,
tracked admission/global resources, independent scientific parity, bounded
correlation, cross-project scientific reuse, per-run/cohort state separation,
prompt engine cancellation, remaining product commands and measured release/scale/
power-loss gates still require their original acceptance evidence.

Verified implementation batch twenty-nine (2026-10-30): independent scientific
references and bounded correlation execution.

The mandatory offline scientific suite now reads pinned independent R output,
checks file/generator/input hashes and compares numerical values, missingness,
rounded counts, identifiers and edge ordering. Owned R generators cover tximport
1.40.0, DESeq2 1.52.0, stats, WGCNA 1.74 and HellCor 1.3. Source artifact hashes,
session versions, regeneration instructions and quantity-specific tolerances are
recorded under `validation/scientific`. The small upstream R fixture is preserved
with its exact revision and MIT attribution. Production/tests require no R runtime.

Import references cover all count modes, offsets, fractional rounding, mapping and
PAR_Y identity, all-zero/zero-abundance rows, zero effective length without mass
and real two-lane pooling. Only zero-length rows with no count/TPM mass are admitted;
positive-mass invalid lengths still fail. Import algorithm revision is 4.
Workflow references cover categorical/numeric designs, numeric/factor interactions,
forward/reverse contrasts, three levels, all-zero and outlier rows, normalized
counts, size factors and VST. The caller now selects the documented local-trend
fallback only for classified parametric-fit failures and records the decision.
Explicit contrast results now receive Cook's filtering and independent filtering
before BH-adjusted values are written; manifest policies capture the decisions.
DE request revision is 7. The rsdeseq2 dependency remains 0.3.1, without a library
fork or upgrade. The user explicitly requested that the new DE library be integrated
at the end; keep the present version during intermediate work and rerun these
reference gates when that final replacement is introduced. Broader engine fixes
belong to the separate library effort.

Correlation separates row preparation, bounded tile computation, heap/stripe output
and payload estimates. Pearson/Spearman/bicor preparation is reused; Kendall and
HellCor use pinned cross-matrix APIs and repeat private preparation per tile.
Top-k retains only B*k candidates plus the input/tile, with exact input-target tie
ordering. Dense and threshold outputs hold one B*G stripe instead of G*G values.
All output modes share constant-row/endpoint policies, and reference edge tests
exercise inclusive thresholds 0, 0.75 and 1 across partial tiles and two threads.
Original R coefficients remain separately available; eight-epsilon endpoint
normalization is explicit in both independent generator and native kernel.
Correlation request revision is 7. Resource allocation stays outside scientific
identity; byte-for-byte allocation tests and CLI historical reuse verify this.

`corr.block_rows`, `max_memory_mb` and `max_output_mb` have validated defaults.
A streaming input scan checks conservative working/output payload estimates before
numeric loading or correlation output. Rejected requests preserve the current
published generation and original files. This is not a global reservation, free-disk
guarantee, OS memory cap or module-enrichment resource budget.

The reproducible process probe records development-profile 2,000/4,000-gene,
16-sample Pearson top-10 runs at 1.498/5.658 seconds and peak RSS
9,125,888/9,846,784 bytes, including process startup/input generation. Both retain
at most 640 edges and 4,096 tile entries, with no top-k dense stripe. The recorded
source/lock hashes match the implementation. These are initial synthetic memory
measurements, not release throughput, production cohort or supported-host evidence.
See `validation/correlation-scale` and `docs/correlation-resources.md`.

The second full verification attempt exposed an intermittent existing concurrent
shared-store initialization failure (`database is locked`). Catalog open now takes
a bounded OS initialization lock before the SQLite connection/WAL selection and
migration transaction, releasing it before blob locks. A deterministic exclusion
test checks that waiting openers cannot touch SQLite before lock release; the
existing eight-opener owner-preservation regression is also repeated separately.
This fixes the journal-setup race rather than accepting an intermittent green rerun.

Final verification passes in the mandatory order: `cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **354 passed, 0 failed, 0 ignored** across
32 test/doc-test suites, terminal exit 0. The authoritative final execution log is
`/tmp/rnaa-scientific-blocks-final-all-checks-third.log`; the first development pass
and second failed attempt are not the final verification. Eighteen tests extend
the 336-test baseline. The shared-store simultaneous initialization/publication
test also passes 20 repeated runs after its fix; all 11 focused store tests pass.
`git diff --check` passes. A durable copy and summary are under `validation/checks`.
The fresh official-registry audit reports 49 packages / 58 direct requirements,
zero failures and one documented Git exception; both lockfiles remain unchanged.

The overall goal remains active. Next correlation work is standardized PCA module
scores with deterministic sign, strict GMT/parameter admission, bounded module
graph/connectivity computation and shared resource limits. Legacy adoption and
unindexed/transfer/support cleanup, cross-project scientific reuse, global resource
reservations, per-run/cohort state separation, engine cancellation, product commands
and full measured release/migration/power-loss gates retain their original scope.

Verified implementation batch thirty (2026-10-30): full-component module PCA,
strict enrichment admission and bounded module execution.

The former monolithic enrichment implementation is replaced by focused graph,
PCA, ranking, admission, output and workflow modules. Mutual positive Spearman
top-20 neighbours are selected in bounded tiles; full connected components and
their positive connectivity sums are retained separately from report truncation.
No module phase allocates a gene-by-gene matrix. Report ordering uses an explicit
1e-10 connectivity grid for roundoff ties, preserving original values in tables.

Module scores now use standardized first-principal-component values from every
member, with unit sample variance, deterministic sign and a canonical leading
eigenspace policy for tied eigenvalues. The old arithmetic mean/top-report-only
score is removed. New canonical `module_membership.tsv` and `module_eigengenes.tsv`
outputs are sealed, registered, required by reuse and retained with their generation.
Module sizes refer to full components even when `--module-top 1` is requested.
Enabling modules no longer overrides the requested ordinary network method:
module graph construction stays Spearman while main-network selection remains explicit.

Strict GMT validation rejects duplicate/empty identities and entries, computes
overlap after excluding constant genes, records filtered/unmatched evidence and
rejects empty eligible pathway populations before graph fitting. Options no longer
silently increase invalid minima or permutation counts. New `--fgsea-min-size` and
`--fgsea-max-size` bind exact pathway bounds to the scientific request. Undefined
constant-row rankings are omitted instead of fabricated zero scores; the existing
finite perfect-correlation t-ranking cap is explicit. Module algorithm revision 2
invalidates old module requests independently of ordinary correlation requests.

The resolved correlation thread budget and block/memory/output limits apply to
module work. An owned local Rayon pool includes rsfgsea's parallel execution.
Payload estimates cover matrix copies, neighbour storage, sample-space PCA,
ranking/permutation scratch and result tables. CLI admission runs before ordinary
correlation work; execution revalidates and snapshots both matrix and GMT before
reading, then checks them again before finalizing. Failed admission preserves the
previous published generation; allocation-only changes reuse the prior producer.

Pinned independent R 4.6.2 / WGCNA 1.74 / fgsea 1.38.0 references cover unequal
scales (explicitly distinguishing PCA from the old mean), zero-average sign,
tied leading components, three samples, full graph/connectivity, report limits,
module scores, gene rankings and seeded simple-permutation enrichment. File and
generator/source hashes are mandatory offline checks. Nondegenerate PCA is checked
against public R `prcomp` and normalized WGCNA outputs; the declared tied-space
policy is independently implemented with R eigenspace projection. Continuous
comparisons retain 1e-9 absolute/relative tolerance.

New failure tests cover ambiguous/no-overlap GMT, invalid report/module/pathway/
permutation options, oversized graph/PCA/permutation/output estimates, nonfinite/
constant/malformed PCA input, bounded neighbour/tile counts and native CLI rollback
and reuse. The R-reference workflow repeats with different report limits, partial
tiles and thread allocations while preserving full membership/scores.

The reusable process probe records complete 2,000/4,000-gene, 16-sample module
workflows in the development profile: 2.827/11.136 seconds, peak RSS
15,155,200/18,137,088 bytes and total outputs 212,547/420,091 bytes. Retained
neighbours grow from 40,000 to 80,000, with at most 4,096 tile entries. Both include
graph, PCA, ranking and 100-permutation enrichment. Source/lock hashes and raw
measurements are saved in `validation/correlation-scale/2026-10-30-modules-debug.json`.
These are synthetic local measurements, not a production/supported-host release
qualification, OS memory cap, shared reservation or prompt-cancellation guarantee.

Final verification passes in the mandatory order: `cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **363 passed, 0 failed, 0 ignored** across
34 test/doc-test suites, terminal exit 0. The authoritative final log is
`/tmp/rnaa-modules-final-all-checks-third.log`; the first two development passes
precede the final reference-tolerance and explicit network-method checks. Nine
tests extend the 354-test baseline: two independent module-reference tests, five
resource/admission/PCA tests and two native CLI tests. `git diff --check` passes.
The durable log and summary are under `validation/checks/2026-10-30-modules.*`.
The fresh registry audit still checks 49 packages / 58 direct requirements with
zero failures and the documented Git exception. Both lockfiles and the production
DESeq2 library remain unchanged; its replacement stays deferred to the end.

The overall goal remains active. Per-run/cohort progress separation, global
CPU/RAM/disk reservations and index sharing, cross-project scientific reuse,
unsealed legacy adoption/unindexed transfer-support cleanup, prompt engine
cancellation, unfinished CLI/product contracts and measured production/release/
migration/power-loss gates retain their original acceptance requirements.

Verified implementation batch thirty-one (2026-10-30): separate run ingest progress
from cohort analysis attempts and outcomes.

Normalization, DE and correlation no longer overwrite member run states on claim,
success, historical selection, failure or interrupted-owner recovery. Their durable
attempts, project-stage records and sealed frozen-cohort requests retain analysis
ownership, inputs and outcomes. Actual invalid-count evidence still invokes the
existing repair path; progress labels do not establish artifact validity. A
validated explicit matrix can be correlated using quantified members without
fabricated DE completion. The production DESeq2 library remains unchanged.

Migration 0021 snapshots all six old DE/correlation run labels with original errors
and update timestamps into immutable receipts, then projects those labels to
`QUANT_DONE`. Existing ingest records are untouched. The transaction creates no
generations, artifact registrations or successful project stages. Old-writer guards
reject subsequent analysis labels in runs. Upgrade holds the actual analysis
writer lock before opening SQLite's migration transaction; an active owner causes
a bounded idle-worker error. Already-upgraded readers can open during analysis.
Competing upgrade callers preserve one original receipt. The read-only shared
ownership reader explicitly admits schema 21, whose ownership columns are unchanged.

Progress totals, observed ingest events and raw bytes use active members only.
Empty cohorts and old analysis run labels cannot claim cohort completion. Mixed
quantification completion keeps all remaining runs in the throughput-based ETA.
The live display distinguishes actual running analysis, failure and idle quantified
cohorts, and counts runs without calling technical replicates biological samples.
Status prints individual cohort-stage outcomes/errors and migration receipt counts.

The final mandatory sequence passes in order: `cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **378 passed, 0 failed, 0 ignored** across
35 test/doc-test suites, terminal exit 0. Fifteen tests extend the 363-test baseline:
five migration/rollback/ownership guards, four progress projection checks and six
native integration checks including the child-worker helper. Native tests compare
entire run records after success/reuse and injected publication failure, narrow the
cohort without losing history, admit an explicit matrix without DE and kill actual
normalize/DE/correlation workers before recovery while independent ingest stays live.

The authoritative test log is `/tmp/rnaa-analysis-progress-final-all-checks-first.log`;
durable output, counts, scope, limitations and changed-source hashes are under
`validation/checks/2026-10-30-analysis-progress.*`. Formatting and Clippy completed
before tests in that same successful shell sequence. The fresh official registry
audit again checks 49 packages / 58 requirements with zero failures and the existing
documented Git exception. Both lockfiles are unchanged from the preceding verified
batch. See `docs/analysis-progress.md` for progress, upgrade and recovery semantics.

The overall goal remains active. Global CPU/RAM/disk reservations and shared index
admission, cross-project scientific reuse, unsealed legacy adoption/unindexed
transfer-support cleanup, prompt engine cancellation, unfinished CLI/product
contracts and measured production/release/migration/power-loss gates retain their
acceptance requirements. The DE engine replacement stays deferred to the end.

Verified implementation batch thirty-two (2026-10-30): shared quantification indexes
and bounded per-run exact-match caches.

The quant scheduler freezes the selected reference index descriptor before claiming
its first quant worker. A lazy group owns one read-only kallistors counting index;
concurrent calls and subsequent worker waves share the same allocation. Mapping
counts, readers, fragment/bias observations, EM state and output remain per run.
Healthy reconciliation/history returns before index loading. The graph is released
when its group and actual workers end; there is no permanent process-global cache.
The existing individual-call quantifier API remains available through a private
group. Every use checks the frozen path/size/SHA-256, first loading checks again
before caching, and computation checks again before manifest creation. Cached
graphs cannot authorize a changed reference. Normal loader errors install no value;
a loader panic fails subsequent group callers explicitly.

The first RSS probe exposed upstream automatic exact-match caches of 512 MiB for
one or two mapping threads, independently multiplied across RNAA runs. RNAA now
passes an explicit allowance from `quant.match_cache_mb`, default 64 MiB across
each run's mapping threads. `quant` and `run` accept `--match-cache-mb`; zero and
unaddressable values fail before allocation. The engine rounds thread tables down
within the allowance. New manifests record requested/resolved cache bytes plus
execution-only index group/loader/load-count observations, and status prints the
configured per-run allowance. Cache allocation/group sharing do not enter scientific
request identity. Strand changes still select distinct results, and historical
selection preserves original manifests and immutable generations.

Five real-index tests establish eight-loader pointer identity, lazy construction,
deallocation, failure/retry behavior, same-size corruption and reference mismatch.
Concurrent shared-index quantification produces exactly equal abundance tables and
run metrics to independent loads across single/paired layouts, all three strand
modes, thread counts and 1/8/512 MiB cache allowances. Two native CLI tests cover
multiple worker waves, resource-only cache/worker changes, historical strand
selection, rejected zero allowance, a failed run and repair while preserving the
other workers' successful generations. One core test covers cache defaults,
legacy configuration deserialization and invalid/overflowing allowances.

The separate-process development-profile probe covers 1,000/2,000 deterministic
600-base transcripts, one/four runs, shared/independent groups and 64/512 MiB caches.
All 16 cases map 64 reads per run and produce identical abundance hashes for each
reference. The 2,000-transcript/four-worker case reaches 2,189,590,528 bytes peak RSS
with independent indexes/512 MiB caches versus 297,009,152 bytes with one index/64 MiB
caches (86.4% lower). This small-reference probe is dominated by the match cache;
serialized index size does not establish resident graph memory. Short-run overlap
varies. The report records source hashes, wall time, toolchain and raw process RSS.
It is not global admission, production sizing or supported-host release evidence.
The earlier upstream-auto diagnostic remains separately marked as historical.

The full mandatory sequence passes on the final implementation:
`cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **386 passed, 0 failed, 0 ignored** across
36 test/doc-test suites, terminal exit 0. Eight tests extend the 378-test baseline.
The authoritative final log is `/tmp/rnaa-index-group-final-all-checks-first.log`;
durable output, scope, source hashes and summary are under
`validation/checks/2026-10-30-quant-indexes.*`. The current scale report and historical
diagnostic are under `validation/quant-index-scale`, and the runtime/reproduction
contract is `docs/quant-index-sharing.md`. The fresh official registry audit again
checks 49 packages / 58 requirements with zero failures and the documented Git
exception. Neither lockfile nor the production DE library changes in this batch.

The overall goal stays active. Shared CPU/RAM/disk reservation/admission, including
decoded multi-lane temporary files, remains the next operational gate; index/cache
sharing does not certify an overall memory ceiling. Cross-project scientific reuse,
unsealed legacy/unindexed transfer-support cleanup, prompt engine cancellation,
remaining product commands and measured production/release/migration/power-loss
gates keep their original acceptance requirements. DE replacement stays deferred.

Verified implementation batch thirty-three (2026-10-30): durable shared resource
reservation backend with cross-process ownership and recovery.

`rnaa_core::resources` now admits CPU threads, working-memory bytes and future disk
writes together in one SQLite transaction. Independent callers using the same
pool see the same reservations. Requests exceeding their own capacities fail
immediately; unavailable capacity returns structured wait reasons. Blocking waits
hold no catalog transaction, check caller ownership/cancellation before each retry
and again before returning admission, and have a caller-defined timeout. Resize
checks the replacement promise atomically and preserves the original on failure.
The backend does not promise FIFO fairness or operating-system scheduling/RSS caps.

The most restrictive active CPU/memory limits and disk headroom apply to new
admissions. A tighter proposal waits if existing usage already exceeds it, including
when requesting zero of that dimension. A new caller cannot silently relax another
live caller's capacity contract. Disk availability is probed through the pinned
`libc` POSIX interface at existing filesystem anchors. All outstanding disk promises
are conservatively summed across destinations, avoiding false independence between
APFS/overlay/subvolume paths with shared backing. Committed bytes remain accounted
by actual observed free space after future-write promises are released. This is
not project quota accounting or protection against unrelated writers.

A catalog-wide kernel lock covers connection/WAL selection, schema checks and
transactions. Every reservation also owns a UUID-named empty lock file. Descriptors
bind device/inode identity; changed, missing, aliased or symlinked ownership evidence
is held explicitly. No timeout or PID observation can release a live owner's
reservation. New admissions and snapshots reclaim dead owners only after acquiring
their exact recorded lock. Explicit recovery also removes empty, unlocked protocol
orphan locks left before failed commits, preserving unknown names, live locks and
nonempty diagnostics. Release refuses to unlink a substituted owner pathname.

Fourteen core tests, including the child helper, cover atomic limits across handles,
memory-constrained contention, restrictive and zero-dimension policies, resize,
future-disk growth/shrink, rollback, orphan maintenance, missing/replaced evidence,
cancellation before/after admission, timeout and unsupported schemas. Six real
child processes compete through multiple waves within one shared memory budget;
another real child is killed while holding capacity before replacement admission.
The final corrected six-process coordinator also passes **20 repeated runs**.
The first full verification exposed a test-harness release bug when overlapping
exits left a final child waiting alone. The coordinator now releases only fresh
owners and handles the final remainder; all mandatory checks were rerun afterward.

The final mandatory sequence passes: `cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features -- -D warnings`, then
`cargo test --workspace --all-features`: **400 passed, 0 failed, 0 ignored** across
36 test/doc-test suites, terminal exit 0. Fourteen tests extend the 386-test baseline.
The authoritative final log is `/tmp/rnaa-resource-pool-final-all-checks-second.log`;
durable output, source hashes, scope and summary are under
`validation/checks/2026-10-30-resource-pool.*`. The fresh official dependency audit
checks **50 registry packages / 59 direct requirements**, zero failures and the
documented Git exception. Its immutable feature snapshot is
`validation/dependencies/2026-10-30-resource-pool.json`.

`libc` 0.2.191 was already locked and is now a direct core dependency. The root
lockfile gains that dependency edge, advancing the existing lock-based scientific
request provenance; older generations remain retained. Other scientific library
versions and the vendor lockfile are unchanged. DE replacement remains deferred.

This is a backend delivery, not the end-to-end admission milestone. Pipeline
producers have not been connected and still use their existing scheduling checks.
Inspection confirms the next required preparation work: current preprocessing
checks each encoded file independently rather than the total decoded mate/lane
group, and native ordered-result maps can retain results beyond channel backlog
behind a slow early batch. Exact decoded-group preflight and a bound on the full
in-flight/reordering window must precede defensible phase claims. Integration must
also preserve shared-index guard lifetime, cheap historical reuse, publication
fences and cancellation. Host policy/default-pool selection, combined-producer
admission, actual usage qualification and production/release/power-loss gates remain
open. See `docs/resource-reservations.md` for the precise API and recovery contract.

The overall objective remains active, including cross-project scientific reuse,
remaining legacy/unindexed transfer-support cleanup, prompt engine cancellation,
unfinished product commands and the final replacement DE engine integration.

Verified implementation batch thirty-four (2026-10-30): exact decoded FASTQ
preflight and bounded native preprocessing preparation.

`FastqPlan` validates and hashes whole mate/lane groups before temporary payload
creation, counts every gzip member, bounds individual decoded lines and separates
wire bytes from normalized LF staging bytes. Materialization is capped to the
planned decoded/staging counts and rejects changed sources or shape; temporary
files are removed on failure. Download validation, preprocessing and quantification
use the configured contract. New settings are `fastq_max_line_mb=1` (positive) and
`fastq_max_decoded_mb=0` (uncapped). The existing preprocessing limit now means
aggregate decoded input, default 1024 MiB, rather than each encoded file separately.
Zero disables its additional cap; a configured global decoded cap still applies.

The vendored native preprocessing entrypoint is extracted into focused producer,
dispatch and window modules. Parallel input uses complete record-balanced pairs
and a credit held from before batch allocation through ordered-result writing,
including pending out-of-order results. Scoped threads join after errors/panic;
blocked credit acquisition is cancelled. The serial path also bounds whole records.
The compatible native API remains, with a configurable-options entrypoint added.
RNAA preprocessing scientific request revision advances to 3 for paired/newline
fixes; quantification revision remains unchanged. Resource-limit changes preserve
healthy cached outputs and avoid heavy preparation on reuse.

All mandatory workspace checks pass: formatting, strict all-target/all-feature
Clippy and **408 tests, 0 failed, 0 ignored across 37 suites**. Eight tests extend
the previous baseline, including two CLI tests for aggregate caps and healthy reuse
after limit changes. The standalone vendor library has **88 passing tests**,
including six native bounded-window/parity/failure tests. The WebAssembly library
check passes with `--lib`; the native CLI requires native features and is outside
that target check. Durable logs, source hashes and completion evidence are under
`validation/checks/2026-10-30-preparation.*`. The fresh dependency audit remains
**50 packages / 59 direct requirements, zero failures**, with the documented Git
exception, saved as `validation/dependencies/2026-10-30-preparation.json`.
Both lockfiles and scientific library versions are unchanged from batch thirty-three.

These preparation bounds do not connect producers to the shared admission backend
or establish an engine RSS ceiling. Common host policy, allocation-lifetime guards,
phase estimates and measured combined-workflow budgets remain the next resource
integration work. Cross-project reuse, legacy/unindexed transfer-support cleanup,
prompt engine cancellation, remaining product commands and supported-host release
qualification remain open. DE replacement stays deferred until the end. The overall
scope-weighted estimate remains about 70%, with a 65–75% range.

Verified implementation batch thirty-five (2026-10-30): common host policy and
actual owned preprocessing resource admission.

`[resources]` now configures a shared local pool, nominal CPU/memory capacity,
disk headroom and bounded capacity waits. Defaults use the same per-user local
pool across projects, available CPU parallelism, 75% of detected memory capacity,
1024 MiB headroom and a 300-second wait. macOS/Linux detection is supported;
standard Linux cgroup limits are intersected when present, with explicit capacity
required for other/nested layouts. These capacities do not measure current free
RAM or control unrelated programs.

Fresh CLI preprocessing admits bounded-line scanning, releases that claim, then
atomically admits engine CPU, structural memory, exact normalized staging and
bounded output promises. Releasing scan claims before engine waits prevents
partial-claim deadlocks. The engine guard outlives temporary payloads and all native
workers and covers its report/manifest writes. Native input additionally enforces
preflight decoded bytes, record counts and minimum/maximum record sizes, so shape
mutation cannot silently undermine those payload estimates. Metadata writes share
a checked byte allowance. Report execution observations include the estimate and
actual reservation descriptor without entering the scientific request identity.

Typed admission failures cannot become a raw-input preprocessing bypass. Healthy
reuse precedes opening the pool, including when its configured path is unavailable.
The ingest driver also cancels and joins remaining sibling workers on an error
instead of detaching them. Active native calls may delay joining until they return;
prompt engine interruption remains open. Native execution and estimate logic are
extracted from the preprocessing root into focused modules with separate tests.

The mandatory formatting, strict all-target/all-feature Clippy and full workspace
test sequence passes with terminal exit 0: **414 passed, 0 failed, 0 ignored across
38 suites**. Six tests extend the 408-test baseline, covering policy defaults and
overflow, typed owned admission/fences, scan-to-engine contention/cancellation,
metadata caps, actual CLI scan/engine rejection with bypass enabled, unavailable-pool
healthy reuse and joined scheduler shutdown. The standalone vendor library passes
**90 tests**, including two new byte/count/record-shape cases through serial and
parallel paths. Its WebAssembly library check also passes. Durable logs, source
hashes and exact scope are in `validation/checks/2026-10-30-preprocessing-reservations.*`.
The fresh dependency audit remains **50 packages / 59 direct requirements, zero
failures and one documented Git exception**; the immutable report is
`validation/dependencies/2026-10-30-preprocessing-reservations.json`. Both lockfiles
and scientific library versions are unchanged. DE replacement remains deferred.

This connects preprocessing, not the complete operational milestone. Quant/index,
downloads/references, analysis and final publication still require admission.
Index loading peak, shared resident lifetime, run-local EC/count/EM state and
packed-read queues must be budgeted without multiplying serialized index size by
an arbitrary factor. Shared-index holders must preserve progress capacity or allow
eviction/retry so independent projects cannot deadlock on working-memory waits.
Supported-host RSS/combined-workflow qualification, cross-project scientific reuse,
remaining cleanup/migration/product/release gates and final DE integration keep
their original acceptance criteria. Overall scope-weighted progress remains about
70% (65–75% range). See `docs/resource-reservations.md` for the exact current contract.

Verified implementation batch thirty-six (2026-10-30): bounded decoded index
layout planning and local loading-peak admission before quantification staging.

`QuantIndexPlan` scans the current counting loader's v13 layout with a fixed 64 KiB
input buffer and bounded bitmap scratch. Physical enclosing sections bound every
serialized count before loops; decoded graphs, transcript names and expanded
bitmaps are not materialized. Observations include decoded DNA/node counts,
minimizer buckets/bitmap values, MPHF representations/hash entries, EC blocks/values
and transcript/onlist sizes. Native/portable arrays, bitsets/runs, Bifrost inline,
tiny and pointer tags and special minimizer counts are supported. Header words
ignored by the engine are retained as observations rather than new format gates.

Initial resident/load estimates account for decoded structures, Vec capacity,
old/new hash tables during growth, serialized/decoded bitmap and node scratch,
transcript/shade metadata and loader stack allowances. Fresh public quantification
rejects an index loading estimate exceeding its declared memory capacity before
FASTQ staging or loader attempts. The plan is retained with the frozen group and
included in execution observations, with no scientific request revision/change.
Healthy cached counts remain reusable after limits are lowered. Owned CLI callbacks
now fence index/FASTQ preflight and staging and run between native numerical phases.
Quantification execution is extracted into `quant/native.rs`, leaving request
identity and reconciliation separate.

Five new layout tests compare real 1/2/8/120-target indexes with independent public
engine metadata/loaded graphs and exercise counts, truncation, overflow, callbacks,
bitmap variants and a valid special minimizer. They also catch the current tiny
RLE decoder's short look-ahead panic before it reaches the loader. A sub-1-MiB
bitmap encoding 2^32 values is counted without expansion, yielding a loading
estimate above 64 GiB and rejection against a 512 MiB capacity. Its expanding
engine decoder is never called. A group test proves memory rejection/cancellation
precede loading, and a CLI test proves rejection, successful retry and healthy reuse
after lowering the limit. Existing layout/strand/thread mapping parity, one-load
worker waves and complete native workflows continue to pass.

The final mandatory formatting, strict all-target/all-feature Clippy and full
workspace test sequence passes with terminal exit 0: **421 passed, 0 failed,
0 ignored across 38 suites**. Seven tests extend the 414-test baseline. An initial
full sequence also passed; the final sequence was rerun after strengthening
temporary hash/EC growth allowances. Durable logs, source hashes and coverage are
in `validation/checks/2026-10-30-index-layout.*`. The fresh dependency audit checks
**50 packages / 59 direct requirements**, zero failures and the documented Git
exception, saved in `validation/dependencies/2026-10-30-index-layout.json`.
Both lockfiles and scientific library versions are unchanged. The prior standalone
preprocessing/WASM evidence remains applicable with all its source hashes
revalidated; that standalone suite was not repeated this batch. DE stays unchanged.

This is local preflight admission, not shared quantification reservation delivery.
Mapping/EC/count/EM state, caches, staging and outputs still need combined requests;
resident index ownership and eviction/progress guarantees remain open. The inspector
bounds/skips opaque DNA/rank/positional fields and special marker widths rather
than proving every engine semantic. Memory constants need supported-host RSS
qualification, and active engines still need prompt interruption. Broader global
resources, cross-project reuse, remaining cleanup/product/release gates and final
deferred DE integration retain their original acceptance criteria. Overall progress
remains about 70% (65–75%). See `docs/quant-index-resources.md` for the precise scope.

Verified implementation batch thirty-seven (2026-10-30): owned shared quantification
scan/index/run admission and resident graph lifecycle.

The resource backend now admits bundles in one catalog transaction while returning
independent lifetime guards. Summed CPU/memory/disk capacity is checked before any
claim is visible. Failed inserts roll back every row without releasing a guard
inside the owning transaction; unlocked empty owner files remain recoverable.
Single-claim admission delegates to the same path. Thread and independent-process
checks verify all-or-none ownership, independent release, resize, cancellation and
the shared budget.

Fresh owned CLI quantification admits bounded layout/FASTQ scanning, releases that
claim, then atomically admits the cold index and first run before staging. Its
combined memory covers the larger of loading plus retained setup and resident
graph plus mapping/EM work. Unused first-run capacity temporarily covers loading;
the index claim shrinks after construction and stays with the graph. Later runs
reserve independent working CPU/memory/staging/output capacity. Healthy reuse
precedes pool opening, and resource settings remain outside scientific identity.

An active `Arc` borrower protects its graph. Idle graph eviction on memory waits
when no CPU owner can progress avoids holding graphs while every group waits for
working memory. The scheduler also releases idle graphs when it has no ready or
active quant work, so preprocessing can feed another wave. CPU/disk contention
alone does not evict useful indexes. Ordinary waves still load once; deliberate
eviction records its counter and can cause a reload. Groups reject mixed pools
and admitted/unadmitted graph use. Waits retain no borrowed graph reference and
check ownership, cancellation and the configured timeout.

`kallistors-mapping-em-v1` accounts for per-run packed read queues, conservative EC
class/member envelopes, sparse EM state, target arrays, lookup caches, worker
stacks and staging/output promises. Paired logical jobs can retain multiple 6,000-
record parser batches while filling an approximately 8 MiB sequence job; minimum
decoded sequence length bounds that capacity. Bounded readers enforce total bytes,
record/line sizes and count, and minimum/maximum sequence lengths before engine
batch allocation. Unsupported 32-bit counter/sequence spans reject explicitly.
Shade references use a target-union envelope. The private default 8 MiB restored
match table is included even with smaller user cache settings. A joined scoped
mapper thread ends its TLS allocations before the run guard drops. Temporary files,
mapping/EM data and engine output/manifest writes also finish under that guard.
Execution manifests record estimates and independently owned index/run descriptors.

Final mandatory root sequence passes formatting, strict all-feature Clippy and
**437 tests, 0 failed, 0 ignored across 39 suites**. Sixteen tests extend the prior
421-test checkpoint: five backend bundle/process cases, nine quant resource/input/
resident cases and two new CLI integrations. The six-process bundled contention
test additionally passes ten repetitions. Two actual CLI processes reach a
confirmed shared-capacity wait; cancelling one leaves no partial claim, and the
surviving project finishes four runs with one resident index claim. A constrained
224 MiB preprocessing/quantification fixture finishes every worker wave. Admitted
paired mapping matches the present engine's abundance TSV exactly. Durable logs,
source hashes and precise limits are in `validation/checks/2026-10-30-quant-admission.*`.

A fresh registry audit passes **50 packages / 59 direct requirements**, no failures
and the same documented Git exception. Its immutable report is
`validation/dependencies/2026-10-30-quant-admission.json`. Both lockfiles and the
production scientific libraries are unchanged. Nine prior vendored preprocessing
source witnesses still match; the standalone native/WASM suite was not rerun.
DESeq2 remains on the present library, with replacement integrated last.

This completes quant scan/index/run integration, not the global operational
milestone. EC ambiguity envelopes can reject jobs with lower actual observed
diversity; constants require supported-host RSS and production-scale qualification.
Compatible library calls without an execution context retain local preflight but
do not join shared run admission. Shared-blob copying, generation sealing and
publication/retention, downloads/references and analysis still need admission.
Prompt native cancellation, cross-project reuse, remaining cleanup/product/release
gates and the deferred final DE integration remain open. Overall progress remains
about **70% (65–75%)**. See `docs/resource-reservations.md` for the current contract.

Verified implementation batch thirty-eight (2026-10-30): owned download admission
from private-cache verification through complete-generation publication.

Fresh downloads and owned input recovery now use the common host pool. Cache
verification, resume planning, streaming transfer, complete FASTQ validation and
publication obtain separate claims, releasing the preceding claim before waiting
for the next. Healthy historical/current reuse still precedes opening the pool.
Resource settings remain outside download request revision 2; scientific libraries,
both lockfiles and the deferred production DESeq2 library are unchanged.

Resume planning uses one CPU, bounded checkpoint metadata and a 16 MiB memory
allowance. It verifies the durable prefix rather than trusting a file's raw length;
ownership checks during hashing cannot discard a prefix on cancellation. Transfer
uses one CPU and a 64 MiB structural memory allowance. Known totals reserve only
remaining writes after the proven prefix. Unknown totals use checked 8 MiB windows
without imposing a total payload-size cap. Capacity is checked before each window
extension; a refused increase saves its proven prefix, returns a typed admission
error and releases claims without a network retry loop. Header size discovery that
needs more capacity unwinds the response and all claims before readmission waits.
An ignored Range reclaims the private old prefix before reserving replacement
writes. Checkpoint progress reduces known future-write promises after synchronization.
Transfer deadlines and current job fences also bound capacity waits.

Private verification/receipt writes have descriptor-dependent metadata bounds.
Completed FASTQ groups use one common scan-request model shared with preprocessing
and quantification. Publication admits inventory-dependent CPU/memory and the full
extra copy plus metadata before starting a generation. Both the generation store
and project catalog destinations are probed, preserving the conservative shared-
backing accounting contract. The copy reader caps writes to the observed source
bytes and checks ownership periodically; source growth/truncation or cancellation
discards temporary output. The publication guard outlives copy/hash, manifest,
seal, registration and transfer-cache removal. Existing source/cohort/attempt
fences and digest checks still determine validity.

The downloader now separates orchestration, private cache handling, resource
estimates, history selection and publication; protocol/admission tests and shared
origin fixtures live in separate files. This keeps resource guards beside their
allocation/write lifetimes instead of adding another monolithic workflow.

Final mandatory root sequence passes formatting, strict all-feature Clippy and
**444 tests, 0 failed, 0 ignored across 40 suites**. Seven tests extend the prior
437-test checkpoint: four HTTP admission/prefix cases, one capped publication-copy
case and two actual CLI integrations. The four HTTP admission cases additionally
pass ten repetitions using the final workspace test binary. Two CLI processes
reach confirmed capacity waits before any HTTP GET; cancellation releases no
partial claim and the survivor publishes under its own reservation. Publication
memory refusal preserves the previous head and the verified private cache; a
larger-capacity retry avoids refetching. Checksum-backed history remains selectable
with an unavailable pool and a lowered limit. The existing seven download history
tests and all native preprocessing/quant/analysis/recovery gates remain green.
Durable logs, source hashes and limits are in
`validation/checks/2026-10-30-download-admission.*`.

The fresh stable dependency audit checks **50 registry packages / 59 requirements**,
with zero failures and the same documented Git exception, in
`validation/dependencies/2026-10-30-download-admission.json`. Nine prior vendored
preprocessing witnesses still match; the standalone native/WASM suite was not rerun.

Reference preparation, analysis and other producer publication/retention admission
remain open. Compatible unowned library download calls and reference HTTP calls
are still unadmitted. Estimates require supported-host RSS/throughput and combined-
workflow qualification; future disk promises are cooperative, not a filesystem
quota. Transfer/support cleanup, cross-project semantic reuse, prompt native
cancellation, remaining CLI/release gates and final DE replacement remain open.
Overall progress stays about **70% (65–75%)**. See `docs/resource-reservations.md`.

Verified implementation batch thirty-nine (2026-10-30): owned reference admission,
bounded external planning and complete multimember gzip preparation.

Fresh references now reserve separate source, planning and build/publication phases
from the common pool. A bounded decoded scanner computes exact clean canonical
k-mer cardinality using sorted external runs and bounded merges; ambiguous windows
retain a conservative upper bound. Line and aggregate decoded limits are enforced
before heavy allocation. Annotation planning and production share normalization
rules. Structural phase estimates cover the native build, annotation, source copies,
decoded staging, output inventory and catalog writes. The native builder runs in a
joined scoped thread, with its actual clipped thread count recorded in provenance.

All planning and decoded staging files belong to the unpublished generation.
Cancellation removes live temporary work; a killed planner leaves tracked work that
existing unpublished-generation trash and expired purge can recover. Independent
references using the same source URL have separate private cache ownership.
Complete gzip decoding consumes all members before the upstream native builder;
reference request revision 2 records that corrected source contract. Healthy history
reuse still precedes pool opening and resource-only changes preserve scientific
reuse. The production scientific libraries and both lockfiles are unchanged.

The final mandatory root sequence passes formatting, strict all-feature Clippy and
**455 tests, 0 failed, 0 ignored across 41 suites**. Eleven tests extend batch
thirty-eight: eight stage/planner cases and three actual CLI resource/cancellation/
crash-cleanup cases. Ten additional repetitions of the two external merge/refusal
tests and three CLI cases pass: **50 extra case executions**. Current release-build
profiles cover 1,000 and 20,000 synthetic transcripts with one/two threads. At 20,000,
measured process RSS is approximately 185/148 MiB and elapsed time is 1.32/1.03 s;
conservative build claims are approximately 2,785/2,789 MiB. Both thread settings
retain matching scientific request/table witnesses, with no remaining resource
claims or workspace files. These synthetic profiles do not qualify production
cohorts, other hosts or combined workflow budgets.

Durable logs and final source hashes are in
`validation/checks/2026-10-30-reference-admission.*`; the current binary-bound profile
is `validation/reference-scale/2026-10-30-local-release.json`. The fresh dependency
audit again covers **50 registry packages / 59 requirements**, zero failures and one
documented Git exception. Nine unchanged vendored source witnesses are revalidated;
the prior standalone native/WASM evidence remains a separately identified checkpoint.

Analysis and publication outside downloads/references, cross-project semantic reuse,
remaining cleanup classes, prompt native cancellation, production qualification and
the deferred DESeq2 replacement remain open. Overall progress remains approximately
**70% (65–75%)**.

Release consolidation (2026-10-30), following the user's request to prepare a release:
new feature expansion is set aside while the current implementation is gathered into
a reviewable preview candidate. `CHANGELOG.md`, `docs/releasing.md`, the offline
installation example and the candidate packager describe scope, upgrades, versioning,
source provenance and distribution gates. Candidate archives and checksums are local
review artifacts. Public release remains gated by distribution licensing (the pinned
mincorr declares GPL-2.0-or-later), supported-host/minimum-toolchain qualification,
the final public version/clean commit and release approval. Production readiness
also requires the outstanding operational/scientific gates; no tag or publication
is implied by preparing the candidate. See `validation/releases/2026-10-30-preview.json`
for the installation and archive checks performed on this candidate.

Release continuation (2026-10-30): the user excluded mincorr licensing review from
the current blockers and explicitly retained rsdeseq2 0.3.1. Version **0.3.0-rc.1**
now identifies the preview's incompatible CLI/catalog changes; the root lockfile
changes only the four workspace versions. No dependency selection or vendor
lockfile changes. Native Ubuntu 24.04.5/Linux ARM64 and macOS ARM64 both pass the
mandatory formatting/strict Clippy/test sequence: **455 tests in 41 suites on each**.
Linux minimum-toolchain compilation and the complete mandatory sequence also pass
on Rust 1.97.1, with the same 455 tests / 41 suites. Runtime qualification and
archive/installation results are retained in
`validation/releases/2026-10-30-linux-rc.json`.

The candidate packager refuses stale executable versions, preserves source hashes
and upstream notices, and records deferred work without treating licensing as a
current technical blocker. Native Linux preprocessing also passes all 90 tests.
The new verifier checks archive/inventory hashes,
executable/source version agreement, installs from extracted source and runs
explicit/automatic-resource examples. Binary bundles now include operator guides
and examples and record the Linux artifact's glibc 2.39 runtime floor.
Corrupt metadata, unsafe checksum paths and
stale executables are refused before execution or candidate creation. A read-only
CI matrix now covers Linux x86-64/ARM64, macOS ARM64 and Linux x86-64 on the minimum
compiler, including packaging and source installation; current stable action
versions and workflow lint have separate evidence. Configuring the matrix does not
claim hosted jobs have executed. Final clean-commit/hosted-CI/publication steps
remain separate, while the larger production roadmap and final DE replacement
remain deferred. No public release has been published.

Follow-up (2026-10-31): `doctor` now inspects a project without initializing it,
opening the catalog read-only for SQLite `quick_check` and the project record,
validating the on-disk configuration and resource policy, and reporting linked
engine versions. `doctor --json` emits a schema-versioned report and returns
failure for catalog, configuration or resource-detection errors. A missing project
or missing configuration is identified explicitly; a missing on-disk config may
fall back to the catalog during normal operation and is not audited by doctor.
This is not a network, full integrity or scientific provenance check. The CLI
integration covers an absent project, healthy project, malformed configuration
and corrupt catalog. SQLite may manage WAL shared-memory files on a live catalog.

The stable-dependency audit now reports 50 registry packages / 59 direct
requirements with zero failures and the same documented Git exception. The vendored
native S3 client is 1.157.0; both lockfiles resolve aws-smithy-schema 1.1.2. The
changed root lock digest can invalidate scientific cache requests despite unchanged
scientific libraries; retained history remains available. Live S3 service behavior
was not tested. rsdeseq2 remains 0.3.1 as requested. On local macOS ARM64 the
mandatory format/strict-Clippy/test sequence passes **456 tests / 42 suites**; all
90 standalone vendor tests and Rust 1.97.1 all-target compilation also pass. See
`validation/checks/2026-10-31-doctor.json`. The overall roadmap remains open.

Subsequent diagnostic follow-up (2026-10-31): when `rnaa.toml` is absent, doctor
now parses and validates the catalog configuration actually used as the startup
fallback, then checks its resource limits. A present file retains precedence.
Integration coverage includes valid and malformed catalog JSON, an unsupported
catalog schema and a valid file overriding a broken catalog copy. JSON engine
versions also include mincorr, its pinned revision and rsfgsea. The earlier doctor
evidence is a historical source snapshot; the current verification is recorded in
`validation/checks/2026-10-31-doctor-fallback.json`.

Structured status follow-up (2026-10-31): `status --json` exposes versioned
project, reference, cohort, per-run eligibility, stage, generation, storage and
progress data without parsing the human display. It performs the same recovery and
reconciliation as text status; the result is observational rather than an atomic
snapshot under concurrent work. A cohort that cannot be verified reports its reason
and null active runs/progress/QC count instead of asserting a zero-run cohort.
Native integration cases cover completed analysis, a narrowed active cohort and a
stale cohort. The mandatory macOS ARM64 sequence passes **456 tests / 42 suites**;
Rust 1.97.1 all-target/all-feature compilation also passes. See
`validation/checks/2026-10-31-status-json.json`. Export, survey, platform CI and
the other outstanding roadmap gates remain open.

Frozen export follow-up (2026-10-31): `export` now selects registered analysis
tables from sealed DE and correlation heads, verifies generation inventories before
and after copying, and writes TSV or CSV with source/output hashes and row/column
counts. Source manifests are retained in the export. Project job ownership, a
generation reader pin and common resource-pool reservations protect the operation;
per-output limits constrain writes. A no-replacement rename and fenced SQLite job
completion publish a complete directory, with rollback on database finalization
failure. Independent tests cover native output equality, quoted fields, capacity
rejection, cancellation, corruption, no-replacement and transaction rollback.

On macOS ARM64 the mandatory format/strict-Clippy/test sequence passes **462 tests
/ 43 suites**. The export integration and three publication cases also pass on
Ubuntu 24.04 ARM64; a fourth conversion-cancellation unit case was added afterward
and passed on macOS. Rust 1.97.1 all-target compilation and the fresh dependency
audit pass. See `validation/checks/2026-10-31-export.json` and `docs/export.md`.
This exports selected frozen generations without proving current scientific
freshness or cross-group consistency. Very conservative write memory estimates,
unbounded inspection-inventory memory, kill/power-loss cleanup and native x86-64
qualification remain open. Survey and library-level export interfaces remain open.

## Remaining implementation sequence and acceptance gates

Linux continuation (2026-10-31): the current uncommitted source passed the full
formatting, strict Clippy and workspace test sequence on Ubuntu 24.04.5 ARM64:
479 tests in 45 suites. The VM used one build job and reduced debug information
to fit its memory and disk. The record includes earlier resource failures and
source witnesses: `validation/checks/2026-10-31-linux-current.json`. This does
not qualify x86-64, production scale, a final commit, or earlier packages.

Packaging continuation (2026-10-31): archive verification now checks every
source-derived bundled file, the complete declared notice inventory, archive roots
and duplicate archive entries. Failure-focused Python tests run in CI. A fresh local
macOS ARM64 candidate from the current source was verified and installed from its
source archive, using an existing build cache; bundled and installed executables
passed both offline resource modes. The report and logs are under
`validation/releases/2026-10-31-macos-current-*`. This is a dirty local snapshot,
not a clean final release or a Linux/x86-64 package qualification.

Candidate source binding (2026-10-31): the packager now builds the executable
from a frozen copy of the source archive inputs and checks that copy for changes
before packaging. An arbitrary prebuilt binary can no longer be substituted based
only on its version string. The frozen-source macOS ARM64 candidate passed archive
verification and source installation with both offline resource modes; compilation
caches were reused. See `validation/releases/2026-10-31-macos-frozen-*`. This local
run does not establish reproducible builds or qualify other architectures. The
same frozen-source packaging and installation checks also passed on Ubuntu
24.04.5 ARM64; Python packaging tests passed there and transferred archive hashes
were rechecked. See `validation/releases/2026-10-31-linux-frozen-*`. The macOS
and Linux source snapshots differ in this plan and validation records only; neither
is a clean final release.

Current module validation (2026-10-31): correlation now records the most recent
module options and resolved pathway path in project configuration and fences
publication against changes to those settings. Strict export and survey rebuild
module requests and verify pathway bytes. Older module results without current
settings still fail closed; independent matrix support is described below.
Running correlation without modules clears the saved module selection.

Current external-matrix validation (2026-10-31): correlation records its selected
external matrix path, and strict export reconstructs the correlation request and
current metadata from that file. Correlation-only export no longer needs a DE
generation for this case. The selected file remains external and its upstream
provenance is not inferred. Older selections lacking configuration may need a
rerun; survey retains its separate DE/VST compatibility requirements.
Metadata preparation now reads only the admitted-size-bounded external matrix
header and uses indexed sample lookup, while correlation execution retains full
matrix validation. The design guard also fences the selected matrix path.
Strict reconstruction and model preparation also use the admitted matrix size
when reopening or hashing the input, rejecting size changes before publication.

Release workflow continuation (2026-10-31): a fresh four-run public yeast project
completed resolution, checksum-verified downloads, quantification, DE, correlation,
strict export and survey on macOS ARM64. Rerunning quantification and analysis
kept all 11 selected generation heads. The tested tutorial and machine-readable
record are in `examples/public_yeast.md` and `validation/real-data-smoke/`.
This is a small execution check, not biological validation or a production-scale
qualification. A source-install attempt exposed stale workspace artifacts when
reusing a build cache; packaging now scopes its cache to the complete source
snapshot, and verification installs with a fresh isolated target.

Artifact witnesses (2026-10-31): snapshots and local manifest verification hash
one open regular file, limit reads to its observed size plus one byte, and reject
observed resize, replacement or metadata changes. This prevents combining a hash
from one observation with a size from another; it is not a lock against external
writers or a guarantee against every possible concurrent edit.

This section is the current execution order. The longer historical sections below
retain the original design discussion; their statements of what already works are
not a substitute for this ledger or executable checks. Each item should land as a
reviewable vertical change with its failure tests, documentation and schema upgrade.

1. **Durable attempts and concurrent ownership.** Replace global interruption
   recovery in `crates/rnaa/src/main.rs` with a `jobs`/`attempts` model in core.
   Claim work in an immediate SQLite transaction, recording stage, subject, request
   digest, owner, lease expiry and a monotonically increasing fencing token.
   Heartbeats renew only the matching live attempt. Completion and artifact
   publication must reject stale fencing tokens. Recover only expired attempts;
   opening another command must not reset another process's active download.
   Add explicit cancellation and structured outcomes. A finite download/run command
   must return a failure status if requested work failed, even when workers continued.
   Acceptance: two processes competing for one run execute it once; independent
   downloader and orchestrator coexist; killing a worker permits expiry/retry;
   a stale worker cannot publish or erase a newer worker's successful result.
   Current delivery includes claims, leases, heartbeats, fenced database mutations,
   cooperative cancellation, ownership-aware recovery and finite failure outcomes.
   Download/preprocessing/quantification/analysis generation publication is delivered.
   Reference preparation also has independent ownership and fenced publication.
   Pending acceptance includes prompt cancellation inside long computations;
   a paused live writer currently retains
   physical exclusion. Per-run/cohort state separation is covered by item 2.

2. **Persist cohorts and separate eligibility from progress.** Store selector
   expressions and resolved inclusion decisions with a cohort revision. Frozen
   analysis inputs must refer to that revision and explicit sample/run membership.
   Keep per-run ingest/quant status separate from cohort DE/correlation attempts;
   avoid overwriting all member runs with project-wide analysis state. QC gate
   decisions should have their own fingerprint and retain enough report evidence to
   reevaluate after raw and trimmed reads have both been removed. Define biological
   sample grouping and aggregation of technical runs before fitting. Add explicit
   numeric/categorical metadata types, missing-value rules and reference levels.
   Acceptance: selector changes exclude previously included runs deterministically;
   changing QC policy after retention updates eligibility without demanding deleted
   reads; two runs for one biological sample do not become two independent replicates;
   rank, replication and contrast failures are reported before expensive work.
   Durable membership, retained QC admission, biological pooling and typed design
   preflight are delivered. Batch thirty-one delivers per-run/cohort progress
   separation, audited legacy projection and interrupted-analysis recovery without
   changing run ingest records. Controlled metadata refresh with retained edit
   history, immutable source witnesses and semantic cache projection is delivered.

3. **Publish immutable generations.** Give each attempt a staging directory and
   produce all outputs there. Verify sizes, hashes, required file sets and manifest
   schemas before exposing a generation. Publish a single generation pointer, then
   register artifacts, attempt outcome, project state and events transactionally.
   Keep prior successful generations available until ownership/retention rules allow
   removal. Startup recovery should either attach a verified completed generation or
   discard an unpublished one. Acceptance: injected termination after each file,
   manifest, rename and database boundary never exposes a mixed set of generations;
   database registration failures do not leave a falsely complete stage; legacy
   hardlinked blobs and the previous valid analysis survive failed replacement.
   Download, preprocessing, quantification, normalization/DE and correlation satisfy
   publication/rollback and pinned-reader checks, including actual crash boundaries.
   Raw/trimmed generation pruning has committed-consumer receipts. Remaining scope
   includes unresolved legacy evidence/unsealed adoption, abandoned transfer/support
   recovery, unindexed ownership and admission reservations.
   Observed logical/inode accounting
   and supported sealed legacy ancestry import are delivered.
   FASTQ/SRA-specific reversible file cleanup is delivered. Native generation dependency/proof reachability, maintenance pins and
   reversible generation GC are delivered. Reference generation publication is delivered.
   Tracked unpublished payloads from terminal attempts support reversible diagnostic
   trash/restore and explicit expired purge without attaching them as results.
   Unindexed directories and transfer/support data remain held for their own
   ownership and recovery protocols.

4. **Complete artifact dependency and reuse semantics.** Store typed requests as
   first-class objects and exact upstream generation edges. Resolve a shared request
   digest to immutable artifacts before recomputation. Validate local views against
   the same generation, including retained/deleted input lineage. Make reference
   requests bind actual source content, annotation-generation policy, index options
   and engine identity, rather than treating an index that loads as proof it matches
   the current reference. Acceptance: two projects with identical semantic inputs
   share one result; replacing input content, adding/removing a lane, changing an
   annotation or scientific option invalidates all dependent outputs; CPU-only
   allocation changes preserve reuse; incomplete or corrupt shared generations are
   never attached as successful results.
   Exact project-scoped historical request lookup/atomic selection is delivered
   for normalization, DE, correlation, preprocessing and native counts, with immutable
   producer history, full registry restoration and current admission fencing. Full
   typed requests and use-time input bindings/retention receipts are cataloged.
   Download/reference history is delivered with current source authority checks.
   QC report reconstruction can recover original inputs and reuse matching counts
   independently of missing/damaged evidence. Cross-project sharing, complete
   legacy ancestry and generalized dependency recovery remain open.

5. **Independent scientific parity gates.** Generate and check in small expected
   outputs from pinned tximport/DESeq2 and independent correlation reference versions,
   recording scripts, versions, source hashes and numerical tolerances. Cover raw,
   scaledTPM and lengthScaledTPM imports, fractional counts and rounding, zero
   abundance/length cases, GENCODE version/PAR_Y identifiers, offsets, numeric and
   categorical designs, interactions, contrast direction and multiple testing.
   For each correlation method cover ties, constant rows, small samples and a known
   reference matrix. Keep the existing synthetic successful workflow as a separate
   integration gate. Acceptance: expected outputs are independent of RNAA's code,
   failures identify the numerical quantity and tolerance, and engine upgrades must
   pass the gate or deliberately bump the scientific algorithm revision with evidence.
   Current delivery includes pinned offline R fixtures, strict provenance/integrity,
   all import modes and technical-run splits, offsets, numeric/categorical designs,
   interactions, both contrast directions, filtering/multiple testing, VST and
   all five correlation methods including independent edge topology references.
   The present rsdeseq2 library remains unchanged; the user has explicitly deferred
   its replacement until the end. Broader DE engine fixes are separate work, and
   the eventual replacement must rerun these gates.

6. **Bound correlation memory and make module semantics explicit.** Replace the
   former full dense allocation for every output mode with
   method-specific reusable row preparation and blockwise pair computation, with
   bounded per-row heaps for top-k and streaming threshold edges. Preserve exact
   tie ordering, diagonal exclusion and edge deduplication semantics. Dense output
   requires an explicit estimated memory/disk admission check. Module enrichment
   formerly built another dense Spearman matrix and called an arithmetic mean an
   eigengene: implement standardized first-principal-component scores with a
   deterministic sign convention or rename/document the existing statistic. Validate
   GMT overlap, minimum size and enrichment parameters before fitting, and apply a
   thread/memory budget to module work. Acceptance: block results match a small dense
   oracle for all supported methods; peak memory for top-k grows with input rows and
   k rather than the square of the gene count; module score tests distinguish PCA
   from a mean; representative scale runs record wall time, RAM and output size.
   Bounded tiles, reusable preparation for Pearson/Spearman/bicor, bounded per-row
   heaps, dense/threshold stripes, conservative payload admission and deterministic
   edge/reference gates are delivered. Kendall/HellCor use pinned cross-matrix
   APIs, with repeated private preparation per tile. Batch thirty implements the
   module PCA/parameter/resource work and independent R PCA/graph/ranking/fgsea
   references. Full component membership is distinct from report truncation; sizes
   and PCA input use every member. Strict pathway overlap/options are checked before
   graph work, undefined constant genes are excluded and perfect-correlation ranking
   policy is explicit. An owned local pool and conservative estimates cover module
   graph/connectivity/PCA/enrichment. Production cohort, supported-host and release
   performance gates remain open; phase estimates are not global reservations.

7. **Admission control and robust transfers.** Reserve CPU, RAM and temporary/final
   disk space before admitting work, including decoded multi-lane FASTQ staging.
   Load one immutable quant index per compatible worker group and share it safely.
   Replace repeated directory scans with tracked reservations and reconcile actual
   usage after completion. For HTTP resumes, validate Content-Range start/total,
   object identity and expected length; safely restart on ignored ranges or changed
   objects. Add bounded retries, jitter, stall detection and cancellation checks.
   Acceptance: local HTTP fixtures cover 200/206/416, changed ETag/length, truncated
   bodies and stalled responses; combined download/quant/preprocess work stays
   within declared budgets; interrupted transfers never become verified artifacts.
   Batch thirty-two delivers one lazy immutable index per compatible in-process
   quant worker group and explicit per-run match-cache allowances, with real-index,
   native reuse/failure and separate-process RSS evidence. Global CPU/RAM/disk
   reservations, decoded staging admission and supported-host combined-budget
   qualification remain open; index file bytes are not a resident-memory estimate.
   Batch thirty-three delivers the shared admission ledger/API and independent
   process/kill/rollback evidence. Batch thirty-four delivers exact decoded mate/lane
   accounting, bounded FASTQ lines, capped normalized staging and the full native
   preprocessing in-flight/result-reordering window. Batch thirty-five connects
   fresh preprocessing scan/engine phases to a common host policy/pool, with exact
   staging and bounded output promises, estimate/claim observations, cancellation
   during waits, non-bypassable admission errors and cache reuse before pool opening.
   Native readers enforce preflight byte/count/record bounds; scheduler errors now
   cancel and join sibling workers. Other producers are not connected yet and
   memory estimates do not establish a process RSS ceiling. Next wire quant/index
   phase requests and shared-index lifetime, then downloads/references, analysis
   and publication, and verify measured combined-workflow budgets.
   Batch thirty-six adds bounded v13 index-layout inspection and initial decoded
   resident/loading estimates, a per-job peak check before FASTQ staging, and owned
   preflight/phase checkpoints. This does not connect shared quant reservations.
   The next index/run integration must admit loading and first-run capacity
   atomically, retain resident ownership, and avoid idle-index/working-memory
   deadlocks while preserving preprocessing progress.
   Batch thirty-seven connects shared scan/index/run admission, atomic cold bundles
   covering sequential loading/mapping peaks, resident lifetime guards and idle
   eviction. Decoded sequence/record bounds cover packed single/paired queues;
   dedicated joined mapper threads release private TLS caches before run claims.
   Independent CLI contention/cancellation and constrained preprocessing feedback
   are covered. Downloads/references, analysis, final shared-blob copy/sealing/
   publication and measured combined/RSS qualification remain open. EC diversity
   envelopes are conservative and can reject workloads with lower observed use.
   Batch thirty-eight connects owned downloads, including recovery, through cache
   verification, resume planning, transfer, FASTQ scan and publication phases.
   Header size discovery releases responses/claims before readmission; unknown
   stream-window growth refuses before unpromised writes and saves its prefix.
   Complete-copy publication is admitted before generation creation and its reader
   caps bytes despite source growth. Resource failure preserves previous heads;
   verified caches and checksum-backed history remain reusable. References,
   analysis and publication outside downloads are the next admission work, with
   supported-host combined/RSS/throughput qualification still required.
   Batch thirty-nine connects reference source acquisition, bounded external planning,
   native building and complete publication to owned phase claims. Full multimember
   gzip staging is byte-bounded, temporary work belongs to unpublished generations,
   and actual process cancellation/kill/recovery tests cover ownership and cleanup.
   Analysis and other publication paths plus measured supported-host combined
   workflow qualification remain open.

   Fresh correlation computation participates in shared CPU admission after
   cache reuse checks. Its network phase now atomically claims CPU, estimated
   memory and primary output space; it releases the earlier CPU claim before
   waiting. Its writer enforces the primary output promise, and the network
   engine checks cancellation periodically between tiles. Cancellation, head
   preservation and contention are tested. Other
   correlation phases, other analysis admission and measured combined-workflow
   qualification remain open.
   Custom matrix sample selection now reads a witnessed header instead of loading
   the numerical matrix before reuse, and an all-gene request skips a redundant
   input load. Filtered selections now validate input row by row and retain
   identifiers and selected positions instead of the full numerical input, with
   bounded top-variance ranking and cancellation checks. Adjustment now streams
   rows under an owned estimated memory and output reservation with input and
   output byte bounds; the model-design sidecar, selection, module work and
   measured combined-workflow qualification remain open.

8. **Shared-store ownership, retention and cleanup.** Introduce project/generation
   references, pins and a transactional GC mark phase in the shared store. Keep
   deletion reversible during the configured trash window and protect active
   attempts and frozen analyses. Implement `TrashCleaner` behavior with an
   inspectable dry run and consistent accounting. Acceptance:
   removing one project's view cannot delete another project's payload; pinned or
   live generations survive GC; interrupted cleanup resumes safely; retention and
   restore tests verify both disk state and database references.
   Shared-blob catalog ownership, recoverable publication handoff, reader/writer
   exclusion, orphan GC preview/apply and recorded-deadline trash/restore/purge are
   delivered. Historical generation owners remain protected until completed retirement.
   Project-generation retirement, binding/pruning reachability, maintenance pins and
   project generation trash/restore/purge are delivered. Count-backed FASTQ/SRA
   trash/restore/purge, automatic `TrashCleaner` and receipt-aware input recovery are
   delivered. Observed logical/inode accounting and tracked unpublished diagnostic
   trash/restore/purge are delivered. Shared schema 2 adds read-only legacy owner
   import with explicit complete project census and accounted hardlinks. Unknown
   projects/evidence and unexplained links remain held; unsealed adoption,
   unindexed/transfer/support cleanup, admission reservations and measured release
   gates remain open.

9. **Extract workflows and finish the CLI contracts.** Move orchestration,
   reconciliation, artifact registration and scientific request construction out of
   the roughly 4,800-line CLI module into testable workflow services. Add validated
   effective configuration, structured command results and real machine-readable
   doctor/status output. Implement export of canonical artifacts with provenance
   and an actual survey workflow in place of both CLI stubs. Choose and document
   canonical matrix storage honestly: current analysis matrices remain TSV, despite
   the historical HDF5-first target. A format migration needs schema/versioned
   readers and validated conversion, not filename changes. Acceptance: each command
   has success/failure integration coverage, JSON is parseable without scraping
   terminal text, exported dimensions/IDs match the frozen analysis, and library
   workflows can run without invoking CLI parsing.

10. **Release evidence and packaging.** Run scientific parity, corruption,
    invalidation, crash, concurrency, migration and scale gates on supported hosts.
    Verify a clean install from the documented source/package and execute the
    quickstart with a small pinned public cohort. Publish dependency and license
    inventories, including the Git-only correlation dependency and vendored patch
    record. Record minimum/toolchain versions, performance limits, unsupported
    modes and a recoverable project upgrade procedure. Acceptance: an operator can
    install, process, restart, export and diagnose the tutorial without editing
    source or silently relying on an undocumented external program. Only after
    every required gate has executable evidence should the overall goal be complete.


> **Status update (2026-10-15):** The production analysis path is Rust-only: reference indexing and quantification use `kallistors`; DE uses the native Rust `tximport`-style adapter and `rsdeseq2`; preprocessing, downloads, normalization, correlation, and enrichment use Rust crates. R, curl, wget, and bioinformatics CLI tools are not runtime requirements. A compact public yeast RNA-seq smoke run is recorded under `validation/real-data-smoke/`. Some legacy integration fixtures still need conversion from mocked R outputs to native abundance fixtures. The detailed stage checklists below are historical planning notes; use the current code and README for delivered behavior.

## Purpose

This document defines the next development stages for RNAA from the current repository state:

- workspace scaffolding exists and the CLI is usable
- SQLite-backed project state, manifests, and artifact tracking exist
- `init`, `add`, `resolve`, `status`, `doctor`, `download`, `refs`, `quant`, `normalize`, `deseq2`, `corr`, and `run` are implemented
- ENA resolution works for runs, studies, projects, and sample selectors and generates a canonical samplesheet
- the downloader is restart-safe, manifest-driven, and independent of later stages
- restart reconciliation exists for quant, normalize/DE, and corr
- `export`, `survey`, cleanup execution, extraction fallback, crash-recovery sweeps, and a real worker-budget scheduler are still incomplete

The objective is not just to "finish features". The objective is to make RNAA credible as a server-first, restart-safe survey engine that can survive long-running public RNA-seq ingestion workflows without corrupting state or locking us into temporary implementations.

## Guiding Strategy

The next stages should follow five rules:

1. Preserve contract-first design.
Every stage must produce stable internal artifacts and manifests so we can replace shell or R implementations later without rewriting the rest of the system.

2. Finish vertical slices, not disconnected helpers.
Each stage should become independently usable from the CLI before moving to the next one. A command that exists but only bails is worse than no command at all.

3. Make failure states explicit.
The pipeline has to be restart-safe under server conditions. That means durable DB state, verifiable manifests, atomic file writes, and predictable retries before optimization.

4. Keep download independence sacred.
`rnaa download --forever` is the operational backbone. Later work must not couple project progress to a single monolithic orchestrator process.

5. Prioritize combined-cohort correctness over feature breadth.
Multi-study harmonization, design validation, metadata discipline, and artifact provenance matter more than adding more aligners quickly.

## Architecture Pivot: HDF5-First + Cross-Project Reuse

This plan now assumes a storage architecture change:

1. Persist quant outputs as HDF5-first artifacts.
- Treat `abundance.h5` as the canonical quant payload.
- Keep `run_info.json` and stage manifests for provenance.
- Do not persist `abundance.tsv` in steady state. If a tool emits TSV by default, treat it as ephemeral and remove after verification.

2. Persist large matrix outputs in HDF5 internally.
- DE counts, normalized counts, transforms, adjusted matrices, and dense correlation outputs should use HDF5 contracts.
- Human-readable TSV/CSV become export-time products, not canonical storage.

3. Introduce a global shared artifact store outside project roots.
- Add a configurable shared root (for example `RNAA_SHARED_ROOT`) with a content-addressed layout by sha256.
- Project directories become views over shared immutable blobs plus project-local manifests/events.
- Prefer hardlinks or reflinks where available; fall back to copy with checksum verification.

4. Make reuse deterministic with cache keys.
- Quant cache key must include:
  - input FASTQ digests
  - reference identity digest
  - quant engine + version
  - quant parameters
  - preprocessing fingerprint
- If the key already exists in shared storage, attach artifacts to the project without recomputation.

5. Add explicit ownership, pinning, and GC.
- Track which projects reference each blob.
- Allow pinning of references/artifacts used by active projects.
- Add garbage collection that removes unreferenced blobs after retention windows.

6. Keep manifests and DB rows project-local, payloads global.
- Project SQLite remains the control plane.
- Payload immutability and dedupe live in the shared artifact store.
- Every project artifact record should resolve to a global blob id plus a materialized path.

7. Support sample-scoped project composition.
- `rnaa add` must support sample-level selectors in addition to project/study/run IDs.
- Accept explicit sample accessions plus metadata predicates for include/exclude filters.
- Persist selectors as first-class inputs so cohort definitions are reproducible and rerunnable.

8. Make FASTQ preprocessing with `fasterp` a first-class pipeline stage before quant.
- `fasterp` is a trimming/QC preprocessor, not an aligner.
- Keep raw downloaded FASTQ immutable; write trimmed FASTQ as derived artifacts.
- Persist preprocessing reports (`json`, optional `html`) as provenance artifacts.
- Include preprocessing tool/version/flags in quant cache identity so reuse is correctness-safe.
- Treat preprocessing as part of the main run graph, not a side feature. The default orchestration path should understand:
  - `VERIFIED -> PREPROCESS_RUNNING -> PREPROCESS_DONE -> QUANT_RUNNING`
  - direct `VERIFIED -> QUANT_RUNNING` only when preprocessing is disabled
- Scheduler decisions must budget preprocessing together with quant because both are run-parallel CPU stages.

## Current State Snapshot

### Working now

- Workspace layout and crate boundaries are in place.
- `rnaa-core` has config, models, DB, state enums, paths, manifests, traits, and samplesheet logic.
- SQLite migrations exist and the state/event/artifact tables are usable.
- ENA resolver fetches `filereport` JSON and expands run metadata.
- Metadata overrides merge into resolved metadata.
- Canonical samplesheet and editable column map are generated.
- Download worker uses `curl`/`wget`, `.part` files, retries, verification, and manifest output.
- Reference preparation is wired through the CLI, uses cached GENCODE presets, and records reference artifacts.
- `fasterp` preprocessing exists and writes derived artifacts before quantification when enabled.
- Quantification runs in-process through kallistors.
- DE and normalization run in-process through the Rust tximport-style adapter and rsdeseq2.
- Correlation works through Rust residualization plus linked `mincorr`; module-level enrichment via `rsfgsea` is present.
- `rnaa run` overlaps download and quant and reconciles existing downstream outputs.
- Tests exist for resolver parsing, metadata merge behavior, state persistence, download templating, normalize workflow, and normalize -> corr -> module enrichment workflow.

### Incomplete or still weak

- Export and survey are still CLI stubs.
- Cleanup policy exists in config and CLI flags, but real move-to-trash plus purge execution is not complete.
- Extraction fallback (`xsra` and `vdb-native`) is still not implemented.
- Crash-recovery is incomplete: there is no sweep that turns stale `*_RUNNING` states into retryable interrupted states after abrupt termination.
- Reference preparation, quantification, DESeq2, and correlation still need deeper integration tests and harder failure-path coverage.
- There is no real worker-budget scheduler yet. `run` overlaps download and quant, but DE/corr scheduling is still simple and stage budgets are not explicit.
- Preprocessing is implemented, but it is not yet modeled as a first-class scheduler stage in the orchestration plan.
- Shared artifact storage is only partially realized; ownership, pinning, and GC remain unfinished.

## Parallelization Snapshot

The current execution model is mixed:

- Download is parallel across runs. `rnaa download --concurrency N` spawns `N` Rust worker threads, and each worker processes one run at a time.
- Preprocessing and quant can run per-run with bounded worker pools, but the orchestration rules are still simple.
- DESeq2 is single-project and single-invocation from RNAA. That is correct at the scheduler level, but RNAA does not yet explicitly cap BLAS/OpenMP threads.
- Correlation is single-project from RNAA. Residualization and `mincorr` execution still need explicit bounded worker ownership.
- `rnaa run` is a pipeline orchestrator, but not yet a real global worker-budget scheduler.

If you give the current code a 100-worker host today, RNAA can overlap I/O and run-parallel compute, but it still does not enforce a single budget across preprocessing, quant, DESeq2, and correlation. That is the remaining scheduler gap.

## Parallelization Strategy

The correct strategy is not "every stage gets all workers". The correct strategy is a stage-aware worker budget with explicit CPU and I/O ownership.

### Rule set

1. Download owns an I/O pool, not the whole machine.
- Default to a small separate I/O budget, for example `download_workers = 1` or `2`.
- Make it configurable, but do not scale it up automatically with host core count.
- Treat download concurrency as an I/O tuning knob, not a CPU saturation target.
- Allow higher manual override on large servers only when the operator explicitly asks for it.
- Keep per-run file handling sequential for now to preserve simple integrity and manifest semantics.

2. Download and processing must overlap.
- The downloader must keep running while quantification is active.
- Verified runs should enter the quant queue immediately; they should not wait for the full project download set to finish.
- The execution model should be pipelined:
  - `RESOLVED -> DOWNLOADING -> VERIFIED` on the I/O side
  - `VERIFIED -> QUANT_RUNNING -> QUANT_DONE` on the CPU side
- `rnaa run` should observe DB state and schedule newly verified runs while the download worker continues independently.
- If an external `rnaa download --forever` worker is already running, `rnaa run` must cooperate with it, not replace it.

3. Quant should saturate the CPU budget by parallelizing across runs first, then within run second.
- Do not give one kallisto job 100 threads by default.
- Introduce:
  - `execution.total_workers`
  - `execution.reserve_workers`
  - `quant.max_concurrent_runs`
  - `quant.threads_per_run`
- Scheduling rule:
  - fill the CPU budget with as many run-level quant jobs as possible
  - only increase `threads_per_run` when there are too few runnable runs to fill the machine
- For a 100-worker host, the default target should look like:
  - many runnable runs: `12 x 8-thread` or `16 x 6-thread` quant jobs
  - few runnable runs: increase per-run threads gradually, but cap a single kallisto job well below 100
- The goal is 100 busy workers in aggregate, not 100 threads per run.

4. Extraction and preprocessing should behave like quant-adjacent run-parallel work.
- `xsra` or future native extraction should use the same CPU pool family as quant.
- `fasterp` must be treated as part of the main run pipeline, not a side option.
- Default preprocessing concurrency should be bounded and coordinated with quant so the two stages do not oversubscribe the machine.
- Default to a smaller concurrency than quant when compression is involved.

5. DESeq2 should remain single-project and non-parallel at the RNAA scheduler level.
- One DESeq2 analysis per project at a time.
- Do not shard one cohort analysis into parallel RNAA jobs.
- Set safe thread caps explicitly through environment variables for the R process:
  - `OMP_NUM_THREADS`
  - `OPENBLAS_NUM_THREADS`
  - `MKL_NUM_THREADS`
- `VECLIB_MAXIMUM_THREADS`
- `NUMEXPR_NUM_THREADS`
- For a 100-worker host, the default DESeq2 cap should be conservative, for example `8` or `16`, not `100`.
- DESeq2 should start only after the cohort barrier is satisfied:
  - either all known project runs are quantified
  - or the user explicitly freezes a partial cohort snapshot

6. Correlation should be single-project but internally parallel with an explicit ceiling.
- Residualization should be parallelized over genes using Rayon.
- `mincorr` should run in a bounded thread pool, not inherit "all cores" implicitly.
- For a 100-worker host, a safe default is a bounded correlation pool such as `48` to `64` workers unless the user overrides it.
- Correlation should not be parallelized "by runs" because it is a cohort matrix stage, not a run stage.

7. `rnaa run` should become a scheduler, not just a sequencer.
- It should compute a worker budget once.
- It should allocate workers per stage.
- It should publish:
  - active stage
  - workers assigned
  - queued runs
  - effective per-run thread count
- Download must remain independently runnable even after this scheduler exists.
- In the steady state, `run` should look like a conveyor:
  - downloader keeps feeding verified runs
  - quant keeps draining verified runs
  - DE/corr wait for the project-level barrier
  - status reports both pipelines at once

### Concrete 100-worker target model

For a 100-worker Linux host, the default strategy should be:

- reserve `4` workers for OS, SQLite, logging, and control-plane work
- allow download to use a small independent I/O pool, for example `1` or `2`, with manual override if needed
- let preprocessing and quant share the CPU pool as run-parallel stages:
  - default `preprocess.max_concurrent_runs = 6`
  - default `quant.threads_per_run = 8`
  - default `quant.max_concurrent_runs = 12`
  - the scheduler should reduce preprocess concurrency when quant is already saturating the machine
- let quant consume the CPU pool aggressively by run-level parallelism:
  - default `threads_per_run = 8`
  - default `max_concurrent_runs = floor(96 / 8) = 12`
  - if fewer than 12 runs are ready, increase per-run threads up to a cap such as `16`
- allow the downloader and quantifier to run at the same time:
  - example steady state: `2` download workers plus a mix of preprocess and quant jobs, bounded to the remaining CPU budget
  - this is valid because download is mostly I/O bound while preprocess and quant are CPU bound
- run DESeq2 as exactly one project job with a hard thread cap such as `8` or `16`
- run correlation as exactly one project job with an explicit pool such as `48` or `64`

That gives the behavior you want:

- alignment uses the whole machine in aggregate
- preprocessing is part of the main data path instead of an afterthought
- downloads keep progressing while alignment is already processing earlier runs
- DESeq2 does not explode into unsafe parallel jobs
- other stages use bounded, stage-appropriate worker counts
- run-parallel work happens where the contracts are naturally per-run

## Priority Gaps

The next work should focus on the things that still block RNAA from being convincingly production-grade.

1. State model correctness
- add project-level stage state instead of overloading run-level state for DE and correlation
- add startup recovery that converts stale `*_RUNNING` states into retryable interrupted failures
- expand shared artifact-store linkage into full ownership, pinning, and GC semantics (basic blob linkage is now in place)

2. Worker scheduling
- implement the worker-budget model described above
- make preprocessing a first-class scheduled stage in `run`
- add explicit thread caps for DESeq2 and correlation
- publish worker assignments and queue state in progress reporting

3. Honest feature surface
- either implement cleanup, export, survey, and extraction, or clearly leave them out of the public promise
- do not keep CLI flags that imply working behavior if the underlying stage is still a placeholder

4. Testing
- add integration tests for:
  - refs prepare
  - preprocess happy path
  - quant wrapper happy path
  - DESeq2 happy path
  - correlation happy path
  - end-to-end `rnaa run`
- add restart/interruption tests for download, quant, and project-level stages

5. Hardening
- add reference-build locking so concurrent processes cannot trample shared cache creation
- add output guardrails for dense correlation outputs
- verify manifest completeness and failure-mode consistency across stages
- add shared-store GC safety rails: retention windows, pinning, and lock-safe sweeping

## Recommended Delivery Sequence

The smart path is:

1. Implement shared artifact store contracts and schema additions.
2. Add project-level stage state and crash-recovery sweeps.
3. Implement a real worker-budget scheduler for download, preprocessing, quant, DE, and corr.
4. Finish reference preparation and locking against shared storage.
5. Finish per-run preprocessing and quant reuse contracts with stable cache keys.
6. Finish project-level DESeq2 and correlation hardening.
7. Finish export, survey, extraction, and cleanup UX.
8. Harden operational behavior, tests, docs, and release packaging.

This sequence is correct because the remaining risk is no longer basic feature existence. The remaining risk is state-model correctness, scheduler correctness, and operational hardening across long-running server workflows.

## Stage 2: Reference Preparation As A Production-Ready Stage

### Goal

Make `rnaa refs prepare` fully operational for human and mouse, deterministic, cacheable, and auditable.

### Current status

The `EnsemblReferenceManager` already contains most of the logic for:

- Ensembl listing resolution
- download/copy of cDNA and GTF
- `tx2gene.tsv` generation
- `gene_annotation.tsv` generation
- kallisto index building
- reference manifest generation

The missing part is integration, testing, and hardening.

### Tasks

- Wire `refs prepare` in the CLI to `EnsemblReferenceManager`.
- Support explicit flags overriding config:
  - `--organism`
  - `--ensembl`
  - `--gtf`
  - `--cdna`
- Persist selected reference info into the DB or a stable reference manifest lookup table.
- Record reference artifacts in the `artifacts` table.
- Materialize reference payloads into the shared artifact store and link them into project views.
- Add checksum verification for downloaded reference files.
- Parse Ensembl `CHECKSUMS` when feasible instead of trusting only local sha256.
- Handle reference reuse:
  - skip rebuild if manifest and artifacts already match
  - skip redownload/reindex if matching blobs already exist in shared storage
  - rebuild if the index is missing or empty
- Add offline unit tests for:
  - listing parsing
  - custom reference mode
  - GTF-derived tx2gene extraction
- Add a live integration test behind an opt-in flag for current Ensembl.

### Design notes

- Reference identity should be stable and compositional: organism + release for presets, content hash for custom inputs.
- The generated `tx2gene.tsv` and `gene_annotation.tsv` are core contracts; later stages should consume them, not re-derive them repeatedly.
- Avoid baking Ensembl-specific assumptions into downstream stages. Downstream stages should consume `ReferenceBundle`, not "human release N" strings.

### Acceptance criteria

- `rnaa refs prepare --organism human --ensembl latest` works end to end.
- Cached reruns are idempotent.
- Two projects using the same reference reuse the same immutable blobs without duplicate payload bytes.
- Reference manifest and artifacts are sufficient for downstream reproducibility.

## Stage 3: Per-Run Quantification (Historical Implementation Plan)

### Goal

Implement the first real quantification engine with clean Rust/R boundaries and safe per-run state transitions.

### Tasks

- Implement transcript indexing and per-run quantification with the linked `kallistors` crate. (Complete in Rust.)
- Use typed Rust configuration and CLI arguments.
- Add a pre-quant preprocessing contract:
  - optional `fasterp` execution per run
  - deterministic output paths for trimmed FASTQ artifacts
  - preprocessing manifest/report capture
- Validate:
  - one FASTQ for single-end
  - two FASTQs for paired-end
  - kallisto index existence
  - output directory writability
- Add tool policy:
  - `fasterp` enabled/disabled via config and CLI
  - strict mode (fail quant on preprocessing failure) vs bypass mode (warn and continue with raw FASTQ)
- Run `fasterp` on FASTQ inputs before quant when enabled.
- Quantify in-process through `kallistors`; do not spawn a separate kallisto executable.
- Capture:
  - tool versions
  - command arguments
  - reference manifest hash
  - input FASTQ checksums
- Capture preprocessing provenance:
  - pre/post read counts
  - trimmed FASTQ checksums
  - `fasterp` version and flags
- Define and persist a deterministic quant cache key from input digests + reference + engine + params + preprocessing fingerprint.
- Reuse existing quant payloads from shared storage when the cache key matches.
- Define success predicate in Rust:
  - output directory exists
  - `abundance.h5` exists and is non-empty
  - `run_info.json` exists and is parseable
- Record quant artifacts into SQLite with global blob linkage.
- Treat `abundance.tsv` as non-canonical; if generated, purge it after validation.
- Advance state:
  - `VERIFIED -> QUANT_RUNNING -> QUANT_DONE`
  - failure should produce `QUANT_FAILED` without damaging downloaded files
- Implement `rnaa quant` filtering logic:
  - run only on `VERIFIED`, `QUANT_FAILED`, or explicitly forced runs
- Implement per-run log files under `logs/quant/` or a structured equivalent.

### Cleanup substage

After quant works, complete the cleanup behavior:

- `--cleanup none|fastq|sra|all`
- `--cleanup-on success|always|never`
- two-step deletion:
  - move into project trash
  - purge after N days
- record cleanup artifacts and events
- ensure repeated cleanup is idempotent

### Why this stage is next

Quantification is the first place where the download layer proves its value. Until quant works per run, the project cannot produce biologically meaningful outputs, and later DE/correlation work would be forced to invent fake inputs.

### Acceptance criteria

- `rnaa quant` works on verified runs with human or mouse references.
- Quant reruns are safe and deterministic.
- Cross-project quant reuse works when inputs/reference/parameters match.
- When preprocessing is enabled, quant consumes trimmed FASTQ and records reproducible preprocessing provenance.
- Cleanup only triggers after verified quant success.

## Stage 4: Project-Level DE And Normalization (Historical Implementation Plan)

### Goal

Turn multiple per-run quant outputs into a coherent project-level cohort analysis.

### Tasks

- Implement tximport-compatible abundance aggregation and differential expression in Rust with `rsdeseq2`.
- Feed it:
  - samplesheet
  - quant directories with canonical `abundance.h5`
  - `tx2gene.tsv`
  - design formula
  - contrasts
  - transform option
- Validate design variables in Rust before fitting with `rsdeseq2`:
  - verify required columns exist
  - reject formulas referencing absent variables
  - provide actionable errors for missing `condition` or `batch`
- Support contrast specification from CLI and persisted DB state.
- Produce:
  - `gene_counts.h5`
  - `gene_norm_counts.h5`
  - `vst.h5`
  - `vst.rds`
  - `de_<contrast>.h5`
  - `size_factors.tsv`
  - optional PCA and sample distance outputs
- Add manifest with:
  - `sessionInfo()`
  - DESeq2 version
  - tximport version
  - design formula
  - contrast definitions
  - sample inclusion/exclusion decisions
- Decide project state semantics carefully:
  - DE is project-level, not per-run
  - do not incorrectly stamp every run as `DE_DONE` without a project-level contract

### Required schema adjustment

The current DB model is run-centric. DE and correlation are cohort-level stages. Before finishing this stage, add a project-stage tracking mechanism, for example:

- `project_stage_runs`
- or `project_artifacts`
- or `project_states`

This prevents awkward overloading of per-run state for cohort-wide outputs.

### Multi-study strategy

- default combined-cohort DE should include `study_accession` as an available design covariate
- provide `--stratify study_accession` later, but combined analysis should remain first-class
- preserve user-visible metadata provenance so harmonization decisions are inspectable

### Acceptance criteria

- `rnaa deseq2 --design "~ batch + condition" --contrast condition A B` produces all expected artifacts.
- Missing design columns fail early with actionable errors.
- Cohort-level provenance is recorded cleanly.
- Internal DE artifacts are HDF5-first; text tables are produced only through export surfaces.

## Stage 5: Adjusted Correlation In Rust + `mincorr`

### Goal

Produce biologically interpretable, metadata-adjusted correlation outputs without routing the core math back through R by default.

### Tasks

- Implement matrix loading from:
  - VST HDF5
  - normalized counts if requested
- Implement geneset grammar parser:
  - `topvar:N`
  - `mean_gt:X`
  - `genes:file.txt`
  - `regex:...`
  - `biotype:...`
- Implement metadata design parsing for adjustment.
- Build residualization in Rust:
  - construct model matrix from samplesheet columns
  - fit per-gene linear model
  - return residual matrix
  - use `f64` internally
- Persist adjusted matrix artifact in HDF5 and manifest linkage into shared storage.
- Replace the current external-binary assumption with linked `mincorr` usage where practical.
- Define output modes:
  - top-k edges per gene
  - thresholded edges
  - dense matrix for small selections only
- Add guardrails:
  - reject dense full-gene outputs above configured size thresholds
  - report filtered gene counts and dropped genes explicitly
- Record method metadata:
  - correlation method
  - adjustment model
  - geneset selection result
  - output mode

### Technical strategy

Use the linked `mincorr` crate as the compute engine where possible, but keep RNAA’s own `Correlator` trait stable. The internal RNAA contract should remain:

- adjusted matrix in
- correlation artifacts out
- manifests capture method and selection

That preserves the ability to swap to direct streaming/top-k or sparse-native implementations later.

### Acceptance criteria

- `rnaa corr --matrix vst --adjust residualize --model "~ batch + condition" --method pearson --geneset topvar:5000 --out edges:topk=50` works end to end.
- Output size guardrails prevent accidental dense explosions.
- Adjustment and correlation are both reproducible and manifest-backed.

## Stage 6: CLI Completion And Operator UX

### Goal

Turn the pipeline from a set of stage commands into an operable product.

### Tasks

- Implement:
  - `rnaa run`
  - `rnaa export`
  - `rnaa survey`
- Keep `run` decoupled from download ownership:
  - it may trigger downloads optionally
  - it must also work when a separate download worker is already active
- Add `--wait-for verified|quant|de` style orchestration options if needed.
- Add status views for:
  - project-stage progress
  - failed runs
  - pending runs
  - artifact disk usage by class
- Improve `doctor`:
  - verify R packages
  - verify kallisto
  - verify ENA/Ensembl connectivity optionally
  - verify writable root and disk space threshold
- Implement `export` formats:
  - TSV/CSV generated from canonical HDF5 artifacts
  - parquet later if value justifies added complexity
- Implement `survey` output layout:
  - metadata summary
  - counts and normalized matrices
  - DE outputs
  - correlation outputs
  - manifests index

### Acceptance criteria

- A new user can run the documented quickstart without reading source code.
- Operators can tell what the system is doing and why from CLI output plus manifests.

## Stage 7: Extraction Support

### Goal

Support `.sra` fallback paths without polluting the main downloader or quantifier contracts.

### Tasks

- Implement `Extractor` trait concretely for `xsra`.
- Detect `xsra` in `doctor`.
- Detect optional `fasterp` in `doctor` for preprocessing-enabled deployments.
- Route SRA-only runs through extractor before quantification.
- Record SRA input artifacts and extracted FASTQ artifacts separately.
- Add feature-flagged `vdb-native` extractor skeleton behind `ncbi-vdb-sys`.

### Strategy note

Do not attempt native VDB extraction until the `xsra` wrapper path is stable. Native VDB is likely the highest operational complexity area in the whole project and should not block a usable public pipeline.

## Stage 8: Production Hardening

### Goal

Make RNAA believable on long-lived Linux servers under systemd.

### Tasks

- Add structured per-stage logging conventions.
- Add log rotation guidance.
- Add a `systemd` unit for the download worker.
- Add a second unit or timer example for orchestrated analysis.
- Add backpressure rules:
  - quant should not start if disk is critically low
  - reference downloads should not duplicate on concurrent invocations
- Add file locking where needed around shared reference/index creation.
- Add file locking around shared artifact-store writes and GC sweeps.
- Add a `rnaa gc` command for shared-store retention and orphan cleanup.
- Add recovery drills:
  - interrupted download
  - interrupted quant
  - interrupted DE
  - partial manifest write
- Add explicit crash-safe temp naming for every writing stage.

### Acceptance criteria

- Stopping and restarting services does not corrupt state.
- Shared references remain reusable under concurrent server processes.
- Shared quant and matrix blobs remain reusable across projects under concurrent server processes.

## Stage 9: Testing Matrix And Validation

### Goal

Move from "compiles and works on one path" to confidence across the main operational scenarios.

### Test layers

- Unit tests
  - resolver parsing
  - metadata merge
  - design parsing
  - geneset grammar
  - residualization math
  - artifact naming
- Offline integration tests
  - small fixtures for resolve -> download planning -> quant wrapper command construction
  - synthetic DESeq2 input/output validation where possible
  - synthetic correlation fixtures
  - shared-store reuse tests across two projects pointing to the same payloads
- Live integration tests
  - opt-in ENA resolution
  - opt-in Ensembl reference acquisition
  - tiny public run download
- Golden output tests
  - manifest shape
  - samplesheet schema
  - reference bundle schema

### Validation philosophy

For the science-facing stages, correctness tests should focus on invariants and reproducibility, not just file existence. Example:

- DE outputs should have expected contrast labels and matrix dimensions.
- Correlation outputs should have symmetric consistency where applicable.
- Residualization should remove modeled intercept/covariate effects in synthetic fixtures.

## Stage 10: Documentation And Release Readiness

### Goal

Finish the product surface so users can install and operate RNAA without reverse-engineering the source.

### Tasks

- Add a real `README.md` quickstart.
- Add example `rnaa.toml`.
- Add systemd examples under `docs/systemd/`.
- Add fixture-driven tutorial project under `examples/`.
- Document:
  - metadata harmonization workflow
  - design formula rules
  - reference caching rules
  - cleanup/trash semantics
  - manifest schema expectations
- Add changelog and versioning policy.

## Cross-Cutting Refactors To Do Opportunistically

These should happen during the main stages, not as a giant isolated rewrite:

- Introduce project-level stage state separate from run-level state.
- Add stronger artifact uniqueness constraints and reference artifact recording.
- Introduce a shared artifact-store abstraction with content addressing and refcount/pinning metadata.
- Add sample-selector input models (sample accession and metadata predicate filters) with deterministic serialization.
- Add a `Preprocessor` trait and artifact kinds for trimmed FASTQ + preprocessing reports.
- Standardize all stage manifest writing through helper functions.
- Add central error types for user-facing failures vs internal failures.
- Add stable command-runner helpers for external tools.
- Consider replacing ad hoc thread spawning in the download worker with a clearer worker abstraction once more stages are active.

## Priority Order Summary

The immediate next coding order should be:

1. repo hygiene and warning cleanup
2. shared artifact-store schema and storage abstraction
3. sample-scoped add/resolve selectors
4. `refs prepare` CLI wiring and tests against shared storage
5. preprocessing contract (`fasterp`) and trimmed FASTQ artifact model
6. quantification via kallistors with HDF5-first persistence
7. cleanup/trash semantics + shared-store pinning/GC policy
8. project-stage schema for DE/corr
9. rsdeseq2 DE implementation with HDF5-first outputs
10. Rust residualization + `mincorr` integration with HDF5 matrix contracts
11. `run`, `export`, `survey`
12. extraction support
13. operational hardening, docs, and release packaging

## Definition Of "Production-Grade" For RNAA

RNAA should only be called production-grade when all of the following are true:

- every implemented stage is restart-safe
- project-wide provenance is manifest-backed
- the downloader can run independently for days
- multi-study metadata harmonization is explicit and inspectable
- reference reuse is deterministic
- quant, DE, and correlation each have clear success predicates
- deletion is reversible during trash retention
- operator diagnostics are good enough to debug failures without reading Rust source
- integration tests cover at least one realistic happy path from public accession to correlation output

Until those conditions are met, RNAA is a promising foundation, not a finished survey system.
