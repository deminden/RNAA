# Count dependency recovery and generation lineage

`rnaa quant` checks completed counts against current scientific options, the
selected source/reference and verified artifacts. Before restaging or fitting,
RNAA looks for a verified historical native count result under the current source,
processing and QC contracts. A changed request with no reusable result can require
reads that retention deliberately removed; RNAA plans and stages those dependencies
before producing replacement counts.

## Recovery choices

The planner chooses the earliest usable input set:

1. Use compatible retained trimmed reads when preprocessing is enabled and current
   QC evidence admits them. Original reads can remain pruned.
2. Use a complete original FASTQ group, preprocessing it when requested.
3. Restage the complete selected FASTQ source through the existing downloader when
   the required inputs are missing, incomplete, invalid or source selection changed.
4. Block with an actionable error when no supported complete source is available
   or verified QC evidence rejects the run.

Historical selection can avoid this reconstruction entirely. It supports restored
scientific settings, retained-trimmed requests, local and managed sources, and
preprocessing on/off round trips. Current input locations are recorded in immutable
[use-time input bindings](scientific-cache.md), while original result files stay
unchanged. Competing commands still use one owned selection attempt.

Preprocessing history can also restore a healthy report for identical source content,
even when its reads are pruned. A new count fit still needs actual processed reads;
an evidence-only preprocessing selection cannot replace their reconstruction.

## QC report reconstruction

`rnaa quant` also plans recovery when the counts still verify but current QC
evidence is missing or damaged and no healthy compatible historical report exists.
`rnaa qc reevaluate` evaluates evidence; it does not download reads or run an engine.

Report reconstruction requires a complete original FASTQ group. Retained originals
can be processed locally. Otherwise the selected supported FASTQ source is staged
through the normal downloader, with its publisher checksum and source/version
checks. Retained trimmed reads cannot recover rejected-read or original-read QC
statistics and are not substituted for originals. A local source with pruned
originals and no downloadable replacement requires restoring/registering those
originals. A blocked report-only plan preserves completed count progress and records
the missing dependency rather than declaring the count payload corrupt.

The operational recovery request (revision 2) records `purpose: "qc_report"` and
the chosen action, prior count generation, target policy and source observations.
The plan is rechecked under run exclusion. Competing commands share the download
and reconstructed report through the usual owned attempts and publication fences.
Original source content and catalog/configuration observations are rechecked before
both fresh preprocessing publication and reuse; allocation-only thread changes do
not change the scientific request.

The new report publishes in a fresh immutable preprocessing generation. Prior
damaged evidence is preserved for diagnosis. Its measured values are evaluated
under the current QC gate, including when engine-error bypass is enabled. Rejected
reads remain ineligible. If the recovered read/processing identities match a healthy
old count request, count selection retains the exact old payload/seal and records
new input bindings; no count engine execution is required. A changed unchecksummed
origin can produce different reads, requiring a new count result after successful
QC. Failed recovery never replaces the earlier count head.

## Input staging

Groups are validated completely, including paired lanes and FASTQ records.
The existence of one file does not establish a usable input set. Disabling
preprocessing requires original reads; retained trimmed reads are not substituted
as original input. Source selection stays independent of run progress and keeps
archive files out of a selected FASTQ group.

Recovery uses the current selected remote descriptors, authoritative checksums
when supplied, HTTP version/resume checks and normal download publication.
It does not invent a source for a locally supplied file or extract an archive.
Restore/register local inputs or resolve a complete downloadable FASTQ source
when required. Initial unresolved downloads remain under `download`/`run`;
`quant` restages stale previous count dependencies rather than downloading every
unrelated newly resolved accession.

The recovery event records the old count generation, desired quantification-policy
fingerprint, source-selection fingerprint, available input snapshots, chosen
action and operational plan digest. The actual count request binds the restored
input bytes after staging. Run state and committed history allow a later command
to reconstruct a pending recovery with current settings after interruption.

## Publication and concurrency

Old generations and their data stay available while replacements are prepared.
A failed transfer or computation cannot replace a complete count head. Recovery
uses the existing run claims, leases, source checks and atomic publishers; two
commands cooperate on one replacement producer. A follower waits for a known
download recovery and cannot announce success while a previous-count dependency
remains unfinished. `run --no-download` reports when required recovery conflicts
with that option.
After historical lookup releases run exclusion, RNAA rechecks the result and
rebuilds any recovery plan under the run lock before publishing it. A competing
command's restored inputs cannot be overwritten by an obsolete download plan.

Count publication rechecks current scientific settings, source selection, layout
and QC admission against the catalog and saved configuration. Thread allocations
remain execution settings. DE/correlation also observe upstream quantification
policy while running. Before consuming the selected DE matrix, correlation
revalidates its current count dependencies after checking its reference binding.
Changing global quantification settings cannot silently reuse an old native count
producer or publish a model fitted under a superseded count policy.

QC admission remains separate from completed count validity. Corrupt or missing
QC evidence blocks analysis and is reported; it does not turn intact count payloads
into corrupt data. Restaging failures and unavailable sources return failing command
statuses while preserving old outputs and durable diagnostics.

## Sealed generation lineage

Migration 0014 indexes exact files owned by committed generations and stores
consumer-to-producer edges. Newly published generation seals bind input roles,
paths, SHA-256/size identities, producer generation IDs and canonical producer-seal
fingerprints. File indexes, edges, head, registry and attempt completion commit
together. Parent catalog changes reject publication; corrupt edge/file indexes or
producer seals reject verification. Managed sources must identify committed
producers through their canonical paths, including when entered through aliases.
Correlation pins a supplied matrix's canonical path before constructing its request.

Preprocessing records its original-read producer. Quantification records available
read producers and the reference index. Fresh DE/correlation seals record their
actual input witnesses, including count, annotation, QC and matrix producers.
Unmanaged inputs retain explicit external witnesses. Internal frozen metadata
belongs to the consumer's own inventory.

Existing seals retain their bytes. The migration backfills file ownership for
committed legacy generations, allowing new consumers to identify those producers.
Old seals without recorded edges report `recorded: false`; their full ancestry is
not reconstructed or claimed complete. For these legacy entries, `external_inputs`
can include witnesses whose producer edges were never stored.

Inspect current or historical upstream evidence:

```sh
rnaa lineage --group quant --subject RUN1
rnaa lineage --group quant --subject RUN1 --history
rnaa lineage --group quant --subject RUN1 --history --verify
rnaa lineage --group de --all-requests --verify
```

The command prints JSON. Without `--verify`, payload files are not rehashed and
`payload_verified` is null. `--all-requests` includes committed request variants
for the selected output subjects; `selected` distinguishes current heads.
Verification can read large files and honors exact
pruning receipts. A failed payload check prints its result and exits unsuccessfully.
Catalog/seal faults stop inspection rather than printing a falsely valid lineage.
Edges describe immutable production history; current requests/settings determine
whether that history can satisfy a new computation.

The JSON separates the original `lineage` from the selected `input_binding` and its
optional verification result. Upstream inspection follows both edge sets. A binding
can authorize retention of newer equivalent input files through its own immutable
receipt, without adding them to an old seal. Only FASTQ absence can be authorized;
altered files, reports and count outputs cannot be hidden by a pruning receipt.
Future GC must retain binding owners and pruning evidence, including historical
bindings needed to explain deletion, in addition to original generation edges.

## Validation and remaining scope

Offline native tests use local HTTP sources and real preprocessing/quantification
to cover all-pruned recovery, retained-trim reuse, preprocessing disablement,
failed restaging/retry, a durable pending state, competing commands and missing
local sources. They verify old generations, exact parent edges, transitive history
and JSON verification failures. Core tests cover index/edge/producer-seal
corruption, rollback, uncommitted/aliased sources and legacy ownership migration.

Analysis stages support [verified historical request selection](scientific-cache.md).
Download and reference history also support verified selection. Restoring a
checksum-identified source can select its old download, preprocessing and counts
without new transfers or engine work. Unchecksummed older download snapshots need
current strong origin evidence; pruned files still need reconstruction for a new fit.
Shared-blob ownership/pins and reversible orphan collection are implemented in the
[shared store](shared-storage.md). This does not close the entire dependency/storage
roadmap. [Project-generation maintenance pins and retirement](generation-cleanup.md)
preserve binding/pruning proofs. Complete legacy
ancestry, immutable external documents and generalized invalidation remain open.
Report reconstruction after pruning can still require explicit operator recovery
when original local inputs have no recoverable source. Missing/damaged report
tests cover complete pruning, retained originals, incomplete evidence with only
trimmed reads, failed transfer/retry, changed source content, QC rejection and
competing commands. Full recovery-specific process-kill, scale and power-loss gates remain
separate work; existing earlier-stage crash tests do not establish those guarantees.
