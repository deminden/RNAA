# Importing ancestry for legacy sealed results

Some older immutable generations predate the producer-edge catalog. Their files
and manifests can still verify, but their missing ancestry prevents dependency
cleanup from proving which history may be collected.

```sh
rnaa storage generations import-lineage
rnaa storage generations import-lineage --apply
rnaa lineage --history --verify
rnaa storage generations gc
```

The import command prints JSON and previews by default. Project opening performs
normal schema upgrades; a preview writes no ancestry receipts or result payloads.
Both preview and apply require project-wide idle maintenance, so live workers and
pinned readers return `busy` rather than racing the import.

## Evidence and publication

Migration 0019 adds a separate ancestry catalog. Import never rewrites an original
seal, manifest, scientific request, result payload, selected head or active artifact
registration. Native generations with recorded ancestry remain unchanged.

For a legacy generation, the importer verifies its existing payload/seal and reads
every registered manifest JSON from that sealed inventory. It reconstructs local
input identities from the original sealed input list, supported stage manifests
and full request descriptors with checksummed physical input witnesses. Registered
manifest coverage, checksums, stages, subjects and outputs are validated. Conflicting
declarations, unsupported schemas and missing checksums hold that generation.

Count manifests must also supply a checksummed reference-index witness, including
the legacy `index`/`index_sha256` fields when present. Download origin URLs and
publisher claims stay in their immutable manifest evidence; ancestry import
performs no network revalidation and grants no new current-origin authority.
Local reference source descriptors must agree with their sealed materialized
transcript/GTF snapshots. Their original source paths join ancestry, protecting
managed producers used to build a custom reference. Missing source-origin evidence
cannot silently turn a materialized reference copy into a source-free result.

The exact consumed input list stays separate from broader ancestry observations.
For example, quantification provenance can name the older raw reads used to create
a reused preprocessing result. Those older producers remain protected, while the
consumed input list continues to describe the selected fit. Adding provenance
must not increase the lane/read set used to admit or reuse a count result.

Managed producer edges require an exact indexed path, SHA-256, length and producer
seal. Missing or ambiguous managed ownership cannot silently become an external
source. Consumer-owned frozen files must appear in the consumer's sealed inventory.
External source witnesses remain historical observations, independently of whether
those mutable sources still match current workflow settings.

Before publication, the complete candidate graph is checked for cycles. New
receipts bind the original catalog inventory, artifact registry, typed seal,
all manifest bytes, consumed inputs and ancestry edges. A second witness is saved
under `.rnaa_support/lineage_imports/IMPORT_ID.json`. One transaction publishes all
verified candidates, their indexed edges and import events. A failed transaction
publishes none of those candidates; a crash before commit can leave diagnostic
witness files, which a retry never adopts merely by filename.

## Inspection, reuse and cleanup

Each report entry identifies its generation, group, subject, action, dependency
count and any held reason. Actions distinguish native ancestry, a proposed import,
an applied import, a previously imported generation and a hold. `complete` describes
ancestry coverage for the inspected catalog; use `lineage --verify` to inspect
current payload integrity as well.

An import can publish verifiable candidates while reporting held generations.
The command returns failure when any generation is held; the JSON records which
candidates were applied. Repeating the command is idempotent. A held legacy
generation continues to trigger conservative history protection.

`lineage` reports `source: native`, `imported` or `unresolved_legacy`. Imported
dependencies participate in downstream-history lookup, generation GC and restore
closure. Their immutable receipts continue to verify after supported trash/purge
operations. Count retention, reversible input trash and recovery read the imported
consumed-input view without editing the producer's original history.

## Limits and verification

This command targets **sealed committed generations**. Unsealed legacy outputs,
old shared-store ownership and unknown files remain in place. Their adoption and
cleanup need their own evidence and migration protocol. Weak historical correlation
manifests that omit adjustment/module inputs remain held unless a complete sealed
input inventory or supported descriptor supplies that evidence. Restore damaged
payloads or recompute unsupported results instead of inventing ancestry.

Core tests cover unchanged seal bytes, exact/extra provenance edges, consumed-input
preservation, cycles, unknown inputs, corruption, pins, transaction rollback,
retirement/restore and real killed processes before/after import publication.
Native CLI tests cover legacy reference/download/preprocessing/count manifests
with and without original sealed input lists, retained count reuse, FASTQ trash and
offline restore, normalization and correlation history.

[Shared legacy ownership import](shared-legacy-import.md) has its own explicit
inventory and read-only project protocol. Unresolved legacy ownership, unsealed
result adoption, unindexed/transfer/support cleanup and measured release/power-loss
gates retain their separate requirements.
[Tracked unpublished cleanup](unpublished-cleanup.md) handles terminal attempts'
diagnostic containers separately.
