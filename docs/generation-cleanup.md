# Project generation retirement

Project generation cleanup is separate from [shared-blob collection](shared-storage.md).
It removes unreachable committed output directories from the active store while
retaining their original catalog metadata, scientific requests and provenance.
The default operation is a preview; collection never replaces the current result.

```sh
rnaa storage generations status
rnaa storage generations gc
rnaa storage generations gc --apply
rnaa storage generations keep GENERATION_ID --label publication
rnaa storage generations unkeep GENERATION_ID
rnaa storage generations restore GENERATION_ID
rnaa storage generations gc --apply --purge-expired
```

These commands use the project selected by `--root`; they do not require a shared
store. `--shared-root` belongs to shared-blob commands and is rejected for project
generation operations. Commands print JSON. A busy project or invalid cleanup
candidate produces a failing command status after the applicable report is printed.

## What remains reachable

Migration 0017 adds project maintenance observations, durable keep labels and
retirement receipts. Cleanup roots include:

- Every selected generation head, including heads for retained inactive runs.
- The selected project reference and active artifact views, including aliases.
- Current QC evidence and explicit `keep` labels.
- Quantification witnesses needed by retained-input records.

Traversal follows original producer edges and every recorded input binding for a
kept consumer, including earlier bindings that are no longer selected. A retained
producer with intentionally pruned reads also needs the exact count consumer or
binding that authorized their absence. Those proof links are followed in the
opposite direction from ordinary producer dependencies. Removing an unused fit
therefore cannot invalidate another result's deletion evidence.

Catalog/seal/index and binding identities are checked before cleanup. An invalid
unreachable payload stays in place and keeps its own required dependencies. An
active alias protects its canonical generation; dangling aliases cannot reactivate
retired artifacts. Unexpected container entries, symlinks, altered files and
ambiguous move destinations stop the corresponding operation.

Current aliases of recorded source inputs also protect their present managed
targets through maintenance-only graph edges. Retargeting an external source
alias therefore cannot make its new managed target collectible. These protection
edges do not rewrite the result's immutable historical ancestry.

Older seals without recorded ancestry, or retained-input records without a sealed
consumer witness, make the project's generation history conservatively held.
Cleanup does not infer unknown historical relationships.
`rnaa storage generations import-lineage` previews reconstruction from supported
checksummed sealed evidence; `--apply` publishes separate verified receipts.
Imported producer edges participate in collection and restore closure while the
original seals remain unchanged. Unsupported histories remain held; see
[legacy ancestry](legacy-ancestry.md).

## Active work and pinned paths

Workers, reconciliation scopes, prepared reuse tokens and generation path snapshots
hold shared project maintenance exclusion. Shared-blob publication also holds a
project pin while reading/copying its source and committing ownership. Database
writes acquire the same boundary before opening their SQLite write transaction.
Collection requires exclusive acquisition and holds it through graph verification,
catalog marking and filesystem operations. A live paused reader or worker cannot
be overridden by lease expiry.

`ProjectPaths::with_generation_heads`, `for_quant_run`, `for_reference`,
`GenerationHead::paths` and `Generation::paths` return pinned snapshots. Retain
the returned `ProjectPaths` for the duration of file consumption; cloning it retains
the same pin. Turning a temporary snapshot into a bare `PathBuf` does not retain
its guard. Library callers holding other raw paths can explicitly retain a
`GenerationStorePin`.

Pins have inspectable observation rows. Insertion is best-effort when a containing
publication already owns SQLite's writer, while OS exclusion remains authoritative.
Observation rows can also outlive a crash. Applied cleanup can clear stale rows
only after acquiring exclusive exclusion. The status list is therefore not a live
worker count. Coordination files under `.rnaa-locks/` are never unlinked.

This is currently project-wide maintenance, rather than per-generation concurrent
collection. A busy collector reports `busy: true`. Other maintenance-bound database
operations/path acquisitions use a bounded wait and then report that a retry is
needed. Large previews can read and hash historical payloads.

## Trash, restore and purge

`gc --apply` marks eligible generations transactionally and moves their complete
directories from `.rnaa_generations/ID` to `trash/generations/ID`. Receipts bind the
original sealed inventory, artifact registry, physical seal witness and exact
observed payload files. Original producer metadata and seal bytes are unchanged.

Retired generations are excluded from historical request lookup and cannot expose
normal pinned paths, become active artifacts/QC/reference evidence or satisfy a
prepared reuse request. `storage generations status` retains their identities and
transition timestamps. The original committed metadata remains available to core
provenance APIs after physical purge.

New trash uses `quant.trash_days`, or `gc --trash-days N`. Each deadline is recorded
at retirement; later options cannot shorten it. Only trash present at the start of
a GC operation can be purged by that invocation. `--purge-expired` is explicit and
still requires `--apply` for deletion. Expired producers remain protected while
another live or restorable generation needs them. Restore is possible after expiry
while the payload still survives and purge has not started.

Restore computes the target's complete dependency/proof closure and verifies all
required trash before changing paths. It marks the closure as one restore operation,
moves the directories back, verifies the restored generations and then makes the
whole restored group available for reuse. It does not change selected heads; repeat
the scientific command to select a restored request. A required purged dependency
prevents restoration instead of admitting a broken result.

Trash and restore intents precede directory renames; purge intent precedes recursive
removal. Interrupted operations resume from their durable state. Partial purge
permits already removed files while requiring every remaining file to match the
recorded inventory. The collector never treats unexpected remaining bytes as
authorized deletion. A completed `PURGED` receipt can release matching ownership in
the shared-blob collector, after receipt and active-reference checks. Project trash
continues to protect its shared copies until that point.

## Verification and remaining scope

Core tests cover transitive producers, keep labels, snapshot clones, prepared tokens,
workers, alias views, damaged outputs/catalogs, legacy holds, recorded deadlines,
restorable dependencies and original seal preservation. Real process-kill tests
cover readers and collectors before/after moves, interrupted restore and partial
purge. Native tests cover preview/keep/trash/restore/reuse, shared-owner release only
after project purge, and surviving input bindings/pruning proofs after another fit
is retired.

This collector targets validated committed generations. A separate
[unpublished cleanup](unpublished-cleanup.md) command handles tracked diagnostic
containers from terminal attempts without publishing them as scientific results.
Unindexed directories and other temporary copies remain diagnostic data. FASTQ/SRA-specific
`CleanupMode`/`TrashCleaner` behavior uses the separate [input cleanup](input-cleanup.md)
workflow. Generation GC does not automatically execute those flags. Restorable file
trash pins its producer/consumer generations; purged input receipts preserve their
count-proof witnesses. [Observed storage usage](storage-usage.md) reports categories
and deduplicated allocated inodes without authorizing collection. Complete legacy
shared ownership migration, unsealed result adoption, resource reservations,
cross-project scientific reuse, scale measurements and full power-loss gates
remain open.
