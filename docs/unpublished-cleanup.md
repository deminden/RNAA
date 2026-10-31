# Reversible cleanup of unpublished generations

Failed or interrupted jobs can leave useful partial output under
`.rnaa_generations/GENERATION_ID/`. These directories are separate from committed
scientific history and from the downloader's retry caches.

```sh
rnaa storage generations unpublished gc
rnaa storage generations unpublished gc --apply
rnaa storage generations unpublished status
rnaa storage generations unpublished restore GENERATION_ID
rnaa storage generations unpublished gc --apply --purge-expired
```

The default operation previews and prints JSON. Commands use the project selected
by `--root`, without requiring a shared store. Live workers and pinned readers
exclude cleanup; a busy command prints `busy: true` and returns failure.

## Eligibility and protection

Cleanup requires an existing generation catalog row in `BUILDING`, `SEALED` or
`ABANDONED`, with its original attempt, project, stage and permanent payload path.
The owning attempt must be terminal. An expired timestamp alone cannot override a
live worker or make a still-running attempt collectible. Normal job recovery must
first recover an expired interrupted attempt; opening the project alone does not
grant deletion authority.

Committed generations are excluded. Published heads, reference selections, keeps,
file indexes, ancestry imports, producer edges, input bindings and pruning proofs
cannot refer to an eligible unpublished generation. The collector checks the
committed graph, registered artifact paths, retained input paths, original and
rebound ancestry, and both persisted and on-disk project settings. Aliases are
resolved, including dangling aliases whose managed target can be determined.
References to either the original container or its trash location hold cleanup.
Unproven sealed ancestry or active unsealed legacy artifacts conservatively hold
unpublished collection until their own ancestry is resolved.

Unindexed directories have no provable attempt ownership. They appear as
`held_unindexed` and remain in place. Unexpected symlinks, special files, altered
receipts, changed catalog identities or ambiguous source/target paths cannot
authorize a move or deletion. A held reference is an ordinary report outcome;
invalid candidates and pending recovery errors also produce a failing exit status.

## Durable trash, diagnostic restore and purge

Migration 0020 records a separate cleanup receipt. It binds the original catalog
identity, original unpublished state, every file's relative path, length and
SHA-256, empty directories, and the recorded trash deadline. A second immutable
witness lives under `.rnaa_support/unpublished_trash/OPERATION_ID.json`.

Apply flushes the observed files, writes the witness and transactionally marks all
new candidates `ABANDONED` with `TRASHING` intents and events. It then renames whole
containers into `trash/unpublished/GENERATION_ID/` and flushes the parent directories.
An interrupted transaction leaves no published cleanup receipts; diagnostic witness
files left before commit are never adopted merely by filename.

`restore` moves the complete verified diagnostic directory back. It preserves its
file and seal bytes, keeps the generation `ABANDONED`, and changes no scientific
request, artifact registration, result head or completed-stage state. Restoring a
failed output does not make it a successful or reusable result. A later scientific
retry publishes its own generation.

New trash uses `quant.trash_days`, or `gc --trash-days N`. The recorded deadline
cannot be shortened by later invocations. Purge requires both `--apply` and
`--purge-expired`; even zero-day trash survives the invocation that first selects
it. Only pre-existing expired trash can be purged, and current references are
checked again. Before partial-purge recovery, every surviving file and directory
must match the receipt; missing entries are allowed, new or modified entries are
held. Receipts and original catalog metadata survive completed purge.

The collector holds exclusive project maintenance throughout each operation.
Intent recovery is performed by applied cleanup and restore commands. Real
process-kill tests exercise unfinished builders, before/after trash moves,
before/after restore moves and interrupted partial purge.

## Accounting and remaining scope

`storage usage` classifies this data as `unpublished_trash` after retirement.
Cleanup reports logical file bytes, including the seal when present. This is not
an estimate of freed space: external hardlinks and copy-on-write sharing can retain
physical storage after a path is purged. Whole-container inventories intentionally
include unsealed files and empty directories, preserving their diagnostic contents.

This workflow handles tracked unpublished generation containers. Unindexed data,
private transfer caches, other support/temporary files and unsealed legacy outputs
still require their own ownership and recovery protocols.
[Shared legacy ownership migration](shared-legacy-import.md) has an explicit
inventory/import workflow; unresolved cases stay held. Admission reservations
remain separate roadmap work. See
[generation cleanup](generation-cleanup.md), [legacy ancestry](legacy-ancestry.md)
and [storage usage](storage-usage.md) for the neighboring contracts.
