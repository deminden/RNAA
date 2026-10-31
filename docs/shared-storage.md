# Shared-blob ownership and collection

Set `storage.shared_root` in `rnaa.toml` or `RNAA_SHARED_ROOT` to enable a shared
store. Configured relative paths are resolved against the project root, including
when `--root` is used from another directory. A `storage --shared-root PATH`
override is interpreted from the invoking directory and can inspect a store
without opening a project database.

Native quantification copies its HDF5 output into content-addressed storage under
`blobs/sha256/AA/BB/DIGEST`. New shared copies own separate inodes from project outputs.
Project generations retain their own complete outputs; the shared store currently
does not replace the scientific count TSV or provide cross-project request reuse.

## Publication handoff

`catalog.sqlite3` is a store-wide SQLite catalog, separate from each project's
database. It records expected blob sizes/states, project database identities,
artifact paths, generation ownership and pins.

Opening the store takes a bounded OS initialization lock before opening SQLite,
selecting WAL mode and applying transactional schema upgrades. This covers the
journal-mode race that SQLite's ordinary busy timeout does not reliably serialize.
The lock is released before blob publication or collection locks; process exit
releases it automatically. A paused initializer can cause a clear busy/retry error
after 15 seconds rather than another opener racing journal setup.

New publication works in this order:

1. Acquire exclusive blob exclusion and persist a `PREPARED` ownership reservation.
2. Copy the payload atomically and verify its actual SHA-256 and size.
3. Commit the project generation, artifact registry, shared reference and successful
   attempt through the normal project transaction.
4. Verify that project publication and finalize store ownership as `COMMITTED`.

This is a recoverable handoff between two databases, rather than a claim that their
commits are one distributed transaction. If the process dies before finalization,
the durable reservation remains. Collection obtains blob exclusion and reads a
consistent, read-only project snapshot. A committed result retains its owner;
an aborted unreferenced publication can become an orphan. Missing, replaced,
unreadable or corrupt project metadata keeps the last known owner protected and
appears in diagnostics. A committed generation whose ownership witness changes
cannot be retired by treating that fault as absence.

Historical generation registrations retain ownership after their active registry
entries are superseded. Managed sources record their generation before copying,
so collection does not need a full history scan for normal committed owners.
Legacy active references remain supported. Existing payloads first encountered
without catalog ownership get a conservative `legacy_hold`; unindexed files are
never direct collection candidates. [Legacy ownership import](shared-legacy-import.md)
can catalog verified old references and release holds with explicit complete
project inventory and accounted hardlinks. The collector cannot infer all older
projects from one project's database.

## Pins and concurrent work

A publication holds exclusive OS exclusion until after the handoff. The library's
`SharedStore::pin` returns a reader that verifies the payload and holds a shared
lock plus an inspectable durable pin until it is dropped. Multiple readers can
coexist. Collection skips busy blobs, including a paused writer or reader; a
heartbeat timeout alone cannot override that exclusion. The implementation uses
Rust's [file locking API](https://doc.rust-lang.org/std/fs/struct.File.html#method.try_lock_shared).

A killed process releases its OS exclusion. Durable reader rows can be retired
only after collection acquires exclusive blob exclusion. Lock files under `.locks/`
are coordination files and are never unlinked. Reader pin rows may therefore
remain visible after a crash until an applied collection reconciles them.

## Operator commands

All storage commands print JSON. Within a configured project:

```sh
rnaa storage status
rnaa storage gc
rnaa storage gc --apply
rnaa storage restore sha256:DIGEST
rnaa storage gc --apply --purge-expired
```

For explicit store access, use `rnaa storage --shared-root /path/to/store status`
and the equivalent `gc` or `restore` subcommand. A blob's bare lowercase SHA-256
also works with `restore`.

`status` reports catalog bytes, state, owner count, prepared owner count, durable
reader count, legacy protection, migration receipt ID and trash timestamps/deadline. These are catalog
quantities, not measured total filesystem usage. Unindexed payloads and abandoned
temporary copies are preserved and are not included as reclaimable bytes.

`gc` defaults to a dry run. It checks ownership and reports the intended action
without reconciling stored owner/pin rows or moving/removing payloads. `--apply`
reconciles ownership and moves verified unowned blobs to `trash/DIGEST.blob`.
Results distinguish `held_owned`, `held_busy`, `held_legacy`, `held_external_links`, `held_error`, `trash`,
`purge`, pending recovery and empty aborted publications. Ownership diagnostics
can keep a blob held even when its project is unavailable. Encountered per-blob
payload verification or store-transition errors retain the file, appear in the
report's errors and cause the CLI to fail after printing that report.

The trash deadline is recorded when a blob is moved. New trash uses
`quant.trash_days` (seven days by default), or `gc --trash-days N`. Later settings
cannot shorten a blob's already recorded deadline. `--purge-expired` is required
to remove expired trash and still defaults to preview unless `--apply` is present.
Restoring verifies the retained payload and moves it back into the live store.

Trash, restore and purge record their intent before filesystem changes. They
flush affected directories and finish their catalog state after moving/removing
the file. An applied retry resumes the same pending operation; ambiguous copies
or unexpected bytes stay protected for diagnosis. A new publication can restore
a matching trashed blob before recording ownership.

## Validation and remaining scope

Core tests cover two-project ownership, orphan handoff, dry-run immutability,
writer/reader exclusion, stale pins, interrupted move/restore intent, retained
deadlines, purge, corrupt payloads, unavailable projects, legacy holds, symlink
boundaries and eight simultaneous publishers. Real process-kill tests cover
uncommitted/committed handoff and reader pin recovery. Native CLI tests cover
generation ownership after head replacement, failed project registration/retry,
concurrent projects, corrupt ownership witnesses, project-relative configuration,
preview, trash and restore.

This collector handles shared blob copies. [Project generation retirement](generation-cleanup.md)
now provides dependency/binding/pruning-proof reachability, maintenance pins,
reversible generation trash and purge. The shared collector releases generation
ownership only after a matching completed purge receipt and active-reference checks;
project trash retains ownership. FASTQ/SRA-specific `TrashCleaner` behavior is
implemented in the separate [input cleanup](input-cleanup.md) workflow. Shared schema 2
adds [legacy ownership import](shared-legacy-import.md), with read-only project
inspection, immutable migration witnesses and link/ownership guards on later GC.
Tracked unpublished diagnostics have their own [cleanup](unpublished-cleanup.md).
Unresolved or unavailable legacy ownership remains held. Resource reservations,
cross-project scientific request reuse, unsealed result adoption, unindexed/transfer/
support cleanup, measured scale and full power-loss gates remain open.
