# Reversible FASTQ and SRA cleanup

`rnaa storage files` manages individual input files after a native count result has
committed. It is separate from shared-blob and whole-generation collection.

```sh
rnaa storage files status
rnaa storage files clean --mode fastq
rnaa storage files clean --mode all --run RUN1 --apply --trash-days 7
rnaa storage files restore TRASH_ID
rnaa storage files gc --purge-expired
rnaa storage files gc --apply --purge-expired
```

Every command prints JSON. `clean` and `gc` preview by default. `--apply` completes
previously authorized transitions before selecting new work. An explicit manual
command controls its own mode and timing; it does not enable automatic cleanup.
Without `--mode`, `clean` uses `quant.cleanup`. `gc` never selects new input files.
`restore` restores the complete run operation containing the requested trash ID.

## Automatic timing and retention

Configure `quant.cleanup = "none" | "fastq" | "sra" | "all"`,
`quant.cleanup_on = "success" | "always" | "never"`, and `quant.trash_days`.
Both `quant` and `run` accept `--cleanup`, `--cleanup-on` and `--trash-days`.

- `none` or `never` disables automatic file trash.
- `success` cleans completed, admitted count results at idle scheduler boundaries,
  and after a successful command. A failed command does not start an additional
  success-only cleanup pass; earlier successful runs may already have been cleaned.
- `always` additionally attempts cleanup after a failed command. It still requires
  a healthy committed count result and exact input proofs. The processing error
  remains the command's primary outcome if cleanup also fails.

Idle cleanup defers while another worker or pinned reader uses the project and
retries at a later scheduler boundary. An idle `download --forever` command releases
its reader pins between attempts, allowing another command to clean completed
runs while that downloader stays alive. A final busy/held cleanup reports an error
without changing count completion. Automatic cleanup never purges expired trash;
permanent removal requires the explicit `gc --apply --purge-expired` command.

When automatic `fastq` or `all` cleanup is enabled, reversible cleanup takes
precedence over `fastq_retention`: selected FASTQs go to trash instead of being
permanently unlinked first. With cleanup disabled, `never`, or `sra` mode, the
existing FASTQ retention policy continues to apply. Whole-generation GC remains
explicit and has its own reachability rules.

## Ownership and count evidence

Cleanup accepts only exclusively registered regular files inside the canonical
project root. External files, hardlink aliases, multiply registered paths, files
under internal trash/support/transfer directories, changed bytes and ambiguous
native producer ownership are held with JSON diagnostics. One held member holds
the entire selected run group. Each selected original/trimmed family must match
the complete consumed count input group, so a missing registry member also holds
its mates and lanes.

For a FASTQ, a receipt binds the exact kind, path, SHA-256 and length to a committed
native count generation and its original input snapshots or immutable use-time
binding. Managed FASTQs also require their download/preprocessing producer and its
sealed artifact registration. The count and producer payloads must verify before
the move is marked. SRA cleanup requires an unused project-owned archive beside
verified native counts for that run. It does not extract SRA reads or establish
support for quantifying an archive directly.

Migration 0018 stores immutable descriptors with a digest in `file_trash`. A second
physical witness lives in `.rnaa_support/file_trash_receipts/TRASH_ID.json`.
The descriptor records producer/consumer seals, any input-binding digest, original
registration, run layout, source selection and the original retention deadline.
Changed catalog/witness data stops cleanup rather than authorizing deletion.

## Moves, restore and interruption

Before moving any file, one database transaction records the complete selected
run group in `MOVING` state and removes its live artifact registrations. Its bytes
move to `trash/files/TRASH_ID/data`. Same-device moves create an exclusive second
name, sync it, then unlink the source. Cross-device moves use a reserved private
copy, verify and sync it, publish without replacing a destination, then unlink the
verified source. Original and destination directories are synced.

States are `MOVING → TRASHED → RESTORING → RESTORED`, or
`TRASHED → PURGING → PURGED`. Durable intents precede physical changes. A retry
finishes pending moves, restores or purges, including a complete copy whose source
was not yet unlinked. Retrying `restore` also succeeds when the preceding recovery
has already completed that operation; it does not reselect a later-stale input set.
A reserved `.rnaa-file-copy-TRASH_ID` may be restarted after
interruption; unrelated files, symlinks, hardlink aliases and conflicting payloads
are preserved and reported. Ingest startup resumes pending operations before
evaluating input availability, subject to project maintenance exclusion.

Restoring first verifies all retained members. Registry restoration happens only
after the full physical group is back, and only if layout, source selection,
producer heads and active input families still agree. Otherwise the bytes return
to their original locations without selecting stale inputs. Conflicting bytes at
an original location are never replaced. A group containing a purged member cannot
be restored as a partial paired/lane set.

Trashed/purged FASTQ receipts supply exact retained-input evidence without
rewriting old manifests or adding permanent-retention receipts. Restore disables
that trash operation's absence permission; an unrelated later loss is damage
again. While trash is restorable, generation GC pins its producers and consumers.
After purge, its count-proof metadata remains protected so missing producer inputs
can still be verified. Original legacy/permanent-retention receipts remain separate.

When a new count fit needs missing reads, recovery prefers a complete verified
trash group whose source and producer selection are current. It restores that
group before downloading. This also works for locally registered reads without a
download source. Report integrity and current QC eligibility remain independent
requirements; count reuse does not bypass a rejected QC decision.

Expiration uses the recorded deadline, not later configuration changes. All
members of an operation must have expired before purge is selected. A newly
trashed group is never purged in the same invocation. Purge verifies retained
bytes, records `PURGING`, unlinks them and retains its receipt/tombstone.

## Verification and remaining limits

Core and native CLI tests cover paired-group holds, pins, hardlink aliases,
corrupted receipts/payloads, persisted deadlines, restore collisions, changing
source selections, retained native counts, historical input bindings, automatic
timing, permanent-retention precedence and offline restoration for a changed fit.
Real killed child processes cover intent-before-move, partial paired moves,
complete/partial copies, partial restores and purge-after-unlink. Reserved copy
recovery is exercised with injected files; actual multiple-device and power-loss
validation remains a release gate.

Reported byte values are validated input sizes. They are not a measured global
physical-storage total: trash still occupies space, hardlink moves do not duplicate
payload bytes, and cross-device copies temporarily need both copies. Complete
legacy ownership import, abandoned/unpublished data cleanup, global accounting,
resource reservations and measured release/scale gates remain open.
