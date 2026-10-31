# Importing legacy shared-blob ownership

Older projects recorded shared blobs only in their own SQLite catalogs. A shared
store cannot discover every such project from one project's references. Import
therefore takes an explicit inventory of project roots:

```sh
rnaa storage --shared-root /data/shared import-legacy \
  --project /data/project-a --project /data/project-b
rnaa storage --shared-root /data/shared import-legacy \
  --project /data/project-a --project /data/project-b --apply
```

Preview is the default. Apply imports verified ownership while retaining legacy
holds unless complete inventory and physical-link coverage are established:

```sh
rnaa storage --shared-root /data/shared import-legacy \
  --project /data/project-a --project /data/project-b \
  --complete-project-inventory --apply
```

`--complete-project-inventory` asserts that every legacy project using this store
is supplied and older clients that bypass the current locking/publication protocol
are stopped. The store checks this assertion against every currently known owner
project and its durable record of previously imported project identities. Missing
known projects prevent release of legacy holds, including holds on otherwise
ownerless unindexed blobs. Unknown projects cannot be inferred from file hashes;
the explicit inventory assertion supplies that part of the migration contract.

## Project and payload evidence

The importer opens each `state.sqlite` read-only and takes exclusive project
maintenance exclusion before inspection. It does not upgrade project schemas,
rewrite project records, create result generations or change selected heads.
Coordination files can be created. Normal shared-store opening upgrades its own
catalog to schema 2, transactionally retaining existing blobs, owners and holds.
Concurrent initializers inspect the version after acquiring SQLite's writer.
Normal project opening/creation takes maintenance exclusion before connection,
schema upgrade and metadata backfill, so those writes cannot race import.
With explicit `--shared-root`, import also bypasses the invoking root's settings
file; old configuration syntax can remain untouched during catalog migration.

Supported project migration histories cover the shared-reference schema through
`0021_run_progress`. This progress-only migration adds audit receipts and guards
without changing ownership columns. Unknown/missing histories stay protected. Future changes
to project reference semantics must explicitly update the shared ownership reader's
supported boundary. Old projects need not gain generation tables to be inspected.

Evidence includes active artifact registrations, local blob descriptions and shared
references, plus retained committed generation registrations and supported retirement
receipts. Shared records must agree on project, blob ID, absolute path, SHA-256 and
length. Generation ownership also requires matching sealed file inventory and seal
witness evidence. Retained generation trash supplies its physical view; completed
project purges are excluded. Pending project maintenance must finish before import.

Existing local views are checked against their declared bytes. A missing unsealed
view can retain its declared shared ownership while its verified shared copy survives;
a missing retained generation view is an error. The shared payload itself must be
a regular file at its exact content-addressed path and match its digest and length.
Scientific result validity remains a separate check: importing ownership never
adopts an unsealed result as a successful scientific generation.

## Complete inventory and hardlinks

On Unix, the importer compares the blob's actual link count with distinct verified
physical view paths sharing its device/inode, including retained project trash.
Aliases naming the same physical path count once. An unexplained hardlink keeps
the legacy hold even when complete project inventory is asserted. Platforms without
the required link metadata can import owners but cannot release legacy holds.

Independent project copies also retain ordinary ownership; matching content does
not imply an unregistered owner. Canonically named unindexed live shared files can be
cataloged after hash verification. Unexpected names, symlinks, temporary files and
other unindexed entries remain in place and appear in the report.
Discovery scans the live `blobs/sha256` namespace, declared project references and
cataloged legacy holds. Unindexed store trash is preserved separately; it is not
enrolled as a live blob by this command.

The report identifies proposed/applied imports, busy blobs, errors, owner counts,
observed/accounted links, and reasons a hold cannot be released. Applied publication
returns `import_id`; shared status exposes each blob's `legacy_import_id` so its
receipt can be located directly. `complete` requires all inspected candidates and
unindexed entries to be resolved. An incomplete import
can add verified protection without clearing holds. Busy projects and evidence
errors produce failing command status after the JSON report; ordinary held reasons
remain inspectable report outcomes. No blob is moved or removed by import.

## Publication, later collection and recovery

An immutable descriptor binds the store, explicit project inventory, original
catalog fingerprints, imported owners, link observations and individual hold
decisions. A second witness lives at `.imports/IMPORT_ID.json`. Shared blobs are
locked in deterministic order, project catalogs and link observations are rechecked,
then one shared transaction publishes the verified candidates, owners, import indexes,
project census and hold changes. A failed transaction publishes none of those
changes. Interrupted pre-commit witness files are never adopted merely by filename.

A subset inventory cannot replace a previously released import proof. Later GC
checks descriptor/catalog/witness identities and still-present imported project
claims; deleting a store ownership row cannot erase an active migrated claim.
Missing, substituted or unsupported project catalogs stay protected. Normal
generation ownership remains until a validated completed project purge.

Once ownership disappears, additional hardlinks still prevent shared trash/purge,
including links added after import or during interrupted cleanup. Durable trash
and purge intents recheck ownership and links before finishing. Reader/writer OS
exclusion remains authoritative over stale observation rows. Completed purge clears
the live blob's import link while retaining immutable import history and project
census, so a later independent publication of the same digest has its own owners.

After migration, use the ordinary shared-store preview, trash, restore and recorded
expiry workflow in [shared storage](shared-storage.md). Import does not remove old
project views or authorize unsealed result cleanup. Unknown/unavailable projects,
unindexed data, transfer/support cleanup, resource reservations, cross-project
scientific request reuse and measured release/scale/power-loss gates retain their
separate plan requirements.

## Verification

Core coverage includes original-schema read-only imports, multiple projects,
unindexed blobs, known project omissions, subset-proof preservation, missing indexes/
owners/witnesses, unsupported schema changes, unknown/new hardlinks, busy projects
and blobs, transaction rollback, concurrent schema initialization and real killed
importers before/after publication. Native CLI coverage preserves original quant
seal bytes and current count reuse, imports historical owners, releases them only
after project purge, and checks standalone old-project inspection without upgrade.
