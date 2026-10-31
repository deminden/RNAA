# Observed project and shared-store usage

```sh
rnaa storage usage
rnaa storage usage --details
rnaa storage usage --include-shared
rnaa storage usage --include-shared --shared-root /data/shared
```

The command prints JSON with scan timestamps, roots, category totals, hardlink
counts, activity flags and errors. `--details` adds individual path observations.
An explicit shared root requires `--include-shared`. Otherwise the shared root
comes from the project configuration/environment using the normal storage rules.
Inspection does not create a missing shared store or delete any data.

## Logical paths and allocated inodes

`logical_bytes` sums the sizes of observed regular-file paths. Category totals use
this quantity and add up to the reported logical total. Two hardlinked paths each
contribute their logical size, because each is a visible file reference.

On Unix, `inode_allocated_bytes` sums allocated blocks once per distinct
`(device,inode)`, including directory and symlink metadata. Hardlinks inside a
project and across included shared storage therefore count once. Sparse files use
their reported allocation rather than their apparent length. Filesystem compression,
copy-on-write clones and shared extents can make reclamation differ from these
inode totals; this value is not a filesystem free-space measurement.

`observed_hardlink_aliases` records extra observed paths to regular-file inodes.
`unobserved_hardlinks` records regular-file link counts beyond the scanned paths;
it does not identify those other owners or authorize collection. Platforms without
the required inode/allocation metadata return null for those quantities.

Canonical overlapping project/shared roots are traversed once. Symlinks are
observed as entries and their targets are not followed. Repeated directory inode
aliases are skipped to prevent traversal cycles, with an explicit counter.

## Categories and consistency

Categories distinguish committed generations, abandoned or unfinished generation
payloads, unindexed generation directories, generation/input/unpublished trash,
transfer caches, support/temporary files, registered legacy files, project metadata, shared blobs,
shared blob trash and other unindexed data. Categories describe catalog/path
associations and are **not collection eligibility decisions**. Active attempts can
own unfinished files; unindexed paths can still be needed by legacy results.

The scan holds a project reader pin so project retirement/trash cannot remove a
generation underneath it. Workers may continue writing. `catalog_changed`,
`shared_catalog_changed` and `active_attempts_observed` qualify that observation;
external filesystem changes also remain possible. This is not an atomic filesystem
snapshot or a CPU/RAM/disk admission reservation.

Traversal errors and inode changes appear in `errors` and set `complete: false`.
The CLI returns failure after printing that incomplete report. Unreadable or
unsupported catalogs can prevent a report from being produced. Concurrent files
may require another observation after work finishes.

Core tests compare hardlinks across separate/nested roots, symlink exclusion,
logical category totals, allocated inode deduplication and unfinished states.
Native CLI coverage measures input trash, retired generations and abandoned/
unindexed payloads while verifying that those bytes and current counts survive.

The existing coarse storage-limit scanner remains separate. Tracked admission
reservations, unresolved legacy ownership, unsealed adoption and ownership-aware
unindexed/transfer/support cleanup are still required before the storage roadmap
is complete. [Shared legacy import](shared-legacy-import.md) and
[tracked unpublished cleanup](unpublished-cleanup.md) are separate
from this observational report.
