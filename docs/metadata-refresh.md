# Editing and refreshing sample metadata

Edit model values in `metadata/samplesheet.tsv`, including `condition`, `batch` and
your declared covariates. Declare their types in `rnaa.toml` as described in
[model designs](model-designs.md). Refresh preserves deliberate edits:

```sh
rnaa metadata refresh
rnaa metadata show
```

`refresh` projects the current stored run metadata and retained read paths; it does
not contact ENA. `show` prints the current publication as JSON and verifies its
retained files against their hashes and committed receipt. Download housekeeping
uses the same merge. It defers refresh while analysis/resolution holds the project
resource, so an independent downloader can still finish successfully. An explicit
refresh requires that resource and records a durable successful or failed attempt.

## Merge rules

RNAA compares the last generated values, your current edits and freshly generated
values, keyed by run accession. It retains an edit when the source value did not
change, takes a source update when you did not edit that value, and accepts a value
both sides agree on. An intentionally blank cell is an edit. New local covariate
columns are retained. Removing a model value does not impute a replacement; model
preflight still checks missingness and replication.

If both sides change a cell differently, refresh stops before publishing a new
metadata or cohort head. The working sheet stays intact. The error lists conflicting
run/column pairs and the directory containing the candidate generated and edited
versions. Inspect them, then explicitly choose a policy for the retry:

```sh
rnaa metadata refresh --conflicts keep-local
rnaa metadata refresh --conflicts use-source
```

Those policies resolve conflicting cells only; unrelated deliberate edits survive.
For a refresh fetched by the resolver, use its corresponding option:

```sh
rnaa resolve --metadata-conflicts keep-local
rnaa resolve --metadata-conflicts use-source
```

`--force` does not discard edits or change the conflict policy. Resolver conflicts
leave the prior cohort, run metadata and metadata publication selected. Candidate
versions remain on disk for review and are not attached as successful results.

## Identities, operational fields and selection

The source/catalog owns project, run, study, sample and experiment identities,
library layout, instrument information and explicit biological identity. Changing
those in the editable sheet fails refresh. For persistent source corrections or
biological grouping changes, register an override with `rnaa add --metadata ...`
and resolve again. Override registration uses unique paths even when filenames
match, and verifies the copy before recording it.

Read paths are generated operational fields; refresh updates them from retained
artifacts. Editing/deleting sample-sheet rows does not change inclusion. Use cohort
selectors and resolve to include/exclude runs. Selectors evaluate resolved source
metadata and registered overrides; local model edits do not silently change them.
The current sheet contains the active cohort only, while retained metadata history
preserves edits for excluded runs. Reincluding a run restores those edits, subject
to the same conflict rules if its source values have changed.

## Retained evidence and recovery

Each publication has a unique directory under `metadata/revisions/` containing
generated, edited and merged tables, generated/merged history, the column map,
resolved run records, copies of registered override sources and a receipt. Original
source locations and hashes remain in provenance; later source edits do not change
the frozen copies. RNAA never rewrites a published revision. Missing/corrupt witnesses
or a receipt that differs from the catalog cause verification to fail.

The metadata head commits in SQLite. Explicit refresh commits its publication event
and successful attempt with that head. Resolution commits the metadata head and
its cohort in one transaction. The editable view installs through an atomic file replacement.
If execution stops between those steps, the next refresh, resolve or analysis
replays the committed view. A later edit is merged against the exact original view;
conflicts stop recovery until explicitly resolved. A missing working sheet can be
restored from the selected publication. Payloads prepared before a failed catalog
commit remain unselected and inspectable.

DE and correlation retain a separate frozen source sheet inside each sealed
analysis generation. Mutable-source observations still fence a running attempt;
later edits do not make a previous result's source witnesses change underneath it.
DE request revision 6 and correlation revision 5 retain source observations
separately from their [scientific metadata cache projection](scientific-cache.md).

Unused operational metadata is excluded from scientific cache identity; current
schema/preflight checks and source observations still protect every attempt.
Automatic reclamation of unselected revisions and process-kill/power-loss recovery release gates remain
separate work. Existing tests cover transaction rollback and both sides of view
installation; they do not establish filesystem power-loss guarantees.
