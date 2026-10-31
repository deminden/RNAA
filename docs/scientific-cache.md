# Scientific metadata and cache reuse

DE request revision 6 and correlation revision 5 separate scientific model values
from complete source/configuration observations. The earlier request form requires
one recomputation. No statistical engine or pooling algorithm changes in this revision.

## What determines reuse

A DE request contains the contributing runs grouped into biological samples and
the actual typed model: referenced factor/numeric values, observed factor levels,
references, sample order, coefficients, matrix, offsets, rank, residual degrees of
freedom and resolved categorical contrasts. Count/reference identities, transform,
import settings, engine versions and QC admission remain dependencies.

Sample-sheet byte order, unrelated annotations, read paths and unused schema
declarations do not change that scientific projection. Finite numeric values such
as `21` and `21.0` share their normalized representation. Factor labels `01` and `1`
remain different categories. A fixed formula uses its referenced variables; `.`
expands the declared predictors present in the sheet, so adding/changing one can
change that model. Literal formula changes remain part of request identity.

Cohort cache identity records resolved run-to-biological-sample membership.
Equivalent selectors or source annotations can reuse the same scientific result
after resolution. Current selectors must still be resolved before processing;
membership, identity, counts or QC changes can invalidate downstream requests.

Correlation records its actual adjustment design. With `adjust = "none"`, model
text and metadata values have no adjustment effect and are excluded from the
scientific projection. Matrix content, method, selected genes, output settings,
annotations and enrichment inputs still determine the correlation request.

Current schema validation and statistical preflight run before accepting a cache
hit. A newly declared missing marker, absent required level, changed covariate type
or insufficient biological replication can therefore reject a request even when
no new fitting has started. Unused values are interpreted when a formula starts
referencing them.

## Producer and consumer metadata

Before consuming the selected DE matrix, correlation reconstructs the current DE
model/membership and compares it with that matrix's recorded scientific contract.
A change to a producer covariate or its effective factor coding requires DE again.

After that check, correlation builds its own adjustment from current grouped
metadata. For example, with DE `~ condition` and correlation `~ age + condition`,
changing `age` permits DE reuse and refits correlation using the new ages. This
also checks agreement across technical contributions for every currently used
adjustment covariate. The original DE sidecar remains the producer's exact record.

A new producer generation with identical matrix bytes can still satisfy an
unchanged downstream correlation request. Scientific identity follows the consumed
values; the original result and each current reuse attempt retain their evidence.

## Observations and retained evidence

`model_design.json` retains the full schema observed during computation.
`configuration_observation` in the sealed manifest/request descriptor records full
settings separately from the scientific request. Raw sheets remain hashed inputs
of the sealed generation through `input_witnesses` and manifest input artifacts.

On reuse, RNAA retains the current attempt's raw tables under
`metadata/analysis_observations/`, addressed by SHA-256, input role and filename.
The reconciliation event records those witnesses and its configuration observation
in the same transaction as successful completion. The original generation remains
unchanged, including its original source sheet and sidecars. Existing observations
are verified and retained without replacement; corrupt evidence stops reuse.

During an attempt, exact source bytes, complete catalog metadata, selectors, the
concrete cohort head and observed schema/settings still fence publication. Even an
unrelated edit during computation causes a retry, preserving a coherent observation
for each published or reused attempt. Hashing/copying occurs before the final writer
transaction; catalog fencing and event/outcome updates occur inside it.

## Returning to an earlier request

Normalization, differential expression and correlation search committed project
history for the exact current request, preferring the selected generation and then
the newest matching result. Migration 0015 adds the scoped request index. The
output group, stage and subject are part of lookup; normalization and contrast
results cannot stand in for each other even though they share the DE output group.

Each candidate must pass its complete seal/inventory and lineage checks, descriptor,
scientific manifest, required output coverage, input hashes and current admission
checks. A corrupt payload is skipped and the request is recomputed if no healthy
candidate remains. Catalog corruption is reported as an error. A metadata edit
followed by restoration can therefore select the original fit rather than producing
a third output set. An adjustment covariate restored to its original values can
similarly select its earlier correlation fit.

Selection restores the original artifact registrations and changes the generation
head in one transaction with current admission/configuration fencing, state/events
and reuse-attempt success. It compares the head observed during verification and
the candidate's catalog/lineage again before selection. Cancellation, stale ownership
or a failing transaction cannot switch the head or leave partial registrations.
The producer's attempt, seal, files and commit timestamp remain unchanged; events
link the new reuse attempt to that original producer and previous selected head.
Readers already pinned to the intervening generation retain their paths.

Use `rnaa lineage --group de --all-requests` or the equivalent `--group corr` to
inspect committed request variants for selected output subjects. JSON includes a
`selected` flag per entry. Add `--history` for upstream ancestry or `--verify` to
rehash payloads. A corrupt historical result makes verification fail even when the
currently selected result is healthy, with the JSON report still printed.

## Ingest requests and input bindings

Migration 0016 stores complete scientific requests as content-addressed catalog
objects. Fresh reference, download, preprocessing, count and analysis generations
register their request in the same transaction as publication. A policy index keeps
scientific options/engines/dependencies while excluding the input multiset; a match
there is candidate discovery, never proof that input data agree. Older native count
requests can be indexed after verifying their full digest against the committed
head, without rewriting their manifests or seals.

Preprocessing and native quantification can select earlier results using current
source/processing identities, QC policy and exact content witnesses. Original
producer paths can differ from the current download's paths. An immutable input
binding records the reuse attempt, original consumer, current input witnesses and
their exact sealed owners. The binding's digest and head commit with selection.
Fresh publication clears the selected binding; historical binding receipts remain
available as durable evidence.

Bindings require the same scientific input identities, including roles, filenames,
hashes, sizes and multiplicity. A count request made from retained trimmed reads may
omit originals. Additional current originals can then serve as source/retention
observations only if their entire identity matches that trimmed producer's original
inputs. Partial groups, altered bytes, different roles and unrelated originals do
not satisfy this proof. Native engine values and request revisions are unchanged.

Count reconciliation compares the selected source and processing content even
after all reads have been pruned; valid old deletion receipts alone cannot make an
old source current. Current report integrity/admission stays separate from count
payload validity. An earlier healthy report can repair QC evidence without reads
or network access when its request/source proof agrees. Reconstructing a report
when no healthy evidence remains uses the
[QC report recovery workflow](dependency-recovery.md#qc-report-reconstruction).
Complete original reads are processed again; matching recovered inputs can reuse
existing counts. Trimmed reads alone cannot reproduce original-read QC evidence.

`lineage` JSON includes the selected `input_binding` and, with `--verify`,
`input_binding_verified`. `--history` traverses both original producer edges and
the selected use-time binding's producer edges. Original seals and sidecars retain
the original production history.

## Download and reference history

Download selection searches exact scoped request history and verifies each sealed
payload, its manifest, complete selected file group, hashes and declared sizes.
Publisher MD5 values are checked against the actual files. A matching publisher
checksum permits selecting earlier history without contacting the source.

Without a publisher checksum, a URL request describes a source selection rather
than a permanent object version. The currently selected verified snapshot keeps
its offline cache behavior. Selecting an older snapshot requires a saved strong
ETag and a current conditional HEAD response that repeats the same strong ETag,
effective URL and compatible representation length/encoding. A bare 304, weak or
changed ETag, unsupported HEAD, missing receipt or unavailable origin does not
prove byte identity. RNAA then fetches a new snapshot; a failed fetch preserves
the selected generation. This follows the strong validator and conditional request
distinction in [HTTP semantics, RFC 9110](https://www.rfc-editor.org/rfc/rfc9110.html#name-validators).

Private file caches surviving a failed multi-file download also require current
origin validation when no publisher checksum exists. A successful individual file
transfer is insufficient evidence to publish stale cached bytes on a later retry.
Reused generation selection records its source-validation evidence with the new
attempt, registrations, run state and generation head in one transaction. Missing
pruned reads are fetched again when actual files are needed; deletion receipts do
not recreate those files.

Reference history uses the scientific source-content, annotation policy, index
settings and engine request. A candidate must satisfy full source SHA-256 checks,
manifest/output coverage and sealed inventory verification. Identical custom
sources moved to new paths can reuse a bundle while recording current source
observations. Changing index settings and restoring them selects the original
healthy generation. Returning to a previously prepared source also updates the
project's reference choice, even when that source's resource head never changed.

Reference choice changes run under an owned reference attempt and commit with the
head, artifact registry, current policy/parent-choice fences and success event.
Unchanged access needs no new attempt. Configuration and original source changes
during preparation reject publication; execution thread changes preserve reuse.
Corrupt historical bundles remain available for diagnosis and are rebuilt into a
new generation. Resolving provider references still uses provider listing/checksum
metadata; local bundle reuse does not imply offline provider resolution.

Inspect variants with `rnaa lineage --group download --subject RUN1 --all-requests`
or `rnaa lineage --group refs --all-requests`. The `selected` flag is per resource
head; the selected project reference choice is a separate relationship. Original
payloads and producer seals remain unchanged by reuse.

This is project-scoped scientific reuse within the current producer/admission
contracts. [Supported pruned count-input recovery and sealed producer edges](dependency-recovery.md)
are implemented. Generalized legacy dependency propagation, cross-project
sharing and observation garbage
collection remain separate work. Synthetic invariants and native workflow tests
do not replace independent scientific reference fixtures or scale/release gates.
