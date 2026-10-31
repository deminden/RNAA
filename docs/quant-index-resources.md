# Quantification index preflight

Fresh native quantification inspects the frozen index before decoded FASTQ staging
or graph loading. `QuantIndexPlan` uses a 64 KiB input buffer and fixed bitmap
scratch space. Serialized counts must fit their physical enclosing sections before
they can drive a loop. Index bytes and graph/node/MPHF section sizes are checked
with overflow-safe arithmetic. No decoded graph, transcript-name table or expanded
bitmap is allocated by this inspector.

The plan covers the current kallistors 0.6.0 counting loader and its kallisto v13
layout, with k/g up to 32. It records graph header words without adding version
restrictions for words the counting loader ignores. It handles native and portable
Roaring arrays, bitsets and runs, Bifrost inline/tiny/pointer tags, special minimizer
counts, EC blocks and transcript/onlist data. Portable container keys, offsets and
actual bitset/run cardinalities are checked. Native array counts remain conservative
upper bounds when the engine sorts/deduplicates their values.

The tiny RLE checks account for the present engine's look-ahead cursor. A short
payload that would cause that decoder to index beyond its input is rejected before
loading. This input guard does not replace or modify the engine library.

## What is counted

| Structure | Observations |
| --- | --- |
| Graph | Unitig/k-mer/abundant nodes, decoded bases and largest unitig |
| Minimizers | Bucket count, decoded graph bitmap values and special markers |
| MPHF | Levels, elements, word/rank arrays, final hash entries and serialized buffer |
| EC data | Block count, decoded values upper bound, largest class/target, subset envelopes and node scratch |
| Transcripts | Count, total/largest name bytes, shade-name detection and decoded onlist values |

Resident and loading estimates include decoded/encoded graph copies, Vec capacity
growth, node metadata and fallback lookups, minimizer positions/directories, EC
directories/values, MPHF representations/hash slack and transcript/shade metadata.
Loading adds serialized-buffer peaks, decoded bitmap/special-marker scratch,
largest node/unitig/name buffers, fixed processing overhead and 2 MiB per requested
loader worker for stack allowance. Sums use wide checked accounting before
conversion to supported byte capacities. Serialized file size is one observation,
rather than a generic multiplier standing in for decoded structures.

The public quantification path rejects a loading estimate exceeding its declared
`[resources].memory_mb` capacity; 0 uses the common policy's detected capacity.
Rejection precedes FASTQ materialization and loader attempts. Increasing the limit
permits a retry, and lowering it does not invalidate healthy cached counts. The
layout plan is retained with the group's frozen index identity and included in
`index_execution.layout_plan` as an execution observation. It does not change the
scientific request key. Ordinary worker waves share one loaded graph; pressure
eviction can require a later reload.

## Shared admission and limits

Owned CLI quantification also uses these observations in the
[shared reservation protocol](resource-reservations.md#current-quantification-integration).
It releases its scan claim before atomically admitting the cold index and first
run. Their combined memory covers the maximum of the loading phase plus setup and
the resident graph plus run working state. The graph retains its own guard after
the loading claim shrinks; later runs each obtain independent CPU/memory/disk claims.
Only idle graphs can be evicted to free capacity. Staging and numerical execution
start after admission, and run guards outlive their buffers, mapper threads and
temporary files. Shared-blob copying and generation publication remain outside this
guard and still need admission.

Working estimates include mapping/EC/count/EM state, caches, packed batches and
staging/output promises. They bound class diversity conservatively using EC subset
envelopes and observed fragments; shade names permit target unions. Such estimates
can reject jobs with much lower actual observed diversity. Numerical constants
require supported-host RSS qualification; estimates are not operating-system
memory ceilings. Library calls without an execution context retain only the local
loading-peak check and do not join shared run admission.

Inspection is not a full semantic validator of every engine structure. Packed DNA,
rank-table contents and positional payloads are bounded/skipped where decoding is
unnecessary for counting allocations. Special marker widths are not reconstructed:
their minimum physical extent and declared graph end bound the owning table.
Producer seals, source hashes and the engine's own semantic checks remain necessary.
Frozen input checks assume cooperative immutable files, rather than an adversarial
writer replacing and restoring a path during inspection.

Owned CLI quantification passes its cancellation/fencing callback into index and
FASTQ preflight/staging and checks it between loading, mapping, EM and output
phases. Active engine calls still run until they return. Numerical execution lives
in `quant/native.rs`; request identity and reconciliation remain separate.

## Evidence

Real generated indexes with 1, 2, 8 and 120 shared targets match the engine's public
metadata and loaded graph. Tests cover huge counts, negative/unsupported fields,
truncated sections, callback cancellation, native/portable bitmap forms and special
minimizer compatibility. A bitmap representing the full 32-bit value space in less
than 1 MiB is counted without expansion and rejected against a 512 MiB capacity;
the expanding engine decoder is never called for that fixture. Native mapping
parity still covers single/paired layouts, strand and thread settings. CLI checks
cover capacity rejection, successful retry and unchanged healthy reuse after lowering
the limit. Further tests verify separately owned index/run claims, atomic cold
admission, cancellation and eviction, same-output paired mapping, competing CLI
processes and constrained preprocessing/quantification progress.
