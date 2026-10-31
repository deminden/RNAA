# Shared quantification indexes

Fresh work now performs a bounded [index layout preflight](quant-index-resources.md)
and checks the estimated loading peak against its declared capacity before FASTQ
staging. The plan is an execution observation shared with the frozen group identity.
Owned CLI work now reserves scan/index/run capacity from the common local pool.
See [resource reservations](resource-reservations.md#current-quantification-integration).

A quantification phase owns one lazy, read-only kallistors counting index for its
compatible worker group. Concurrent workers and later waves share that allocation.
Each worker still owns its FASTQ readers, pseudoalignment counts, fragment/bias
observations, EM state and output files. Threads, strand and layout options remain
per run. The loader and scientific algorithms are unchanged.

The scheduler freezes the selected reference index descriptor before the first
quant worker claim. Constructing the group does not load an index. Historical
selection and reconciliation return before index loading, so healthy cached count
outputs do not require a loaded graph. The first actual quantification loads it;
other workers wait for that load and then borrow the same immutable allocation.
The group retains it for subsequent waves while capacity permits and releases it
when the group and its workers finish. Its resident memory claim lives with the
graph. Cold loading and first-run working capacity are admitted atomically; later
runs reserve only their own working state. There is no permanent process-global cache.

An idle graph may be evicted on a memory wait when the pool has no active CPU
owners, or by the scheduler when no quant work is ready or active. Active runs
protect their graph references. Eviction allows other projects and upstream
preprocessing to make progress and may require a later reload. Separate processes
coordinate capacity but still load independent graphs. A group rejects mixed
resource pools or admitted/unadmitted graph use.

The descriptor binds the reference path, SHA-256 and serialized byte count. Every
use verifies the frozen file, and an initial load is checked again before it is
installed. Quantification checks it once more before writing its manifest. A cached
graph cannot authorize a changed file or a different selected reference. Existing
reference-generation, input and publication fences remain in force. A normal
loading error installs no cached value and leaves the group available for retry;
a loader panic poisons the group and fails subsequent callers explicitly.

The public `QuantIndexGroup` implements the existing `Quantifier` trait. Library
callers can create a group from a frozen `REFERENCE_INDEX` manifest artifact and
share it across their own worker calls. The existing `RustKallistoQuantifier`
entrypoint remains available for an individual call, using a private group.
Callers that publish results remain responsible for the surrounding input,
reference and transaction guards provided by RNAA's command workflows. Compatible
trait calls do not supply shared run admission. Library callers can opt into
`quantify_with_resources` with an owned `ExecutionResources` context, as CLI work does.

## Execution observations

`quant.match_cache_mb` sets a separate exact-match lookup-cache allowance across
the mapping threads of each run, defaulting to 64 MiB. `quant` and `run` accept
`--match-cache-mb N`. Zero and unaddressable sizes fail validation. The library
rounds individual thread tables down to powers of two, so the actual aggregate
active lookup cache is at most the requested allowance. Manifests record both requested and
allocated bytes; status prints the configured per-run MiB. This allowance covers
only the exact-match cache, not the graph, other lookup structures, reads or EM.
The engine restores its default 8 MiB thread-local table on scope exit; run memory
accounting includes that transient capacity. Mapping executes in a joined scoped
thread so private thread-local allocations end before its resource guard drops.

Upstream's automatic policy assigns 512 MiB for one or two mapping threads and
larger tiers for more threads. Applying that policy independently to many RNAA
runs multiplied memory use even for tiny inputs. RNAA now passes an explicit
per-run allowance. Cache size is an execution setting; numerical equivalence and
resource-only reuse checks include different cache allowances.

New quantification manifests include `parameters.index_execution` with a process
group label, loader/engine identity, index digest, serialized byte count and load
attempt/success/eviction counters. `parameters.resource_execution` adds the working
estimate and independently owned index/run reservation descriptors for admitted
work. These observations do not enter the scientific request.
Changing worker allocation or using a shared index does not invalidate healthy
count generations or rewrite their original manifests. Scientific strand changes
still select or compute the appropriate distinct generation.

`index_bytes` describes the serialized file, not resident memory. Sharing an
index does not impose a memory ceiling. Shared CLI reservations account for its
estimated resident lifetime and per-run work. Downloads also reserve their transfer
and publication phases. References admit decoded planning and native build/publication.
Analysis, quant publication admission and
measured combined-workflow qualification remain open.

## Verification and scale probe

Real-index tests check pointer identity across eight competing loaders, lazy
construction, deallocation, failed loads, same-size corruption and reference
mismatch. Concurrent quantification matches independently loaded abundance tables
and run metrics for single/paired layouts, all three strand modes and different
thread allocations. Native CLI checks exercise multiple worker waves, resource-only
reuse, scientific strand changes, historical selection and one failed run followed
by repair while preserving the other successful generations.

To measure loading separately from reference construction:

```sh
cargo build --release -p rnaa-stages --example quant_index_scale
python3 scripts/benchmark_quant_indexes.py \
  --transcripts 1000 2000 --length 600 --workers 1 4 --match-cache-mb 64 512 \
  --report target/release-evidence/quant-index-scale.json
```

The probe builds each deterministic reference in a separate process. Each measured
process then loads and quantifies 64 mapped reads per worker, retaining all loaded
groups until the workers finish. Independent mode gives each worker its own group;
shared mode gives all workers one group. With one mapping thread per run, 512 MiB
matches upstream's former automatic allowance. Reports include observed load counts,
identical abundance hashes, elapsed time, process peak RSS, platform/toolchain and
source hashes. Peak RSS includes validation, mapping, EM and output, not just the
index. These synthetic measurements are not production throughput or supported-host
release qualification.
