# Correlation computation and limits

Correlation prepares input rows once for Pearson, Spearman and bicor, then
computes bounded tiles. Spearman ranks the original values with average ties;
bicor uses WGCNA's default biweight preparation and individual Pearson fallback
for zero MAD. Kendall and HellCor use the pinned mincorr cross-matrix APIs; their
private row preparations are currently repeated per tile. All five methods pass
the checked-in [independent R references](../validation/scientific/README.md).

Top-k keeps at most k candidates per row in the current stripe and never builds
the full gene-by-gene matrix. For G genes, S samples and fixed block size B, its
numerical storage grows with `G*S + B*k + B*B`, plus bounded method scratch.
Dense and threshold outputs keep one stripe, bounded by `B*G`, and write rows in
input order. Dense disk output still grows quadratically. Threshold output may
also contain every distinct pair when the threshold is zero.

Top-k edges are directed, exclude self/undefined pairs and sort by descending
absolute coefficient, then original target order for exact ties. Threshold edges
use the inclusive absolute threshold, emit only the upper triangle and preserve
source/target input order. Constant-row pairs remain undefined; dense TSV writes
`NaN` for them and edge outputs omit them. Variable-row self-correlations are one.
Coefficients within eight machine epsilons of -1/1 become exact endpoints, with
raw independent estimates retained in the validation corpus. Algorithm revision 7
invalidates old correlation requests for this policy/kernel change.

The example configuration includes:

```toml
[corr]
block_rows = 64
max_memory_mb = 1024
max_output_mb = 65536
```

Both limits use MiB. Before loading the numerical matrix, a streaming shape scan
estimates the numerical working payload, growable buffers, method scratch and
worst-case output size. Invalid limits/modes or estimates over either cap fail
before creating correlation output. Output estimates allow unusually long decimal
representations and quoted identifiers; they are deliberately conservative.
Reduce genes/k/block size or raise the declared limit when appropriate. Resource
settings do not enter the scientific request identity. Changes to block size or
thread count preserve exact edge ordering in the allocation regression gate.

These payload admission estimates are not process RSS or OS memory caps. For fresh
CLI work, adjustment and network computation also make cooperative shared resource
reservations; the network reserves estimated memory and primary-output disk space
before computation. Other phases remain outside complete shared admission, and
runtime/allocator overhead and unrelated processes can exceed estimates. See
[resource reservations](resource-reservations.md#correlation-resource-coordination)
for the precise scope. Module enrichment uses
[bounded graph/PCA/enrichment work](module-enrichment.md) with the same limits
and an owned thread pool. Source input is
expected to remain immutable during computation; owned CLI generations provide
the existing publication/content guards.

Manifests record configured limits, estimated bytes, evaluated/undefined entries,
maximum tile/stripe size and maximum retained edge count. A process-level scale
probe is available:

```sh
cargo build -p rnaa-stages --example correlation_scale --release
/usr/bin/time -l target/release/examples/correlation_scale 4000 16 10 64 2 pearson
```

On Linux use `/usr/bin/time -v`. Arguments are genes, samples, k, block rows,
threads and method. The probe generates deterministic input in a temporary project,
reports correlation wall time/output bytes/counters as JSON and removes its data.
The external measurement includes process startup and input generation. Release
profiles, representative public cohorts and supported-host measurements remain
required before setting production performance guarantees.

The initial [development-profile measurements](../validation/correlation-scale/README.md)
record 2,000/4,000-gene top-k runs with raw process RSS, output sizes and source hashes.
