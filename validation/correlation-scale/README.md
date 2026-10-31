# Correlation scale evidence

`2026-10-31-network-cancellation-debug.json` records the current network tile
execution on the same synthetic top-k probe: 2,000 genes took 1.532 s and 4,000
genes took 5.864 s in separate macOS ARM64 debug processes. This uses the library's
no-op cancellation callback, so it does not measure the CLI's database checks or
establish production throughput. Source hashes and process RSS are in the report.

`2026-10-31-selection.json` records a separate selection-only comparison on a
synthetic 50,000-gene, 500-sample matrix, selecting 1,000 genes with equal
variance. One macOS ARM64 debug run measured 222,838,784 bytes peak RSS for a
reconstruction of the previous full-matrix selection and 10,256,384 bytes for
the current streaming selection; output checksums matched. The record includes source hashes and raw measurements. The original local probe
and standalone manifest were removed from tracked release source because they
contained machine-specific paths. This single
observation excludes downstream analysis, does not control cold caches, and is
not a production throughput or whole-pipeline memory guarantee.

`2026-10-30-modules-debug.json` adds the bounded module graph, full-component PCA,
gene ranking and simple-permutation enrichment workload. It uses the same
deterministic input generator, 16 samples, 64-row tiles, two threads, report-top=10,
100 permutations and seed 42. Each size forms one full module in this probe.

| Genes | Module workflow time | Process peak RSS | All output bytes | Retained neighbours |
| --- | ---: | ---: | ---: | ---: |
| 2,000 | 2.827 s | 15,155,200 | 212,547 | 40,000 |
| 4,000 | 11.136 s | 18,137,088 | 420,091 | 80,000 |

Both use at most 4,096 tile entries. Doubling genes adds 2,981,888 RSS bytes while
one hypothetical dense matrix grows from 32,000,000 to 128,000,000 bytes. The
report includes recursive module source hashes, declared estimates, PCA/module
diagnostics and raw process measurements. The older top-k report below describes
the batch-twenty-nine implementation; its source hashes are historical evidence.

```sh
cargo build -p rnaa-stages --example module_scale
python3 scripts/benchmark_correlation.py --binary target/debug/examples/module_scale --genes 2000 4000 --method spearman --report /tmp/module-scale.json
```

These are development-profile synthetic measurements, not supported-host release
qualification. Allocation estimates do not cap OS RSS or preceding CLI selection
memory, and this probe does not exercise multiple concurrent module jobs.

`2026-10-30-debug.json` records separate-process synthetic Pearson top-k runs on
macOS ARM64 with Rust 1.99.0, 16 samples, k=10, 64-row tiles and two threads.
Source/lock hashes, raw `/usr/bin/time -l` measurements, exact output sizes and
allocation counters are included in the report.

| Genes | Correlation wall time | Process peak RSS | One hypothetical dense f64 matrix |
| --- | ---: | ---: | ---: |
| 2,000 | 1.498 s | 9,125,888 bytes | 32,000,000 bytes |
| 4,000 | 5.658 s | 9,846,784 bytes | 128,000,000 bytes |

Both runs retained at most 640 edges and 4,096 tile entries at once; neither
allocated a dense output stripe for top-k. Doubling the gene count increased
measured RSS by 720,896 bytes. The dense column is the size of one f64 matrix,
not a measurement of the previous implementation or an end-to-end speed comparison.

These measurements use the development profile. RSS includes process startup and
input generation, while the JSON correlation time starts after input generation.
CPU work still grows quadratically with the number of genes. The data demonstrate
bounded top-k working buffers and observed memory growth for this workload; they
do not establish production throughput, every method's scaling, module-enrichment
memory, global reservations or supported-host release readiness.

Reproduce from repository root:

```sh
cargo build -p rnaa-stages --example correlation_scale
python3 scripts/benchmark_correlation.py --binary target/debug/examples/correlation_scale --genes 2000 4000 --report /tmp/correlation-scale.json
```

For release evidence, build with `--release` and select the release binary. Each
report records its actual build profile; compare source hashes before interpreting
numbers as measurements of another implementation. See the [resource contract](../../docs/correlation-resources.md).
