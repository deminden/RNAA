# RNAA maintenance of fasterp 0.2.1

RNAA uses this vendored copy through the root Cargo patch. Upstream remains at
0.2.1. Keep Cargo.toml and Cargo.toml.orig dependency requirements synchronized;
Cargo.toml is the build manifest. Cargo.lock here governs standalone vendor checks;
the root Cargo.lock governs RNAA builds.

On 2026-10-31 the optional native S3 client was updated to aws-sdk-s3 1.157.0;
both lockfiles also resolve aws-smithy-schema 1.1.2. S3 network integration needs
separate live-service verification; the offline vendor tests do not cover it.

Besides the existing workspace Clippy fixes, local changes include:

- Stable dependency updates, including gzp 2 and ureq 3. Optional WebAssembly
  dependencies were refreshed to js-sys 0.3.107 / wasm-bindgen 0.2.130 on 2026-10-30. The parallel compressor
  now names its owned File writer and static lifetime; HTTP downloads use the
  response-body reader and explicitly enable Rustls.
- Multi-member gzip decoding on all native input paths. A single FASTQ file may
  concatenate multiple gzip members, and all reads must reach preprocessing.
- Local HTTP and parallel gzip round-trip tests for the migrated APIs. RNAA also
  checks multi-member input through its actual preprocessing stage.
- The normalized manifest keeps the pure Rust deflate backend used by RNAA. The
  original upstream manifest's native zlib-ng setting remains historical context.
- RNAA's native entrypoints live in `src/native/`. The compatible `process_paths`
  API remains; `process_paths_with_options` accepts explicit batch, complete-record
  and in-flight limits. A record-balanced paired reader replaces independently
  queued mate chunks on this path. CRLF input is normalized.
- A batch credit is acquired before allocation and follows queued input, workers
  and out-of-order results until the ordered merger writes or drops the payload.
  Scoped threads join on producer/writer errors and panic; cancellation wakes
  blocked credit acquisition. Other upstream CLI/Python producers retain their
  existing behavior. The one-thread native path has a complete-record read cap.
- Six native tests cover serial/parallel output and QC parity, unequal mate lengths,
  CRLF and gzip members, invalid options, reordered-result credits and bounded
  failure shutdown. The standalone library suite now has 88 passing tests. RNAA
  preprocessing request revision 3 invalidates older results for these input fixes.
- Native options additionally enforce optional decoded byte/count and minimum
  record promises. Both serial and parallel readers enforce them before accepting
  excess input. Bounded readers return at record boundaries so short reads do not
  repeatedly scan a whole decoder buffer. Two further tests exercise byte growth
  and changed record shape; the standalone library suite now has 90 tests.

Run standalone library checks without the upstream integration suite's external
fastp binary/data requirements:

```sh
cargo test --manifest-path vendor/fasterp/Cargo.toml --lib --locked
cargo check --manifest-path vendor/fasterp/Cargo.toml --lib --target wasm32-unknown-unknown --no-default-features --features wasm --locked
```

RNAA's full workspace verification is still required separately. Python 3.11+
`scripts/check_dependencies.py` checks stable registry versions and both lockfiles.
