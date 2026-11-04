# RNAA maintenance of fasterp

This vendored source incorporates GitHub `drbh/fasterp` main revision
`dcca6ad25e238ecfb4c2372ae5d4c83b6cc68281` (2026-11-02). Upstream's latest
published crate remains 0.2.1; RNAA reports `0.2.1+git.dcca6ad25e23` so this
unpublished source revision cannot be confused with the crates.io release.
`LICENSE.fastp` preserves the upstream notice for its fastp-derived code.

RNAA maintains additional integration changes:

- `src/native/` provides `process_paths_with_options` with explicit decoded-input,
  record and in-flight batch limits. Paired reads remain record-balanced, CRLF is
  normalized and gzip input handles concatenated members.
- Scoped workers and batch credits bound queued and reordered input, join on errors,
  and unblock on cancellation. The upstream-compatible CLI engine is included but
  RNAA continues using its bounded native path rather than silently changing engines.
- The newer `gzp` 2 and `ureq` 3 APIs, current stable dependency requirements and
  Rust deflate backend are retained. HTTP input streams response bodies; the upstream
  parallel libdeflate gzip encoder was adapted to `gzp` 2's writer type.
- The upstream randomized overlap test uses the current `rand` API. Root and vendor
  lockfiles are maintained separately and both must be checked after updates.

Run the local vendor gates:

```sh
cargo test --manifest-path vendor/fasterp/Cargo.toml --lib --test compatibility_v140 --locked
cargo check --manifest-path vendor/fasterp/Cargo.toml --lib --target wasm32-unknown-unknown --no-default-features --features wasm --locked
```

The upstream differential tests against an external fastp 1.4.0 executable require
that exact external tool and are separate from the self-contained compatibility suite.
Run the full RNAA workspace checks and independent scientific references as well.
