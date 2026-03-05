# Architecture

RNAA is organized as a contract-first Rust workspace.

## Crates

- `crates/rnaa`: CLI entrypoint and orchestration
- `crates/rnaa-core`: config, schemas, SQLite DB, manifests, state machine, traits
- `crates/rnaa-stages`: resolver, downloader, references, quantification, DESeq2, correlation
- `crates/rnaa-formats`: matrix readers and writers

## Design intent

- stage contracts live in Rust traits
- artifacts are written to stable per-stage locations
- every stage writes or records enough provenance to be rerun safely
- shell and R integrations are implementations, not the architecture

## Current caveat

Project-level stages such as DESeq2 and correlation still reuse run-level state transitions. That is functional today, but a dedicated project-stage state model is still planned.
