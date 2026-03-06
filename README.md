# RNAA

Rust Nucleotide Alignment and Analytics.

RNAA is an RNA-seq survey pipeline for public cohorts. It resolves ENA accessions, downloads raw data with durable SQLite state, preprocesses reads, quantifies with kallisto, runs DESeq2, and computes adjusted correlation outputs.

## Version

Current workspace version: `0.2.0`.

## Recent Changes (Since Last Commit)

- Added sample-scoped cohort composition in `rnaa add`:
  - `--sample`
  - `--include-sample`
  - `--exclude-sample`
- Added metadata predicate selectors in `rnaa add`:
  - `--where field=value`
  - `--exclude-where field!=value`
- Added optional in-process preprocessing with `fasterp` (Rust integration, not CLI shell-out), with strict/bypass modes and reports.
- Switched quant canonical payload to `abundance.h5` (HDF5-first) and updated DESeq2 input flow accordingly.
- Added shared artifact linkage fields and shared-blob tables in SQLite to support cross-project payload reuse foundations.
- Added shared storage config support via `[storage].shared_root` and `RNAA_SHARED_ROOT`.

## Status

Implemented now:

- project init, add, resolve, status, doctor
- ENA resolution for run/study/project/sample accessions
- sample include/exclude selectors and metadata include/exclude predicates
- restart-safe download worker with verification and manifests
- reference preparation with cached Ensembl/custom bundles
- optional in-process `fasterp` preprocessing before quant
- quantification through `Rscript` plus kallisto (HDF5 canonical output)
- DESeq2 through `Rscript` plus tximport/DESeq2
- adjusted correlation through Rust residualization plus `mincorr`
- sequential orchestration via `rnaa run`

Still incomplete:

- `export` and `survey` are stubs
- cleanup policy flags exist, but full trash/purge execution is not finished
- `.sra` extraction fallback is not finished
- project-level stage state still needs to be separated from run-level state
- scheduler-level overlap of download plus quant is planned, not implemented yet

## Requirements

- Rust stable toolchain
- `Rscript`
- `kallisto`
- `curl` or `wget`
- R packages: `jsonlite`, `tximport`, `DESeq2`, `rhdf5`

Optional:

- `xsra`
- `mincorr` CLI is not required; RNAA links the `mincorr` crate directly

## Quickstart

```bash
git clone https://github.com/deminden/RNAA
cd RNAA
cargo build --release

./target/release/rnaa init --root /data/my-rnaa-project
./target/release/rnaa add --root /data/my-rnaa-project --id SRP114962
./target/release/rnaa add --root /data/my-rnaa-project --include-sample SAMN09909193
./target/release/rnaa add --root /data/my-rnaa-project --where instrument_platform=ILLUMINA
./target/release/rnaa resolve --root /data/my-rnaa-project
./target/release/rnaa refs prepare --root /data/my-rnaa-project --organism human --ensembl latest
./target/release/rnaa download --root /data/my-rnaa-project --forever
./target/release/rnaa quant --root /data/my-rnaa-project --preprocess
./target/release/rnaa deseq2 --root /data/my-rnaa-project --design "~ batch + condition" --contrast condition A B
./target/release/rnaa corr --root /data/my-rnaa-project --model "~ batch + condition" --geneset topvar:5000 --out edges:topk=50
```

Single-command orchestration:

```bash
./target/release/rnaa run --root /data/my-rnaa-project
```

## Configuration

See [`configs/rnaa.example.toml`](configs/rnaa.example.toml).

Shared artifact storage can be enabled by setting one of:

- `storage.shared_root` in project config
- `RNAA_SHARED_ROOT` environment variable

## Repository Layout

- `crates/rnaa`: CLI binary
- `crates/rnaa-core`: schemas, config, DB, traits, state, manifests
- `crates/rnaa-stages`: resolver, downloader, refs, preprocess, quant, DESeq2, correlation
- `crates/rnaa-formats`: matrix readers and writers
- `r/`: R entrypoints used by quant and DESeq2 stages
- `fixtures/`: offline test fixtures
- `docs/systemd/`: example Linux services

## Contributing

Run before submitting:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

## License

MIT.
