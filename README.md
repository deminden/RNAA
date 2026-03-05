# RNAA

Rust Nucleotide Alignment and Analytics.

RNAA is a RNA-seq survey pipeline for public cohorts. It resolves public accessions, downloads raw data with durable state in SQLite, quantifies with kallisto through an R runner, performs DESeq2 analysis and normalisation, and computes adjusted correlation outputs through Rust plus `mincorr`.

## Status

Implemented now:

- project init, add, resolve, status, doctor
- ENA metadata resolution for run/study/project accessions
- restart-safe download worker with verification and manifests
- reference preparation with cached Ensembl/custom bundles
- quantification through `Rscript` plus kallisto
- DESeq2 through `Rscript` plus tximport/DESeq2
- adjusted correlation through Rust residualization plus `mincorr`
- orchestration via `rnaa run`

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
- R packages: `jsonlite`, `tximport`, `DESeq2`

Optional:

- `xsra`
- `mincorr` CLI is not required; RNAA links the `mincorr` crate directly

## Quickstart

```bash
cargo build --release

./target/release/rnaa init --root /data/my-rnaa-project
./target/release/rnaa add --root /data/my-rnaa-project --id SRP114962
./target/release/rnaa resolve --root /data/my-rnaa-project
./target/release/rnaa refs prepare --root /data/my-rnaa-project --organism human --ensembl latest
./target/release/rnaa download --root /data/my-rnaa-project --forever
./target/release/rnaa quant --root /data/my-rnaa-project
./target/release/rnaa deseq2 --root /data/my-rnaa-project --design "~ batch + condition" --contrast condition A B
./target/release/rnaa corr --root /data/my-rnaa-project --model "~ batch + condition" --geneset topvar:5000 --out edges:topk=50
```

Single-command orchestration:

```bash
./target/release/rnaa run --root /data/my-rnaa-project
```

## Configuration

See [`configs/rnaa.example.toml`](configs/rnaa.example.toml).

The example keeps download concurrency intentionally small by default. Download is an I/O pool, not a CPU saturation target.

## Repository Layout

- `crates/rnaa`: CLI binary
- `crates/rnaa-core`: schemas, config, DB, traits, state, manifests
- `crates/rnaa-stages`: resolver, downloader, refs, quant, DESeq2, correlation
- `crates/rnaa-formats`: matrix readers and writers
- `r/`: R entrypoints used by quant and DESeq2 stages
- `fixtures/`: offline test fixtures
- `docs/systemd/`: example Linux services

## Development


## Contributing

Contributions are welcome and encouraged! If you’d like to help improve `RNAA`, feel free to open an issue for bugs, feature requests, performance ideas.

Pull requests are especially appreciated for:
- performance improvements
- additional methods or variants
- more or better tests
- docs, examples, and benchmarks

Plese run, before submitting:
```bash
cargo fmt --all -- --chec
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
```

## License

MIT. See [`LICENSE`](LICENSE).
