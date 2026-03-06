# RNAA

Rust Nucleotide Alignment and Analytics.

RNAA is an RNA-seq survey pipeline for public cohorts. It resolves ENA accessions, downloads raw data with durable SQLite state, prepares transcriptome references, quantifies with kallisto, runs DESeq2 normalization/DE, and computes adjusted correlation outputs.

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
- reference preparation with cached GENCODE preset/custom bundles
  default human and mouse presets use the matching GENCODE `transcripts.fa.gz` plus `annotation.gtf.gz`
  reference caches are provider/versioned so stale reference layouts are not silently reused
- optional in-process `fasterp` preprocessing before quant
- quantification through `Rscript` plus kallisto (HDF5 canonical output)
- restart-safe quant reconciliation from existing artifacts and matching reference manifests
- normalization-only stage via `rnaa normalize` (counts, normalized counts, VST)
- DESeq2 through `Rscript` plus tximport/DESeq2
  transcript/gene IDs are normalized across quant and tx2gene inputs, and annotated outputs are written alongside raw ENSG outputs
- adjusted correlation through Rust residualization plus `mincorr`
- restart-safe DE/corr reconciliation from request manifests and existing outputs
- optional module-level enrichment from Spearman modules via Rust `rsfgsea`
- orchestration via `rnaa run` with download plus quant overlap
- quant run-level parallelism (`[quant].workers` or `rnaa quant --workers`)
- operator-facing run/status summaries with stage-aware project progress

Still incomplete:

- `export` and `survey` are stubs
- cleanup policy flags exist, but full trash/purge execution is not finished
- `.sra` extraction fallback is not finished
- project-level stage state still needs to be separated from run-level state
- build output still includes warnings from vendored `fasterp`

## Requirements

- Rust stable toolchain
- `Rscript`
- `kallisto`
- `curl` or `wget`
- R packages: `jsonlite`, `tximport`, `DESeq2`, `rhdf5`

Optional:

- `xsra`
- `mincorr` CLI is not required; RNAA links the `mincorr` crate directly
- `rsfgsea` is linked directly for module enrichment

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
./target/release/rnaa normalize --root /data/my-rnaa-project --design "~ batch + condition"
./target/release/rnaa deseq2 --root /data/my-rnaa-project --design "~ batch + condition" --contrast condition A B
./target/release/rnaa corr --root /data/my-rnaa-project --model "~ batch + condition" --geneset topvar:5000 --out edges:topk=50
./target/release/rnaa corr --root /data/my-rnaa-project --module-gmt /data/pathways.gmt --module-top 500 --fgsea-permutations 2000
```

Single-command orchestration:

```bash
./target/release/rnaa run --root /data/my-rnaa-project --log-file /data/my-rnaa-project/logs/run-live.log
```

`rnaa run` prints compact stage summaries such as:

```text
Project  <project_id>
Elapsed  00:02:14
Stage    quantification | download 6/6 | quant 2/6 | normalize 0/1 | corr 0/1
ETA      00:18:00 total | 00:00:00 download | 00:18:00 processing | limiting processing
Active   quantifying 2 run(s); 2/6 complete
Paths    /data/my-rnaa-project/metadata/samplesheet.tsv | /data/my-rnaa-project/logs
```

`rnaa status` uses the same project-level stage model, so `normalize` and `corr` are reported as `0/1` or `1/1` instead of fake per-run progress.

Annotated sidecar outputs are written where applicable, for example:

- `gene_counts_annotated.tsv`
- `gene_norm_counts_annotated.tsv`
- `vst_annotated.tsv`
- `edges_topk_annotated.tsv`
- `modules_top_genes_annotated.tsv`

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
- `crates/rnaa/tests`: integration tests for normalize and normalize->corr module workflows
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
