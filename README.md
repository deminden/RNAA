# RNAA

Rust Nucleotide Alignment and Analytics.

RNAA is an RNA-seq survey pipeline for public cohorts. It resolves ENA accessions, downloads raw data with durable SQLite state, prepares transcriptome references, quantifies with kallisto, runs DESeq2 normalization/DE, and computes adjusted correlation outputs.

## Version

Current workspace version: `0.2.1`.

## Capabilities

- project init, add, resolve, status, doctor
- ENA resolution for run, study, project, and sample accessions
- sample include/exclude selectors and metadata include/exclude predicates
- restart-safe download worker with verification and manifests
- cached GENCODE preset references for human and mouse
  default presets use matching `transcripts.fa.gz` plus `annotation.gtf.gz`
- optional in-process `fasterp` preprocessing before quantification
- FASTQ retention after successful quantification: `both`, `original`, `trimmed`, or `none`
- threaded `fasterp` preprocessing with configurable per-run thread count
- normalized preprocess QC storage with conservative default gating before quant/downstream stages
- quantification through `Rscript` plus kallisto
- restart-safe quant reconciliation from existing artifacts and matching reference manifests
- normalization-only stage via `rnaa normalize`
- DESeq2 through `Rscript` plus tximport/DESeq2
- adjusted correlation through Rust residualization plus `mincorr`
- restart-safe DE/corr reconciliation from request manifests and existing outputs
- optional module-level enrichment from Spearman modules via Rust `rsfgsea`
- orchestration via `rnaa run` with download plus quant overlap
- quant run-level parallelism via `[quant].workers` or `rnaa quant --workers`

## Limitations

- `export` and `survey` are still stubs
- cleanup policy flags exist, but full trash/purge execution is not finished
- `.sra` extraction fallback is not finished
- project-level stage state is still inferred from run-level state
- vendored `fasterp` remains patched for workspace quality gates

## Requirements

- Rust stable toolchain
- `Rscript`
- `kallisto`
- native HTTP downloads via `reqwest` are built in
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
./target/release/rnaa quant --root /data/my-rnaa-project --preprocess --preprocess-threads 16 --fastq-retention both
./target/release/rnaa qc reevaluate --root /data/my-rnaa-project
./target/release/rnaa normalize --root /data/my-rnaa-project --design "~ batch + condition"
./target/release/rnaa deseq2 --root /data/my-rnaa-project --design "~ batch + condition" --contrast condition A B
./target/release/rnaa corr --root /data/my-rnaa-project --model "~ batch + condition" --geneset topvar:5000 --out edges:topk=50
./target/release/rnaa corr --root /data/my-rnaa-project --module-gmt /data/pathways.gmt --module-top 500 --fgsea-permutations 2000
```

Single-command orchestration:

```bash
./target/release/rnaa run --root /data/my-rnaa-project --preprocess-threads 16 --fastq-retention trimmed --log-file /data/my-rnaa-project/logs/run-live.log
```

## Runtime Behavior

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

Downloaded raw data, prepared references, quant outputs, and DE/corr outputs are persisted on disk and tracked in SQLite. Re-running commands reconciles existing artifacts where possible instead of redoing work blindly.

When preprocessing is enabled, RNAA can prune FASTQ inputs after a successful quantification:

- `--fastq-retention both` keeps original and trimmed FASTQ files
- `--fastq-retention original` keeps only original FASTQ files
- `--fastq-retention trimmed` keeps only trimmed FASTQ files
- `--fastq-retention none` removes both original and trimmed FASTQ files and keeps only downstream quant outputs
- `--preprocess-threads N` sets per-run `fasterp` threads when preprocessing is enabled

Removed FASTQ files are permanently deleted after successful quantification and removed from the active artifact set.

Preprocess QC is parsed into a normalized SQLite table and gated by default before quant/downstream analysis. The built-in defaults are conservative fail-fast thresholds intended to reject catastrophically bad files:

```toml
[quant.qc_gate]
min_pass_rate = 0.60
max_low_quality_rate = 0.30
max_too_many_n_rate = 0.10
max_too_short_rate = 0.40
```

You can override them in `rnaa.toml`, for example:

```toml
[quant.qc_gate]
min_pass_rate = 0.8
max_low_quality_rate = 0.2
max_too_short_rate = 0.2
```

Runs that fail the configured QC gate are recorded in `preprocess_qc` with a gate reason and are excluded from quant/downstream processing.

If you change the gate after preprocessing has already been completed, rerun the gate without rerunning `fasterp`:

```bash
./target/release/rnaa qc reevaluate --root /data/my-rnaa-project
```

This recomputes QC verdicts from stored preprocess reports, updates run eligibility, and marks project-level normalize/corr outputs as pending when cohort membership changes.

`rnaa status` includes the current skipped count as `runs_skipped_qc`.

For DE and correlation outputs, ENSG-keyed files remain the primary computational outputs. Annotated sidecar outputs are also written where applicable:

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

Project configuration currently uses the `[refs].ensembl` key for the preset release selector. In the current implementation this selects the matching GENCODE preset release, with `"latest"` mapping to the latest supported preset for the chosen organism.

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

The workspace currently includes a vendored `fasterp` patch to keep the documented clippy gate passing under the workspace settings.

## License

MIT.
