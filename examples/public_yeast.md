# Public yeast workflow

This small example processes four public *Saccharomyces cerevisiae* RNA-seq runs
from [PRJNA675088](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA675088), using
two controls and two treated samples. It exercises the installed CLI, ENA
resolution and download, reference preparation, quantification, DE, correlation,
restart, export and diagnosis. It is an execution check, not a biological
conclusion or a DESeq2 parity validation.

You need a built or installed `rnaa`, network access, Python 3 and `curl` for the
reference downloads, and several gigabytes of free local disk. The four compressed
read files total about 72 MiB; decoded and generated files need additional space.
Commands below use a new temporary directory and leave it in place for inspection.
Run from a workspace root after `cargo build --release --locked -p rnaa`, or replace
`rnaa_bin` with the absolute path of an installed executable.

```bash
set -eu
rnaa_bin="$PWD/target/release/rnaa"
work_dir="$(mktemp -d "${TMPDIR:-/tmp}/rnaa-yeast.XXXXXXXX")"
project_dir="$work_dir/project"
mkdir -p "$work_dir/reference"
curl --fail --location \
  https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/fungi/fasta/saccharomyces_cerevisiae/cdna/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa.gz \
  --output "$work_dir/reference/yeast.fa.gz"
curl --fail --location \
  https://ftp.ebi.ac.uk/ensemblgenomes/pub/release-63/fungi/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.63.gtf.gz \
  --output "$work_dir/reference/yeast.gtf.gz"
python3 - "$work_dir/reference" <<'PY'
import hashlib, pathlib, sys
root = pathlib.Path(sys.argv[1])
for name, expected in {
    "yeast.fa.gz": "3cc5b6bb7ff166cffac5c45317dad7052d06431f43f709ed14b33779702dc6a4",
    "yeast.gtf.gz": "253bd4013785b3c85a6a488bd541c12d8264aa4989046c3530d0d9ea2ecc931e",
}.items():
    actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
    if actual != expected:
        raise SystemExit(f"Reference checksum mismatch: {name}")
PY
cat > "$work_dir/metadata.tsv" <<'TSV'
run_accession	condition
SRR12999707	DMSO
SRR12999708	DMSO
SRR12999711	Cordycepin
SRR12999712	Cordycepin
TSV
"$rnaa_bin" --root "$project_dir" init
"$rnaa_bin" --root "$project_dir" add \
  --id SRR12999707 --id SRR12999708 --id SRR12999711 --id SRR12999712 \
  --metadata "$work_dir/metadata.tsv"
"$rnaa_bin" --root "$project_dir" resolve
"$rnaa_bin" --root "$project_dir" refs prepare \
  --organism saccharomyces_cerevisiae --ensembl 63 \
  --cdna "$work_dir/reference/yeast.fa.gz" --gtf "$work_dir/reference/yeast.gtf.gz"
"$rnaa_bin" --root "$project_dir" download --concurrency 2
"$rnaa_bin" --root "$project_dir" quant --worker-budget 2 --threads 2 --workers 1 --no-preprocess
"$rnaa_bin" --root "$project_dir" deseq2 --design '~ condition' --contrast condition DMSO Cordycepin
"$rnaa_bin" --root "$project_dir" corr --method spearman --adjust residualize \
  --model '~ condition' --geneset topvar:1000 --out edges:topk=10
"$rnaa_bin" --root "$project_dir" export --what all --format tsv --require-current
"$rnaa_bin" --root "$project_dir" survey --stratify condition --require-current
"$rnaa_bin" --root "$project_dir" status --json
"$rnaa_bin" --root "$project_dir" doctor --json
printf 'Project: %s\n' "$project_dir"
```

The survey summary should report four biological samples, with two in each
condition. Repeating the quantification, DE and correlation commands should reuse
their current generations when inputs and settings are unchanged. Input checksums
and the observed run are recorded in the [validation report](../validation/real-data-smoke/2026-10-31-yeast-release.md).
Review exported metadata and local paths before sharing a bundle.
