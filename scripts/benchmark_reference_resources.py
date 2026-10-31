#!/usr/bin/env python3
"""Profile complete admitted reference preparation in independent CLI processes."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import sqlite3
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fixture(root, transcripts, length):
    fasta, gtf = root / "source.fa", root / "source.gtf"
    random = 17
    sequence = ""
    with fasta.open("w") as fa, gtf.open("w") as annotation:
        for transcript in range(transcripts):
            gene, isoform = divmod(transcript, 4)
            if isoform == 0:
                bases = []
                for _ in range(length + 120):
                    random = (random * 6364136223846793005 + 1) & ((1 << 64) - 1)
                    bases.append("ACGT"[(random >> 32) & 3])
                sequence = "".join(bases)
            read = sequence[isoform * 40:isoform * 40 + length]
            fa.write(f">tx{transcript}.1\n")
            for start in range(0, len(read), 60):
                fa.write(read[start:start + 60] + "\n")
            annotation.write(
                f'chr1\tsynthetic\ttranscript\t1\t{length}\t.\t+\t.\t'
                f'gene_id "G{gene}"; transcript_id "tx{transcript}.1"; gene_name "Gene{gene}";\n')
    return fasta, gtf


def configure(path, settings):
    section, observed, output = "", set(), []
    for line in path.read_text().splitlines():
        if line.startswith("["):
            section = line.strip("[]")
        key = line.split("=", 1)[0].strip()
        if (section, key) in settings:
            line = f"{key} = {json.dumps(settings[section, key])}"
            observed.add((section, key))
        output.append(line)
    if observed != set(settings):
        raise ValueError(f"Unrecognized configuration fields: {set(settings) - observed}")
    path.write_text("\n".join(output) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/rnaa")
    parser.add_argument("--transcripts", type=int, nargs="+", default=[1000, 20000])
    parser.add_argument("--length", type=int, default=600)
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 2])
    parser.add_argument("--memory-mb", type=int, default=4096)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if platform.system() not in {"Darwin", "Linux"}:
        raise ValueError("Requires macOS time -l or Linux time -v")
    if args.length < 100 or min(args.transcripts + args.threads + [args.memory_mb]) < 1:
        raise ValueError("Require positive counts and sequence length >= 100")
    binary, rows = str(args.binary.resolve()), []
    with tempfile.TemporaryDirectory(prefix="rnaa-reference-profile-") as temporary:
        work = Path(temporary)
        for transcripts in args.transcripts:
            source = work / f"source-{transcripts}"
            source.mkdir()
            fasta, gtf = fixture(source, transcripts, args.length)
            science = set()
            for threads in args.threads:
                project = work / f"project-{transcripts}-{threads}"
                project.mkdir()
                subprocess.run([binary, "--root", str(project), "init"], cwd=ROOT,
                               check=True, capture_output=True, text=True)
                configure(project / "rnaa.toml", {
                    ("refs", "custom_cdna"): str(fasta), ("refs", "custom_gtf"): str(gtf),
                    ("refs", "index_threads"): threads,
                    ("resources", "pool_root"): str(project / "pool"),
                    ("resources", "cpu_threads"): threads,
                    ("resources", "memory_mb"): args.memory_mb,
                    ("resources", "disk_headroom_mb"): 0,
                })
                scratch = project / "scratch"
                scratch.mkdir()
                command = ["/usr/bin/time", "-l" if platform.system() == "Darwin" else "-v",
                           binary, "--root", str(project), "refs", "prepare"]
                started = time.monotonic()
                run = subprocess.run(command, cwd=ROOT, env={**os.environ, "TMPDIR": str(scratch)},
                                     check=True, capture_output=True, text=True, timeout=1800)
                elapsed = time.monotonic() - started
                manifest_path = next(line.split("\t", 1)[1] for line in run.stdout.splitlines()
                                     if line.startswith("manifest\t"))
                manifest = json.loads(Path(manifest_path).read_text())
                match = re.search(r"(\d+)\s+maximum resident set size", run.stderr) if platform.system() == "Darwin" else re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", run.stderr)
                if match is None:
                    raise ValueError(f"Missing process RSS: {run.stderr}")
                rss = int(match.group(1)) * (1 if platform.system() == "Darwin" else 1024)
                parameters = manifest["parameters"]
                resource = parameters["resource_execution"]
                claim = resource["build_reservation"]["requested"]
                outputs = {a["kind"]: a for a in manifest["output_artifacts"]}
                if parameters["mapping"]["mapped_transcripts"] != transcripts:
                    raise ValueError("Unexpected mapped transcript count")
                science.add((parameters["request_sha256"], outputs["REFERENCE_TX2GENE"]["checksum"],
                             outputs["REFERENCE_GENE_ANNOTATION"]["checksum"], resource["build_report"]["kmers"]))
                pool_catalog = (project / "pool" / "reservations.sqlite").resolve()
                if not pool_catalog.is_file():
                    raise ValueError(f"Missing private pool catalog {pool_catalog}; observed {sorted(str(p.relative_to(project)) for p in project.rglob('*.sqlite'))}")
                try:
                    with sqlite3.connect(pool_catalog.as_uri() + "?mode=ro", uri=True) as db:
                        if db.execute("SELECT COUNT(*) FROM reservations").fetchone()[0] != 0:
                            raise ValueError("Reference process left resource claims")
                except sqlite3.OperationalError as error:
                    raise RuntimeError(f"Cannot read the finished pool with SQLite {sqlite3.sqlite_version}; retry with a Python 3.11+ runtime containing a current SQLite library") from error
                if list(scratch.iterdir()):
                    raise ValueError("Reference process left planning/staging files")
                if list(project.glob(".rnaa_generations/*/payload/.rnaa_reference_work")):
                    raise ValueError("Reference process retained its private generation workspace")
                row = dict(transcripts=transcripts, sequence_length=args.length, requested_threads=threads,
                           wall_seconds=elapsed, process_peak_rss_bytes=rss,
                           build_memory_claim_bytes=claim["memory_bytes"],
                           process_peak_within_build_claim=rss <= claim["memory_bytes"],
                           source_sha256={"fasta": sha(fasta), "gtf": sha(gtf)},
                           request_sha256=parameters["request_sha256"], estimate=resource["estimate"],
                           build_report=resource["build_report"],
                           remaining_resource_rows=0, remaining_scratch_files=0,
                           output_bytes={k: v["bytes"] for k, v in outputs.items()}, time_stderr=run.stderr)
                rows.append(row)
                print(f"{transcripts} transcripts, {threads} threads: {elapsed:.3f}s, peak RSS {rss} bytes, build claim {claim['memory_bytes']} bytes", flush=True)
            if len(science) != 1:
                raise ValueError("Scientific requests, annotation tables or distinct graph counts changed with threads")
    sources = [ROOT / "Cargo.lock", ROOT / "crates/rnaa-core/src/config.rs", ROOT / "crates/rnaa-core/src/util.rs",
               *sorted((ROOT / "crates/rnaa-stages/src/refs").rglob("*.rs")), Path(__file__).resolve()]
    report = dict(schema_version=1, recorded_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  platform=platform.platform(), architecture=platform.machine(),
                  rustc=subprocess.check_output(["rustc", "-Vv"], text=True), python=sys.version,
                  sqlite=sqlite3.sqlite_version, binary_sha256=sha(binary),
                  measurement="Separate CLI process and empty private pool per case. Source generation and init excluded; RSS includes discovery, acquisition, external planning, native building, validation and publication.",
                  scope="Deterministic overlapping synthetic references. Single local host; not production-cohort or supported-host release qualification. Process RSS includes baseline/coordinator overhead; cooperative phase claims are not OS memory ceilings.",
                  source_sha256={str(p.relative_to(ROOT)): sha(p) for p in sources}, runs=rows)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
