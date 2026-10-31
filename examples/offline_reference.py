#!/usr/bin/env python3
"""Exercise an installed RNAA without network access or existing project data."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import tomllib


def exercise(binary, root, automatic_resources=False):
    def run(*arguments):
        result = subprocess.run([str(binary), "--root", str(root), *arguments],
                                check=True, capture_output=True, text=True, timeout=120)
        return result.stdout

    version = run("--version").strip()
    run("init")
    fasta, gtf = root / "source.fa", root / "source.gtf"
    fasta.write_text(">tx1.1\n" + "ACGT" * 24 + "\n>tx2.1\n" + "TGCA" * 24 + "\n")
    gtf.write_text(
        'chr1\tsynthetic\ttranscript\t1\t96\t.\t+\t.\tgene_id "G1"; transcript_id "tx1.1"; gene_name "one";\n'
        'chr1\tsynthetic\ttranscript\t1\t96\t.\t+\t.\tgene_id "G2"; transcript_id "tx2.1"; gene_name "two";\n')
    settings = {("refs", "custom_cdna"): str(fasta),
                ("refs", "custom_gtf"): str(gtf), ("refs", "index_threads"): 1,
                ("resources", "pool_root"): str(root / "private-pool"),
                ("resources", "cpu_threads"): 0 if automatic_resources else 1,
                ("resources", "memory_mb"): 0 if automatic_resources else 256,
                ("resources", "disk_headroom_mb"): 0}
    config = root / "rnaa.toml"
    section, observed, lines = "", set(), []
    for line in config.read_text().splitlines():
        if line.startswith("["):
            section = line.strip("[]")
        key = section, line.split("=", 1)[0].strip()
        if key in settings:
            line = f"{key[1]} = {json.dumps(settings[key])}"
            observed.add(key)
        lines.append(line)
    if observed != set(settings):
        raise ValueError("Generated configuration no longer matches the example")
    config.write_text("\n".join(lines) + "\n")
    tomllib.loads(config.read_text())
    first = run("refs", "prepare")
    manifest = Path(next(line.split("\t", 1)[1] for line in first.splitlines()
                         if line.startswith("manifest\t")))
    limits = json.loads(manifest.read_text())["parameters"]["resource_execution"]["build_reservation"]["limits"]
    if min(limits["cpu_threads"], limits["memory_bytes"]) <= 0:
        raise ValueError("Invalid detected resource capacity")
    digest = hashlib.sha256(manifest.read_bytes()).hexdigest()
    second = run("refs", "prepare")
    repeated = Path(next(line.split("\t", 1)[1] for line in second.splitlines()
                        if line.startswith("manifest\t")))
    if repeated != manifest or hashlib.sha256(repeated.read_bytes()).hexdigest() != digest:
        raise ValueError("Repeated preparation did not reuse the reference")
    if "reference_ready\tyes" not in run("status"):
        raise ValueError("Prepared reference is not ready")
    doctor = run("doctor")
    if "analysis_engines\tOK\t" not in doctor:
        raise ValueError("Installed engines were not reported")
    mode = "automatic" if automatic_resources else "explicit"
    print(f"PASS {version}: initialize, prepare two transcripts, reuse, status, doctor; {mode} resources")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=Path("target/release/rnaa"))
    parser.add_argument("--automatic-resources", action="store_true",
                        help="Exercise host CPU/memory detection with a private resource pool")
    args = parser.parse_args()
    binary = args.binary.resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="rnaa-offline-example-") as temporary:
        exercise(binary, Path(temporary), args.automatic_resources)


if __name__ == "__main__":
    main()
