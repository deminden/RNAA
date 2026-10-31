#!/usr/bin/env python3
"""Seal independently generated R output. Regeneration is an explicit review step."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PINS = {"R": "4.6.2", "tximport": "1.40.0", "DESeq2": "1.52.0",
        "WGCNA": "1.74", "HellCor": "1.3", "fgsea": "1.38.0"}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("family", choices=["tximport", "correlation", "de_workflows", "modules"])
    args = parser.parse_args()
    directory = ROOT / args.family
    with (directory / "versions.tsv").open() as stream:
        versions = {row["package"]: row["version"] for row in csv.DictReader(stream, delimiter="\t")}
    required = {"tximport": {"R", "tximport", "DESeq2"}, "correlation": {"R", "WGCNA", "HellCor"},
                "de_workflows": {"R", "DESeq2"}, "modules": {"R", "WGCNA", "fgsea"}}[args.family]
    if set(versions) != required:
        raise ValueError(f"Expected exactly these pinned oracle packages: {sorted(required)}")
    for name, version in versions.items():
        if PINS.get(name) != version:
            raise ValueError(f"Unexpected oracle version {name}={version}")
    sources = json.loads((ROOT / "oracle_sources.json").read_text())
    generator = ROOT / {"tximport": "generate_import.R", "correlation": "generate_correlation.R",
                        "de_workflows": "generate_de.R", "modules": "generate_modules.R"}[args.family]
    files = {}
    for path in sorted(directory.rglob("*")):
        if path.is_file() and path.name != "provenance.json":
            files[path.relative_to(directory).as_posix()] = {
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    report = {"schema_version": 1, "oracle": "independent R package execution",
              "versions": versions,
              "generator": {"path": generator.name,
                            "sha256": hashlib.sha256(generator.read_bytes()).hexdigest()},
              "sources": {name: sources[name] for name in versions}, "files": files,
              "tolerances": {"continuous": {"absolute": 1e-9, "relative": 1e-9},
                             "wald_iterative": {"absolute": 1e-7, "relative": 1e-7},
                             "rounded_counts": "exact", "undefined": "matching missingness"}}
    if args.family == "de_workflows":
        source = ROOT.parents[1] / "fixtures/native_workflow_counts.tsv"
        report["generator"]["inputs"] = {"fixtures/native_workflow_counts.tsv": {
            "sha256": hashlib.sha256(source.read_bytes()).hexdigest()}}
    (directory / "provenance.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Sealed {len(files)} files for {args.family}; versions {versions}")

if __name__ == "__main__":
    main()
