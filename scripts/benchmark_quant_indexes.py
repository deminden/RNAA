#!/usr/bin/env python3
"""Measure index sharing in separate processes, with reference build RSS excluded."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/examples/quant_index_scale")
    parser.add_argument("--transcripts", type=int, nargs="+", default=[1000, 2000])
    parser.add_argument("--length", type=int, default=600)
    parser.add_argument("--workers", type=int, nargs="+", default=[1, 4])
    parser.add_argument("--match-cache-mb", type=int, nargs="+", default=[64, 512])
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if platform.system() not in {"Darwin", "Linux"}:
        raise ValueError("RSS measurement requires macOS time -l or Linux time -v")
    binary = str(args.binary.resolve())
    rows = []
    with tempfile.TemporaryDirectory(prefix="rnaa-index-scale-") as work:
        for transcripts in args.transcripts:
            reference = str(Path(work) / str(transcripts))
            subprocess.run([binary, "build", reference, str(transcripts), str(args.length)],
                           cwd=ROOT, check=True, capture_output=True, text=True)
            hashes = set()
            for workers in args.workers:
                for cache_mb in args.match_cache_mb:
                    for mode in ["independent", "shared"]:
                        command = ["/usr/bin/time", "-l" if platform.system() == "Darwin" else "-v",
                                   binary, mode, reference, str(workers), str(cache_mb)]
                        run = subprocess.run(command, cwd=ROOT, check=True, capture_output=True, text=True)
                        result = json.loads(run.stdout)
                        if platform.system() == "Darwin":
                            match = re.search(r"(\d+)\s+maximum resident set size", run.stderr)
                            rss = int(match.group(1)) if match else None
                        else:
                            match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", run.stderr)
                            rss = int(match.group(1)) * 1024 if match else None
                        if rss is None:
                            raise ValueError(f"Missing process RSS measurement: {run.stderr}")
                        expected = 1 if mode == "shared" else workers
                        if result["successful_index_loads"] != expected or result["loaded_index_groups"] != expected:
                            raise ValueError(f"Unexpected index load count: {result}")
                        hashes.add(result["abundance_sha256"])
                        result.update(transcripts=transcripts, transcript_length=args.length,
                                      process_max_rss_bytes=rss, time_stderr=run.stderr)
                        rows.append(result)
                        print(f"{transcripts} transcripts, {workers} workers, {cache_mb} MiB/run, {mode}: RSS {rss} bytes", flush=True)
            if len(hashes) != 1:
                raise ValueError("Abundance outputs differ between sharing modes, cache allowances or worker counts")
    sources = [ROOT / "Cargo.lock", ROOT / "crates/rnaa-core/src/config.rs", ROOT / "crates/rnaa-stages/src/quant.rs",
               ROOT / "crates/rnaa-stages/src/quant/index_group.rs",
               ROOT / "crates/rnaa-stages/examples/quant_index_scale.rs", Path(__file__).resolve()]
    report = dict(schema_version=1, recorded_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  platform=platform.platform(), architecture=platform.machine(),
                  rustc=subprocess.check_output(["rustc", "-Vv"], text=True),
                  measurement="Separate process per case; index build excluded; loaded groups retained until all workers finish; RSS includes loading, validation, mapping, EM and outputs.",
                  scope="Synthetic index-sharing probe; serialized index bytes are not resident-memory estimates; not global admission or production/release qualification.",
                  source_sha256={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
                  runs=rows)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
