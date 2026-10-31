#!/usr/bin/env python3
"""Record deterministic correlation scale probes with process RSS and source hashes."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, default=ROOT / "target/release/examples/correlation_scale")
    parser.add_argument("--genes", type=int, nargs="+", default=[2000, 4000])
    parser.add_argument("--samples", type=int, default=16)
    parser.add_argument("--k", type=int, default=10)
    parser.add_argument("--block-rows", type=int, default=64)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--method", default="pearson")
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    if platform.system() not in {"Darwin", "Linux"}:
        raise ValueError("RSS measurement requires macOS time -l or Linux time -v")
    rows = []
    for genes in args.genes:
        command = ["/usr/bin/time", "-l" if platform.system() == "Darwin" else "-v", str(args.binary.resolve()),
                   str(genes), str(args.samples), str(args.k), str(args.block_rows), str(args.threads), args.method]
        run = subprocess.run(command, cwd=ROOT, text=True, capture_output=True, check=True)
        result = json.loads(run.stdout)
        if platform.system() == "Darwin":
            match = re.search(r"(\d+)\s+maximum resident set size", run.stderr)
            rss = int(match.group(1)) if match else None
        else:
            match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", run.stderr)
            rss = int(match.group(1)) * 1024 if match else None
        if rss is None:
            raise ValueError(f"Missing process RSS measurement: {run.stderr}")
        result["process_max_rss_bytes"] = rss
        result["time_stderr"] = run.stderr
        result["hypothetical_single_dense_f64_bytes"] = genes * genes * 8
        rows.append(result)
        seconds = result.get("correlation_wall_seconds", result.get("module_wall_seconds"))
        print(f"{genes} genes: {seconds:.3f}s, RSS {rss} bytes", flush=True)
    sources = [ROOT / "Cargo.lock", ROOT / "crates/rnaa-core/src/config.rs",
               ROOT / "crates/rnaa-stages/examples/correlation_scale.rs",
               ROOT / "crates/rnaa-stages/examples/module_scale.rs",
               ROOT / "crates/rnaa-stages/examples/support/mod.rs"]
    sources.extend(sorted((ROOT / "crates/rnaa-stages/src/corr").rglob("*.rs")))
    report = {"schema_version": 1, "recorded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "platform": platform.platform(), "architecture": platform.machine(),
              "rustc": subprocess.check_output(["rustc", "-Vv"], text=True),
              "measurement": "separate process per size; RSS includes startup and input generation",
              "scope": "deterministic synthetic network/module probe; not a production cohort or supported-host release qualification",
              "source_sha256": {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
              "runs": rows}
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
