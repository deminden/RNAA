#!/usr/bin/env python3
"""Assemble verified platform candidates into uniquely named release assets."""
import hashlib
import json
from pathlib import Path
import sys
import tarfile
import tomllib

PLATFORMS = {
    "ubuntu-24.04": "x86_64-unknown-linux-gnu",
    "ubuntu-24.04-arm": "aarch64-unknown-linux-gnu",
    "macos-15": "aarch64-apple-darwin",
}


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def assemble(inputs, output, tag, commit, workspace_version):
    if tag not in (workspace_version, f"v{workspace_version}"):
        raise ValueError("Release tag must match the workspace version, optionally prefixed with v")
    expected_dirs = {f"rnaa-candidate-{os}-rust-stable" for os in PLATFORMS}
    actual_dirs = {p.name for p in inputs.iterdir() if p.is_dir()}
    if actual_dirs != expected_dirs or any(not p.is_dir() for p in inputs.iterdir()):
        raise ValueError("Release candidate platforms are missing or unexpected")
    candidates = []
    snapshots = set()
    for os_name, target in PLATFORMS.items():
        candidate = inputs / f"rnaa-candidate-{os_name}-rust-stable"
        manifest = json.loads((candidate / "RELEASE-MANIFEST.json").read_text())
        if (manifest["target"] != target or manifest["git_head"] != commit
                or manifest["dirty"] or manifest["executable_version"] != f"rnaa {workspace_version}"):
            raise ValueError(f"Candidate identity does not match the release: {os_name}")
        checksums = {}
        for line in (candidate / "SHA256SUMS").read_text().splitlines():
            expected, name = line.split("  ", 1)
            if Path(name).name != name or name in checksums:
                raise ValueError(f"Invalid candidate checksum entry: {os_name}")
            checksums[name] = expected
        required = {"RELEASE-MANIFEST.json", "DEPENDENCIES.json", "source.tar.gz", f"{target}.tar.gz"}
        if set(checksums) != required or {p.name for p in candidate.iterdir()} != required | {"SHA256SUMS"}:
            raise ValueError(f"Invalid candidate inventory: {os_name}")
        if any(digest(candidate / name) != expected for name, expected in checksums.items()):
            raise ValueError(f"Candidate checksum mismatch: {os_name}")
        snapshots.add(manifest["source_snapshot_sha256"])
        candidates.append((target, candidate))
    if len(snapshots) != 1:
        raise ValueError("Platform candidates do not share the same source snapshot")
    output.mkdir(parents=True, exist_ok=False)
    assets = []
    for target, candidate in candidates:
        asset = output / f"rnaa-{tag}-{target}.tar.gz"
        with tarfile.open(asset, "w:gz") as archive:
            archive.add(candidate, arcname=f"rnaa-{tag}-{target}")
        assets.append(asset)
    (output / "SHA256SUMS").write_text("".join(f"{digest(asset)}  {asset.name}\n" for asset in assets))


def main():
    if len(sys.argv) != 5:
        raise SystemExit("usage: assemble_release.py INPUTS OUTPUT TAG COMMIT")
    version = tomllib.loads(Path("Cargo.toml").read_text())["workspace"]["package"]["version"]
    assemble(Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3], sys.argv[4], version)


if __name__ == "__main__":
    main()
