#!/usr/bin/env python3
"""Verify candidate archives and exercise both bundled and source-installed RNAA."""
import argparse
from contextlib import ExitStack
import datetime
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import tarfile
import tempfile
import tomllib

PRIVATE_PATH = re.compile(rb"/(?:Users|home)/[A-Za-z0-9._-]+/|/(?:private)/var/folders/")


def digest(path):
    hasher = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def checked_path(root, name):
    relative = Path(name)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError(f"Unsafe candidate path: {name}")
    return root / relative


def extract(archive, destination, expected_root):
    with tarfile.open(archive) as source:
        names = set()
        for member in source.getmembers():
            checked_path(destination, member.name)
            parts = Path(member.name).parts
            if not parts or parts[0] != expected_root:
                raise ValueError(f"Unexpected archive root: {member.name}")
            normalized = str(Path(member.name))
            if normalized in names:
                raise ValueError(f"Duplicate archive entry: {member.name}")
            names.add(normalized)
            if not member.isfile() and not member.isdir():
                raise ValueError(f"Non-regular archive entry: {member.name}")
        source.extractall(destination, filter="data")


def verify_binary_inventory(binary_root, expected_files, dependencies):
    bundled_source = {
        name for name in expected_files
        if name in ("LICENSE", "README.md", "CHANGELOG.md", "vendor/fasterp/RNAA_PATCHES.md")
        or name.startswith(("docs/", "configs/", "examples/", "validation/"))
    }
    notices = set()
    for component in dependencies["components"]:
        for name in component["license_files"]:
            if not name.startswith("licenses/"):
                raise ValueError(f"Unsafe bundled notice path: {name}")
            checked_path(binary_root, name)
            notices.add(name)
    expected = bundled_source | notices | {
        "bin/rnaa", "RELEASE-MANIFEST.json", "DEPENDENCIES.json"
    }
    actual = {str(path.relative_to(binary_root)) for path in binary_root.rglob("*")
              if path.is_file()}
    if actual != expected:
        raise ValueError("Bundled file inventory differs from source and dependency manifests")
    for name in bundled_source:
        if digest(checked_path(binary_root, name)) != expected_files[name]:
            raise ValueError(f"Bundled source file hash mismatch: {name}")


def fresh_install_target(parent, stack):
    if parent is not None:
        parent.mkdir(parents=True, exist_ok=True)
    return Path(stack.enter_context(tempfile.TemporaryDirectory(
        prefix="rnaa-install-target-", dir=parent)))


def source_snapshot_key(root, expected_files):
    inventory = {name: {"sha256": expected, "mode": checked_path(root, name).stat().st_mode & 0o777}
                 for name, expected in expected_files.items()}
    encoded = (json.dumps(inventory, indent=2, sort_keys=True) + "\n").encode()
    return hashlib.sha256(encoded).hexdigest()


def verify(candidate, install_target):
    manifest = json.loads((candidate / "RELEASE-MANIFEST.json").read_text())
    if manifest["schema_version"] != 1:
        raise ValueError("Unsupported release manifest schema")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.-]{0,79}", manifest["label"]):
        raise ValueError("Invalid candidate label")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", manifest["target"]):
        raise ValueError("Invalid candidate target")
    checksums = {}
    for line in (candidate / "SHA256SUMS").read_text().splitlines():
        expected, name = line.split("  ", 1)
        if name in checksums or digest(checked_path(candidate, name)) != expected:
            raise ValueError(f"Checksum mismatch or duplicate entry: {name}")
        checksums[name] = expected
    required = {"RELEASE-MANIFEST.json", "DEPENDENCIES.json", "source.tar.gz",
                f'{manifest["target"]}.tar.gz'}
    if set(checksums) != required:
        raise ValueError("Candidate checksum inventory is incomplete or unexpected")
    with tempfile.TemporaryDirectory(prefix="rnaa-candidate-check-") as temporary, ExitStack() as stack:
        work = Path(temporary)
        source_name = f'rnaa-{manifest["label"]}-source'
        binary_name = f'rnaa-{manifest["label"]}-{manifest["target"]}'
        extract(candidate / "source.tar.gz", work, source_name)
        extract(candidate / f'{manifest["target"]}.tar.gz', work, binary_name)
        source = work / source_name
        binary_root = work / binary_name
        binary = binary_root / "bin/rnaa"
        expected_files = manifest["source_sha256"]
        actual_files = {str(path.relative_to(source)) for path in source.rglob("*") if path.is_file()}
        if actual_files != set(expected_files):
            raise ValueError("Archived source file inventory differs from release manifest")
        for name, expected in expected_files.items():
            if digest(checked_path(source, name)) != expected:
                raise ValueError(f"Archived source hash mismatch: {name}")
        if "source_snapshot_sha256" in manifest and source_snapshot_key(source, expected_files) != manifest["source_snapshot_sha256"]:
            raise ValueError("Archived source snapshot identity mismatch")
        for name in expected_files:
            if PRIVATE_PATH.search(checked_path(source, name).read_bytes()):
                raise ValueError(f"Archived source contains a machine-specific path: {name}")
        for name in manifest.get("omitted_historical_validation_files", []):
            if not name.startswith("validation/") or name in expected_files:
                raise ValueError(f"Invalid omitted historical validation file: {name}")
        if PRIVATE_PATH.search(binary.read_bytes()):
            raise ValueError("Archived executable contains a machine-specific home or temporary path")
        if digest(binary) != manifest["binary_sha256"]:
            raise ValueError("Archived executable hash mismatch")
        for name in ("RELEASE-MANIFEST.json", "DEPENDENCIES.json"):
            if (binary_root / name).read_bytes() != (candidate / name).read_bytes():
                raise ValueError(f"Bundled {name} differs from candidate")
        dependencies = json.loads((candidate / "DEPENDENCIES.json").read_text())
        verify_binary_inventory(binary_root, expected_files, dependencies)
        if digest(source / "Cargo.lock") != manifest["lockfile_sha256"]:
            raise ValueError("Archived lockfile hash mismatch")
        version = tomllib.loads((source / "Cargo.toml").read_text())["workspace"]["package"]["version"]
        expected_version = f"rnaa {version}"
        if manifest["executable_version"] != expected_version:
            raise ValueError("Manifest and source workspace versions disagree")
        runtime = manifest["runtime_requirements"]
        if runtime["os"] != platform.system():
            raise ValueError(f'Run this candidate on {runtime["os"]} to verify it')
        if runtime.get("libc") == "glibc":
            name, version = platform.libc_ver()
            required_glibc = tuple(map(int, runtime["minimum_glibc"].split(".")))
            if name != "glibc" or tuple(map(int, version.split("."))) < required_glibc:
                raise ValueError(f'Candidate requires glibc {runtime["minimum_glibc"]} or newer')
        example = source / "examples/offline_reference.py"
        build_target = fresh_install_target(install_target, stack)
        for label, executable in [("bundled", binary), ("installed", work / "install/bin/rnaa")]:
            if label == "installed":
                command = ["cargo", "install", "--locked", "--all-features", "--path", str(source / "crates/rnaa"),
                           "--root", str(work / "install")]
                command.extend(["--target-dir", str(build_target)])
                subprocess.run(command, check=True)
            output = subprocess.check_output([str(executable), "--version"], text=True).strip()
            if output != expected_version:
                raise ValueError(f"{label} executable version mismatch")
            subprocess.run([sys.executable, str(example), "--binary", str(executable)], check=True)
            subprocess.run([sys.executable, str(example), "--binary", str(executable),
                            "--automatic-resources"], check=True)
        installed_hash = digest(work / "install/bin/rnaa")
    return {"schema_version": 1,
            "recorded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            "platform": platform.platform(), "machine": platform.machine(),
            "label": manifest["label"], "target": manifest["target"],
            "executable_version": expected_version, "source_files_verified": len(expected_files),
            "source_snapshot_identity_verified": "source_snapshot_sha256" in manifest,
            "local_path_scan_passed": True,
            "checksums": checksums, "bundled_binary_sha256": manifest["binary_sha256"],
            "runtime_requirements": runtime, "bundled_operator_guides_verified": True,
            "bundled_source_inventory_verified": True,
            "source_installed_binary_sha256": installed_hash,
            "source_install_target": (f"fresh isolated target under {install_target}" if install_target
                                      else "fresh temporary target"),
            "bundled_smoke_passed": True, "source_install_passed": True,
            "source_installed_smoke_passed": True, "automatic_resources_smoke_passed": True,
            "passed": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--install-target-dir", type=Path,
                        help="Parent directory for an isolated fresh build target")
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    target = args.install_target_dir.resolve() if args.install_target_dir else None
    report = verify(args.candidate.resolve(strict=True), target)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(f'PASS {report["target"]}: archive hashes, {report["source_files_verified"]} source files, '
          'bundled executable, source installation and both offline checks')


if __name__ == "__main__":
    main()
