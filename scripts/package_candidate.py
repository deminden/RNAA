#!/usr/bin/env python3
"""Prepare local review archives; never tag, publish or change source versions."""
import argparse
import datetime
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import tarfile
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parents[1]
PRIVATE_PATH = re.compile(rb"/(?:Users|home)/[A-Za-z0-9._-]+/|/(?:private)/var/folders/")


def command(*args, cwd=ROOT):
    return subprocess.check_output(args, cwd=cwd, text=True)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def json_bytes(document):
    return (json.dumps(document, indent=2, sort_keys=True) + "\n").encode()


def snapshot_key(sources):
    inventory = {name: {"sha256": sha(data), "mode": mode}
                 for name, (data, mode) in sources.items()}
    return sha(json_bytes(inventory))


def runtime_requirements(binary):
    requirements = {"os": platform.system()}
    if platform.system() == "Linux":
        symbols = command("readelf", "--version-info", str(binary))
        versions = re.findall(r"Name: GLIBC_([0-9.]+)", symbols)
        if not versions:
            raise ValueError("Expected glibc symbol versions for the native Linux GNU executable")
        requirements["libc"] = "glibc"
        requirements["minimum_glibc"] = max(versions, key=lambda v: tuple(map(int, v.split("."))))
        requirements["dynamic_dependencies"] = sorted({line.strip().split()[0]
                                                        for line in command("ldd", str(binary)).splitlines()
                                                        if line.strip()})
    return requirements


# Keep scientific fixtures and the small set of evidence cited by the guides.
# Complete historical logs and prior candidate records remain in the repository.
VALIDATION_GUIDES = {
    "validation/correlation-scale/README.md",
    "validation/correlation-scale/2026-10-30-debug.json",
    "validation/correlation-scale/2026-10-30-modules-debug.json",
    "validation/correlation-scale/2026-10-31-network-cancellation-debug.json",
    "validation/correlation-scale/2026-10-31-selection.json",
    "validation/dependencies/2026-10-31-release-wrap.json",
    "validation/dependencies/2026-11-04-0.3.1.json",
    "validation/checks/2026-10-31-release-wrap.json",
    "validation/real-data-smoke/2026-10-31-yeast-release.md",
    "validation/real-data-smoke/2026-10-31-yeast-release.json",
    "validation/real-data-smoke/2026-10-31-yeast-linux.json",
    "validation/real-data-smoke/manifest.json",
}


def source_snapshot(omitted=None):
    names = subprocess.check_output(
        ["git", "ls-files", "-z", "--cached", "--others", "--exclude-standard"], cwd=ROOT)
    snapshot = {}
    for name in sorted(set(names.decode().split("\0")) - {""}):
        if name.startswith("validation/") and not name.startswith(("validation/scientific/", "validation/native-fixture/")) and name not in VALIDATION_GUIDES:
            if omitted is not None:
                omitted.append(name)
            continue
        path = ROOT / name
        if path.is_symlink():
            raise ValueError(f"Review source symlink before packaging: {name}")
        if path.is_file():
            data = path.read_bytes()
            if PRIVATE_PATH.search(data):
                raise ValueError(f"Source contains a machine-specific path: {name}")
            snapshot[name] = (data, 0o755 if path.stat().st_mode & 0o111 else 0o644)
    for required in ["Cargo.lock", "vendor/fasterp/Cargo.toml", "LICENSE"]:
        if required not in snapshot:
            raise ValueError(f"Source snapshot is missing {required}")
    return snapshot


def inventory(host, root=ROOT):
    metadata = json.loads(command("cargo", "metadata", "--locked", "--all-features",
                                  "--format-version", "1", "--filter-platform", host, cwd=root))
    packages = {p["id"]: p for p in metadata["packages"]}
    nodes = {p["id"]: p for p in metadata["resolve"]["nodes"]}
    pending = [p["id"] for p in packages.values() if p["name"] == "rnaa"]
    seen, components, notices = set(), [], {}
    while pending:
        identity = pending.pop()
        if identity in seen:
            continue
        seen.add(identity)
        package = packages[identity]
        base = Path(package["manifest_path"]).parent
        files = [p for p in base.iterdir() if p.is_file() and
                 p.name.lower().startswith(("license", "copying", "notice", "copyright"))]
        if package["license_file"]:
            files.append(base / package["license_file"])
        if identity in metadata["workspace_members"]:
            files.append(root / "LICENSE")
        # The pinned Git dependency keeps its license at the repository root.
        if package["name"] == "mincorr":
            files.append(base.parent.parent / "LICENSE")
        stem = f'{package["name"]}-{package["version"]}'
        notice_names = []
        for path in sorted(set(files)):
            if not path.is_file():
                raise ValueError(f"Missing declared license file: {path}")
            name = f"licenses/{stem}/{path.name}"
            notices[name] = (path.read_bytes(), 0o644)
            notice_names.append(name)
        components.append({"name": package["name"], "version": package["version"],
                           "source": package["source"] or "workspace/local",
                           "license": package["license"], "repository": package["repository"],
                           "license_files": sorted(set(notice_names))})
        pending.extend(dep["pkg"] for dep in nodes[identity]["deps"]
                       if any(kind["kind"] in (None, "build") for kind in dep["dep_kinds"]))
    return sorted(components, key=lambda p: (p["name"], p["version"])), notices


def archive(path, prefix, entries):
    with path.open("xb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as gz:
        with tarfile.open(fileobj=gz, mode="w|") as tar:
            for name, (data, mode) in sorted(entries.items()):
                relative = Path(name)
                if relative.is_absolute() or ".." in relative.parts:
                    raise ValueError(f"Unsafe archive entry: {name}")
                entry = tarfile.TarInfo(f"{prefix}/{name}")
                entry.size, entry.mode, entry.mtime = len(data), mode, 0
                tar.addfile(entry, io.BytesIO(data))


def materialize_source(root, sources):
    for name, (data, mode) in sources.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        path.chmod(mode)


def verify_frozen_source(root, sources):
    actual = {str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()}
    if actual != set(sources):
        raise ValueError("Frozen source inventory changed during candidate build")
    for name, (data, mode) in sources.items():
        path = root / name
        if path.is_symlink() or path.read_bytes() != data or (path.stat().st_mode & 0o777) != mode:
            raise ValueError(f"Frozen source changed during candidate build: {name}")


def build_from_snapshot(root, sources, host, target_dir):
    materialize_source(root, sources)
    env = os.environ.copy()
    flags = env.get("CARGO_ENCODED_RUSTFLAGS", "").split("\x1f") if env.get("CARGO_ENCODED_RUSTFLAGS") else []
    for source, replacement in [(str(Path.home()), "/rnaa-build-user"),
                                (str(Path(tempfile.gettempdir()).resolve()), "/rnaa-build-tmp"),
                                (str(root), "/rnaa-build-source")]:
        flags.append(f"--remap-path-prefix={source}={replacement}")
    env["CARGO_ENCODED_RUSTFLAGS"] = "\x1f".join(flags)
    subprocess.run(["cargo", "build", "--release", "--all-features", "--locked", "-p", "rnaa",
                    "--target", host, "--target-dir", str(target_dir)], cwd=root, check=True, env=env)
    verify_frozen_source(root, sources)
    return target_dir / host / "release" / ("rnaa.exe" if os.name == "nt" else "rnaa")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--label", required=True)
    parser.add_argument("--build-target-dir", type=Path, default=ROOT / "target/candidate-build",
                        help="Compiler cache for the build from frozen source")
    parser.add_argument("--allow-dirty", action="store_true")
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.-]{0,79}", args.label):
        raise ValueError("Label must contain at most 80 letters, digits, dots or hyphens")
    status = command("git", "status", "--porcelain=v1")
    if status and not args.allow_dirty:
        raise ValueError("Review the dirty source tree and explicitly use --allow-dirty")
    host = next(line.split(": ", 1)[1] for line in command("rustc", "-vV").splitlines()
                if line.startswith("host: "))
    omitted = []
    sources = source_snapshot(omitted)
    source_key = snapshot_key(sources)
    expected = tomllib.loads(sources["Cargo.toml"][0].decode())["workspace"]["package"]["version"]
    with tempfile.TemporaryDirectory(prefix="rnaa-candidate-build-") as temporary:
        frozen = Path(temporary) / "source"
        frozen.mkdir()
        binary = build_from_snapshot(frozen, sources, host,
                                     args.build_target_dir.resolve() / f"snapshot-{source_key}")
        version = command(str(binary), "--version").strip()
        if version != f"rnaa {expected}":
            raise ValueError(f"Executable version {version!r} does not match workspace {expected}")
        components, notices = inventory(host, frozen)
        verify_frozen_source(frozen, sources)
        binary_data = binary.read_bytes()
        if PRIVATE_PATH.search(binary_data):
            raise ValueError("Built executable contains a machine-specific home or temporary path")
        runtime = runtime_requirements(binary)
    source_hashes = {name: sha(value[0]) for name, value in sources.items()}
    report = {"schema_version": 1, "prepared_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "label": args.label, "executable_version": version, "target": host,
              "git_head": command("git", "rev-parse", "HEAD").strip(), "dirty": bool(status),
              "source_sha256": source_hashes, "binary_sha256": sha(binary_data),
              "source_snapshot_sha256": source_key,
              "omitted_historical_validation_files": omitted,
              "binary_build": "cargo build --release --all-features --locked from frozen source",
              "lockfile_sha256": sha(sources["Cargo.lock"][0]),
              "rustc": command("rustc", "--version").strip(),
              "runtime_requirements": runtime,
              "distribution_approved": False,
              "deferred_work": ["mincorr licensing review (deferred by user)",
                                "Production-scale scientific and operational qualification remains incomplete"],
              "public_release_blockers": ["Verify platform and installation evidence for this exact source snapshot",
                                          "Verify a clean final commit and approve publication"],
              "components_without_root_license_files": [f'{p["name"]} {p["version"]}'
                                                        for p in components if not p["license_files"]]}
    output = ROOT / "target/release-candidates" / args.label
    output.mkdir(parents=True, exist_ok=False)
    manifest = json_bytes(report)
    dependencies = json_bytes({"scope": "Host normal/build dependency closure; not a compatibility verdict",
                               "target": host, "components": components})
    bundle = {name: value for name, value in sources.items()
              if name in ("LICENSE", "README.md", "CHANGELOG.md", "vendor/fasterp/RNAA_PATCHES.md")
              or name.startswith(("docs/", "configs/", "examples/", "validation/"))}
    bundle.update(notices)
    bundle.update({"bin/rnaa": (binary_data, 0o755),
                   "RELEASE-MANIFEST.json": (manifest, 0o644),
                   "DEPENDENCIES.json": (dependencies, 0o644)})
    archive(output / "source.tar.gz", f"rnaa-{args.label}-source", sources)
    archive(output / f"{host}.tar.gz", f"rnaa-{args.label}-{host}", bundle)
    (output / "RELEASE-MANIFEST.json").write_bytes(manifest)
    (output / "DEPENDENCIES.json").write_bytes(dependencies)
    sums = [f"{sha(path.read_bytes())}  {path.name}" for path in sorted(output.iterdir())]
    (output / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    print(f"Prepared local candidate: {output}")
    print(f"{len(sources)} source files; {len(components)} dependency components; distribution_approved=false")


if __name__ == "__main__":
    main()
