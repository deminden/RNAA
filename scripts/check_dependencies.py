#!/usr/bin/env python3
"""Audit direct registry dependencies against non-yanked stable crates.io releases.

Python 3.11+, standard library only. Reads both maintained lockfiles; never edits
manifests or resolves dependencies. Network failures are failures, not fresh results.
"""

import argparse
import concurrent.futures
import datetime
import hashlib
import json
from pathlib import Path
import sys
import tomllib
import urllib.request


ROOT = Path(__file__).resolve().parents[1]
USER_AGENT = "RNAA-dependency-audit (https://github.com/deminden/RNAA)"


def read_json(url):
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def stable_release(name):
    url = f"https://crates.io/api/v1/crates/{name}"
    try:
        response = read_json(url)
        stable = response["crate"]["max_stable_version"]
        release = next(item for item in response["versions"] if item["num"] == stable)
        if release["yanked"] or "-" in stable.split("+")[0]:
            raise ValueError("registry did not return a non-yanked stable release")
        return name, {
            "stable": stable,
            "rust_version": release["rust_version"],
            "license": release["license"],
            "source": url,
        }
    except Exception as error:
        return name, {"error": str(error), "source": url}


def version_core(version):
    return version.lstrip("=^~").split("+")[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", type=Path)
    args = parser.parse_args()
    specifications = []
    git_dependencies = []
    lock_hashes = {}
    for manifest, lock, sections in [
        ("Cargo.toml", "Cargo.lock", [("workspace", "dependencies")]),
        (
            "vendor/fasterp/Cargo.toml",
            "vendor/fasterp/Cargo.lock",
            [("dependencies",), ("dev-dependencies",)],
        ),
    ]:
        document = tomllib.loads((ROOT / manifest).read_text())
        raw_lock = (ROOT / lock).read_bytes()
        lock_hashes[lock] = hashlib.sha256(raw_lock).hexdigest()
        locked = tomllib.loads(raw_lock.decode())["package"]
        for section in sections:
            entries = document
            for component in section:
                entries = entries[component]
            for alias, specification in entries.items():
                entry = {"version": specification} if isinstance(specification, str) else specification
                name = entry.get("package", alias)
                if "git" in entry:
                    git_dependencies.append({"name": name, "manifest": manifest, **entry})
                    continue
                if "version" not in entry:
                    continue
                versions = [
                    package["version"]
                    for package in locked
                    if package["name"] == name
                    and package.get("source", "").startswith("registry+")
                ]
                # fasterp itself is patched to the vendored copy in the workspace.
                if name == "fasterp" and manifest == "Cargo.toml":
                    versions = [
                        tomllib.loads((ROOT / "vendor/fasterp/Cargo.toml").read_text())["package"]["version"]
                    ]
                specifications.append({
                    "name": name,
                    "manifest": manifest,
                    "requirement": entry["version"],
                    "locked": sorted(versions),
                })
    names = sorted({entry["name"] for entry in specifications})
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        releases = dict(pool.map(stable_release, names))
    failures = []
    for entry in specifications:
        latest = releases[entry["name"]]
        stable = latest.get("stable")
        entry["latest_stable"] = stable
        entry["current"] = bool(
            stable
            and version_core(entry["requirement"]) == version_core(stable)
            and any(version_core(version) == version_core(stable) for version in entry["locked"])
        )
        if not entry["current"]:
            failures.append(entry)

    git_results = []
    for entry in git_dependencies:
        result = dict(entry)
        if entry["name"] == "mincorr" and entry["git"] == "https://github.com/deminden/min-rust-corr":
            url = "https://api.github.com/repos/deminden/min-rust-corr/releases"
            result["release_source"] = url
            try:
                releases_on_github = read_json(url)
                stable = [release["tag_name"] for release in releases_on_github
                          if not release["draft"] and not release["prerelease"]]
                result["published_stable_releases"] = stable
                result["status"] = "needs_review" if stable else "documented_exception"
                result["reason"] = (
                    "No published stable release. Retain the exact Git revision with "
                    "Kendall tie-handling fixes absent from the older v0.2.2-gpl tag."
                )
                if stable:
                    failures.append(result)
            except Exception as error:
                result["error"] = str(error)
                failures.append(result)
        else:
            result["status"] = "unreviewed_git_dependency"
            failures.append(result)
        git_results.append(result)

    report = {
        "checked_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "scope": "Direct registry dependencies, including vendored optional/dev dependencies; compatible transitive updates use cargo update.",
        "lockfile_sha256": lock_hashes,
        "registry_releases": releases,
        "dependencies": specifications,
        "git_dependencies": git_results,
        "passed": not failures,
    }
    if args.write_report:
        args.write_report.parent.mkdir(parents=True, exist_ok=True)
        args.write_report.write_text(json.dumps(report, indent=2) + "\n")
    for failure in failures:
        print(json.dumps(failure), file=sys.stderr)
    print(f"Checked {len(names)} registry packages across {len(specifications)} direct requirements; "
          f"{len(failures)} failures, {len(git_results)} documented Git dependency decisions.")
    return int(bool(failures))


if __name__ == "__main__":
    sys.exit(main())
