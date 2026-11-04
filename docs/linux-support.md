# Linux release qualification

A 0.3.0-rc.1 release-preparation snapshot passed formatting, Clippy and
**508 tests in 45 suites** on native **Linux ARM64** in an Ubuntu 24.04.5 LTS VM
(Linux 6.8, six CPUs, 12 GiB memory). An all-target minimum-toolchain check also
passed. These results execute Linux binaries against a Linux kernel and filesystem;
they are not a macOS cross-build. See the [source-bound results](../validation/checks/2026-10-31-release-wrap.json).
Subsequent source and documentation edits are not covered by those results; rerun
qualification and rebuild for the final release source.

Local qualification establishes this Ubuntu/ARM64 target. Linux x86-64 is included
in the CI matrix for both stable Rust and the declared minimum, but hosted runs
against the final commit remain a separate gate. Other distributions, older glibc,
musl, Windows and other architectures have not been qualified by these results.
No supported-host performance claim follows from a VM smoke check.

## Build and installation

On an Ubuntu 24.04 machine with a complete workspace checkout/source archive:

```bash
sudo apt-get update
sudo apt-get install build-essential pkg-config python3 git curl ca-certificates time
# Install Rust 1.97.1 or newer using the official rustup installer.
cargo build --release --all-features --locked -p rnaa
python3 examples/offline_reference.py --binary target/release/rnaa
python3 examples/offline_reference.py --binary target/release/rnaa --automatic-resources
cargo install --locked --path crates/rnaa --root "$PWD/target/install"
python3 examples/offline_reference.py --binary target/install/bin/rnaa
```

Python 3.12 is used on Ubuntu 24.04. The executable itself does not require Python,
R or external bioinformatics tools. Building needs a C toolchain for native
dependencies, including bundled SQLite and compression/cryptography libraries.
The release report records the executable's actual dynamic dependencies; no
static-linking or older-distribution compatibility claim is made.
On a memory- or disk-limited build host, setting `CARGO_BUILD_JOBS=1` and
`CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0` reduces the resources needed
for development checks. The full current Linux test run used those settings after
parallel linking exhausted VM memory and old build caches filled its disk.
The locally built ARM64 bundle uses **GLIBC_2.39** symbols and requires glibc
2.39 or newer. Its release manifest records that floor. Older Linux distributions
need a separately qualified build; the Ubuntu 24.04 artifact is not advertised as
portable to them.

The automatic-resources example exercises Linux CPU/memory detection with a
private reservation pool. Pool files and catalogs must reside on a local writable
filesystem supporting locks and atomic renames. Distributed filesystem behavior
has not been qualified. Use explicit `[resources]` limits when automatic detection
does not describe the workload's assigned capacity.

## Archive and CI gates

```bash
python3 scripts/package_candidate.py --label 0.3.1 --allow-dirty
python3 scripts/verify_candidate.py target/release-candidates/0.3.1 --report target/release-evidence/linux-install.json
```

The verifier checks archive checksums, exact file inventories, source/executable
version agreement, source hashes and machine-specific home/temporary paths.
Unneeded historical validation files are excluded and recorded in the manifest. It installs from extracted source and runs
both executables with explicit and automatically detected resources. Passing
`--install-target-dir` chooses a parent for the fresh isolated build directory. The small
offline example validates installation and reference reuse; the workspace suite
provides the broader scientific and recovery checks.

CI runs the workspace checks on Ubuntu 24.04 x86-64/ARM64 and macOS ARM64,
with a separate Ubuntu x86-64 minimum-Rust job. A pushed matching version tag also
builds and verifies the stable-Rust platform candidates; after all jobs pass, the
release job publishes them as GitHub Release assets. See [release preparation](releasing.md)
for the tag and source-identity checks. Passing local tests does not imply that an
unexecuted hosted workflow has passed.
