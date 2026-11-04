# Preparing an RNAA release

The current workspace is 0.3.1, an unreleased update after the published 0.3.0 release.
Release readiness and completion of the improvement roadmap are separate decisions.

## Versioning and distribution

Use `0.MINOR.PATCH` while the CLI, catalog and scientific contracts are evolving.
Use a new minor version for incompatible CLI/configuration changes or catalog
upgrades that prevent downgrade. Use patch versions for compatible fixes. Public
previews use `-rc.N`; update the workspace version and lockfile together only when
the target version is chosen. Record scientific algorithm revisions separately
from the package version and rerun the independent scientific gates for engine
changes. Keep completed release entries immutable in `CHANGELOG.md`.

Distribute the CLI and the complete workspace source. A standalone `cargo publish`
package is not prepared: the workspace has local path dependencies, a patched
vendored library and a build script that reads the workspace lockfile. Build or
install from the complete checkout/source archive:

```bash
cargo build --release --locked -p rnaa
cargo install --locked --path crates/rnaa --root /your/install/prefix
/your/install/prefix/bin/rnaa --version
```

Rust 1.97.1 or newer is declared. The retained [release-preparation record](../validation/checks/2026-10-31-release-wrap.json)
describes an earlier snapshot; its historical log digests refer to local evidence
that is no longer tracked. Rebuild and qualify the exact source intended for release.

## Candidate gates

Run from the workspace root:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace --all-features
cargo test --manifest-path vendor/fasterp/Cargo.toml --lib --test compatibility_v140 --locked
cargo build --release --all-features --locked -p rnaa
python3 -m unittest discover -s scripts/tests
python3 scripts/check_dependencies.py --write-report target/release-evidence/dependency-audit.json
python3 examples/offline_reference.py --binary target/release/rnaa
python3 scripts/package_candidate.py --label 0.3.1 --allow-dirty
python3 scripts/verify_candidate.py target/release-candidates/0.3.1 --report target/release-evidence/candidate-install.json
```

Python helpers require Python 3.11+. The packager records a dirty working tree and
includes non-ignored untracked source files, so a candidate includes the actual
reviewed work. Unneeded historical validation files are omitted and listed by relative name in
the manifest; scientific reference files and selected evidence cited by guides are
retained and checked. The candidate source is therefore a release snapshot, not
a complete archive of every local validation log. The packager rejects local
home and temporary paths in other source files and in the compiled executable. Candidate bundles go under ignored `target/release-candidates/`.
The packager builds its executable from a temporary frozen copy of the exact
archived source and rejects changes to that copy during the build. It uses
`target/candidate-build` for compilation artifacts by default; pass
`--build-target-dir` to choose a different cache root. Cached artifacts are scoped
to the complete source snapshot and file modes. A separately built executable
cannot be supplied to the packager.
The bundles contain a source snapshot, a local executable bundle, dependency metadata,
source hashes, operator guides, examples, the selected validation evidence, runtime
requirements and checksums.
They are local review artifacts, not permission to
publish. The dependency inventory includes the selected host's normal/build
dependency closure; it is not a security scan or a license compatibility verdict.

The verifier checks the checksum inventory, rejects duplicate or misplaced archive
entries, compares every archived source hash and every source-derived bundled file,
checks the bundled dependency-notice inventory and executable/workspace version
agreement, and
installs from the extracted source before running both executables through the
offline example. Installation always uses a fresh isolated build directory;
`--install-target-dir target` chooses its parent but does not reuse old build
artifacts. This prevents stale compiled workspace libraries from being mistaken
for a successful source installation. Save the exact commands and outcomes in
the chosen release-evidence location.

Native Linux execution and installation scope are described in
[Linux qualification](linux-support.md). The CI matrix covers Ubuntu 24.04
x86-64 and ARM64 plus macOS ARM64, and also runs
the declared minimum Rust toolchain on Ubuntu x86-64. Every job runs the mandatory
checks and vendored preprocessing tests. On a pushed version tag, the three stable-Rust
platform jobs also build candidates and verify their archives and source installation.
Runner labels follow the [official GitHub runner list](https://docs.github.com/en/actions/reference/runners/github-hosted-runners).

After all four jobs succeed, the release job checks that the tagged commit is in
`main`, that the tag exactly matches the workspace version, with an optional `v` prefix, and that the
three candidates identify the same source snapshot and pass checksum validation.
It creates a draft, uploads a self-contained candidate archive for each platform
plus a `SHA256SUMS` file, and publishes the release. Prerelease version tags such
as `v0.3.0-rc.1` are marked as prereleases. Each platform archive contains its
binary and source archives, manifests and their original checksums. GitHub Actions
artifacts transfer candidates between jobs; the durable downloads are release assets.
No release is created by a normal branch push or pull request.

A configured job is not evidence of a successful run. Before a commit and again
before a push, rerun the mandatory root sequence in `AGENTS.md`. Pushing a matching
version tag triggers publication; create and push a tag only when the exact commit
is approved for release.

## Current public-release blockers

| Gate | Current evidence / required action |
| --- | --- |
| Workspace correctness | The macOS and Ubuntu Linux ARM64 release-preparation sequences each passed 508 tests after the readiness-marker fix. Rerun the full sequence on the final clean release source. |
| Scientific engine scope | Validate `rsdeseq2 0.4.0` and the updated `fasterp` revision against the independent references and preprocessing regression tests; numerical and behavior changes require review. |
| Distribution licensing | Review deferred by explicit user instruction. Dependency license declarations and bundled notices remain intact; no license is changed. This is outside the current technical release blockers. |
| Platform / installation support | Run native Linux checks, optimized build and packaged-source installation; retain exact host/architecture/compiler evidence. Hosted x86-64/ARM64/macOS matrix must pass on the final commit before advertising those targets. |
| Operational qualification | A four-run public yeast workflow passed on macOS and Linux ARM64; see `validation/real-data-smoke/2026-10-31-yeast-release.md`. On macOS, streamed adjustment reproduced the earlier candidate's adjusted matrix byte for byte; this is one small dataset. Broader scale, migration, interruption and combined-workflow qualification remain necessary for production claims. |
| Product scope | Extraction remains incomplete; the roadmap also retains broader freshness coverage, analysis/publication admission, cross-project semantic reuse and remaining cleanup/cancellation work. |
| Release identity | Version 0.3.1 is selected. Create a reviewed clean commit and push its matching version tag only when ready to publish. |

## Upgrade and rollback

Stop all writers using the project and shared store. Back up the entire project
and relevant shared-store payloads/catalogs together; copying only the SQLite main
file while a writer has an active WAL is not a complete backup. Preserve the old
executable and configuration alongside the backup.

Opening a project applies catalog migrations. Run `rnaa status`, review active
cohort/reference selection and inspect retained results before resuming work.
Reference request revision 2 can require a rebuild because the candidate now
consumes complete multimember gzip input. Verify resource limits and disk headroom
before that rebuild. Reuse and cleanup behavior are described in the linked
operator guides in the README.

Roll back by restoring the complete idle backup and its matching old executable.
Do not point an older executable at an upgraded catalog. Existing integration
tests cover migration/recovery cases; a rehearsal with representative retained
projects is still required for the public release.
