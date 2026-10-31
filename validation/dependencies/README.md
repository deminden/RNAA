# Dependency evidence

[The 2026-10-31 audit](2026-10-31-release-wrap.json) records the direct registry
dependency review for the workspace and vendored preprocessing library. It is a
historical record, not a guarantee that versions remain current or a security audit.

Run a fresh audit after changing either lockfile, using Python 3.11 or newer:

```sh
python3 scripts/check_dependencies.py --write-report target/dependency-audit.json
```

The audit checks stable registry releases, direct requirements and resolved versions.
Network failures are reported as failures, not treated as current evidence. Upstream
constraints can retain older transitive versions. Run the workspace checks in the
[release guide](../../docs/releasing.md) after dependency updates.
