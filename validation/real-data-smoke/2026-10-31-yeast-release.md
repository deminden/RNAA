# Public yeast release smoke (2026-10-31)

The commands in [the public yeast example](../../examples/public_yeast.md) were
run in a fresh temporary project on macOS ARM64 with the bundled CLI from local
candidate `0.3.0-rc.1-macos-exact`. Its runtime source files match the workspace
used for this run; the candidate was a dirty local snapshot and predates later
packaging-script and documentation changes. It was not published.

RNAA resolved all four runs, downloaded **75,523,534 compressed bytes**, and
quantified all four. Independently calculated SHA-256 values matched the pinned
read checksums in the [earlier dataset manifest](manifest.json); the two reference
files also matched their pinned hashes. The project explicitly selected
`saccharomyces_cerevisiae` and release `63` for the custom reference.

DE, residualized Spearman correlation, strict export and strict survey completed.
The survey reported **4 samples**, two per condition, **6,612 genes** and 1,000
correlation genes. The unannotated top-k edge table contained 10,000 rows. On
rerunning quantification, DE and correlation with the same inputs, all 11 selected
generation heads remained unchanged. `doctor --json` reported success.

The same tutorial also completed on Ubuntu 24.04.5 ARM64 with the Linux
candidate, including read checksum verification and unchanged heads on repeat.
The counts, normalized counts, VST and edge tables had identical SHA-256 values
across these two runs. The DE table differed in 13 numeric cells;
the largest absolute difference was approximately `3.11e-15`. This single
dataset does not establish general cross-platform scientific parity. See the
[Linux record](2026-10-31-yeast-linux.json). The original tutorial log is excluded from the repository because it contains
machine-specific paths.

The [machine-readable record](2026-10-31-yeast-release.json) contains input and
output hashes, engine versions, tutorial hash and log hashes. The full tutorial and reuse logs are excluded from the tracked release source. This shallow four-sample workflow is an execution check,
not a biological conclusion, DESeq2 parity test, production-scale qualification
or a substitute for Linux x86-64 and final clean-release checks. DESeq2 remains
at the user-requested existing version.
