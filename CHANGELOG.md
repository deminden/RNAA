# Changelog

## 0.3.1 — unreleased

- Update `rsdeseq2` from 0.3.1 to 0.4.0, including upstream model and numerical
  corrections. Independent pinned R-reference tests remain the acceptance gate.
- Integrate `fasterp` GitHub main revision `dcca6ad25e238ecfb4c2372ae5d4c83b6cc68281`
  (upstream package version 0.2.1), preserving RNAA's bounded native API. RNAA
  continues to use its bounded native path; upstream's newer compatibility engine
  is available in the vendor library but is not silently substituted for it.
- Bump preprocessing request and processing identity to revision 4. Earlier
  preprocessing results cannot be mistaken for outputs from the updated engine.
  Reprocessing may require recoverable original inputs; retained historical
  results are preserved. Dependency lock changes also affect cache identities.
- Update direct registry dependencies and compatible lockfile selections; see the
  [dependency audit](validation/dependencies/2026-11-04-0.3.1.json).
- Add upstream compatibility and WebAssembly checks to CI, alongside independent
  scientific references and bounded processing regression tests.

## 0.3.0 — 2026-10-31

This release consolidated the development work since 0.2.1, including incompatible
catalog and CLI changes. See the gates in [the release guide](docs/releasing.md).

### Added

- Durable attempts, leases, cancellation and fenced publication of complete output
  generations, with recovery and retained producer history.
- Cohort membership, biological-sample pooling, typed designs, controlled metadata
  refresh, retained QC decisions and separate ingest/analysis progress.
- Content-addressed shared storage, dependency lineage and reversible cleanup of
  eligible inputs, generations and tracked unpublished work.
- Independent, pinned scientific fixtures for import, normalization, differential
  expression, five correlation methods, PCA and module enrichment.
- A common local resource pool for downloads, reference preparation, preprocessing
  and quantification, plus bounded input/index preflight and shared resident indexes.
  Correlation adjustment and network computation also reserve estimated resources;
  other analysis phases are not yet fully admitted.
- Operator guides, dependency audits, scale probes and source-bound check evidence.
- Native Ubuntu ARM64 qualification and a CI matrix covering Linux x86-64/ARM64,
  macOS ARM64 and the declared minimum Rust compiler.
- Candidate archive/source-install verification, automatic-resource smoke checks
  and refusal to package an executable with the wrong workspace version.

### Changed

- Doctor provides read-only catalog/configuration checks and a JSON report. The
  vendored S3 client is updated to 1.157.0. Its lockfile change alters the recorded
  dependency fingerprint, so previously cached references, preprocessing,
  quantification and analysis may need recomputation; retained history is preserved.
- `status --json` reports structured active-cohort progress, stage outcomes,
  reference selection, storage, generation heads and diagnostic errors. A stale
  cohort produces an explicit reason and null active progress.
- `export` writes TSV or CSV tables from verified sealed analysis generations,
  verifies the bytes actually consumed, and records provenance and output hashes
  with bounded, cancellable manifests and atomic no-replacement publication. Exported
  generations are frozen snapshots and are not revalidated against current settings.
- `survey` combines compatible frozen analysis and correlation generations with
  sample metadata summaries, provenance and checks on sample and gene identities.
- Interrupted tracked exports and surveys are reconciled on the next export or
  survey after the prior lease expires, preserving staging for inspection.
- Export and survey can optionally require current active-cohort membership;
  other scientific inputs and settings are not revalidated by this check.
- An additional admission check compares QC eligibility, quantified abundances
  and the selected reference without refreshing QC decisions.
- An optional design check reconstructs the recorded model from current used
  sample metadata; other scientific settings remain outside its scope.
- `--require-current` reconstructs supported current scientific requests and
  checks saved module options, pathway contents and selected external matrices;
  it fails closed when the relevant current inputs cannot be established.
- Correlation computes tiled outputs with bounded top-k state; module scores use
  standardized first principal components with deterministic signs. Gene selection
  and matrix adjustment process rows incrementally; adjustment and primary network
  output enforce their reserved byte limits, and network work checks cancellation
  between groups of tiles.
- Reference preparation consumes every gzip member, plans decoded input and graph
  shape before building, and keeps temporary work inside recoverable generations.
- Reference requests use algorithm revision 2. Older references may rebuild once;
  retained generations remain available. Resource-only changes preserve valid reuse.
- Large workflow files are split around ownership, input planning, publication,
  resource guards and recovery responsibilities.
- Direct registry dependencies use the audited stable releases. `mincorr` remains
  pinned to a Git revision; vendored `fasterp` retains documented RNAA patches.

### Upgrade notes

- Back up an idle project, including its SQLite catalog, configuration, generations
  and any shared store, before opening it with this release. Schema upgrades run
  when the catalog opens. Downgrading an upgraded catalog is not supported; restore
  the backup with its matching executable instead.
- All commands sharing a resource pool must use the same pool configuration.
  Analysis and some publication paths still fall outside this pool.
- Consult [FASTQ limits](docs/fastq-preparation.md) and
  [reference/resource limits](docs/resource-reservations.md) before large jobs.
- DESeq2 replacement is deferred. Production continues to use `rsdeseq2 0.3.1`;
  this release does not claim to fix the separately tracked engine issues.

### Known limits

SRA extraction is incomplete. Export and survey use frozen results; current-settings validation is optional
and limited to supported inputs. Cross-project
semantic result reuse, some cleanup/recovery classes, prompt cancellation inside
native computations and supported-host production qualification remain open.
The release inventory retains dependency licenses; mincorr licensing review is
deferred by user instruction. Production-scale qualification remains incomplete.
