# RNAA Development Plan

## Purpose

This document defines the next development stages for RNAA from the current repository state:

- workspace scaffolding exists and the CLI is usable
- SQLite-backed project state, manifests, and artifact tracking exist
- `init`, `add`, `resolve`, `status`, `doctor`, `download`, `refs`, `quant`, `normalize`, `deseq2`, `corr`, and `run` are implemented
- ENA resolution works for runs, studies, projects, and sample selectors and generates a canonical samplesheet
- the downloader is restart-safe, manifest-driven, and independent of later stages
- restart reconciliation exists for quant, normalize/DE, and corr
- `export`, `survey`, cleanup execution, extraction fallback, crash-recovery sweeps, and a real worker-budget scheduler are still incomplete

The objective is not just to "finish features". The objective is to make RNAA credible as a server-first, restart-safe survey engine that can survive long-running public RNA-seq ingestion workflows without corrupting state or locking us into temporary implementations.

## Guiding Strategy

The next stages should follow five rules:

1. Preserve contract-first design.
Every stage must produce stable internal artifacts and manifests so we can replace shell or R implementations later without rewriting the rest of the system.

2. Finish vertical slices, not disconnected helpers.
Each stage should become independently usable from the CLI before moving to the next one. A command that exists but only bails is worse than no command at all.

3. Make failure states explicit.
The pipeline has to be restart-safe under server conditions. That means durable DB state, verifiable manifests, atomic file writes, and predictable retries before optimization.

4. Keep download independence sacred.
`rnaa download --forever` is the operational backbone. Later work must not couple project progress to a single monolithic orchestrator process.

5. Prioritize combined-cohort correctness over feature breadth.
Multi-study harmonization, design validation, metadata discipline, and artifact provenance matter more than adding more aligners quickly.

## Architecture Pivot: HDF5-First + Cross-Project Reuse

This plan now assumes a storage architecture change:

1. Persist quant outputs as HDF5-first artifacts.
- Treat `abundance.h5` as the canonical quant payload.
- Keep `run_info.json` and stage manifests for provenance.
- Do not persist `abundance.tsv` in steady state. If a tool emits TSV by default, treat it as ephemeral and remove after verification.

2. Persist large matrix outputs in HDF5 internally.
- DE counts, normalized counts, transforms, adjusted matrices, and dense correlation outputs should use HDF5 contracts.
- Human-readable TSV/CSV become export-time products, not canonical storage.

3. Introduce a global shared artifact store outside project roots.
- Add a configurable shared root (for example `RNAA_SHARED_ROOT`) with a content-addressed layout by sha256.
- Project directories become views over shared immutable blobs plus project-local manifests/events.
- Prefer hardlinks or reflinks where available; fall back to copy with checksum verification.

4. Make reuse deterministic with cache keys.
- Quant cache key must include:
  - input FASTQ digests
  - reference identity digest
  - quant engine + version
  - quant parameters
  - preprocessing fingerprint
- If the key already exists in shared storage, attach artifacts to the project without recomputation.

5. Add explicit ownership, pinning, and GC.
- Track which projects reference each blob.
- Allow pinning of references/artifacts used by active projects.
- Add garbage collection that removes unreferenced blobs after retention windows.

6. Keep manifests and DB rows project-local, payloads global.
- Project SQLite remains the control plane.
- Payload immutability and dedupe live in the shared artifact store.
- Every project artifact record should resolve to a global blob id plus a materialized path.

7. Support sample-scoped project composition.
- `rnaa add` must support sample-level selectors in addition to project/study/run IDs.
- Accept explicit sample accessions plus metadata predicates for include/exclude filters.
- Persist selectors as first-class inputs so cohort definitions are reproducible and rerunnable.

8. Make FASTQ preprocessing with `fasterp` a first-class pipeline stage before quant.
- `fasterp` is a trimming/QC preprocessor, not an aligner.
- Keep raw downloaded FASTQ immutable; write trimmed FASTQ as derived artifacts.
- Persist preprocessing reports (`json`, optional `html`) as provenance artifacts.
- Include preprocessing tool/version/flags in quant cache identity so reuse is correctness-safe.
- Treat preprocessing as part of the main run graph, not a side feature. The default orchestration path should understand:
  - `VERIFIED -> PREPROCESS_RUNNING -> PREPROCESS_DONE -> QUANT_RUNNING`
  - direct `VERIFIED -> QUANT_RUNNING` only when preprocessing is disabled
- Scheduler decisions must budget preprocessing together with quant because both are run-parallel CPU stages.

## Current State Snapshot

### Working now

- Workspace layout and crate boundaries are in place.
- `rnaa-core` has config, models, DB, state enums, paths, manifests, traits, and samplesheet logic.
- SQLite migrations exist and the state/event/artifact tables are usable.
- ENA resolver fetches `filereport` JSON and expands run metadata.
- Metadata overrides merge into resolved metadata.
- Canonical samplesheet and editable column map are generated.
- Download worker uses `curl`/`wget`, `.part` files, retries, verification, and manifest output.
- Reference preparation is wired through the CLI, uses cached GENCODE presets, and records reference artifacts.
- `fasterp` preprocessing exists and writes derived artifacts before quantification when enabled.
- Quantification works through the Rust wrapper plus `Rscript` plus kallisto.
- Normalize and DESeq2 work through the Rust wrapper plus `Rscript` plus tximport/DESeq2.
- Correlation works through Rust residualization plus linked `mincorr`; module-level enrichment via `rsfgsea` is present.
- `rnaa run` overlaps download and quant and reconciles existing downstream outputs.
- Tests exist for resolver parsing, metadata merge behavior, state persistence, download templating, normalize workflow, and normalize -> corr -> module enrichment workflow.

### Incomplete or still weak

- Export and survey are still CLI stubs.
- Cleanup policy exists in config and CLI flags, but real move-to-trash plus purge execution is not complete.
- Extraction fallback (`xsra` and `vdb-native`) is still not implemented.
- Project-stage state is still overloaded onto per-run state for normalize, DE, and correlation.
- Crash-recovery is incomplete: there is no sweep that turns stale `*_RUNNING` states into retryable interrupted states after abrupt termination.
- Reference preparation, quantification, DESeq2, and correlation still need deeper integration tests and harder failure-path coverage.
- There is no real worker-budget scheduler yet. `run` overlaps download and quant, but DE/corr scheduling is still simple and stage budgets are not explicit.
- Preprocessing is implemented, but it is not yet modeled as a first-class scheduler stage in the orchestration plan.
- Shared artifact storage is only partially realized; ownership, pinning, and GC remain unfinished.

## Parallelization Snapshot

The current execution model is mixed:

- Download is parallel across runs. `rnaa download --concurrency N` spawns `N` Rust worker threads, and each worker processes one run at a time.
- Preprocessing and quant can run per-run with bounded worker pools, but the orchestration rules are still simple.
- DESeq2 is single-project and single-invocation from RNAA. That is correct at the scheduler level, but RNAA does not yet explicitly cap BLAS/OpenMP threads.
- Correlation is single-project from RNAA. Residualization and `mincorr` execution still need explicit bounded worker ownership.
- `rnaa run` is a pipeline orchestrator, but not yet a real global worker-budget scheduler.

If you give the current code a 100-worker host today, RNAA can overlap I/O and run-parallel compute, but it still does not enforce a single budget across preprocessing, quant, DESeq2, and correlation. That is the remaining scheduler gap.

## Parallelization Strategy

The correct strategy is not "every stage gets all workers". The correct strategy is a stage-aware worker budget with explicit CPU and I/O ownership.

### Rule set

1. Download owns an I/O pool, not the whole machine.
- Default to a small separate I/O budget, for example `download_workers = 1` or `2`.
- Make it configurable, but do not scale it up automatically with host core count.
- Treat download concurrency as an I/O tuning knob, not a CPU saturation target.
- Allow higher manual override on large servers only when the operator explicitly asks for it.
- Keep per-run file handling sequential for now to preserve simple integrity and manifest semantics.

2. Download and processing must overlap.
- The downloader must keep running while quantification is active.
- Verified runs should enter the quant queue immediately; they should not wait for the full project download set to finish.
- The execution model should be pipelined:
  - `RESOLVED -> DOWNLOADING -> VERIFIED` on the I/O side
  - `VERIFIED -> QUANT_RUNNING -> QUANT_DONE` on the CPU side
- `rnaa run` should observe DB state and schedule newly verified runs while the download worker continues independently.
- If an external `rnaa download --forever` worker is already running, `rnaa run` must cooperate with it, not replace it.

3. Quant should saturate the CPU budget by parallelizing across runs first, then within run second.
- Do not give one kallisto job 100 threads by default.
- Introduce:
  - `execution.total_workers`
  - `execution.reserve_workers`
  - `quant.max_concurrent_runs`
  - `quant.threads_per_run`
- Scheduling rule:
  - fill the CPU budget with as many run-level quant jobs as possible
  - only increase `threads_per_run` when there are too few runnable runs to fill the machine
- For a 100-worker host, the default target should look like:
  - many runnable runs: `12 x 8-thread` or `16 x 6-thread` quant jobs
  - few runnable runs: increase per-run threads gradually, but cap a single kallisto job well below 100
- The goal is 100 busy workers in aggregate, not 100 threads per run.

4. Extraction and preprocessing should behave like quant-adjacent run-parallel work.
- `xsra` or future native extraction should use the same CPU pool family as quant.
- `fasterp` must be treated as part of the main run pipeline, not a side option.
- Default preprocessing concurrency should be bounded and coordinated with quant so the two stages do not oversubscribe the machine.
- Default to a smaller concurrency than quant when compression is involved.

5. DESeq2 should remain single-project and non-parallel at the RNAA scheduler level.
- One DESeq2 analysis per project at a time.
- Do not shard one cohort analysis into parallel RNAA jobs.
- Set safe thread caps explicitly through environment variables for the R process:
  - `OMP_NUM_THREADS`
  - `OPENBLAS_NUM_THREADS`
  - `MKL_NUM_THREADS`
- `VECLIB_MAXIMUM_THREADS`
- `NUMEXPR_NUM_THREADS`
- For a 100-worker host, the default DESeq2 cap should be conservative, for example `8` or `16`, not `100`.
- DESeq2 should start only after the cohort barrier is satisfied:
  - either all known project runs are quantified
  - or the user explicitly freezes a partial cohort snapshot

6. Correlation should be single-project but internally parallel with an explicit ceiling.
- Residualization should be parallelized over genes using Rayon.
- `mincorr` should run in a bounded thread pool, not inherit "all cores" implicitly.
- For a 100-worker host, a safe default is a bounded correlation pool such as `48` to `64` workers unless the user overrides it.
- Correlation should not be parallelized "by runs" because it is a cohort matrix stage, not a run stage.

7. `rnaa run` should become a scheduler, not just a sequencer.
- It should compute a worker budget once.
- It should allocate workers per stage.
- It should publish:
  - active stage
  - workers assigned
  - queued runs
  - effective per-run thread count
- Download must remain independently runnable even after this scheduler exists.
- In the steady state, `run` should look like a conveyor:
  - downloader keeps feeding verified runs
  - quant keeps draining verified runs
  - DE/corr wait for the project-level barrier
  - status reports both pipelines at once

### Concrete 100-worker target model

For a 100-worker Linux host, the default strategy should be:

- reserve `4` workers for OS, SQLite, logging, and control-plane work
- allow download to use a small independent I/O pool, for example `1` or `2`, with manual override if needed
- let preprocessing and quant share the CPU pool as run-parallel stages:
  - default `preprocess.max_concurrent_runs = 6`
  - default `quant.threads_per_run = 8`
  - default `quant.max_concurrent_runs = 12`
  - the scheduler should reduce preprocess concurrency when quant is already saturating the machine
- let quant consume the CPU pool aggressively by run-level parallelism:
  - default `threads_per_run = 8`
  - default `max_concurrent_runs = floor(96 / 8) = 12`
  - if fewer than 12 runs are ready, increase per-run threads up to a cap such as `16`
- allow the downloader and quantifier to run at the same time:
  - example steady state: `2` download workers plus a mix of preprocess and quant jobs, bounded to the remaining CPU budget
  - this is valid because download is mostly I/O bound while preprocess and quant are CPU bound
- run DESeq2 as exactly one project job with a hard thread cap such as `8` or `16`
- run correlation as exactly one project job with an explicit pool such as `48` or `64`

That gives the behavior you want:

- alignment uses the whole machine in aggregate
- preprocessing is part of the main data path instead of an afterthought
- downloads keep progressing while alignment is already processing earlier runs
- DESeq2 does not explode into unsafe parallel jobs
- other stages use bounded, stage-appropriate worker counts
- run-parallel work happens where the contracts are naturally per-run

## Priority Gaps

The next work should focus on the things that still block RNAA from being convincingly production-grade.

1. State model correctness
- add project-level stage state instead of overloading run-level state for DE and correlation
- add startup recovery that converts stale `*_RUNNING` states into retryable interrupted failures
- expand shared artifact-store linkage into full ownership, pinning, and GC semantics (basic blob linkage is now in place)

2. Worker scheduling
- implement the worker-budget model described above
- make preprocessing a first-class scheduled stage in `run`
- add explicit thread caps for DESeq2 and correlation
- publish worker assignments and queue state in progress reporting

3. Honest feature surface
- either implement cleanup, export, survey, and extraction, or clearly leave them out of the public promise
- do not keep CLI flags that imply working behavior if the underlying stage is still a placeholder

4. Testing
- add integration tests for:
  - refs prepare
  - preprocess happy path
  - quant wrapper happy path
  - DESeq2 happy path
  - correlation happy path
  - end-to-end `rnaa run`
- add restart/interruption tests for download, quant, and project-level stages

5. Hardening
- add reference-build locking so concurrent processes cannot trample shared cache creation
- add output guardrails for dense correlation outputs
- verify manifest completeness and failure-mode consistency across stages
- add shared-store GC safety rails: retention windows, pinning, and lock-safe sweeping

## Recommended Delivery Sequence

The smart path is:

1. Implement shared artifact store contracts and schema additions.
2. Add project-level stage state and crash-recovery sweeps.
3. Implement a real worker-budget scheduler for download, preprocessing, quant, DE, and corr.
4. Finish reference preparation and locking against shared storage.
5. Finish per-run preprocessing and quant reuse contracts with stable cache keys.
6. Finish project-level DESeq2 and correlation hardening.
7. Finish export, survey, extraction, and cleanup UX.
8. Harden operational behavior, tests, docs, and release packaging.

This sequence is correct because the remaining risk is no longer basic feature existence. The remaining risk is state-model correctness, scheduler correctness, and operational hardening across long-running server workflows.

## Stage 2: Reference Preparation As A Production-Ready Stage

### Goal

Make `rnaa refs prepare` fully operational for human and mouse, deterministic, cacheable, and auditable.

### Current status

The `EnsemblReferenceManager` already contains most of the logic for:

- Ensembl listing resolution
- download/copy of cDNA and GTF
- `tx2gene.tsv` generation
- `gene_annotation.tsv` generation
- kallisto index building
- reference manifest generation

The missing part is integration, testing, and hardening.

### Tasks

- Wire `refs prepare` in the CLI to `EnsemblReferenceManager`.
- Support explicit flags overriding config:
  - `--organism`
  - `--ensembl`
  - `--gtf`
  - `--cdna`
- Persist selected reference info into the DB or a stable reference manifest lookup table.
- Record reference artifacts in the `artifacts` table.
- Materialize reference payloads into the shared artifact store and link them into project views.
- Add checksum verification for downloaded reference files.
- Parse Ensembl `CHECKSUMS` when feasible instead of trusting only local sha256.
- Handle reference reuse:
  - skip rebuild if manifest and artifacts already match
  - skip redownload/reindex if matching blobs already exist in shared storage
  - rebuild if the index is missing or empty
- Add offline unit tests for:
  - listing parsing
  - custom reference mode
  - GTF-derived tx2gene extraction
- Add a live integration test behind an opt-in flag for current Ensembl.

### Design notes

- Reference identity should be stable and compositional: organism + release for presets, content hash for custom inputs.
- The generated `tx2gene.tsv` and `gene_annotation.tsv` are core contracts; later stages should consume them, not re-derive them repeatedly.
- Avoid baking Ensembl-specific assumptions into downstream stages. Downstream stages should consume `ReferenceBundle`, not "human release N" strings.

### Acceptance criteria

- `rnaa refs prepare --organism human --ensembl latest` works end to end.
- Cached reruns are idempotent.
- Two projects using the same reference reuse the same immutable blobs without duplicate payload bytes.
- Reference manifest and artifacts are sufficient for downstream reproducibility.

## Stage 3: Per-Run Quantification Via R + Kallisto

### Goal

Implement the first real quantification engine with clean Rust/R boundaries and safe per-run state transitions.

### Tasks

- Create `r/quant_kallisto.R`.
- Add robust argument parsing in R.
- Add a pre-quant preprocessing contract:
  - optional `fasterp` execution per run
  - deterministic output paths for trimmed FASTQ artifacts
  - preprocessing manifest/report capture
- Validate:
  - one FASTQ for single-end
  - two FASTQs for paired-end
  - kallisto index existence
  - output directory writability
- Add tool policy:
  - `fasterp` enabled/disabled via config and CLI
  - strict mode (fail quant on preprocessing failure) vs bypass mode (warn and continue with raw FASTQ)
- Run `fasterp` on FASTQ inputs before quant when enabled.
- Run `kallisto quant` via `system2`, never shell concatenation.
- Capture:
  - tool versions
  - command arguments
  - reference manifest hash
  - input FASTQ checksums
- Capture preprocessing provenance:
  - pre/post read counts
  - trimmed FASTQ checksums
  - `fasterp` version and flags
- Define and persist a deterministic quant cache key from input digests + reference + engine + params + preprocessing fingerprint.
- Reuse existing quant payloads from shared storage when the cache key matches.
- Define success predicate in Rust:
  - output directory exists
  - `abundance.h5` exists and is non-empty
  - `run_info.json` exists and is parseable
- Record quant artifacts into SQLite with global blob linkage.
- Treat `abundance.tsv` as non-canonical; if generated, purge it after validation.
- Advance state:
  - `VERIFIED -> QUANT_RUNNING -> QUANT_DONE`
  - failure should produce `QUANT_FAILED` without damaging downloaded files
- Implement `rnaa quant` filtering logic:
  - run only on `VERIFIED`, `QUANT_FAILED`, or explicitly forced runs
- Implement per-run log files under `logs/quant/` or a structured equivalent.

### Cleanup substage

After quant works, complete the cleanup behavior:

- `--cleanup none|fastq|sra|all`
- `--cleanup-on success|always|never`
- two-step deletion:
  - move into project trash
  - purge after N days
- record cleanup artifacts and events
- ensure repeated cleanup is idempotent

### Why this stage is next

Quantification is the first place where the download layer proves its value. Until quant works per run, the project cannot produce biologically meaningful outputs, and later DE/correlation work would be forced to invent fake inputs.

### Acceptance criteria

- `rnaa quant` works on verified runs with human or mouse references.
- Quant reruns are safe and deterministic.
- Cross-project quant reuse works when inputs/reference/parameters match.
- When preprocessing is enabled, quant consumes trimmed FASTQ and records reproducible preprocessing provenance.
- Cleanup only triggers after verified quant success.

## Stage 4: Project-Level DESeq2 And Normalization

### Goal

Turn multiple per-run quant outputs into a coherent project-level cohort analysis.

### Tasks

- Create `r/tximport_deseq2.R`.
- Feed it:
  - samplesheet
  - quant directories with canonical `abundance.h5`
  - `tx2gene.tsv`
  - design formula
  - contrasts
  - transform option
- Build design validation in Rust before spawning R:
  - verify required columns exist
  - reject formulas referencing absent variables
  - provide actionable errors for missing `condition` or `batch`
- Support contrast specification from CLI and persisted DB state.
- Produce:
  - `gene_counts.h5`
  - `gene_norm_counts.h5`
  - `vst.h5`
  - `vst.rds`
  - `de_<contrast>.h5`
  - `size_factors.tsv`
  - optional PCA and sample distance outputs
- Add manifest with:
  - `sessionInfo()`
  - DESeq2 version
  - tximport version
  - design formula
  - contrast definitions
  - sample inclusion/exclusion decisions
- Decide project state semantics carefully:
  - DE is project-level, not per-run
  - do not incorrectly stamp every run as `DE_DONE` without a project-level contract

### Required schema adjustment

The current DB model is run-centric. DE and correlation are cohort-level stages. Before finishing this stage, add a project-stage tracking mechanism, for example:

- `project_stage_runs`
- or `project_artifacts`
- or `project_states`

This prevents awkward overloading of per-run state for cohort-wide outputs.

### Multi-study strategy

- default combined-cohort DE should include `study_accession` as an available design covariate
- provide `--stratify study_accession` later, but combined analysis should remain first-class
- preserve user-visible metadata provenance so harmonization decisions are inspectable

### Acceptance criteria

- `rnaa deseq2 --design "~ batch + condition" --contrast condition A B` produces all expected artifacts.
- Missing design columns fail early with actionable errors.
- Cohort-level provenance is recorded cleanly.
- Internal DE artifacts are HDF5-first; text tables are produced only through export surfaces.

## Stage 5: Adjusted Correlation In Rust + `mincorr`

### Goal

Produce biologically interpretable, metadata-adjusted correlation outputs without routing the core math back through R by default.

### Tasks

- Implement matrix loading from:
  - VST HDF5
  - normalized counts if requested
- Implement geneset grammar parser:
  - `topvar:N`
  - `mean_gt:X`
  - `genes:file.txt`
  - `regex:...`
  - `biotype:...`
- Implement metadata design parsing for adjustment.
- Build residualization in Rust:
  - construct model matrix from samplesheet columns
  - fit per-gene linear model
  - return residual matrix
  - use `f64` internally
- Persist adjusted matrix artifact in HDF5 and manifest linkage into shared storage.
- Replace the current external-binary assumption with linked `mincorr` usage where practical.
- Define output modes:
  - top-k edges per gene
  - thresholded edges
  - dense matrix for small selections only
- Add guardrails:
  - reject dense full-gene outputs above configured size thresholds
  - report filtered gene counts and dropped genes explicitly
- Record method metadata:
  - correlation method
  - adjustment model
  - geneset selection result
  - output mode

### Technical strategy

Use the linked `mincorr` crate as the compute engine where possible, but keep RNAA’s own `Correlator` trait stable. The internal RNAA contract should remain:

- adjusted matrix in
- correlation artifacts out
- manifests capture method and selection

That preserves the ability to swap to direct streaming/top-k or sparse-native implementations later.

### Acceptance criteria

- `rnaa corr --matrix vst --adjust residualize --model "~ batch + condition" --method pearson --geneset topvar:5000 --out edges:topk=50` works end to end.
- Output size guardrails prevent accidental dense explosions.
- Adjustment and correlation are both reproducible and manifest-backed.

## Stage 6: CLI Completion And Operator UX

### Goal

Turn the pipeline from a set of stage commands into an operable product.

### Tasks

- Implement:
  - `rnaa run`
  - `rnaa export`
  - `rnaa survey`
- Keep `run` decoupled from download ownership:
  - it may trigger downloads optionally
  - it must also work when a separate download worker is already active
- Add `--wait-for verified|quant|de` style orchestration options if needed.
- Add status views for:
  - project-stage progress
  - failed runs
  - pending runs
  - artifact disk usage by class
- Improve `doctor`:
  - verify R packages
  - verify kallisto
  - verify ENA/Ensembl connectivity optionally
  - verify writable root and disk space threshold
- Implement `export` formats:
  - TSV/CSV generated from canonical HDF5 artifacts
  - parquet later if value justifies added complexity
- Implement `survey` output layout:
  - metadata summary
  - counts and normalized matrices
  - DE outputs
  - correlation outputs
  - manifests index

### Acceptance criteria

- A new user can run the documented quickstart without reading source code.
- Operators can tell what the system is doing and why from CLI output plus manifests.

## Stage 7: Extraction Support

### Goal

Support `.sra` fallback paths without polluting the main downloader or quantifier contracts.

### Tasks

- Implement `Extractor` trait concretely for `xsra`.
- Detect `xsra` in `doctor`.
- Detect optional `fasterp` in `doctor` for preprocessing-enabled deployments.
- Route SRA-only runs through extractor before quantification.
- Record SRA input artifacts and extracted FASTQ artifacts separately.
- Add feature-flagged `vdb-native` extractor skeleton behind `ncbi-vdb-sys`.

### Strategy note

Do not attempt native VDB extraction until the `xsra` wrapper path is stable. Native VDB is likely the highest operational complexity area in the whole project and should not block a usable public pipeline.

## Stage 8: Production Hardening

### Goal

Make RNAA believable on long-lived Linux servers under systemd.

### Tasks

- Add structured per-stage logging conventions.
- Add log rotation guidance.
- Add a `systemd` unit for the download worker.
- Add a second unit or timer example for orchestrated analysis.
- Add backpressure rules:
  - quant should not start if disk is critically low
  - reference downloads should not duplicate on concurrent invocations
- Add file locking where needed around shared reference/index creation.
- Add file locking around shared artifact-store writes and GC sweeps.
- Add a `rnaa gc` command for shared-store retention and orphan cleanup.
- Add recovery drills:
  - interrupted download
  - interrupted quant
  - interrupted DE
  - partial manifest write
- Add explicit crash-safe temp naming for every writing stage.

### Acceptance criteria

- Stopping and restarting services does not corrupt state.
- Shared references remain reusable under concurrent server processes.
- Shared quant and matrix blobs remain reusable across projects under concurrent server processes.

## Stage 9: Testing Matrix And Validation

### Goal

Move from "compiles and works on one path" to confidence across the main operational scenarios.

### Test layers

- Unit tests
  - resolver parsing
  - metadata merge
  - design parsing
  - geneset grammar
  - residualization math
  - artifact naming
- Offline integration tests
  - small fixtures for resolve -> download planning -> quant wrapper command construction
  - synthetic DESeq2 input/output validation where possible
  - synthetic correlation fixtures
  - shared-store reuse tests across two projects pointing to the same payloads
- Live integration tests
  - opt-in ENA resolution
  - opt-in Ensembl reference acquisition
  - tiny public run download
- Golden output tests
  - manifest shape
  - samplesheet schema
  - reference bundle schema

### Validation philosophy

For the science-facing stages, correctness tests should focus on invariants and reproducibility, not just file existence. Example:

- DE outputs should have expected contrast labels and matrix dimensions.
- Correlation outputs should have symmetric consistency where applicable.
- Residualization should remove modeled intercept/covariate effects in synthetic fixtures.

## Stage 10: Documentation And Release Readiness

### Goal

Finish the product surface so users can install and operate RNAA without reverse-engineering the source.

### Tasks

- Add a real `README.md` quickstart.
- Add example `rnaa.toml`.
- Add systemd examples under `docs/systemd/`.
- Add fixture-driven tutorial project under `examples/`.
- Document:
  - metadata harmonization workflow
  - design formula rules
  - reference caching rules
  - cleanup/trash semantics
  - manifest schema expectations
- Add changelog and versioning policy.

## Cross-Cutting Refactors To Do Opportunistically

These should happen during the main stages, not as a giant isolated rewrite:

- Introduce project-level stage state separate from run-level state.
- Add stronger artifact uniqueness constraints and reference artifact recording.
- Introduce a shared artifact-store abstraction with content addressing and refcount/pinning metadata.
- Add sample-selector input models (sample accession and metadata predicate filters) with deterministic serialization.
- Add a `Preprocessor` trait and artifact kinds for trimmed FASTQ + preprocessing reports.
- Standardize all stage manifest writing through helper functions.
- Add central error types for user-facing failures vs internal failures.
- Add stable command-runner helpers for external tools.
- Consider replacing ad hoc thread spawning in the download worker with a clearer worker abstraction once more stages are active.

## Priority Order Summary

The immediate next coding order should be:

1. repo hygiene and warning cleanup
2. shared artifact-store schema and storage abstraction
3. sample-scoped add/resolve selectors
4. `refs prepare` CLI wiring and tests against shared storage
5. preprocessing contract (`fasterp`) and trimmed FASTQ artifact model
6. quantification via R + kallisto with HDF5-first persistence
7. cleanup/trash semantics + shared-store pinning/GC policy
8. project-stage schema for DE/corr
9. DESeq2 implementation with HDF5-first outputs
10. Rust residualization + `mincorr` integration with HDF5 matrix contracts
11. `run`, `export`, `survey`
12. extraction support
13. operational hardening, docs, and release packaging

## Definition Of "Production-Grade" For RNAA

RNAA should only be called production-grade when all of the following are true:

- every implemented stage is restart-safe
- project-wide provenance is manifest-backed
- the downloader can run independently for days
- multi-study metadata harmonization is explicit and inspectable
- reference reuse is deterministic
- quant, DE, and correlation each have clear success predicates
- deletion is reversible during trash retention
- operator diagnostics are good enough to debug failures without reading Rust source
- integration tests cover at least one realistic happy path from public accession to correlation output

Until those conditions are met, RNAA is a promising foundation, not a finished survey system.
