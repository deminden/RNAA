# RNAA Development Plan

## Purpose

This document defines the next development stages for RNAA after the current foundation work:

- workspace scaffolding exists
- SQLite-backed project state exists
- `init`, `add`, `resolve`, `status`, `doctor`, and `download` are functional
- ENA resolution works for runs/studies/projects and generates a canonical samplesheet
- the downloader is restart-safe, manifest-driven, and independent of later stages
- `refs`, `quant`, `deseq2`, `corr`, and `run` are now implemented, but not yet hardened enough for production use
- `export`, `survey`, cleanup execution, extraction fallback, and worker-budget scheduling are still incomplete

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

## Current State Snapshot

### Working now

- Workspace layout and crate boundaries are in place.
- `rnaa-core` has config, models, DB, state enums, paths, manifests, traits, and samplesheet logic.
- SQLite migrations exist and the state/event/artifact tables are usable.
- ENA resolver fetches `filereport` JSON and expands run metadata.
- Metadata overrides merge into resolved metadata.
- Canonical samplesheet and editable column map are generated.
- Download worker uses `curl`/`wget`, `.part` files, retries, verification, and manifest output.
- Reference preparation is wired through the CLI and records reference artifacts.
- Quantification works through the Rust wrapper plus `Rscript` plus kallisto.
- DESeq2 works through the Rust wrapper plus `Rscript` plus tximport/DESeq2.
- Correlation works through Rust residualization plus linked `mincorr`.
- `rnaa run` orchestrates the implemented stages sequentially.
- Tests exist for resolver parsing, metadata merge behavior, state persistence, and download command templating.

### Incomplete or still weak

- Export and survey are still CLI stubs.
- Cleanup policy exists in config and CLI flags, but real move-to-trash plus purge execution is not complete.
- Extraction fallback (`xsra` and `vdb-native`) is still not implemented.
- Project-stage state is still overloaded onto per-run state for DE and correlation.
- Reference preparation, quantification, DESeq2, and correlation still need live integration tests and harder failure-path coverage.
- There is no crash-recovery sweep yet for stale `*_RUNNING` states after abrupt process death.
- There is no worker-budget scheduler yet; only download is truly parallel at the RNAA scheduler level.
- The repo contains nested `.git` directories under crate folders from `cargo new`; that should be cleaned before real repository hygiene work proceeds.
- There is still no real `README.md` or `docs/systemd/` operational surface for a first public commit.

## Parallelization Snapshot

The current execution model is mixed:

- Download is parallel across runs. `rnaa download --concurrency N` spawns `N` Rust worker threads, and each worker processes one run at a time.
- Quant is not parallel across runs yet. The RNAA CLI loops through runs serially, and each run gets one `Rscript`/kallisto process with `config.quant.threads` threads.
- DESeq2 is single-project and single-invocation from RNAA. That is correct at the scheduler level, but RNAA does not yet explicitly cap BLAS/OpenMP threads.
- Correlation is single-project from RNAA. Residualization is currently serial over genes, while `mincorr` itself uses internal parallelism.
- `rnaa run` is a sequential stage orchestrator, not a global worker-budget scheduler.

If you give the current code a 100-worker host today, only download can use 100 RNAA-managed workers. Quant can use 100 threads only by giving a single run all 100 kallisto threads, which is not the right long-term strategy. DESeq2 and correlation each run as one project job.

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

4. Extraction and FASTQ-prep should behave like quant-adjacent run-parallel work.
- `xsra` or future native extraction should use the same CPU pool family as quant.
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
- let quant consume the CPU pool aggressively by run-level parallelism:
  - default `threads_per_run = 8`
  - default `max_concurrent_runs = floor(96 / 8) = 12`
  - if fewer than 12 runs are ready, increase per-run threads up to a cap such as `16`
- allow the downloader and quantifier to run at the same time:
  - example steady state: `2` download workers plus `12 x 8-thread` quant jobs
  - this is valid because download is mostly I/O bound and quant is CPU bound
- run DESeq2 as exactly one project job with a hard thread cap such as `8` or `16`
- run correlation as exactly one project job with an explicit pool such as `48` or `64`

That gives the behavior you want:

- alignment uses the whole machine in aggregate
- downloads keep progressing while alignment is already processing earlier runs
- DESeq2 does not explode into unsafe parallel jobs
- other stages use bounded, stage-appropriate worker counts
- run-parallel work happens where the contracts are naturally per-run

## First Commit Readiness Gaps

The first commit should not present RNAA as production-ready yet. The minimum things still not prepared for a credible first public commit are:

1. Repository hygiene
- remove nested crate-local `.git` directories under `crates/`
- decide whether this is the actual top-level repo root and normalize git state

2. Operator surface
- add a real `README.md`
- add a real quickstart that matches the implemented commands
- add `docs/systemd/` examples for the download worker and orchestration

3. State model correctness
- add project-level stage state instead of overloading run-level state for DE and correlation
- add startup recovery that converts stale `*_RUNNING` states into retryable interrupted failures

4. Worker scheduling
- implement the worker-budget model described above
- make download and quant overlap by default instead of serializing the project behind a stage boundary
- add explicit thread caps for DESeq2 and correlation
- make quant parallel across runs instead of serial-at-the-RNAA-layer

5. Honest feature surface
- either implement cleanup, export, survey, and extraction, or clearly leave them out of the first public promise
- do not keep CLI flags that imply working behavior if the underlying stage is still a placeholder

6. Testing
- add integration tests for:
  - refs prepare
  - quant wrapper happy path
  - DESeq2 happy path
  - correlation happy path
  - end-to-end `rnaa run`
- add restart/interruption tests for download, quant, and project-level stages

7. Hardening
- add reference-build locking so concurrent processes cannot trample shared cache creation
- add output guardrails for dense correlation outputs
- verify manifest completeness and failure-mode consistency across stages

## Recommended Delivery Sequence

The smart path is:

1. Finish reference preparation and CLI wiring.
2. Finish per-run quantification via R + kallisto.
3. Finish project-level DESeq2.
4. Finish adjusted correlation in Rust + `mincorr`.
5. Finish orchestration, export, and survey UX.
6. Add extraction and cleanup policy completion.
7. Harden operational behavior, tests, docs, and release packaging.

This sequence is correct because each later stage depends on validated outputs from the earlier one, and because reference/quant/DE/correlation form the real scientific pipeline backbone.

## Stage 1: Foundation Hardening And Repo Hygiene

### Goal

Turn the current scaffold into a clean base before deeper pipeline work multiplies complexity.

### Tasks

- Remove nested `.git` directories under `crates/` created by `cargo new`.
- Add a real top-level repository if one is intended here.
- Normalize root docs structure:
  - `README.md`
  - `docs/architecture.md`
  - `docs/state-machine.md`
  - `docs/artifact-contracts.md`
- Add CLI smoke tests for:
  - `init`
  - `add`
  - `resolve`
  - `status`
  - `doctor`
  - `download`
- Tighten `doctor` output into machine-readable and human-readable sections.
- Add a stable version field to manifests and start documenting compatibility expectations.
- Add helper utilities for atomic write patterns used by future stages.

### Why this matters

Without hygiene now, every next stage will repeat avoidable mistakes:

- inconsistent manifest writing
- inconsistent error semantics
- ad hoc path handling
- poor operator experience on Linux servers

### Acceptance criteria

- `cargo check` and `cargo test` are clean with no avoidable warnings.
- Core CLI commands have regression coverage.
- Repository layout no longer contains accidental nested git repositories.

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
- Add checksum verification for downloaded reference files.
- Parse Ensembl `CHECKSUMS` when feasible instead of trusting only local sha256.
- Handle reference reuse:
  - skip rebuild if manifest and artifacts already match
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
- Reference manifest and artifacts are sufficient for downstream reproducibility.

## Stage 3: Per-Run Quantification Via R + Kallisto

### Goal

Implement the first real quantification engine with clean Rust/R boundaries and safe per-run state transitions.

### Tasks

- Create `r/quant_kallisto.R`.
- Add robust argument parsing in R.
- Validate:
  - one FASTQ for single-end
  - two FASTQs for paired-end
  - kallisto index existence
  - output directory writability
- Run `kallisto quant` via `system2`, never shell concatenation.
- Capture:
  - tool versions
  - command arguments
  - reference manifest hash
  - input FASTQ checksums
- Define success predicate in Rust:
  - output directory exists
  - `abundance.tsv` exists and is non-empty
  - `run_info.json` exists and is parseable
- Record quant artifacts into SQLite.
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
- Cleanup only triggers after verified quant success.

## Stage 4: Project-Level DESeq2 And Normalization

### Goal

Turn multiple per-run quant outputs into a coherent project-level cohort analysis.

### Tasks

- Create `r/tximport_deseq2.R`.
- Feed it:
  - samplesheet
  - quant directories
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
  - `gene_counts.tsv`
  - `gene_norm_counts.tsv`
  - `vst.tsv`
  - `vst.rds`
  - `de_<contrast>.tsv`
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

## Stage 5: Adjusted Correlation In Rust + `mincorr`

### Goal

Produce biologically interpretable, metadata-adjusted correlation outputs without routing the core math back through R by default.

### Tasks

- Implement matrix loading from:
  - VST TSV
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
- Persist adjusted matrix artifact and manifest.
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
  - TSV/CSV first
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
- Add recovery drills:
  - interrupted download
  - interrupted quant
  - interrupted DE
  - partial manifest write
- Add explicit crash-safe temp naming for every writing stage.

### Acceptance criteria

- Stopping and restarting services does not corrupt state.
- Shared references remain reusable under concurrent server processes.

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
- Standardize all stage manifest writing through helper functions.
- Add central error types for user-facing failures vs internal failures.
- Add stable command-runner helpers for external tools.
- Consider replacing ad hoc thread spawning in the download worker with a clearer worker abstraction once more stages are active.

## Priority Order Summary

The immediate next coding order should be:

1. repo hygiene and warning cleanup
2. `refs prepare` CLI wiring and tests
3. quantification via R + kallisto
4. cleanup/trash semantics
5. project-stage schema for DE/corr
6. DESeq2 implementation
7. Rust residualization + `mincorr` integration
8. `run`, `export`, `survey`
9. extraction support
10. operational hardening, docs, and release packaging

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
