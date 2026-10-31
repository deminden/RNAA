# Shared resource reservation protocol

`rnaa_core::resources` provides a persistent, cooperative CPU/memory/disk admission
ledger for independent processes using the same pool directory. The host-local
backend requires reliable SQLite and kernel-lock semantics; it is not a distributed
scheduler. Fresh owned downloads, preprocessing and quantification use phase
requests, including download publication copies and shared resident indexes.
Reference acquisition, planning and complete build/publication now use the same
pool. Analysis and publication outside downloads/references remain open.

## Reservation contract

Each request names its project, attempt and stage, optional run, CPU threads,
working-memory bytes and future disk writes. A single catalog transaction admits
all dimensions together or leaves the request unadmitted. Requests exceeding their
own CPU/memory capacity fail immediately. Capacity contention returns a structured
wait reason; the blocking API checks the supplied ownership/cancellation fence
before each retry and again before returning an admitted guard. Waiting holds no
catalog transaction. The capacity wait has a caller-specified timeout; catalog
exclusion has a separate bounded 15-second wait. Admission has no FIFO guarantee.

Callers supply capacities. Among live reservations, the smallest CPU and memory
limits and largest disk headroom apply to new admissions. A caller cannot relax
another live caller's contract. A tighter proposal waits if existing usage already
exceeds it, even when the proposed request uses zero of that resource. Limits are
not permanent host configuration: once reservations end, later callers can use
their own limits. Producers must choose a common pool and consistent host policy;
the backend does not discover host limits or select a default pool location.
The `ExecutionResources` policy layer supplies that default for references, downloads, preprocessing and quantification;
see the integration contract below.

The returned guard holds the reservation until explicit release or drop. Resize
checks the entire replacement request atomically. A failed growth leaves the old
reservation unchanged. CPU and memory must remain reserved until their allocations
and workers have ended; a shared index needs a guard that lives with its group,
independently of individual runs. These are accounting promises, not operating-system
memory ceilings or proof that an engine's allocations fit a supplied estimate.

`try_reserve_bundle` and `reserve_bundle` admit multiple claims in one transaction.
Their summed CPU/memory/disk promises must fit before any row is inserted. Every
claim receives an independent owner guard after commit, allowing lifetimes to
diverge. A failed insert rolls back the whole bundle; unlocked empty owner files
are recovered by the existing orphan protocol. Rollback never invokes a releasing
guard while the admission transaction holds the catalog lock.

Disk bytes mean additional writes still promised, not existing occupied storage.
Reduce that promise only after corresponding writes are physically visible, so
free-space observations account for them. Releasing the promise after completion
does not manufacture free disk: later admissions inspect current free space.
Future output paths resolve to their nearest existing filesystem ancestor.

Disk admission sums outstanding promises across the pool and subtracts them plus
headroom from available space at every new destination. Separate paths and device
identifiers can share backing storage, including APFS volumes and overlays. The
conservative sum prevents treating such paths as independent copies of free bytes;
it can also underuse genuinely independent filesystems. This does not implement
project active-size quotas, reserve space against unrelated applications or replace
the filesystem's own full-disk errors. Consumers must reserve before writes and
bound staging/output growth. Free-space probing currently supports Unix hosts and
uses `statvfs` available blocks, respecting filesystem availability to the caller.

## Ownership and recovery

`catalog.lock` serializes catalog initialization, WAL selection, admission, resize,
release and recovery. SQLite commits with full synchronization. Each admitted row
also has its own UUID-named, empty owner lock file, held exclusively by the guard.
The descriptor binds that file's device and inode. Project/run input or scientific
catalogs are not changed by resource admission.

A live, paused or slow owner retains its kernel lock. No clock, PID lookup or lease
expiration can reclaim that reservation. Opening a second handle probes the exact
recorded file identity and tries to acquire its lock. Only an available lock proves
the original owner is gone. The row is then retired transactionally before the old
empty lock is removed. New admissions and snapshots recover such abandoned rows.
Process termination releases locks even when Rust destructors do not run.

Missing, replaced, aliased or symlinked ownership evidence is held as an error;
absence does not establish an owner's death. Changed descriptors and unsupported
catalog schemas also fail explicitly. Owner release can retire its own unchanged
record but does not unlink a replacement file at the old pathname.

`ResourcePool::recover` additionally inspects the private owner directory for
empty, unlocked, canonical UUID lock files with no catalog record, such as files
left before a failed catalog commit. It preserves active locks, unknown names,
nonempty diagnostics and unrecognized evidence. Recovery operates only on pool
coordination metadata and never deletes project payloads. This local crash protocol
does not establish power-loss qualification.

## Evidence and remaining integration

Tests cover atomic multi-resource admission, independent handles, restrictive
policies, zero-dimension policy changes, resizing, disk growth, rollback, orphan
maintenance, missing/replaced ownership proof, cancellation before/after admission,
timeouts and unsupported schemas. Eight competing threads and six separate child
processes pass multiple waves through the same memory budget. A real child is killed
while holding capacity; a live child remains protected and the dead owner's row and
lock are recovered before replacement admission.

The remaining producer changes must reserve from this backend before heavy preparation,
not merely before the final numerical call. They must cover raw validation and
decoded temporary FASTQs, reference/index ownership, worker-local caches and counts,
normalization/design preparation, matrix selection/adjustment, network/module
phases and final publication. Reuse should avoid unnecessary heavy claims. Actual
supported-host measurements remain necessary to qualify those estimates and the
combined workflow. The complete operational resource milestone remains open.

Decoded-group preflight and bounded native preprocessing are now delivered.
`FastqPlan` counts all decoded mates/lanes, caps lines and reports exact normalized
staging bytes before temporary payload creation. Materialization cannot exceed its
preflight promise. Native preprocessing acquires a credit before producing each
batch and holds it through queued work, execution and ordered-result writing;
paired input is read together without independent unbounded carryover. See
[FASTQ preparation](fastq-preparation.md) for the precise limits and evidence.
Fresh preprocessing now uses these bounds to obtain actual shared reservations,
as described below. Quantification also uses actual scan/index/run claims. Remaining producers require phase estimates and guards
covering actual allocation lifetime. No measured RSS ceiling is established.

## Current preprocessing integration

`[resources]` defaults to the common per-user directory
`$HOME/.cache/rnaa/resources/v1`, available CPU parallelism, 75% of detected memory
capacity, 1024 MiB disk headroom and a 300-second capacity wait. A configured pool
must be absolute. Use a local directory shared by the participating processes;
different pool paths do not coordinate. This is per user, not cross-user admission.
The home directory default assumes local storage. CPU worker budgets continue to
limit each command separately.

Automatic memory detection supports macOS and Linux. Linux also intersects the
standard cgroup v2 `memory.max` and v1 `memory.limit_in_bytes` locations when present.
Set `memory_mb` explicitly for other/nested mount layouts or a different operating
system. Capacity is not currently free RAM; unrelated applications and unconnected
RNAA producers are outside this ledger. Configured limits are cooperative contracts,
not operating-system controls.

The production preprocessing path first admits one nominal CPU and bounded-line
scan memory. It releases that claim after preflight before waiting for the complete
engine request. Waiting while retaining a smaller claim could deadlock concurrent
workers. The engine request includes nominal native worker threads, a structural
memory estimate, exact normalized temporary bytes and bounded future output writes.
Its guard lives through staging, native thread joins, source verification and the
engine's report/manifest writes. Temporary payloads drop before the guard does.
Generation sealing/publication and other pipeline stages still need admission.

The `fasterp-native-defaults-v1` estimate includes Vec spare capacity and
original/trimmed buffers, record offsets/result pieces across the full in-flight
window, before/after position arrays and merger state, report trees/serialization,
fixed decoder/writer overhead and stack allowances. It uses decoded shape rather
than a compressed-file multiplier. Native input enforces preflight total bytes,
record count and conservative minimum/maximum record sizes; position/report
allowances conservatively use the maximum record size. Report/manifest writes
also share a checked metadata byte cap. Current trimming defaults preserve or
shrink normalized FASTQ bytes, with UMI/header expansion disabled. Changing these
defaults requires revisiting the estimate and write contract.

Memory constants remain conservative estimates requiring supported-host RSS
qualification, including allocator overhead. They do not establish a process ceiling.
Disk promises remain conservative until the engine returns; visible writes are
not progressively subtracted. Healthy historical/current reuse precedes admission,
including when the pool is unavailable. Typed admission errors cannot become an
engine-error bypass. Ownership checks cancel capacity waits and periodic decoded
scans; stopping an active numerical/native engine promptly remains separate work.

The ingest scheduler now cancels and joins remaining workers when its driver exits
with an error. It does not detach them and return while their allocations or owned
attempts remain live. A native computation may still delay that join until it ends.

## Current quantification integration

Owned CLI quantification performs healthy historical/current reuse before opening
the pool. Fresh work admits bounded-line scan memory and one nominal CPU before
index/FASTQ inspection. It releases that claim before waiting for heavy capacity.
Mapping threads are clipped to the shared policy's CPU capacity before estimating
memory and cache tables. Admission happens before temporary FASTQs or output
directories are materialized.

Cold loading atomically admits two claims: the index and the first run. The total
memory promise covers the larger of loading plus retained setup, or resident graph
plus run state. The first run's allowance remains unused during loading and covers
part of that transient peak. Once loading ends, the index claim shrinks to its
resident estimate. This avoids holding an index while waiting to admit its first
run. Later runs reserve only their own working CPU/memory/staging/output allowance.

The graph owns its memory guard and is freed before the guard drops. Active runs
retain the graph through an `Arc`; each has an independent run guard. The index
claim may name an already finished first attempt: its kernel lock, rather than
the attempt's status or a clock, establishes allocation ownership. A group cannot
mix different pools or admitted and unadmitted loads.

On a memory wait, a group may evict its own idle graph when the pool has no CPU
owners left to make progress. It does not evict active graphs or evict solely for
CPU/disk contention. Waiting retains no borrowed graph reference. The ingest
scheduler also releases idle graphs when it has neither ready nor active quant
work, allowing upstream preprocessing to feed later waves. Ordinary waves retain
one graph when capacity permits; intentional eviction can require another load.
Waits check cancellation/ownership, report their first capacity reason and honor
the configured timeout. There is no cross-process shared graph allocation.

The `kallistors-mapping-em-v1` model includes packed-read queues, EC keys/counts,
sparse EM state, target arrays, lookup caches, worker stacks, decoded staging and
output promises. Paired parallel mapping groups parser batches into logical jobs
of approximately 8 MiB of sequence data; the minimum decoded sequence length bounds
how many preallocated batches can coexist. Readers enforce preflight minimum and
maximum sequence lengths, line/record/count and total-byte bounds before engine
batch allocation. Unsupported 32-bit counter/span requests fail explicitly.

The engine restores an 8 MiB default thread-local match table on scope exit even
when a smaller budget was requested. The estimate includes that restored capacity;
mapping runs in a joined scoped thread so those private caches die before the run
guard releases. Temporary FASTQs, mapping/EM allocations and engine output writes
also finish under that guard. `parameters.resource_execution` records the estimate
and index/run reservation descriptors without changing scientific request identity.

EC class/member estimates use conservative subset envelopes of decoded index
classes and observed fragments. Shade references use a target-union envelope.
Complex ambiguity can therefore reject a workload whose actual observed EC
diversity would fit. These constants need supported-host RSS and production-scale
qualification; they do not establish a hard process memory ceiling. The compatible
library API without an `ExecutionResources` context retains per-job preflight and
does not enforce shared run admission. Production owned CLI work supplies context.

Real-index tests cover shared ownership, active-run protection, cancellation during
loading, idle eviction and mapping parity. Six independent processes contend with
two-claim bundles without exposing partial admission. Two actual CLI processes
share one pool: cancellation during a confirmed capacity wait releases no partial
claim and the surviving project finishes. A constrained preprocessing/quantification
workflow finishes all worker waves. Shared-blob copying, generation sealing,
quant publication/retention and analysis still need their own
admission. Prompt interruption inside active engine calls remains separate work.

## Current download integration

Fresh owned downloads, including `recover_inputs`, use the same host policy and
pool as preprocessing/quantification. Historical/current reuse precedes opening
the pool; verified checksum-backed history remains reusable even with an
unavailable pool or a newly lowered memory limit. The compatible `Downloader`
trait method without owned job context retains its existing unadmitted behavior.
Reference-source HTTP now uses this admitted transfer protocol as part of the
reference workflow described below.

Phases release their claims before waiting for the next phase:

1. Private-cache verification reserves one CPU, descriptor-dependent bounded
   receipt memory and metadata writes. Receipts are read through a byte cap and
   still need their original source/digest/version evidence.
2. Resume planning reserves one CPU, 16 MiB and 4 MiB of future coordination
   writes. It checks prefix identity/hash and only trusts the durable checkpoint,
   not a partial file's raw length. Prefix hashing checks ownership/cancellation
   every approximately 8 MiB; cancellation cannot reset a valid prefix.
3. Transfer reserves one CPU, 64 MiB for bounded streaming/client/hash state and
   future file writes plus 4 MiB of metadata. A known total reserves only bytes
   after the proven prefix. Unknown totals start with an 8 MiB write window and
   check capacity before extending it; this is not a total download-size limit.
4. Completed private files are verified and receipted under a fresh verification
   claim. A complete FASTQ group then obtains the same bounded-line scan request
   used by preprocessing and quantification. SRA payloads do not claim FASTQ
   validation or extraction.
5. Publication reserves one CPU, inventory-dependent memory and the full extra
   copy plus metadata bytes before beginning a generation. Disk requests name
   both the generation directory and the project catalog, so each destination's
   available space is checked even when they reside on different filesystems.
   The backend's conservative shared-backing sum still applies. Its guard outlives
   copying, hashing, manifest writes, sealing, registration and transfer-cache
   removal. The data-copy reader cannot write beyond its preflight byte promise;
   truncation/growth or cancellation discards its temporary output. Source hashes
   and the existing source/cohort/attempt fences remain required before commit.

When headers reveal a total that cannot fit the existing window, the response and
all CPU/memory/disk guards unwind before readmission waits for the complete size.
The next request runs under that full allowance. If a server ignores Range, the
private prefix is reset and its real storage reclaimed before reserving the full
replacement. Known-size promises shrink after flushed checkpoint progress.
Unknown-size growth uses a nonblocking resize: refusal checkpoints the proven
prefix and exits with a typed admission error, releasing all claims. It does not
wait while retaining a response or treat capacity refusal as another network retry.
Retrying later uses the original validator-aware resume protocol.

Request/retry buffers finish before a capacity wait or network backoff. Capacity
waits honor both the configured pool timeout and the transfer deadline, as well as
the current job fence. The generic blocking admission path logs its first wait
reason. HTTP checkpoint JSON is capped at 1 MiB; invalid resume evidence cannot
authorize reuse. A premature EOF, source change or failed admission cannot publish
a partial generation.

`parameters.resource_execution` records download publication estimates/ownership
and per-file transfer observations. Resource settings remain outside download
request revision 2. Publication capacity failure preserves the existing head and
the private verified cache; a larger-capacity retry can publish without refetching.
These are conservative allocation/write models needing supported-host RSS and
combined-workflow qualification. Download copies still coexist with their transfer
cache until commit, and their disk promise is conservatively retained for the whole
publication phase. Unrelated processes and unadmitted producers remain outside the
ledger; complete analysis admission and other producer publication remain open.

## Correlation resource coordination

After checking for a reusable result, fresh correlation work reserves its selected
worker count in the shared pool. The count is capped by the configured shared CPU
capacity. After selection, adjustment releases that initial claim and requests
CPU, estimated memory and space for the adjusted matrix together. It rechecks the
estimate after admission, processes the input a row at a time, bounds input and
output bytes, and releases the claim after writing. The output estimate assumes
worst-case floating-point formatting, so it can reject work whose actual output
would fit. The separately written model-design JSON is not included in this disk
promise.

Before network computation it releases the intervening CPU claim and requests
CPU, the engine's estimated working memory plus 64 MiB, and the engine's estimated
primary output size together. The engine repeats preflight and refuses a changed
estimate after admission. Its output writer refuses to exceed the promised byte
count. Capacity waits observe the job cancellation fence and do not hold CPU
while waiting for memory. Network execution also checks cancellation between
input-preparation steps and periodically between tiles; an individual read,
preparation step or tile is not interruptible. Cancellation preserves the prior
published generation. The primary output promise is released after it is written;
subsequent disk checks see the space actually consumed.

These are cooperative estimates, not operating-system limits. Input discovery,
cache checks, selection, model-design JSON, module prevalidation and enrichment,
annotation, manifests and publication do not yet have complete shared memory or
disk claims; the network output promise covers only the primary network table.
Some preflight and design construction reads inputs before the relevant claim.
DE analysis and other analysis publication paths remain outside this coordination.

For an external matrix, sample selection reads its header rather than loading the
numeric matrix before reuse. A size-bounded file witness is checked against the
scientific request and again before reuse or publication. These observations detect
ordinary input changes; they do not lock out an external writer. With an `all`
gene set, adjustment performs the first full matrix load and validation, avoiding
an earlier redundant load. Filtered selections scan and validate one row at a
time, retaining gene identifiers and selected row positions rather than a full
numerical matrix. Top-variance selection keeps at most the requested number of
ranked positions; it preserves stable ties and output order. Adjustment streams the selected rows; later correlation loads the selected
matrix, and the identifier set grows with input rows. Selection still precedes the adjustment claim; later network
and module work have their own allocation behavior.

## Current reference integration

Fresh owned reference work acquires source verification/transfer, decoded planning
and native build/publication claims. Healthy current/historical reuse happens
before opening the pool. Bounded provider discovery and local source hashing still
precede reference identity/ownership; provider listings/checksum catalogs have a
1 MiB decoded cap. This metadata work does not establish a process-wide resource
ceiling. The preparation callback checks both the parent analysis owner, when
present, and the independent reference attempt.

`refs.max_source_line_mb` defaults to 1 MiB per decoded FASTA/GTF line.
`refs.max_decoded_mb` defaults to 0, with no aggregate cap; a positive value applies
across the full decoded FASTA and GTF group. New work requires UTF-8 headers/GTF and
ASCII sequence symbols. Long unwrapped sequences need an appropriate line setting.
The scan counts sequences, bases, names, k-mer positions, annotation rows/strings,
decoded hashes and largest lines without retaining all sequences or annotations.

FASTA planning counts distinct clean canonical k-mers exactly using bounded sorted
runs. Case and U/T normalization match the native builder; windows containing other
symbols contribute a conservative additional cardinality allowance. Poly-A trimming
can lower the builder's actual count. The planner retains a 4 MiB key buffer,
merges at most 32 runs at once and compacts levels online. It reserves one CPU,
64 MiB plus sixteen line allowances, and checked 16 MiB future-write windows plus
4 MiB metadata on the generation filesystem. Buffered bytes flush before renewing a
window. Refused scratch growth fails without retaining CPU in a capacity wait;
all owned temporary files and buffers drop before the guard. Final graph/cardinality
observations escape, not the sorted runs. GTF planning uses a separate bounded claim.

A small admitted setup step creates the generation coordinator and private
`.rnaa_reference_work` directory before scanning. Sort runs and decoded staging
live inside that owned unpublished container. Normal completion/cancellation
removes them; after a killed process, kernel-lock recovery releases its claims and
the existing unpublished-generation preview/trash/purge workflow can recover the
physical files. They are not untracked system-temporary data. Successful sealing
requires removal of the work directory.

Only after planning ends does the workflow admit the native phase. Its memory is
the maximum of decoded graph/occurrence/serialization and annotation/validation
allowances, plus metadata; serialized input file size is not a graph multiplier.
Distinct clean keys reduce excessive estimates for overlapping transcripts. Bounds
still use conservative graph fragmentation and occurrence envelopes and may reject
workloads with lower actual use. Requested index threads are clipped to pool CPU
capacity; manifests record requested and effective values. Constants and planner
overhead require supported-host RSS and production throughput qualification.

One native-phase guard covers bounded original-source copies, annotation tables,
decoded FASTA staging, index construction/inspection, mapping validation, manifests,
generation sealing/registration and successful transfer-cache removal. Its future
disk promises cover the generation store, project catalog and decoded staging
filesystem before source copies, staging or native output are materialized. The
small coordinator setup has its own earlier claim. Memory/file promises
are retained conservatively through the phase. Index report dimensions, serialized
index bytes and table sizes must fit the plan before publication. A joined scoped
index thread ends before the guard releases. Active numerical calls still need
prompt cancellation support; waits, planning, copies and annotation work use fences.

Compressed FASTA is decoded completely into bounded temporary plain FASTA before
calling the present builder, which otherwise reads only its first gzip member.
Staging must match the planned decoded byte count and hash. GTF also reads every
gzip member. Reference request revision 2 records
`complete_multimember_fasta_gtf_v1`; old revision-1 bundles rebuild once, retaining
their original files/history. Byte-identical rebuilt indexes preserve count-binding
identity. Limits, thread allocation and reservations remain execution observations.

Transfer caches include the reference subject in their ownership key, so one
reference cannot remove another reference's identical-URL source during cleanup.
Private copies and original content hashes remain checked across preparation.
Real CLI tests confirm two projects wait before heavy work, cancellation leaves no
partial claim, the surviving project publishes, and cancellation during external
planning removes temporary runs. Capacity rejection preserves prior bundles;
healthy historical selection still works with a lowered limit and unavailable pool.
An actual killed planner is also recovered through unpublished trash and expired
purge without publishing a reference or deleting its original source.

To repeat the local synthetic profile with Python 3.11+ and a current SQLite runtime:

```sh
cargo build --release --all-features -p rnaa
python3 scripts/benchmark_reference_resources.py \
  --transcripts 1000 20000 --threads 1 2 \
  --report target/release-evidence/reference-scale.json
```

Each case has a fresh process, project, source fixture and private pool. The report
records full-process peak RSS, elapsed time, declared build memory, actual distinct
keys/node counts, output sizes and source/binary/toolchain witnesses. It checks
scientific request/mapping/graph-count consistency across thread allocations and
that no resource rows or temporary files remain. Input generation and project init
are excluded from timing; acquisition, planning, native work and publication are
included. These overlapping synthetic references do not establish production or
cross-platform release qualification, and process RSS includes baseline overhead
outside the individual cooperative phase claims.
