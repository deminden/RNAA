# FASTQ preparation and native preprocessing bounds

Download validation, preprocessing and quantification use the same `FastqPlan`
input contract. It validates complete single/paired lane groups, rejects duplicate
paths and mixed raw/trimmed inputs, and checks records, mate counts and identifiers.
Every gzip member is decoded, including footer validation. Inputs must be regular
files. Source content is hashed before inspection and checked afterward and before
materialization; stages also check frozen inputs before publication. These are
content-change checks under cooperative immutable input use, not filesystem snapshots
against an adversarial writer.

## Decoded limits and exact staging

All size settings below are MiB in `[quant]`:

| Setting | Default | Contract |
| --- | --- | --- |
| `fastq_max_line_mb` | 1 | Positive limit for each decoded line, including its terminator. |
| `fastq_max_decoded_mb` | 0 | Total decoded wire bytes across all mates/lanes; 0 disables this cap. |
| `preprocess_max_input_mb` | 1024 | Additional aggregate decoded cap for preprocessing; 0 disables this component. |

`--preprocess-max-mb` sets the last value. Its aggregate decoded meaning replaces
the previous per-file encoded-size check. If both aggregate caps are configured,
preprocessing applies their minimum. Checked conversions reject overflow. Line
reading is bounded before allocation; a huge unterminated line cannot bypass the cap.

Preflight scans the whole group without creating temporary payloads. Its summary
distinguishes decoded wire bytes, which include CRLF, from normalized LF output
bytes. Multiple lanes, mismatched compression naming, CRLF or a missing final
newline require normalized staging. A compatible single-file group can remain
at its original paths, reporting zero temporary staging bytes.

`FastqPlan::summary().staging_bytes()` provides an exact temporary byte promise.
Materialization checks sources first, then creates a private temporary directory
and scans again while writing. It caps decoded and normalized writes to the
preflight counts regardless of optional user caps. Changed counts/content fail;
temporary files are removed on errors and on normal drop. This adds scan/hash I/O,
especially for compressed multi-lane input. Ownership/cancellation checks are
periodic during decoding rather than a database operation on every line.

Healthy cache reuse precedes heavy preparation, so lowering a resource limit does
not invalidate an existing scientific result. Preprocessing's scientific request
revision advances to 3 for the paired reading and newline fixes. Quantification's
request revision stays unchanged; observed execution summaries and resource limits
do not enter its scientific identity.

## Native work window

RNAA uses `fasterp::process_paths_with_options`. Parallel native preprocessing
defaults to 32 MiB batches and at most `threads + 1` in-flight batches. An input
batch obtains a credit before allocation and holds it through queueing, worker
execution and ordered-result merging. An early slow batch therefore cannot allow
the pending-results map to grow beyond the same window. Payloads drop before
their credit is returned. One paired fragment is read together, so unequal mate
lengths do not accumulate independent decoder carryover.

Each batch buffer can exceed its target by one bounded complete record. RNAA
derives the native record bound from preflight, with checked newline allowance.
Both paths enforce the planned total decoded bytes, record count and conservative
minimum/maximum record sizes. The reader returns at record boundaries to avoid
rescanning an entire decoder buffer for each short FASTQ line. Options are validated before
output creation. Scoped producer/worker/merger threads are joined; errors and
panics cancel blocked credit acquisition. Execution reports expose options and
peak/active credits. Successful completion has zero active credits.

The window bounds payload lifetime, not all engine allocations or process RSS.
Filter state, output/compression buffers, index allocations and other producers
need their own measured phase estimates. Fresh CLI preprocessing now reserves its
scan and engine phases from the shared pool. Quantification, index loading and
other producers remain unconnected, so no pipeline-wide CPU/RAM/disk guarantee
follows. See [resource reservation integration](resource-reservations.md).

## Verification

FASTQ unit tests cover exact wire/normalized counts, aggregate compressed lanes,
line/staging limits, source mutation before/during preparation, cancellation,
corrupt gzip footers and temporary cleanup. Native tests compare serial/parallel
single/paired output and QC, exercise CRLF/multiple gzip members, hold a reordered
result's credit, and bound producer/writer failure shutdown in child processes.
CLI tests reject aggregate groups, allow them after raising limits, and preserve
healthy cached heads/attempts after limits are lowered again.
