# Biological samples and technical runs

RNAA fits one model row and emits one matrix column per biological sample.
The identity is explicit `biological_sample_id` metadata when supplied, otherwise
the sample accession. A missing identity falls back to the run accession and is
labeled `run_fallback`; that label does not establish biological independence.

The DE request records a typed `sample_plan`: each sample ID, identity source,
contributing run IDs and model metadata. Every run must belong to one sample.
The exact abundance/provenance bytes for every contribution remain request inputs.
Adding, replacing, excluding or rejecting one contribution changes the request.
Output column IDs now use biological identities; run IDs remain quantification IDs.

`sample_metadata.tsv` is saved with each DE generation. Its `run_accessions` column
is a JSON array, so IDs containing punctuation remain unambiguous. The original
filtered run-level metadata is also frozen in that generation. Model and contrast
variables must have one identical, nonmissing value across a sample's runs.
Conflicting values fail before fitting. Other varying columns are left blank in
the sample table; their original values remain in the run-level snapshot.

## Pooling policy

Import algorithm revision 3 uses RNAA's `fractional_transcript_mass_v1` policy for
samples with multiple runs. For transcript `t`, estimated fragment count `c`,
effective length `l`, and contributing runs `r`:

```
C[t] = sum_r c[t,r]
M[t] = sum_r c[t,r] / l[t,r]
L[t] = C[t] / M[t]
TPM[t] = 1,000,000 * M[t] / sum_t M[t]
```

All-zero targets use the mean of their positive run effective lengths, with zero
counts and TPM. Entirely zero pooled libraries have zero TPM. Inputs must have
matching transcript IDs and lengths, positive finite effective lengths, and finite
nonnegative counts/TPM. Duplicate normalized transcript IDs and reused run IDs fail.
Traversal is sorted by sample/run ID. A sample with one run preserves that run's
original estimates.

Gene summarization, length correction, scaling and integer rounding occur after
pooling. Supported `no`, `scaledTPM` and `lengthScaledTPM` modes retain their import
contracts. In particular, mean gene lengths are calculated across biological
samples, so additional sequencing runs do not give a sample extra weight in that
mean. Fractional counts are rounded only after transcript and gene aggregation.

This is a defined pooling of separately estimated transcript tables. It does not
rerun quantification jointly on combined FASTQs, combine bootstrap replicates or
claim equivalence to jointly rerunning the quantifier. Closed-form numerical and
split-library invariance tests cover the implementation; independent scientific
oracle fixtures and joint-quantification comparisons remain release work.

DESeq2 describes summing technical-run counts and explicitly leaves combination
of other assays to the caller because it is ambiguous. RNAA therefore records its
length/abundance pooling rule explicitly rather than treating it as a supplied
DESeq2 algorithm. See the [DESeq2 technical-replicate implementation](https://github.com/thelovelab/DESeq2/blob/devel/R/helper.R)
and [tximport's gene summarization](https://github.com/thelovelab/tximport/blob/devel/R/summarizeToGene.R).

## Metadata and correlation

Normalization/DE capture the canonical samplesheet as their metadata source,
preserve the selected rows, and check its hash before publication. Registered
biological identities must agree with any identity explicitly present in that
sheet. Metadata edits now affect the fitted model; the source hash is a request
dependency alongside the frozen model rows.

Default correlation compares the current producer model/membership with the frozen
DE contract, then builds adjustment from current grouped metadata. Changes to used
DE covariates require normalization/DE again; a covariate used only by adjustment
can refit correlation while reusing DE. [Unrelated metadata](scientific-cache.md)
does not change scientific cache identity. A supplied custom
matrix using known biological sample IDs gets grouped current metadata; legacy
run-ID metadata is still supported. Known inactive or QC-ineligible contributions,
ambiguous aliases, and multiple columns representing the same known biological
sample are rejected. External sample identities remain caller-provided.

Factor/numeric roles, missing rules and references are configured explicitly; see
[typed metadata and model designs](model-designs.md). [Controlled refresh](metadata-refresh.md)
preserves local model edits, including across cohort exclusion/reinclusion, and
retains immutable source witnesses. Conflicting source updates fail until resolved
with an explicit policy. Identity and selector corrections use registered overrides
and resolution.
