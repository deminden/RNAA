# Independent scientific references

These expected values come from executing R package APIs. Mandatory offline Rust
tests compare actual RNAA workflows with the saved values; R is not a production
or test-runtime dependency. Missing fixtures, mismatched file/generator/input
hashes, missing numerical values and differences outside the declared tolerances
fail the tests. Expected values are never derived from RNAA output or refreshed
automatically in CI.

The reference families are:

| Family | Independent source | Coverage |
| --- | --- | --- |
| `deseq2` | Upstream R 4.6.1 / DESeq2 1.52.0 fixture at the exact revision in its provenance | Default four-gene fit, parametric failure/local fallback, normalized counts, size factors, Wald results |
| `tximport` | R 4.6.2, tximport 1.40.0, DESeq2 1.52.0 | All three import modes, fractional counts and half-to-even rounding, zero abundance, zero effective length without mass, unmapped mass, versioned GENCODE/PAR_Y IDs, gene lengths and complete offsets; single-run and split-lane native import routes |
| `de_workflows` | R 4.6.2, DESeq2 1.52.0 | Categorical/numeric designs, numeric and factor interactions, both contrast directions, three levels, all-zero/low-expression/outlier rows, Cook's filtering, independent filtering/BH, normalized counts, size factors and VST |
| `correlation` | R 4.6.2 stats, WGCNA 1.74, HellCor 1.3 | All five methods; ties, outliers, zero MAD, constants, two/three samples, near-perfect pairs; dense values, directed top-2 and inclusive thresholds 0, 0.75 and 1 |
| `modules` | R 4.6.2 stats, WGCNA 1.74, fgsea 1.38.0 | Standardized PCA with unequal scales, zero average, tied leading space and three samples; full graph membership/connectivity, report limits, score/rank values and seeded simple-permutation enrichment |

Every family has a `provenance.json` with file hashes, versions, generator identity
and tolerances. `oracle_sources.json` records the source/runtime artifact URLs and
SHA-256 values. The original upstream DE fixture and its MIT attribution are
preserved byte-for-byte, including its session information and generator URL.
The larger workflow generator's input file is also hashed.

Continuous import, normalization, VST and correlation comparisons generally allow
absolute plus relative error of `1e-9` each; import lengths/offsets and size factors
use tighter checks. Larger iterative Wald fits allow `1e-7` absolute plus relative
error: independently computed coefficients agree very closely, while a few
interaction standard errors/statistics differ by roughly `2.3e-8`–`6.5e-8` between
the numerical solvers. The four-gene reference retains its tighter `1e-8` gate.
Count rounding, gene/sample/edge identity, ordering and missingness are exact.
Native DE TSV empty cells correspond to R `NA`.

Correlation `*_raw.tsv` files retain unadapted R estimates. The expected matrix
contract declares every pair involving a constant row undefined, uses a unit
diagonal for variable rows, and snaps coefficients within eight machine epsilons
of either endpoint to exactly -1 or 1. This prevents a floating-point representation
of a perfect relationship from changing top-k ties or threshold=1 decisions.
The same policy is applied independently in the R generator and native kernel;
near-perfect relationships outside that tolerance remain unchanged. HellCor's
raw scalar self-estimator is not necessarily one, and the raw files preserve that
distinction. WGCNA uses `maxPOutliers=1`, individual Pearson fallback and one thread;
HellCor uses Kmax=Lmax=20 and alpha=6 without permutation p-values.

To regenerate the owned references, install the exact pinned R/packages in a
separate development environment and run from repository root:

```sh
Rscript validation/scientific/generate_import.R validation/scientific/tximport
Rscript validation/scientific/generate_de.R validation/scientific/de_workflows fixtures/native_workflow_counts.tsv
Rscript validation/scientific/generate_correlation.R validation/scientific/correlation
Rscript validation/scientific/generate_modules.R validation/scientific/modules
python3 validation/scientific/freeze_references.py tximport
python3 validation/scientific/freeze_references.py de_workflows
python3 validation/scientific/freeze_references.py correlation
python3 validation/scientific/freeze_references.py modules
cargo test -p rnaa-stages --test scientific_reference
```

The freezer is an explicit review step. Review scientific and missingness changes
before resealing; changing hashes alone is not evidence of correctness. Packages
used for the current generation were installed in an isolated temporary library
with official CRAN/Bioconductor artifacts; nothing was installed into RNAA or a
global R environment. Session files record platform and transitive R versions.

The current production DE library stays at rsdeseq2 0.3.1. Its replacement is
deferred until the end at the user's request, when it must pass these gates.
These fixtures do not establish equivalence for every possible DESeq2 setting,
automatic count replacement with seven or more replicates, rlog, every offset
formula or every singular/near-singular dataset. Further DE-engine work belongs
to that separate library effort. Broad scale/release evidence and process-wide
resource reservations remain separate roadmap gates. Module PCA and enrichment
now have an independent gate and [declared policies](../../docs/module-enrichment.md).

The DESeq2 reference session information has its original host-specific BLAS path
redacted. Its provenance records both the upstream hash and the redacted-file hash;
the numerical reference fixtures are unchanged.
