# Native workflow fixture

`generate.py` creates `fixtures/native_workflow_counts.tsv` using Python's standard
library, seed 20261029, and a gamma–Poisson model. There are 500 genes and eight
independent samples, with four samples in each condition. Genes G0000–G0049 have
fourfold higher expected counts in B; G0050–G0099 have fourfold lower expected
counts. The remaining genes have equal expected counts across conditions.

Regenerate from the repository root:

```sh
python3 validation/native-fixture/generate.py
```

`native_analysis_integration.rs` creates real kallistors HDF5/TSV abundance files
from these counts and runs native DE, residualization, correlation and module
enrichment. It checks the expected effect directions, successful enrichment,
unchanged-request reuse, and invalidation after changes to abundance, metadata,
gene-list content, pathway content and output files. Run it with:

```sh
cargo test -p rnaa --test native_analysis_integration
```

The same suite also exercises successful normalization -> correlation, complete
generation reuse, preservation of the previous successful generation after partial
DE/correlation failure, and rejection/repair of a damaged upstream matrix. Core
generation tests separately kill real child workers around the publication boundary
and check rollback of head, artifact registry, state, events and job outcome.

These are synthetic workflow and regression checks. They do not independently
validate rsdeseq2 statistics, kallistors quantification or parity with R/DESeq2.
Pinned independent scientific reference fixtures remain a separate release gate.
