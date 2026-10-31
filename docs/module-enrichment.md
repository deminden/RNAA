# Module scores and enrichment

`rnaa corr --module-gmt pathways.gmt` uses a mutual positive Spearman neighbour
graph. Each variable gene selects its strongest 20 positive neighbours (or fewer
when unavailable). Edges require mutual selection; components meeting
`--module-min-size` are the full modules. This is RNAA's graph algorithm and does
not claim to reproduce WGCNA clustering.

Module construction always uses Spearman. The ordinary network independently
honors `corr.method` or `--method`; enabling modules no longer silently overrides it.

Scores are the first principal component of **all genes in the module**, after
centering each gene and dividing by its sample standard deviation. Scores have
mean zero and sample variance one. Their sign aligns with standardized average
expression. If that direction has zero projection, the first maximal sample
coordinate within numerical roundoff is positive. A tied leading eigenspace uses
its first nonzero canonical sample-axis projection before alignment. Numerical
ties use the recorded 64-machine-epsilon/sample-count rule. PCA failures fail the
attempt; there is no arithmetic-mean fallback.

The score convention is checked against independent R `prcomp` and normalized
WGCNA module-eigengene outputs. See the [reference corpus](../validation/scientific/README.md)
and [R's PCA documentation](https://search.r-project.org/R/refmans/stats/html/prcomp.html).
Tied eigenspaces are deliberately canonicalized instead of treating an arbitrary
rotation within the leading space as a uniquely defined score.

`--module-top` limits reporting only: membership, connectivity and PCA use the full
component. Modules sort by full size, then first input gene. Connectivity sums
positive Spearman coefficients to every other component member. Report ordering
rounds connectivity to an absolute `1e-10` grid, then uses input gene order, so
roundoff in mathematical ties cannot shuffle reporting. Files retain the original
unrounded connectivity values.

| File | Contents |
| --- | --- |
| `module_membership.tsv` | All component members, full size, connectivity and report rank |
| `modules_top_genes.tsv` | At most `--module-top` members per module; size still means full size |
| `module_eigengenes.tsv` | Full-component PCA scores with biological sample columns |
| `module_t_values.tsv` | Gene/module Spearman coefficients and finite t-derived ranking scores |
| `module_fgsea.tsv` | Simple-permutation enrichment and BH adjustment per module |
| `module_fgsea_manifest.json` | Policies, diagnostics, exclusions, resource estimates and hashes |

Constant genes have undefined module correlations and are omitted from the ranking
universe and pathway overlap. Other undefined correlations fail the attempt.
Ranking converts Spearman coefficients using `samples - 2` degrees of freedom,
capping absolute correlation at `0.9999999` first, including perfect relationships.
These are ranking scores, not separately reported correlation-test p-values.

GMT requires unique nonempty pathway names, a description field (possibly empty)
and unique nonempty gene IDs matching the computational matrix exactly. Unknown
and constant genes are removed from overlap; unmatched entries and filtered
pathways are counted in the manifest. Ambiguous GMT or no remaining pathway within
the size bounds fails before graph fitting.

```sh
rnaa corr --adjust none --geneset topvar:5000 --out edges:topk=50 \
  --module-gmt pathways.gmt --module-min-size 25 --module-top 100 \
  --fgsea-permutations 1000 --fgsea-seed 42 --fgsea-min-size 10 --fgsea-max-size 500
```

Minimum module size must be at least two, report size positive and permutations
at least 100. Pathway bounds must satisfy
`1 <= min_size <= max_size < nonconstant_gene_count`; omitting the maximum uses
universe size minus one. Values are validated as supplied rather than silently
clamped. The simple-permutation mode is independently checked against R fgsea
1.38.0 with fixed seeds and matching universes.

The graph stores at most 20 neighbours per gene. Neighbours/connectivity use
prepared Spearman rows and bounded tiles; PCA uses a samples-by-samples Gram
matrix. Neither phase allocates a gene-by-gene matrix. Memory/output estimates
cover input copies, graph storage, PCA, ranking/permutation scratch and worst-case
tables. `corr.block_rows`, `max_memory_mb`, `max_output_mb` and resolved thread
allocation apply to module work. CLI admission runs before ordinary correlation
and before module execution. An owned local pool includes the linked enrichment
library's parallel work.

These phase-specific payload estimates do not cap OS RSS, reserve shared disk/CPU,
bound preceding selection/adjustment memory or provide prompt cancellation.
Global reservations and supported-host production validation remain roadmap work.
Initial [process measurements](../validation/correlation-scale/README.md) include
graph, PCA, ranking and enrichment on deterministic synthetic input.

Module algorithm revision 2 is included in scientific requests. Previous module
requests recompute while their original generations remain available. Allocation
changes preserve identity; reporting/pathway/score parameters change it. The CLI
seals and registers both new tables with the generation, and failed admission
preserves the previous published result.
