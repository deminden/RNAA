"""Deterministic synthetic NB workflow fixture, not a reference-engine parity oracle.
Seed 20261029, 500 genes, 8 independent samples, gamma-Poisson NB counts.
The first 50 genes increase 4x in B; the next 50 decrease 4x.
Python standard library only. Run from the repository root.
"""
import math
import pathlib
import random
rng = random.Random(20261029)
def poisson(mean):
    # Add independent Poisson blocks to avoid underflow for large means.
    total = 0
    while mean > 0:
        block = min(mean, 20.0)
        product, count = 1.0, -1
        limit = math.exp(-block)
        while product > limit:
            product *= rng.random()
            count += 1
        total += count
        mean -= block
    return total
out = pathlib.Path('fixtures/native_workflow_counts.tsv')
with out.open('w') as handle:
    handle.write('gene_id\t' + '\t'.join(f'RUN{i+1}' for i in range(8)) + '\n')
    for gene in range(500):
        base = 20.0 + gene % 200 * 4
        dispersion = 0.03 + 10.0/base
        counts = []
        for sample in range(8):
            effect = (4 if gene < 50 else .25 if gene < 100 else 1) if sample >= 4 else 1
            mean = base * effect
            counts.append(poisson(rng.gammavariate(1/dispersion, mean*dispersion)))
        handle.write(f'G{gene:04}\t' + '\t'.join(map(str,counts)) + '\n')
