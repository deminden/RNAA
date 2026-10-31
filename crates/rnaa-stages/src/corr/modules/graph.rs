//! Mutual positive Spearman neighbours and full-component connectivity in tiles.
use super::super::{
    prepared::Prepared,
    stream::{Edge, offer},
};
use serde::Serialize;
use std::collections::BinaryHeap;

pub(super) struct Module {
    pub members: Vec<usize>,
    pub ranked: Vec<(usize, f64)>,
}
#[derive(Default, Serialize)]
pub(super) struct GraphStatistics {
    pub peak_tile_entries: usize,
    pub retained_neighbors: usize,
    pub components: usize,
}

pub(super) fn detect(
    prepared: &Prepared,
    min_size: usize,
    block_rows: usize,
) -> (Vec<Module>, GraphStatistics) {
    let n = prepared.constant.len();
    let k = 20.min(n - 1);
    let mut neighbors = vec![Vec::<Edge>::new(); n];
    let mut statistics = GraphStatistics::default();
    for start in (0..n).step_by(block_rows) {
        let end = (start + block_rows).min(n);
        let mut heaps = (start..end).map(|_| BinaryHeap::new()).collect::<Vec<_>>();
        for target in (0..n).step_by(block_rows) {
            let tile = prepared.block(start..end, target..(target + block_rows).min(n));
            statistics.peak_tile_entries = statistics.peak_tile_entries.max(tile.len());
            for ((row, column), &value) in tile.indexed_iter() {
                let other = target + column;
                if start + row != other && value.is_finite() && value > 0.0 {
                    offer(
                        &mut heaps[row],
                        Edge {
                            target: other,
                            value,
                        },
                        k,
                    );
                }
            }
        }
        for (index, heap) in heaps.into_iter().enumerate() {
            neighbors[start + index] = heap.into_vec();
            neighbors[start + index].sort();
        }
    }
    statistics.retained_neighbors = neighbors.iter().map(Vec::len).sum();
    let mut adjacency = vec![Vec::new(); n];
    for (i, list) in neighbors.iter().enumerate() {
        for edge in list {
            let j = edge.target;
            if j > i && neighbors[j].iter().any(|candidate| candidate.target == i) {
                adjacency[i].push(j);
                adjacency[j].push(i);
            }
        }
    }
    let mut visited = vec![false; n];
    let mut modules = Vec::new();
    for i in 0..n {
        if visited[i] {
            continue;
        }
        visited[i] = true;
        let mut stack = vec![i];
        let mut members = Vec::new();
        while let Some(node) = stack.pop() {
            members.push(node);
            for &other in &adjacency[node] {
                if !visited[other] {
                    visited[other] = true;
                    stack.push(other);
                }
            }
        }
        if members.len() >= min_size {
            members.sort_unstable();
            modules.push(Module {
                members,
                ranked: vec![],
            });
        }
    }
    modules.sort_by(|a, b| {
        b.members
            .len()
            .cmp(&a.members.len())
            .then(a.members[0].cmp(&b.members[0]))
    });
    let mut membership = vec![None; n];
    for (index, module) in modules.iter().enumerate() {
        for &gene in &module.members {
            membership[gene] = Some(index);
        }
    }
    let mut connectivity = vec![0.0; n];
    // Accumulate targets in input order, independent of tile size/thread count.
    for start in (0..n).step_by(block_rows) {
        let end = (start + block_rows).min(n);
        for target in (0..n).step_by(block_rows) {
            let tile = prepared.block(start..end, target..(target + block_rows).min(n));
            for ((row, column), &value) in tile.indexed_iter() {
                let (i, j) = (start + row, target + column);
                if i != j
                    && membership[i].is_some()
                    && membership[i] == membership[j]
                    && value.is_finite()
                    && value > 0.0
                {
                    connectivity[i] += value;
                }
            }
        }
    }
    for module in &mut modules {
        module.ranked = module
            .members
            .iter()
            .map(|&gene| (gene, connectivity[gene]))
            .collect();
        // Connectivity is a sum: changing which self entry is omitted can move
        // mathematical ties by a few ulps. Use an explicit reporting-order grid;
        // retain the original connectivity values in all output tables.
        module.ranked.sort_by(|a, b| {
            (b.1 * 1e10)
                .round_ties_even()
                .total_cmp(&(a.1 * 1e10).round_ties_even())
                .then(a.0.cmp(&b.0))
        });
    }
    statistics.components = modules.len();
    (modules, statistics)
}
