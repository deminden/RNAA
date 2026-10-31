//! Validate parameters, declared payload budgets and effective GMT overlap first.
use super::super::limits::{CorrelationLimits, inspect_shape};
use super::ModuleGseaParams;
use anyhow::{Result, bail};
use rnaa_formats::matrix::{MatrixData, read_matrix_tsv};
use rsfgsea::prelude::Pathway;
use serde::Serialize;
use std::{
    collections::BTreeSet,
    fs,
    io::{BufRead, BufReader},
    path::Path,
};

#[derive(Serialize)]
pub(super) struct Admission {
    pub working_bytes: u128,
    pub output_bytes: u128,
    pub eligible_genes: usize,
    pub excluded_constant_genes: Vec<String>,
    pub pathway_min_size: usize,
    pub pathway_max_size: usize,
    pub input_pathways: usize,
    pub filtered_pathways: usize,
    pub unmatched_gene_entries: usize,
}

impl ModuleGseaParams {
    pub fn validate(&self) -> Result<()> {
        if self.top_genes == 0 || self.min_module_size < 2 {
            bail!("module-top must be positive and module-min-size must be at least 2");
        }
        if self.permutations < 100 {
            bail!("fgsea-permutations must be at least 100; values are not silently increased");
        }
        if self.pathway_min_size == 0
            || self
                .pathway_max_size
                .is_some_and(|max| max < self.pathway_min_size)
        {
            bail!("fgsea pathway sizes must be positive with max-size >= min-size");
        }
        Ok(())
    }
}

pub(super) fn load(
    path: &Path,
    params: &ModuleGseaParams,
    limits: CorrelationLimits,
    threads: usize,
) -> Result<(MatrixData, Vec<Pathway>, Admission)> {
    params.validate()?;
    if limits.block_rows == 0
        || limits.block_rows > 1024
        || limits.max_memory_mb == 0
        || limits.max_output_mb == 0
    {
        bail!("invalid module resource limits");
    }
    let shape = inspect_shape(path)?;
    if shape.samples < 3 {
        bail!("module enrichment requires at least three samples");
    }
    if params.min_module_size > shape.genes {
        bail!("module-min-size exceeds the number of input genes");
    }
    let gmt_bytes = fs::metadata(&params.gmt_path)?.len() as u128;
    let (n, s, b) = (
        shape.genes as u128,
        shape.samples as u128,
        limits.block_rows.min(shape.genes) as u128,
    );
    let max_modules = n / params.min_module_size as u128;
    let base_working = n * s * 64
        + s * s * 128
        + n * 2048
        + b * b * 16
        + shape.id_bytes * 8
        + shape.header_bytes * 3
        + shape.max_record_bytes * 4
        + gmt_bytes * 12;
    let memory_limit = limits.max_memory_mb as u128 * 1048576;
    if base_working > memory_limit {
        bail!(
            "module working estimate {base_working} bytes exceeds max_memory_mb={}",
            limits.max_memory_mb
        );
    }
    // Parse strictly: duplicate pathway identities or gene entries are ambiguous.
    let mut pathways = Vec::new();
    let mut names = BTreeSet::new();
    let mut max_pathway_name = 0usize;
    for (index, line) in BufReader::new(fs::File::open(&params.gmt_path)?)
        .lines()
        .enumerate()
    {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let parts = line.split('\t').collect::<Vec<_>>();
        if parts.len() < 3
            || parts[0].trim().is_empty()
            || !names.insert(parts[0].trim().to_string())
        {
            bail!("invalid or duplicate GMT pathway on line {}", index + 1);
        }
        let mut unique = BTreeSet::new();
        let mut genes = Vec::new();
        for gene in &parts[2..] {
            let gene = gene.trim();
            if gene.is_empty() || !unique.insert(gene) {
                bail!("empty or duplicate GMT gene on line {}", index + 1);
            }
            genes.push(gene.to_string());
        }
        max_pathway_name = max_pathway_name.max(parts[0].len());
        pathways.push(Pathway {
            name: parts[0].trim().into(),
            description: (!parts[1].is_empty()).then(|| parts[1].into()),
            genes,
        });
    }
    if pathways.is_empty() {
        bail!("module GMT has no pathways");
    }
    let working_bytes = base_working
        + n * threads.max(1) as u128 * 256
        + params.permutations as u128 * pathways.len() as u128 * 32;
    // Reserve all module rankings and the worst-case pathway result fields on disk.
    let output_bytes = max_modules
        * (n * (shape.max_id as u128 * 2 + 384)
            + s * 384
            + pathways.len() as u128 * (max_pathway_name as u128 * 2 + 2048))
        + n * (shape.max_id as u128 * 4 + 768)
        + 16384;
    if working_bytes > memory_limit {
        bail!(
            "module working estimate {working_bytes} bytes exceeds max_memory_mb={}",
            limits.max_memory_mb
        );
    }
    if output_bytes > limits.max_output_mb as u128 * 1048576 {
        bail!(
            "module output estimate {output_bytes} bytes exceeds max_output_mb={}",
            limits.max_output_mb
        );
    }
    let matrix = read_matrix_tsv(path)?;
    if matrix.row_ids.len() != shape.genes || matrix.col_ids.len() != shape.samples {
        bail!("module input changed after resource preflight");
    }
    let mut eligible = BTreeSet::new();
    let mut excluded_constant_genes = Vec::new();
    for (gene, row) in matrix.row_ids.iter().zip(&matrix.values) {
        if row.iter().all(|value| *value == row[0]) {
            excluded_constant_genes.push(gene.clone());
        } else {
            eligible.insert(gene.clone());
        }
    }
    if eligible.len() < params.min_module_size || eligible.len() < 3 {
        bail!("too few nonconstant genes for module enrichment");
    }
    let max_size = params.pathway_max_size.unwrap_or(eligible.len() - 1);
    if max_size >= eligible.len() || params.pathway_min_size > max_size {
        bail!(
            "fgsea pathway bounds must satisfy min-size <= max-size < the nonconstant gene universe ({})",
            eligible.len()
        );
    }
    let input_pathways = pathways.len();
    let mut unmatched_gene_entries = 0;
    for pathway in &mut pathways {
        let before = pathway.genes.len();
        pathway.genes.retain(|gene| eligible.contains(gene));
        unmatched_gene_entries += before - pathway.genes.len();
    }
    pathways.retain(|pathway| (params.pathway_min_size..=max_size).contains(&pathway.genes.len()));
    if pathways.is_empty() {
        bail!("no GMT pathways have eligible gene overlap within the requested size bounds");
    }
    let report = Admission {
        working_bytes,
        output_bytes,
        eligible_genes: eligible.len(),
        excluded_constant_genes,
        pathway_min_size: params.pathway_min_size,
        pathway_max_size: max_size,
        input_pathways,
        filtered_pathways: input_pathways - pathways.len(),
        unmatched_gene_entries,
    };
    Ok((matrix, pathways, report))
}
