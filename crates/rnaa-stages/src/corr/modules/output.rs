//! Module tables expose full membership and scores independently of report limits.
use super::{graph::Module, pca::ModuleScore};
use anyhow::Result;
use csv::{Writer, WriterBuilder};
use rnaa_formats::matrix::MatrixData;
use rsfgsea::prelude::EnrichmentResult;
use std::{
    fs::File,
    path::{Path, PathBuf},
};

pub(super) struct Tables {
    top: Writer<File>,
    membership: Writer<File>,
    scores: Writer<File>,
    ranking: Writer<File>,
    enrichment: Writer<File>,
    pub paths: Vec<PathBuf>,
}
impl Tables {
    pub fn open(directory: &Path, samples: &[String]) -> Result<Self> {
        let names = [
            "modules_top_genes.tsv",
            "module_membership.tsv",
            "module_eigengenes.tsv",
            "module_t_values.tsv",
            "module_fgsea.tsv",
        ];
        let paths = names
            .iter()
            .map(|name| directory.join(name))
            .collect::<Vec<_>>();
        let mut top = WriterBuilder::new().delimiter(b'\t').from_path(&paths[0])?;
        let mut membership = WriterBuilder::new().delimiter(b'\t').from_path(&paths[1])?;
        let mut scores = WriterBuilder::new().delimiter(b'\t').from_path(&paths[2])?;
        let mut ranking = WriterBuilder::new().delimiter(b'\t').from_path(&paths[3])?;
        let mut enrichment = WriterBuilder::new().delimiter(b'\t').from_path(&paths[4])?;
        top.write_record([
            "module_id",
            "module_size",
            "gene",
            "rank_in_module",
            "connectivity",
        ])?;
        membership.write_record([
            "module_id",
            "module_size",
            "gene",
            "rank_in_module",
            "connectivity",
        ])?;
        scores
            .write_record(std::iter::once("module_id").chain(samples.iter().map(String::as_str)))?;
        ranking.write_record(["module_id", "gene", "spearman_r", "t_value"])?;
        enrichment.write_record([
            "module_id",
            "module_size",
            "pathway",
            "pathway_size",
            "es",
            "nes",
            "p_value",
            "padj",
            "log2err",
            "leading_edge_size",
        ])?;
        Ok(Self {
            top,
            membership,
            scores,
            ranking,
            enrichment,
            paths,
        })
    }
    pub fn module(
        &mut self,
        id: &str,
        module: &Module,
        matrix: &MatrixData,
        score: &ModuleScore,
        report_top: usize,
    ) -> Result<()> {
        for (rank, &(gene, connectivity)) in module.ranked.iter().enumerate() {
            let record = [
                id.to_string(),
                module.members.len().to_string(),
                matrix.row_ids[gene].clone(),
                (rank + 1).to_string(),
                connectivity.to_string(),
            ];
            self.membership.write_record(&record)?;
            if rank < report_top {
                self.top.write_record(&record)?;
            }
        }
        self.scores.write_field(id)?;
        for value in &score.values {
            self.scores.write_field(value.to_string())?;
        }
        self.scores.write_record(std::iter::empty::<&str>())?;
        Ok(())
    }
    pub fn rank(&mut self, id: &str, gene: &str, r: f64, score: f64) -> Result<()> {
        self.ranking
            .write_record([id, gene, &r.to_string(), &score.to_string()])?;
        Ok(())
    }
    pub fn enrichment(
        &mut self,
        id: &str,
        size: usize,
        results: &[EnrichmentResult],
    ) -> Result<()> {
        for item in results {
            self.enrichment.write_record([
                id,
                &size.to_string(),
                item.pathway_name.as_str(),
                &item.size.to_string(),
                &item.es.to_string(),
                &item.nes.unwrap_or(f64::NAN).to_string(),
                &item.p_value.to_string(),
                &item.padj.unwrap_or(f64::NAN).to_string(),
                &item.log2err.unwrap_or(f64::NAN).to_string(),
                &item.leading_edge.len().to_string(),
            ])?;
        }
        Ok(())
    }
    pub fn finish(mut self) -> Result<Vec<PathBuf>> {
        self.top.flush()?;
        self.membership.flush()?;
        self.scores.flush()?;
        self.ranking.flush()?;
        self.enrichment.flush()?;
        Ok(self.paths)
    }
}
