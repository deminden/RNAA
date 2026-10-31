//! Bounded layout inspection of the current counting loader's kallisto v13 input.
mod bitmap;
mod reader;
#[cfg(test)]
mod tests;
use anyhow::{Context, Result, bail};
use reader::{Reader, add};
use std::path::Path;

#[derive(Debug, Default, Clone, serde::Serialize)]
pub struct QuantIndexPlan {
    pub model: &'static str,
    pub file_bytes: u64,
    pub index_version: u64,
    pub k: u32,
    pub g: u32,
    pub graph_header: u64,
    pub metadata_header: u64,
    pub unitigs: u64,
    pub kmer_nodes: u64,
    pub abundant_nodes: u64,
    pub decoded_graph_bases: u64,
    pub max_unitig_bases: u64,
    pub minimizer_buckets: u64,
    pub graph_bitmap_values: u64,
    pub special_minimizers: u64,
    pub special_bitmap_values: u64,
    pub mphf_words: u64,
    pub mphf_levels: u64,
    pub mphf_elements: u64,
    pub mphf_serialized_bytes: u64,
    pub mphf_final_entries: u64,
    pub dlist_entries: u64,
    pub nodes: u64,
    pub ec_blocks: u64,
    pub ec_values_upper: u64,
    pub max_ec_cardinality: u64,
    pub non_single_subset_classes_upper: u64,
    pub non_single_subset_members_upper: u64,
    pub max_ec_target: Option<u32>,
    pub max_node_payload_bytes: u64,
    pub max_bitmap_payload_bytes: u64,
    pub transcripts: u64,
    pub transcript_name_bytes: u64,
    pub max_transcript_name_bytes: u64,
    pub has_shade_names: bool,
    pub onlist_values: u64,
    pub resident_estimate_bytes: u64,
    pub load_peak_estimate_bytes: u64,
}
impl QuantIndexPlan {
    /// Include loader worker stack allowances. This checks one job's capacity;
    /// it does not reserve resources against concurrent producers.
    pub fn load_memory_bytes(&self, threads: usize) -> Result<u64> {
        (threads.max(1) as u64)
            .checked_mul(2 * 1024 * 1024)
            .and_then(|stacks| self.load_peak_estimate_bytes.checked_add(stacks))
            .context("index loader memory allowance overflows")
    }
    pub fn check_memory_limit(&self, capacity: u64, threads: usize) -> Result<()> {
        let requested = self.load_memory_bytes(threads)?;
        if requested > capacity {
            return Err(rnaa_core::resources::ResourceAdmissionError::new(anyhow::anyhow!(
                "quant index loading estimate ({requested} bytes) exceeds declared memory capacity ({capacity} bytes); reduce loader threads or increase resources.memory_mb"
            )).into());
        }
        Ok(())
    }
    pub fn inspect(path: &Path, check: &impl Fn() -> Result<()>) -> Result<Self> {
        let mut r = Reader::open(path, check)?;
        let mut s = Self {
            model: "kallistors-counting-v1-layout",
            file_bytes: r.end,
            ..Self::default()
        };
        s.index_version = r.u64()?;
        if s.index_version != 13 {
            bail!(
                "index preflight supports kallisto v13, found {}",
                s.index_version
            );
        }
        let graph = r.u64()?;
        r.section(graph, |r| s.graph(r))?;
        let mphf = r.u64()?;
        s.mphf_serialized_bytes = mphf;
        r.section(mphf, |r| s.mphf(r))?;
        s.dlist_entries = r.u64()?;
        r.u64()?; // D-list overhang
        r.skip(r.fits(s.dlist_entries, 8)?)?;
        s.nodes = r.u64()?;
        if s.nodes
            != s.unitigs
                .checked_add(s.kmer_nodes)
                .and_then(|n| n.checked_add(s.abundant_nodes))
                .context("index node count overflows")?
        {
            bail!("index graph and node counts differ");
        }
        r.fits(s.nodes, u64::from(s.k) + 4)?;
        for _ in 0..s.nodes {
            r.skip(u64::from(s.k))?;
            let size = u64::from(r.u32()?);
            s.max_node_payload_bytes = s.max_node_payload_bytes.max(size);
            r.section(size, |r| s.node(r))?;
        }
        let transcripts = r.u32()?;
        if transcripts == 0 || transcripts > i32::MAX as u32 {
            bail!("invalid index transcript count");
        }
        s.transcripts = u64::from(transcripts);
        r.fits(s.transcripts, 4)?;
        for _ in 0..s.transcripts {
            if r.u32()? > i32::MAX as u32 {
                bail!("negative index transcript length");
            }
        }
        r.fits(s.transcripts, 8)?;
        for _ in 0..s.transcripts {
            let bytes = r.u64()?;
            s.has_shade_names |= r.contains(bytes, b"_shade_")?;
            add(&mut s.transcript_name_bytes, bytes)?;
            s.max_transcript_name_bytes = s.max_transcript_name_bytes.max(bytes);
        }
        let onlist = r.u64()?;
        s.max_bitmap_payload_bytes = s.max_bitmap_payload_bytes.max(onlist);
        if onlist > 0 {
            s.onlist_values = bitmap::compatible(&mut r, onlist)?.values;
        }
        if s.max_ec_target
            .is_some_and(|id| u64::from(id) >= s.transcripts)
        {
            bail!("index equivalence-class target exceeds its transcript table");
        }
        s.estimate()?;
        check()?;
        Ok(s)
    }
    fn graph(&mut self, r: &mut Reader<'_>) -> Result<()> {
        self.graph_header = r.u64()?;
        self.k = r.u32()?;
        self.g = r.u32()?;
        if self.k == 0 || self.k > 32 || self.g == 0 || self.g > self.k {
            bail!("invalid index k/g");
        }
        self.unitigs = r.u64()?;
        r.fits(self.unitigs, 8)?;
        for _ in 0..self.unitigs {
            let bases = r.u64()?;
            if bases < u64::from(self.k) {
                bail!("index unitig is shorter than k");
            }
            r.skip(bases.div_ceil(4))?;
            add(&mut self.decoded_graph_bases, bases)?;
            self.max_unitig_bases = self.max_unitig_bases.max(bases);
        }
        self.kmer_nodes = r.u64()?;
        r.skip(r.fits(self.kmer_nodes, 8)?)?;
        self.abundant_nodes = r.u64()?;
        r.skip(r.fits(self.abundant_nodes, 8)?)?;
        let singleton = self
            .kmer_nodes
            .checked_add(self.abundant_nodes)
            .and_then(|n| n.checked_mul(u64::from(self.k)))
            .context("index decoded graph size overflows")?;
        add(&mut self.decoded_graph_bases, singleton)?;
        self.metadata_header = r.u64()?;
        r.u64()?; // checksum
        if r.u64()? != self.unitigs
            || r.u64()? != self.kmer_nodes
            || r.u64()? != self.abundant_nodes
        {
            bail!("index graph metadata counts differ");
        }
        self.minimizer_buckets = r.u64()?;
        for _ in 0..2 {
            let bitmaps = r.u64()?;
            r.fits(bitmaps, 8)?;
            for _ in 0..bitmaps {
                let values = self.graph_bitmap(r)?;
                add(&mut self.graph_bitmap_values, values)?;
            }
        }
        self.special_minimizers = r.u64()?;
        let pages = (self.special_minimizers >> 32) + 1;
        r.fits(pages, 16)?;
        for _ in 0..pages {
            let values = self.graph_bitmap(r)?;
            add(&mut self.special_bitmap_values, values)?;
        }
        for _ in 0..pages {
            let values = self.graph_bitmap(r)?;
            add(&mut self.special_bitmap_values, values)?;
        }
        // Every special entry owns one (minimizer, position) pair after loading.
        // Its serialized minimizer alone takes 8 bytes. Optional 0/4/8-byte
        // marker fields do not change allocation size. The declared graph end
        // locates the next section without materializing membership sets.
        r.fits(self.special_minimizers, 8)?;
        r.skip(r.end - r.pos)?;
        Ok(())
    }
    fn graph_bitmap(&mut self, r: &mut Reader<'_>) -> Result<u64> {
        let start = r.pos;
        let values = bitmap::bitcontainer(r)?.values;
        self.max_bitmap_payload_bytes = self.max_bitmap_payload_bytes.max(r.pos - start);
        Ok(values)
    }
    fn mphf(&mut self, r: &mut Reader<'_>) -> Result<()> {
        r.u64()?; // gamma f64 bits
        let levels = r.u32()?;
        if levels > i32::MAX as u32 {
            bail!("negative index MPHF level count");
        }
        self.mphf_levels = u64::from(levels);
        r.u64()?; // last bitset rank
        self.mphf_elements = r.u64()?;
        r.fits(self.mphf_levels, 24)?;
        for _ in 0..levels {
            r.u64()?; // bit domain
            let words = r.u64()?;
            r.skip(r.fits(words, 8)?)?;
            add(&mut self.mphf_words, words)?;
            let ranks = r.u64()?;
            r.skip(r.fits(ranks, 8)?)?;
            add(&mut self.mphf_words, ranks)?;
        }
        self.mphf_final_entries = r.u64()?;
        r.skip(r.fits(self.mphf_final_entries, 16)?)?;
        Ok(())
    }
    fn node(&mut self, r: &mut Reader<'_>) -> Result<()> {
        r.u32()?; // node identity
        let blocks = match r.u8()? {
            0 => 0,
            1 => 1,
            2 => r.u64()?,
            _ => bail!("invalid index EC block flag"),
        };
        r.fits(blocks, 24)?;
        add(&mut self.ec_blocks, blocks)?;
        for _ in 0..blocks {
            let lower = r.u32()?;
            let upper = r.u32()?;
            if upper < lower {
                bail!("index EC block bounds are reversed");
            }
            let bytes = r.u64()?;
            let ec = bitmap::native(r, bytes)?;
            add(&mut self.ec_values_upper, ec.values)?;
            self.max_ec_cardinality = self.max_ec_cardinality.max(ec.values);
            let variants = if ec.values >= 64 {
                u64::MAX
            } else {
                (1_u64 << ec.values).saturating_sub(ec.values + 1)
            };
            let members = if ec.values >= 64 {
                u64::MAX
            } else if ec.values > 1 {
                ec.values.saturating_mul((1_u64 << (ec.values - 1)) - 1)
            } else {
                0
            };
            self.non_single_subset_classes_upper = self
                .non_single_subset_classes_upper
                .saturating_add(variants);
            self.non_single_subset_members_upper =
                self.non_single_subset_members_upper.saturating_add(members);
            if let Some(id) = ec.maximum {
                self.max_ec_target = Some(self.max_ec_target.map_or(id, |old| old.max(id)));
            }
            let positions = r.u64()?;
            r.fits(positions, 8)?;
            for _ in 0..positions {
                let size = r.u64()?;
                r.skip(size)?;
            }
        }
        Ok(())
    }
    fn estimate(&mut self) -> Result<()> {
        // Explicit per-structure allowances for the current ASCII counting
        // loader, including Vec growth and hash table slack. Not an RSS ceiling.
        let dna = u128::from(self.decoded_graph_bases);
        let nodes = u128::from(self.nodes);
        let positions = u128::from(self.graph_bitmap_values) + u128::from(self.special_minimizers);
        let resident = 32 * 1024 * 1024_u128
            + dna * 2
            + nodes * 128
            + positions * 32
            + u128::from(self.minimizer_buckets) * 32
            + u128::from(self.ec_blocks) * 64
            + u128::from(self.ec_values_upper) * 8
            + u128::from(self.mphf_words) * 48
            + (u128::from(self.mphf_final_entries) + u128::from(self.dlist_entries)) * 64
            + u128::from(self.transcript_name_bytes) * 6
            + u128::from(self.transcripts) * 128;
        let temporary = 64 * 1024 * 1024_u128
            + u128::from(self.mphf_serialized_bytes) * 2
            + u128::from(self.max_bitmap_payload_bytes) * 2
            + (u128::from(self.graph_bitmap_values)
                + u128::from(self.onlist_values))
                * 16
            // Special markers are inserted into growing hash sets while their
            // decoded vectors still live. Include old/new tables during growth.
            + u128::from(self.special_bitmap_values)*48
            // EC arena growth and the optional shade-name hash table may also
            // overlap their old allocation before it is released.
            + u128::from(self.ec_values_upper)*4
            + u128::from(self.transcripts)*128
            + u128::from(self.max_node_payload_bytes) * 2
            + u128::from(self.max_unitig_bases) * 2
            + u128::from(self.special_minimizers) * 80
            + u128::from(self.max_transcript_name_bytes) * 6;
        self.resident_estimate_bytes =
            u64::try_from(resident).context("index resident memory estimate overflows")?;
        self.load_peak_estimate_bytes =
            u64::try_from(resident + temporary).context("index load memory estimate overflows")?;
        Ok(())
    }
}
