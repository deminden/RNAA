//! Shared deterministic input for network and module process-level scale probes.
use anyhow::Result;
use rnaa_core::model::AdjustedMatrix;
use std::{fs, io::Write, path::Path};
pub fn input(root: &Path, genes: usize, samples: usize) -> Result<AdjustedMatrix> {
    let path = root.join("input.tsv");
    let mut writer = std::io::BufWriter::new(fs::File::create(&path)?);
    write!(writer, "gene_id")?;
    for sample in 0..samples {
        write!(writer, "\tsample{sample}")?;
    }
    writeln!(writer)?;
    for gene in 0..genes {
        write!(writer, "gene{gene:07}")?;
        for sample in 0..samples {
            let value = ((gene as u64 + 1)
                .wrapping_mul(6364136223846793005)
                .wrapping_add((sample as u64).wrapping_mul(1442695040888963407)))
            .rotate_left((sample % 63) as u32)
                % 100000;
            write!(writer, "\t{value}")?;
        }
        writeln!(writer)?;
    }
    writer.flush()?;
    drop(writer);
    Ok(AdjustedMatrix {
        matrix_path: path.clone(),
        genes: vec![],
        samples: vec![],
    })
}
