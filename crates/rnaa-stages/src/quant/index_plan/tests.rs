use super::*;
use std::{
    fs,
    sync::atomic::{AtomicUsize, Ordering},
};

fn index(root: &Path, targets: usize) -> Result<std::path::PathBuf> {
    let fasta = root.join("transcripts.fa");
    let mut random = 91_u64;
    let sequence = (0..1200)
        .map(|_| {
            random = random.wrapping_mul(6364136223846793005).wrapping_add(1);
            b"ACGT"[((random >> 32) & 3) as usize] as char
        })
        .collect::<String>();
    let mut text = String::new();
    for target in 0..targets {
        text.push_str(&format!(">tx{target}\n{sequence}\n"));
    }
    fs::write(&fasta, text)?;
    let path = root.join("reference.index");
    kallistors::index::build_index(
        &path,
        &[fasta],
        kallistors::index::IndexBuildOptions {
            k: 21,
            ..Default::default()
        },
    )
    .map_err(|e| anyhow::anyhow!("{e}"))?;
    Ok(path)
}

#[test]
fn layouts_match_independent_public_engine_metadata_and_loaded_graphs() -> Result<()> {
    for targets in [1, 2, 8, 120] {
        let root = tempfile::tempdir()?;
        let path = index(root.path(), targets)?;
        let plan = QuantIndexPlan::inspect(&path, &|| Ok(()))?;
        let metadata = kallistors::index::Index::load(&path).map_err(|e| anyhow::anyhow!("{e}"))?;
        let loaded = kallistors::pseudoalign::build_bifrost_index_quant_threaded(&path, false, 2)
            .map_err(|e| anyhow::anyhow!("{e}"))?;
        assert_eq!(plan.k, metadata.k);
        assert_eq!(plan.transcripts, metadata.transcripts.len() as u64);
        assert_eq!(
            plan.transcript_name_bytes,
            metadata
                .transcripts
                .iter()
                .map(|t| t.name.len() as u64)
                .sum::<u64>()
        );
        assert_eq!(plan.unitigs, loaded.unitigs.len() as u64);
        assert_eq!(plan.kmer_nodes, loaded.km_unitigs.len() as u64);
        assert_eq!(plan.abundant_nodes, loaded.h_kmers.len() as u64);
        assert_eq!(
            plan.decoded_graph_bases,
            loaded
                .unitigs
                .iter()
                .chain(&loaded.km_unitigs)
                .chain(&loaded.h_kmers)
                .map(|s| s.len() as u64)
                .sum::<u64>()
        );
        assert_eq!(plan.max_ec_cardinality, targets as u64);
        assert_eq!(plan.max_ec_target, Some(targets as u32 - 1));
        assert_eq!(plan.onlist_values, targets as u64);
        assert!(plan.load_peak_estimate_bytes > plan.resident_estimate_bytes);
    }
    Ok(())
}

#[test]
fn impossible_counts_and_truncation_fail_without_engine_allocations() -> Result<()> {
    let root = tempfile::tempdir()?;
    let path = index(root.path(), 2)?;
    let original = fs::read(&path)?;
    let graph_end = 16 + u64::from_le_bytes(original[8..16].try_into()?) as usize;
    for at in [0, 8, 32, 40, graph_end, graph_end + 8 + 36] {
        let mut bytes = original.clone();
        bytes[at..at + 8].copy_from_slice(&u64::MAX.to_le_bytes());
        fs::write(&path, bytes)?;
        assert!(
            QuantIndexPlan::inspect(&path, &|| Ok(())).is_err(),
            "field {at}"
        );
    }
    for end in [0, 7, 16, original.len() / 2, original.len() - 1] {
        fs::write(&path, &original[..end])?;
        assert!(
            QuantIndexPlan::inspect(&path, &|| Ok(())).is_err(),
            "length {end}"
        );
    }
    fs::write(&path, &original)?;
    let calls = AtomicUsize::new(0);
    let result = QuantIndexPlan::inspect(&path, &|| {
        if calls.fetch_add(1, Ordering::SeqCst) > 0 {
            anyhow::bail!("cancelled index preflight");
        }
        Ok(())
    });
    assert!(format!("{:#}", result.err().unwrap()).contains("cancelled index preflight"));
    Ok(())
}

#[test]
fn compressed_full_u32_onlist_is_budgeted_without_materializing_its_values() -> Result<()> {
    let root = tempfile::tempdir()?;
    let path = index(root.path(), 1)?;
    let mut original = fs::read(&path)?;
    // One full run in each 16-bit key container represents 2^32 values in less
    // than 1 MiB. Never hand this bitmap to a decoder that expands it into Vec.
    let containers = 65536_u32;
    let mut encoded = (12347_u32 | ((containers - 1) << 16))
        .to_le_bytes()
        .to_vec();
    encoded.extend_from_slice(&[255_u8; 8192]);
    for key in 0..containers {
        encoded.extend_from_slice(&(key as u16).to_le_bytes());
        encoded.extend_from_slice(&65535_u16.to_le_bytes());
    }
    let first = 4 + 8192 + containers * 8;
    for key in 0..containers {
        encoded.extend_from_slice(&(first + key * 6).to_le_bytes());
    }
    for _ in 0..containers {
        encoded.extend_from_slice(&1_u16.to_le_bytes());
        encoded.extend_from_slice(&0_u16.to_le_bytes());
        encoded.extend_from_slice(&65535_u16.to_le_bytes());
    }
    assert!(encoded.len() < 1024 * 1024);
    // The independent writer emits an 18-byte portable singleton onlist.
    let at = original.len() - 26;
    assert_eq!(u64::from_le_bytes(original[at..at + 8].try_into()?), 18);
    original.truncate(at);
    original.extend_from_slice(&(encoded.len() as u64).to_le_bytes());
    original.extend_from_slice(&encoded);
    fs::write(&path, original)?;
    let plan = QuantIndexPlan::inspect(&path, &|| Ok(()))?;
    assert_eq!(plan.onlist_values, 1_u64 << 32);
    assert!(plan.file_bytes < 2 * 1024 * 1024);
    assert!(plan.load_memory_bytes(1)? > 64 * 1024 * 1024 * 1024);
    assert!(plan.check_memory_limit(512 * 1024 * 1024, 1).is_err());
    Ok(())
}

fn bitmap(bytes: &[u8], native: bool) -> Result<bitmap::Bitmap> {
    let root = tempfile::tempdir()?;
    let path = root.path().join("bitmap");
    fs::write(&path, bytes)?;
    let mut reader = Reader::open(&path, &|| Ok(()))?;
    if native {
        bitmap::native(&mut reader, bytes.len() as u64)
    } else {
        bitmap::compatible(&mut reader, bytes.len() as u64)
    }
}

#[test]
fn bitmap_arrays_runs_and_bitsets_are_counted_and_malformed_cards_are_rejected() -> Result<()> {
    let mut array = vec![1];
    array.extend_from_slice(&3_u32.to_le_bytes());
    for value in [9_u32, 2, 9] {
        array.extend_from_slice(&value.to_le_bytes());
    }
    let plan = bitmap(&array, true)?;
    assert_eq!(plan.values, 3); // conservative bound before engine deduplication
    assert_eq!(plan.maximum, Some(9));
    let mut runs = 12347_u32.to_le_bytes().to_vec();
    runs.push(1); // one run container, no offsets below four containers
    runs.extend_from_slice(&0_u16.to_le_bytes());
    runs.extend_from_slice(&65535_u16.to_le_bytes());
    runs.extend_from_slice(&1_u16.to_le_bytes());
    runs.extend_from_slice(&0_u16.to_le_bytes());
    runs.extend_from_slice(&65535_u16.to_le_bytes());
    assert_eq!(bitmap(&runs, false)?.values, 65536);
    let mut bitset = 12346_u32.to_le_bytes().to_vec();
    bitset.extend_from_slice(&1_u32.to_le_bytes());
    bitset.extend_from_slice(&0_u16.to_le_bytes());
    bitset.extend_from_slice(&8191_u16.to_le_bytes());
    bitset.extend_from_slice(&16_u32.to_le_bytes());
    for _ in 0..128 {
        bitset.extend_from_slice(&u64::MAX.to_le_bytes());
    }
    bitset.resize(16 + 8192, 0);
    assert_eq!(bitmap(&bitset, false)?.values, 8192);
    bitset[12..16].copy_from_slice(&17_u32.to_le_bytes());
    assert!(bitmap(&bitset, false).is_err());
    bitset[12..16].copy_from_slice(&16_u32.to_le_bytes());
    bitset[10..12].copy_from_slice(&8190_u16.to_le_bytes());
    assert!(bitmap(&bitset, false).is_err());
    runs.pop();
    assert!(bitmap(&runs, false).is_err());
    Ok(())
}

#[test]
fn local_bitcontainers_and_special_entries_keep_counting_loader_compatibility() -> Result<()> {
    let root = tempfile::tempdir()?;
    for (tag, payload, card) in [
        (1_u64 | (1 << 3) | (1 << 63), Vec::new(), 2_u64),
        ((123_u64 << 3) | 2, Vec::new(), 1),
        (0, vec![((6_u16 << 3) | 2), 3, 0, 7, 11, 20], 3),
        (0, vec![((7_u16 << 3) | 4), 5, 0, 5, 9, 0, 0], 5),
        (0, vec![(4_u16 << 3), 3, 0, 0b101010], 3),
    ] {
        let mut bytes = tag.to_le_bytes().to_vec();
        for value in payload {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        let path = root.path().join("bitmap");
        fs::write(&path, bytes)?;
        let mut reader = Reader::open(&path, &|| Ok(()))?;
        let actual = bitmap::bitcontainer(&mut reader)?;
        assert_eq!(actual.values, card);
        let expected =
            kallistors::index::bifrost::read_bitcontainer_values(&mut std::fs::File::open(&path)?)
                .map_err(|e| anyhow::anyhow!("{e}"))?;
        assert_eq!(expected.len() as u64, card);
    }
    let mut malformed = 0_u64.to_le_bytes().to_vec();
    for value in [(5_u16 << 3) | 4, 5, 0, 5, 9] {
        malformed.extend_from_slice(&value.to_le_bytes());
    }
    let invalid = root.path().join("bad-tiny");
    fs::write(&invalid, malformed)?;
    let mut reader = Reader::open(&invalid, &|| Ok(()))?;
    assert!(bitmap::bitcontainer(&mut reader).is_err());
    let path = index(root.path(), 2)?;
    let mut bytes = fs::read(&path)?;
    let graph = u64::from_le_bytes(bytes[8..16].try_into()?);
    let end = 16 + graph as usize;
    // The current writer ends graph metadata with the special count and two
    // empty 16-byte pointer bitmaps. Add one normal special minimizer entry.
    assert_eq!(&bytes[end - 40..end - 32], &0_u64.to_le_bytes());
    bytes[end - 40..end - 32].copy_from_slice(&1_u64.to_le_bytes());
    bytes.splice(end..end, [0_u8; 16]);
    bytes[8..16].copy_from_slice(&(graph + 16).to_le_bytes());
    fs::write(&path, bytes)?;
    let plan = QuantIndexPlan::inspect(&path, &|| Ok(()))?;
    assert_eq!(plan.special_minimizers, 1);
    assert!(kallistors::pseudoalign::build_bifrost_index_quant_threaded(&path, false, 1).is_ok());
    Ok(())
}
