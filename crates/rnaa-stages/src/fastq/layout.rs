use anyhow::{Context, Result, bail};
use rnaa_core::{
    LibraryLayout,
    model::{ArtifactKind, VerifiedFile},
};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

#[derive(Debug)]
pub(super) struct Lane {
    pub first: PathBuf,
    pub second: Option<PathBuf>,
}
#[derive(Debug)]
pub(super) struct Layout {
    pub lanes: Vec<Lane>,
    pub first_kind: ArtifactKind,
    pub staging: bool,
}
pub(super) fn inspect(layout: LibraryLayout, files: &[VerifiedFile]) -> Result<Layout> {
    let (mut singles, mut first, mut second) = (Vec::new(), Vec::new(), Vec::new());
    let mut paths = BTreeSet::new();
    let mut trimmed = BTreeSet::new();
    let mut staging = false;
    for file in files {
        let path = fs::canonicalize(&file.path)
            .with_context(|| format!("missing FASTQ {}", file.path.display()))?;
        if !fs::metadata(&path)?.is_file() {
            bail!(
                "FASTQ inputs must be regular files: {}",
                file.path.display()
            );
        }
        if !paths.insert(path) {
            bail!("duplicate FASTQ input {}", file.path.display());
        }
        match file.artifact_kind {
            ArtifactKind::FastqSingle => {
                singles.push(file);
                trimmed.insert(false);
            }
            ArtifactKind::FastqTrimmedSingle => {
                singles.push(file);
                trimmed.insert(true);
            }
            ArtifactKind::FastqR1 => {
                first.push(file);
                trimmed.insert(false);
            }
            ArtifactKind::FastqTrimmedR1 => {
                first.push(file);
                trimmed.insert(true);
            }
            ArtifactKind::FastqR2 => {
                second.push(file);
                trimmed.insert(false);
            }
            ArtifactKind::FastqTrimmedR2 => {
                second.push(file);
                trimmed.insert(true);
            }
            other => bail!("unsupported FASTQ artifact {}", other.as_str()),
        }
        let mut reader = BufReader::new(File::open(&file.path)?);
        let gzip = reader.fill_buf()?.starts_with(&[0x1f, 0x8b]);
        staging |= gzip
            != (file
                .path
                .extension()
                .and_then(|extension| extension.to_str())
                == Some("gz"));
    }
    if trimmed.len() > 1 {
        bail!("cannot mix original and preprocessed FASTQs");
    }
    let pairs = match layout {
        LibraryLayout::Single if !singles.is_empty() && first.is_empty() && second.is_empty() => {
            singles.sort_by_key(|file| &file.path);
            singles
                .into_iter()
                .map(|file| (file, None))
                .collect::<Vec<_>>()
        }
        LibraryLayout::Paired
            if singles.is_empty() && !first.is_empty() && first.len() == second.len() =>
        {
            if first.len() == 1 {
                vec![(first[0], Some(second[0]))]
            } else {
                let mut mates = BTreeMap::new();
                for file in second {
                    if mates.insert(lane_key(&file.path, 2)?, file).is_some() {
                        bail!("ambiguous R2 lane names");
                    }
                }
                first.sort_by_key(|file| &file.path);
                let mut paired = Vec::new();
                for file in first {
                    let key = lane_key(&file.path, 1)?;
                    let mate = mates
                        .remove(&key)
                        .with_context(|| format!("missing R2 for lane {key}"))?;
                    paired.push((file, Some(mate)));
                }
                paired
            }
        }
        _ => bail!(
            "FASTQ files do not match declared {} layout: require single-end files or complete R1/R2 lane pairs",
            layout.as_str()
        ),
    };
    staging |= pairs.len() > 1;
    Ok(Layout {
        first_kind: pairs[0].0.artifact_kind,
        lanes: pairs
            .into_iter()
            .map(|(first, second)| Lane {
                first: first.path.clone(),
                second: second.map(|file| file.path.clone()),
            })
            .collect(),
        staging,
    })
}
fn lane_key(path: &Path, mate: u8) -> Result<String> {
    let name = path
        .file_name()
        .and_then(|name| name.to_str())
        .context("FASTQ basename must be UTF-8")?;
    let expression = regex::Regex::new(&format!(r"(?i)(^|[._-])R?{mate}([._-]|$)"))?;
    if expression.find_iter(name).count() != 1 {
        bail!("ambiguous mate marker in multi-lane FASTQ {name}");
    }
    Ok(expression.replace(name, "${1}MATE${2}").into_owned())
}
