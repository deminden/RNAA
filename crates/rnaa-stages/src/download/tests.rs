use super::*;
use rnaa_core::LibraryLayout;
use rnaa_core::config::{DownloadMethod, ProjectConfig};

#[test]
fn publication_copy_enforces_its_byte_promise_and_keeps_prior_outputs_on_change_or_cancellation()
-> Result<()> {
    let root = tempfile::tempdir()?;
    let source = root.path().join("cache");
    let target = root.path().join("published");
    fs::write(&target, b"prior")?;
    for changed in [b"abc".as_slice(), b"abcdef".as_slice()] {
        fs::write(&source, changed)?;
        let result = publish_file(&target, |temporary| {
            let result = resources::copy_bounded(&source, temporary, 5, &|| Ok(()));
            assert!(fs::metadata(temporary)?.len() <= 5);
            result
        });
        assert!(result.is_err());
        assert_eq!(fs::read(&target)?, b"prior");
    }
    let data = vec![b'x'; 9 * 1024 * 1024];
    fs::write(&source, &data)?;
    let calls = AtomicUsize::new(0);
    assert!(
        publish_file(&target, |temporary| {
            resources::copy_bounded(&source, temporary, data.len() as u64, &|| {
                if calls.fetch_add(1, Ordering::SeqCst) > 0 {
                    bail!("cancelled copy");
                }
                Ok(())
            })
        })
        .is_err()
    );
    assert_eq!(fs::read(&target)?, b"prior");
    assert_eq!(fs::read_dir(root.path())?.count(), 2);
    Ok(())
}

#[test]
fn remote_selection_prefers_fastq() {
    let run = RunRecord {
        project_id: "p".to_string(),
        study_accession: "SRP1".to_string(),
        source_study_accession: None,
        run_accession: "SRR1".to_string(),
        sample_accession: None,
        experiment_accession: None,
        library_layout: LibraryLayout::Paired,
        instrument_platform: None,
        instrument_model: None,
        remote_files: vec![
            RemoteFile {
                url: "https://example.com/a_1.fastq.gz".to_string(),
                md5: None,
                bytes: None,
                kind: RemoteFileKind::FastqR1,
                basename: "a_1.fastq.gz".to_string(),
            },
            RemoteFile {
                url: "https://example.com/a_2.fastq.gz".to_string(),
                md5: None,
                bytes: None,
                kind: RemoteFileKind::FastqR2,
                basename: "a_2.fastq.gz".to_string(),
            },
            RemoteFile {
                url: "https://example.com/a.sra".to_string(),
                md5: None,
                bytes: None,
                kind: RemoteFileKind::Sra,
                basename: "a.sra".to_string(),
            },
        ],
        metadata: BTreeMap::new(),
        state: RunState::Resolved,
        last_error: None,
        updated_at: now_rfc3339(),
    };
    let selected = ReqwestDownloader.select_remote_files(&run, DownloadPreference::Fastq);
    assert_eq!(selected.len(), 2);
    assert_eq!(selected[0].kind, RemoteFileKind::FastqR1);
    let expected = ReqwestDownloader
        .source_request_digest(&run, DownloadPreference::Fastq)
        .unwrap();
    for state in [
        RunState::Resolved,
        RunState::Downloading,
        RunState::Verified,
        RunState::QuantDone,
        RunState::DeRunning,
        RunState::DeDone,
        RunState::CorrDone,
        RunState::Cleaned,
    ] {
        let mut observed = run.clone();
        observed.state = state;
        let files = ReqwestDownloader.select_remote_files(&observed, DownloadPreference::Fastq);
        assert_eq!(files.len(), 2);
        assert!(files.iter().all(|file| file.kind != RemoteFileKind::Sra));
        assert_eq!(
            ReqwestDownloader
                .source_request_digest(&observed, DownloadPreference::Fastq)
                .unwrap(),
            expected
        );
        assert!(ReqwestDownloader.fastq_recovery_available(&observed, DownloadPreference::Fastq));
        assert!(!ReqwestDownloader.fastq_recovery_available(&observed, DownloadPreference::Sra));
    }
}

#[test]
fn config_selects_native_reqwest_downloads() {
    let config: ProjectConfig = toml::from_str(
        r#"
        [project]
        name = "test"
        schema_version = 1

        [download]
        method = "reqwest"
        "#,
    )
    .expect("config should parse");
    assert_eq!(config.download.method, DownloadMethod::Reqwest);
}
