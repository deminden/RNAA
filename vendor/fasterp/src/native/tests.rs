use super::*;
use anyhow::Result;
use std::{
    fs,
    io::Write,
    path::Path,
    sync::{mpsc, Arc},
    time::Duration,
};

fn options(threads: usize) -> NativeProcessOptions {
    NativeProcessOptions {
        threads,
        batch_bytes: 512,
        max_in_flight_batches: 2,
        max_record_bytes: 8192,
        input_byte_limits: None,
        min_record_bytes: 0,
        max_input_records: None,
    }
}
fn fastq(reads: usize, length: usize, crlf: bool, mate: usize) -> String {
    let newline = if crlf { "\r\n" } else { "\n" };
    let seq = "ACGT".repeat(length.div_ceil(4));
    let seq = &seq[..length];
    (0..reads)
        .map(|read| {
            format!(
                "@read{read}/{mate}{newline}{seq}{newline}+{newline}{}{newline}",
                "I".repeat(length)
            )
        })
        .collect()
}
fn metrics(result: &NativeProcessResult) -> serde_json::Value {
    let mut value: serde_json::Value = serde_json::from_str(&result.json_report).unwrap();
    value.as_object_mut().unwrap().remove("threads");
    value.as_object_mut().unwrap().remove("execution");
    value
}
fn run(
    root: &Path,
    input: &Path,
    input2: Option<&Path>,
    label: &str,
    options: NativeProcessOptions,
) -> Result<(NativeProcessResult, Vec<u8>, Option<Vec<u8>>)> {
    let out = root.join(format!("{label}-1.fastq"));
    let out2 = root.join(format!("{label}-2.fastq"));
    let result = process_paths_with_options(
        input.to_str().unwrap(),
        out.to_str().unwrap(),
        input2.map(|path| path.to_str().unwrap()),
        input2.map(|_| out2.to_str().unwrap()),
        options,
    )?;
    let report: serde_json::Value = serde_json::from_str(&result.json_report)?;
    assert!(
        report["execution"]["batch_window"]["peak"]
            .as_u64()
            .unwrap()
            <= options.max_in_flight_batches as u64
    );
    assert_eq!(report["execution"]["batch_window"]["active"], 0);
    Ok((
        result,
        fs::read(out)?,
        input2.map(|_| fs::read(out2)).transpose()?,
    ))
}

#[test]
fn bounded_single_and_paired_processing_match_serial_with_unequal_mate_lengths() -> Result<()> {
    let root = tempfile::tempdir()?;
    let first = root.path().join("r1.fastq");
    let second = root.path().join("r2.fastq");
    fs::write(&first, fastq(100, 1200, false, 1))?;
    fs::write(&second, fastq(100, 60, false, 2))?;
    for paired in [false, true] {
        let mate = paired.then_some(second.as_path());
        let baseline = run(root.path(), &first, mate, "serial", options(1))?;
        for threads in [2, 4] {
            let result = run(
                root.path(),
                &first,
                mate,
                &format!("parallel-{threads}"),
                options(threads),
            )?;
            assert_eq!(baseline.1, result.1);
            assert_eq!(baseline.2, result.2);
            assert_eq!(metrics(&baseline.0), metrics(&result.0));
        }
    }
    Ok(())
}

#[test]
fn crlf_and_multimember_gzip_match_normalized_plain_input() -> Result<()> {
    let root = tempfile::tempdir()?;
    let plain = root.path().join("plain.fastq");
    let gzip = root.path().join("members.fastq.gz");
    fs::write(&plain, fastq(20, 100, false, 1))?;
    let mut compressed = Vec::new();
    for start in [0, 10] {
        let text = fastq(10, 100, true, 1).replace("@read", &format!("@part{start}-read"));
        let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(text.as_bytes())?;
        compressed.extend(encoder.finish()?);
    }
    fs::write(&gzip, compressed)?;
    let serial = run(root.path(), &plain, None, "plain", options(1))?;
    let crlf = run(root.path(), &gzip, None, "gzip", options(3))?;
    assert_eq!(metrics(&serial.0), metrics(&crlf.0));
    let crlf_plain = root.path().join("crlf.fastq");
    fs::write(&crlf_plain, fastq(20, 100, true, 1))?;
    assert_eq!(
        serial.1,
        run(root.path(), &crlf_plain, None, "crlf", options(3))?.1
    );
    Ok(())
}

#[test]
fn total_decoded_byte_promises_reject_growth_in_serial_and_parallel_paths() -> Result<()> {
    let root = tempfile::tempdir()?;
    let input = root.path().join("input.fastq");
    let mate = root.path().join("mate.fastq");
    let first = fastq(4, 100, false, 1);
    let second = fastq(4, 60, false, 2);
    fs::write(&input, &first)?;
    fs::write(&mate, &second)?;
    for threads in [1, 3] {
        for paired in [false, true] {
            let mut limits = options(threads);
            limits.input_byte_limits = Some([first.len() as u64, second.len() as u64]);
            let name = format!("bounded-{threads}-{paired}");
            let _ = run(
                root.path(),
                &input,
                paired.then_some(mate.as_path()),
                &name,
                limits,
            )?;
            limits.input_byte_limits.as_mut().unwrap()[usize::from(paired)] -= 1;
            let error = run(
                root.path(),
                &input,
                paired.then_some(mate.as_path()),
                &format!("short-{threads}-{paired}"),
                limits,
            )
            .unwrap_err();
            assert!(format!("{error:#}").contains("preflight decoded byte promise"));
        }
    }
    Ok(())
}

#[test]
fn record_shape_promises_reject_shorter_or_additional_records() -> Result<()> {
    let root = tempfile::tempdir()?;
    let path = root.path().join("reads.fastq");
    let text = fastq(4, 100, false, 1);
    fs::write(&path, &text)?;
    for threads in [1, 3] {
        let mut limits = options(threads);
        limits.max_input_records = Some(4);
        limits.min_record_bytes = 200;
        let _ = run(
            root.path(),
            &path,
            None,
            &format!("shape-ok-{threads}"),
            limits,
        )?;
        limits.max_input_records = Some(3);
        let error = run(
            root.path(),
            &path,
            None,
            &format!("count-fail-{threads}"),
            limits,
        )
        .err()
        .unwrap();
        assert!(format!("{error:#}").contains("preflight record count"));
        limits.max_input_records = Some(4);
        limits.min_record_bytes = 300;
        let error = run(
            root.path(),
            &path,
            None,
            &format!("minimum-fail-{threads}"),
            limits,
        )
        .err()
        .unwrap();
        assert!(format!("{error:#}").contains("preflight minimum"));
    }
    Ok(())
}

#[test]
fn invalid_options_and_incomplete_pair_arguments_fail_before_output_creation() -> Result<()> {
    let root = tempfile::tempdir()?;
    let input = root.path().join("in.fastq");
    let output = root.path().join("out.fastq");
    fs::write(&input, fastq(1, 40, false, 1))?;
    let mut bad = options(2);
    bad.max_in_flight_batches = 0;
    assert!(process_paths_with_options(
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        None,
        None,
        bad
    )
    .is_err());
    assert!(process_paths_with_options(
        input.to_str().unwrap(),
        output.to_str().unwrap(),
        Some(input.to_str().unwrap()),
        None,
        options(2)
    )
    .is_err());
    assert!(!output.exists());
    Ok(())
}

fn result(id: u64, permit: Arc<BatchPermit>) -> crate::pipeline::WorkerResult {
    use crossbeam_channel::{bounded, unbounded};
    let text = format!("@read{id}\n{}\n+\n{}\n", "ACGT".repeat(10), "I".repeat(40));
    let header = text.find('\n').unwrap() + 1;
    let (tx, rx) = bounded(2);
    let (out_tx, out_rx) = unbounded();
    tx.send(Some(crate::pipeline::Batch {
        id,
        buf: Arc::new(text.into_bytes()),
        recs: vec![[0, header, header + 41, header + 43]],
        permit: Some(permit),
    }))
    .unwrap();
    tx.send(None).unwrap();
    drop(tx);
    crate::pipeline::worker_thread(
        rx,
        out_tx,
        15,
        5,
        15,
        40,
        0,
        false,
        30,
        false,
        default_trimming(),
    );
    out_rx.recv().unwrap().unwrap()
}

#[test]
fn an_out_of_order_result_holds_its_credit_until_the_merger_writes_it() -> Result<()> {
    let root = tempfile::tempdir()?;
    let output = root.path().join("ordered.fastq");
    let window = window::BatchWindow::new(2);
    let zero = result(0, window.acquire()?);
    let one = result(1, window.acquire()?);
    let (tx, rx) = crossbeam_channel::bounded(2);
    let merger = std::thread::spawn({
        let output = output.clone();
        move || {
            crate::pipeline::merger_thread(
                rx,
                output.display().to_string(),
                2,
                None,
                crate::split::SplitConfig::default(),
                false,
            )
        }
    });
    tx.send(Some(one))?;
    let (started_tx, started_rx) = mpsc::channel();
    let (admitted_tx, admitted_rx) = mpsc::channel();
    let producer = std::thread::spawn({
        let window = Arc::clone(&window);
        move || {
            started_tx.send(()).unwrap();
            let _permit = window.acquire().unwrap();
            admitted_tx.send(()).unwrap();
        }
    });
    started_rx.recv()?;
    assert!(admitted_rx.recv_timeout(Duration::from_millis(50)).is_err());
    tx.send(Some(zero))?;
    admitted_rx.recv_timeout(Duration::from_secs(2))?;
    tx.send(None)?;
    tx.send(None)?;
    drop(tx);
    producer.join().unwrap();
    merger.join().unwrap()?;
    let text = fs::read_to_string(output)?;
    assert!(text.find("@read0").unwrap() < text.find("@read1").unwrap());
    assert_eq!(window.stats()?.active, 0);
    assert_eq!(window.stats()?.peak, 2);
    Ok(())
}

#[test]
fn native_error_child() -> Result<()> {
    let Some(root) = std::env::var_os("RNAA_NATIVE_ERROR_TEST") else {
        return Ok(());
    };
    let root = Path::new(&root);
    let input = root.join("input.fastq");
    fs::write(&input, fastq(1000, 100, false, 1))?;
    assert!(process_paths_with_options(
        input.to_str().unwrap(),
        root.to_str().unwrap(),
        None,
        None,
        options(3)
    )
    .is_err());
    let bad = root.join("bad.fastq");
    fs::write(&bad, "@broken\nACGT\n+\nI\n")?;
    assert!(process_paths_with_options(
        bad.to_str().unwrap(),
        root.join("bad-output").to_str().unwrap(),
        None,
        None,
        options(3)
    )
    .is_err());
    for threads in [1, 3] {
        let mut opts = options(threads);
        opts.max_record_bytes = 64;
        assert!(process_paths_with_options(
            input.to_str().unwrap(),
            root.join("limited-output").to_str().unwrap(),
            None,
            None,
            opts
        )
        .is_err());
    }
    Ok(())
}

#[test]
fn producer_and_writer_errors_stop_every_scoped_thread() -> Result<()> {
    let root = tempfile::tempdir()?;
    let mut child = TestChild(
        std::process::Command::new(std::env::current_exe()?)
            .args([
                "--exact",
                "native::tests::native_error_child",
                "--nocapture",
            ])
            .env("RNAA_NATIVE_ERROR_TEST", root.path())
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .spawn()?,
    );
    let deadline = std::time::Instant::now() + Duration::from_secs(10);
    loop {
        if let Some(status) = child.0.try_wait()? {
            assert!(status.success());
            break;
        }
        if std::time::Instant::now() > deadline {
            panic!("native pipeline did not stop after failure");
        }
        std::thread::sleep(Duration::from_millis(10));
    }
    Ok(())
}
struct TestChild(std::process::Child);
impl Drop for TestChild {
    fn drop(&mut self) {
        let _ = self.0.kill();
        let _ = self.0.wait();
    }
}
