//! Advisory current-version checks for a historical snapshot without a checksum.
use super::{Checkpoint, Client, DownloadConfig, Result};
use anyhow::Context;
use reqwest::{
    StatusCode,
    header::{ACCEPT_ENCODING, CONTENT_ENCODING, CONTENT_LENGTH, ETAG, IF_NONE_MATCH},
};
use serde_json::{Value, json};
use std::time::Duration;

fn strong_etag(tag: &str) -> bool {
    let bytes = tag.as_bytes();
    bytes.len() >= 2
        && bytes.first() == Some(&b'"')
        && bytes.last() == Some(&b'"')
        && bytes[1..bytes.len() - 1]
            .iter()
            .all(|byte| *byte == 0x21 || (0x23..=0x7e).contains(byte))
}

pub(crate) fn validate_saved_version(
    url: &str,
    sha256: &str,
    bytes: u64,
    expected_bytes: Option<u64>,
    receipt: &Value,
    config: &DownloadConfig,
    check: &impl Fn() -> Result<()>,
) -> Result<Option<Value>> {
    check()?;
    let Ok(saved) = serde_json::from_value::<Checkpoint>(receipt.clone()) else {
        return Ok(None);
    };
    let Some(etag) = saved.etag.as_deref().filter(|tag| strong_etag(tag)) else {
        return Ok(None);
    };
    if saved.schema_version != 1
        || serde_json::to_value(&saved)? != *receipt
        || saved.url != url
        || saved.prefix_sha256 != sha256
        || saved.committed_bytes != bytes
        || saved.expected_md5.is_some()
        || saved.expected_bytes != expected_bytes
        || saved.total_bytes.is_some_and(|total| total != bytes)
    {
        return Ok(None);
    }
    config.validate()?;
    let timeout = Duration::from_secs(config.read_timeout_seconds.min(config.max_transfer_seconds));
    let _ = std::time::Instant::now()
        .checked_add(timeout)
        .context("HTTP origin validation timeout is too large")?;
    let client = Client::builder()
        .timeout(timeout)
        .connect_timeout(Duration::from_secs(config.connect_timeout_seconds).min(timeout))
        .no_gzip()
        .no_brotli()
        .no_deflate()
        .no_zstd()
        .build()?;
    let response = client
        .head(url)
        .header(ACCEPT_ENCODING, "identity")
        .header(IF_NONE_MATCH, etag)
        .send();
    check()?;
    let Ok(response) = response else {
        return Ok(None);
    };
    if !matches!(response.status(), StatusCode::OK | StatusCode::NOT_MODIFIED)
        || response.url().as_str() != saved.effective_url
        || response
            .headers()
            .get(ETAG)
            .and_then(|value| value.to_str().ok())
            != Some(etag)
        || response
            .headers()
            .get(CONTENT_ENCODING)
            .is_some_and(|value| value.as_bytes() != b"identity")
    {
        return Ok(None);
    }
    if let Some(length) = response.headers().get(CONTENT_LENGTH)
        && length
            .to_str()
            .ok()
            .and_then(|value| value.parse::<u64>().ok())
            != Some(bytes)
    {
        return Ok(None);
    }
    Ok(Some(
        json!({"schema_version":1,"kind":"http_strong_etag","method":"HEAD","url":url,
        "effective_url":response.url().as_str(),"etag":etag,"bytes":bytes,"sha256":sha256,
        "status":response.status().as_u16(),"checked_at":rnaa_core::util::now_rfc3339()}),
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        io::{BufRead, BufReader, Write},
        net::TcpListener,
        thread,
    };

    fn serve(headers: String) -> (String, thread::JoinHandle<()>) {
        let listener = TcpListener::bind("127.0.0.1:0").unwrap();
        let url = format!("http://{}/reads.fastq", listener.local_addr().unwrap());
        let handle = thread::spawn(move || {
            let (mut stream, _) = listener.accept().unwrap();
            stream
                .set_read_timeout(Some(Duration::from_secs(5)))
                .unwrap();
            let mut reader = BufReader::new(stream.try_clone().unwrap());
            let mut line = String::new();
            reader.read_line(&mut line).unwrap();
            assert!(line.starts_with("HEAD /reads.fastq "), "{line}");
            let mut conditional = false;
            loop {
                line.clear();
                assert!(reader.read_line(&mut line).unwrap() > 0);
                if line == "\r\n" {
                    break;
                }
                conditional |= line
                    .to_ascii_lowercase()
                    .starts_with("if-none-match: \"v1\"");
            }
            assert!(conditional);
            stream.write_all(headers.as_bytes()).unwrap();
        });
        (url, handle)
    }
    fn saved(url: &str) -> Value {
        json!({"schema_version":1,"url":url,"expected_md5":null,"expected_bytes":4,"effective_url":url,
            "etag":"\"v1\"","last_modified":null,"total_bytes":4,"committed_bytes":4,"prefix_sha256":"a".repeat(64)})
    }
    #[test]
    fn strong_origin_checks_require_matching_representation_not_just_success_or_304() -> Result<()>
    {
        for (status, etag, length, encoding, accepted) in [
            ("200 OK", Some("\"v1\""), Some("4"), None, true),
            ("304 Not Modified", Some("\"v1\""), None, None, true),
            ("304 Not Modified", None, None, None, false),
            ("304 Not Modified", Some("W/\"v1\""), None, None, false),
            ("200 OK", Some("\"v2\""), Some("4"), None, false),
            ("200 OK", Some("\"v1\""), Some("5"), None, false),
            ("200 OK", Some("\"v1\""), Some("4"), Some("gzip"), false),
            (
                "405 Method Not Allowed",
                Some("\"v1\""),
                Some("4"),
                None,
                false,
            ),
        ] {
            let mut headers = format!("HTTP/1.1 {status}\r\nConnection: close\r\n");
            if let Some(etag) = etag {
                headers.push_str(&format!("ETag: {etag}\r\n"));
            }
            if let Some(length) = length {
                headers.push_str(&format!("Content-Length: {length}\r\n"));
            }
            if let Some(encoding) = encoding {
                headers.push_str(&format!("Content-Encoding: {encoding}\r\n"));
            }
            headers.push_str("\r\n");
            let (url, handle) = serve(headers);
            let result = validate_saved_version(
                &url,
                &"a".repeat(64),
                4,
                Some(4),
                &saved(&url),
                &DownloadConfig::default(),
                &|| Ok(()),
            )?;
            handle.join().unwrap();
            assert_eq!(result.is_some(), accepted, "status={status} etag={etag:?}");
        }
        Ok(())
    }
    #[test]
    fn invalid_saved_validators_and_mismatched_witnesses_cannot_prove_current_content() -> Result<()>
    {
        for tag in ["W/\"v1\"", "v1", "\"a\",\"b\"", "\"space tag\"", "\"a\n\""] {
            assert!(!strong_etag(tag));
        }
        assert!(strong_etag("\"v1\""));
        assert!(strong_etag("\"a,b\""));
        let url = "http://127.0.0.1:1/reads.fastq";
        for (field, value) in [
            ("prefix_sha256", json!("b".repeat(64))),
            ("committed_bytes", json!(3)),
            ("expected_bytes", json!(5)),
            ("unknown", json!(true)),
        ] {
            let mut receipt = saved(url);
            receipt[field] = value;
            assert!(
                validate_saved_version(
                    url,
                    &"a".repeat(64),
                    4,
                    Some(4),
                    &receipt,
                    &DownloadConfig::default(),
                    &|| Ok(())
                )?
                .is_none()
            );
        }
        let mut receipt = saved(url);
        receipt["etag"] = json!("W/\"v1\"");
        assert!(
            validate_saved_version(
                url,
                &"a".repeat(64),
                4,
                Some(4),
                &receipt,
                &DownloadConfig::default(),
                &|| Ok(())
            )?
            .is_none()
        );
        Ok(())
    }
    #[test]
    fn cancellation_after_origin_validation_propagates_instead_of_becoming_a_cache_miss() {
        let (url, handle) =
            serve("HTTP/1.1 304 Not Modified\r\nETag: \"v1\"\r\nConnection: close\r\n\r\n".into());
        let calls = std::sync::atomic::AtomicUsize::new(0);
        let result = validate_saved_version(
            &url,
            &"a".repeat(64),
            4,
            Some(4),
            &saved(&url),
            &DownloadConfig::default(),
            &|| {
                if calls.fetch_add(1, std::sync::atomic::Ordering::SeqCst) > 0 {
                    anyhow::bail!("cancelled");
                }
                Ok(())
            },
        );
        handle.join().unwrap();
        assert!(result.unwrap_err().to_string().contains("cancelled"));
    }
}
