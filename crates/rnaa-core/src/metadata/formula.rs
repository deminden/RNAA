//! Discover metadata references; rsdeseq2 remains the formula syntax validator.
use super::{ColumnType, MetadataConfig};
use anyhow::{Result, bail};
use std::collections::BTreeSet;

pub fn variables(
    formula: &str,
    headers: &[String],
    schema: &MetadataConfig,
) -> Result<BTreeSet<String>> {
    let (lhs, rhs) = formula
        .split_once('~')
        .ok_or_else(|| anyhow::anyhow!("design formula must start with '~'"))?;
    if !lhs.trim().is_empty() || rhs.trim().is_empty() {
        bail!("design must have a nonempty right-hand side and no response variable");
    }
    let bytes = rhs.as_bytes();
    let mut i = 0;
    let mut result = BTreeSet::new();
    while i < bytes.len() {
        match bytes[i] {
            b'\'' | b'"' => {
                let quote = bytes[i];
                i += 1;
                while i < bytes.len() && bytes[i] != quote {
                    if bytes[i] == b'\\' {
                        i += 1;
                    }
                    i += 1;
                }
                if i >= bytes.len() {
                    bail!("unterminated string in design formula");
                }
                i += 1;
            }
            b'%' if rhs[i..].starts_with("%in%") => i += 4,
            b'0'..=b'9' => i = number_end(bytes, i),
            b'A'..=b'Z' | b'a'..=b'z' | b'_' => {
                let begin = i;
                i += 1;
                while i < bytes.len()
                    && (bytes[i].is_ascii_alphanumeric() || matches!(bytes[i], b'_' | b'.'))
                {
                    i += 1;
                }
                let name = &rhs[begin..i];
                let mut next = i;
                while next < bytes.len() && bytes[next].is_ascii_whitespace() {
                    next += 1;
                }
                if next < bytes.len() && matches!(bytes[next], b'(' | b'=') {
                    continue;
                }
                if !matches!(name, "TRUE" | "FALSE") {
                    result.insert(name.to_string());
                }
            }
            b'.' if i + 1 < bytes.len() && bytes[i + 1].is_ascii_digit() => {
                i = number_end(bytes, i)
            }
            b'~' => bail!("design formula must contain exactly one unquoted '~'"),
            b'.' => {
                for name in headers {
                    if schema.column(name).is_some_and(|column| {
                        matches!(column.kind, ColumnType::Factor | ColumnType::Numeric)
                    }) {
                        result.insert(name.clone());
                    }
                }
                i += 1;
            }
            byte if !byte.is_ascii() => {
                bail!("design identifiers must use ASCII letters, digits and underscores")
            }
            _ => i += 1,
        }
    }
    let missing = result
        .iter()
        .filter(|name| !headers.contains(name))
        .cloned()
        .collect::<Vec<_>>();
    if !missing.is_empty() {
        bail!(
            "design formula '{formula}' references missing samplesheet columns: {}",
            missing.join(", ")
        );
    }
    for name in &result {
        let column = schema.column(name).ok_or_else(|| anyhow::anyhow!("undeclared covariate '{name}'; declare metadata.columns.{name}.type as factor or numeric"))?;
        if !matches!(column.kind, ColumnType::Factor | ColumnType::Numeric) {
            bail!("identifier/ignored metadata '{name}' cannot be a model covariate");
        }
    }
    Ok(result)
}

fn number_end(bytes: &[u8], mut index: usize) -> usize {
    while index < bytes.len() && (bytes[index].is_ascii_digit() || bytes[index] == b'.') {
        index += 1;
    }
    if index < bytes.len() && matches!(bytes[index], b'e' | b'E') {
        index += 1;
        if index < bytes.len() && matches!(bytes[index], b'+' | b'-') {
            index += 1;
        }
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
    }
    index
}
