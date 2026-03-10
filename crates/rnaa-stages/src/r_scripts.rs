use std::fs;
use std::path::PathBuf;

use anyhow::{Context, Result};
use rnaa_core::paths::ProjectPaths;

const QUANT_KALLISTO_R: &str = include_str!("../../../r/quant_kallisto.R");
const TXIMPORT_DESEQ2_R: &str = include_str!("../../../r/tximport_deseq2.R");

pub fn resolve_r_script(paths: &ProjectPaths, script_name: &str) -> Result<PathBuf> {
    if let Ok(dir) = std::env::var("RNAA_R_SCRIPTS_DIR") {
        let candidate = PathBuf::from(dir).join(script_name);
        if candidate.exists() {
            return Ok(candidate);
        }
    }
    materialize_embedded_r_script(paths, script_name)
}

fn materialize_embedded_r_script(paths: &ProjectPaths, script_name: &str) -> Result<PathBuf> {
    let payload = embedded_r_script(script_name).with_context(|| {
        format!(
            "no embedded payload available for required R script {}",
            script_name
        )
    })?;
    let target = paths.support_r_script_path(script_name);
    fs::create_dir_all(
        target
            .parent()
            .ok_or_else(|| anyhow::anyhow!("invalid support path {}", target.display()))?,
    )
    .with_context(|| format!("failed to create support dir for {}", target.display()))?;
    let should_write = match fs::read_to_string(&target) {
        Ok(existing) => existing != payload,
        Err(_) => true,
    };
    if should_write {
        fs::write(&target, payload)
            .with_context(|| format!("failed to write embedded R script {}", target.display()))?;
    }
    Ok(target)
}

fn embedded_r_script(script_name: &str) -> Result<&'static str> {
    match script_name {
        "quant_kallisto.R" => Ok(QUANT_KALLISTO_R),
        "tximport_deseq2.R" => Ok(TXIMPORT_DESEQ2_R),
        _ => anyhow::bail!("unknown embedded R script: {script_name}"),
    }
}

#[cfg(test)]
mod tests {
    use super::materialize_embedded_r_script;
    use rnaa_core::paths::ProjectPaths;

    #[test]
    fn materializes_embedded_quant_script_when_no_external_copy_exists() {
        let temp = tempfile::tempdir().expect("failed to create tempdir");
        let paths = ProjectPaths::new(temp.path());
        paths
            .ensure_layout()
            .expect("failed to create project layout");
        let resolved = materialize_embedded_r_script(&paths, "quant_kallisto.R")
            .expect("failed to materialize embedded script");
        assert_eq!(resolved, paths.support_r_script_path("quant_kallisto.R"));
        assert!(resolved.exists(), "embedded script should be materialized");
    }
}
