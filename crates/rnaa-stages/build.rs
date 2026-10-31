use sha2::{Digest, Sha256};
use std::{env, fs, path::Path};

fn main() {
    let manifest = env::var("CARGO_MANIFEST_DIR").expect("Cargo supplies the manifest directory");
    let lock = Path::new(&manifest).join("../../Cargo.lock");
    println!("cargo:rerun-if-changed={}", lock.display());
    let bytes = fs::read(&lock).expect("the workspace lockfile must be available");
    let digest = hex::encode(Sha256::digest(&bytes));
    println!("cargo:rustc-env=RNAA_DEPENDENCY_LOCK_SHA256={digest}");
    let document: toml::Value =
        toml::from_str(std::str::from_utf8(&bytes).expect("Cargo.lock is UTF-8"))
            .expect("Cargo.lock is valid TOML");
    let packages = document["package"].as_array().expect("lockfile packages");
    for (package, variable) in [
        ("kallistors", "RNAA_KALLISTORS_VERSION"),
        ("rsdeseq2", "RNAA_RSDESEQ2_VERSION"),
        ("rsfgsea", "RNAA_RSFGSEA_VERSION"),
        ("mincorr", "RNAA_MINCORR_VERSION"),
    ] {
        let matches: Vec<_> = packages
            .iter()
            .filter(|entry| entry["name"].as_str() == Some(package))
            .collect();
        assert_eq!(matches.len(), 1, "expected one locked {package} engine");
        let entry = matches[0];
        let version = entry["version"].as_str().expect("engine version");
        println!("cargo:rustc-env={variable}={version}");
        if package == "mincorr" {
            let source = entry["source"].as_str().expect("mincorr Git source");
            let (_, revision) = source.rsplit_once('#').expect("locked Git revision");
            println!("cargo:rustc-env=RNAA_MINCORR_REVISION={revision}");
        }
    }
}
