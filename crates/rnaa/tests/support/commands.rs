use std::{
    path::Path,
    process::{Command, Output},
};

pub fn command(root: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_rnaa"))
        .arg("--root")
        .arg(root)
        .args(args)
        .output()
        .unwrap()
}
pub fn success(root: &Path, args: &[&str]) {
    let result = command(root, args);
    assert!(
        result.status.success(),
        "{args:?}: {}\n{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
}
