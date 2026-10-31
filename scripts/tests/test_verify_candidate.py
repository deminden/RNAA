import hashlib
from contextlib import ExitStack
import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest import mock


spec = importlib.util.spec_from_file_location(
    "verify_candidate", Path(__file__).resolve().parents[1] / "verify_candidate.py"
)
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)
package_spec = importlib.util.spec_from_file_location(
    "package_candidate", Path(__file__).resolve().parents[1] / "package_candidate.py"
)
packager = importlib.util.module_from_spec(package_spec)
package_spec.loader.exec_module(packager)


class CandidateVerificationTests(unittest.TestCase):
    def test_build_caches_are_scoped_to_complete_snapshot_and_mode(self):
        sources = {"script": (b"original", 0o644)}
        key = packager.snapshot_key(sources)
        self.assertNotEqual(key, packager.snapshot_key({"script": (b"changed!", 0o644)}))
        self.assertNotEqual(key, packager.snapshot_key({"script": (b"original", 0o755)}))
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            packager.materialize_source(root, sources)
            hashes = {"script": hashlib.sha256(b"original").hexdigest()}
            self.assertEqual(key, verifier.source_snapshot_key(root, hashes))

    def test_install_target_is_fresh_and_preserves_existing_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            parent = Path(directory)
            sentinel = parent / "previous-artifact"
            sentinel.write_bytes(b"keep")
            with ExitStack() as stack:
                first = verifier.fresh_install_target(parent, stack)
                second = verifier.fresh_install_target(parent, stack)
                self.assertNotEqual(first, second)
                self.assertFalse(list(first.iterdir()))
                self.assertTrue(sentinel.exists())
            self.assertFalse(first.exists())
            self.assertEqual(sentinel.read_bytes(), b"keep")

    def test_candidate_build_uses_and_checks_the_frozen_source(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "source"
            root.mkdir()
            target = Path(directory) / "target"
            sources = {"Cargo.toml": (b"frozen", 0o644)}

            def fake_build(command, cwd, check, env):
                self.assertEqual(cwd, root)
                self.assertTrue(check)
                self.assertIn("--all-features", command)
                self.assertIn("--locked", command)
                self.assertIn("--remap-path-prefix=", env["CARGO_ENCODED_RUSTFLAGS"])
                self.assertEqual((root / "Cargo.toml").read_bytes(), b"frozen")
                binary = target / "test-host/release/rnaa"
                binary.parent.mkdir(parents=True)
                binary.write_bytes(b"built")

            with mock.patch.object(packager.subprocess, "run", side_effect=fake_build):
                binary = packager.build_from_snapshot(root, sources, "test-host", target)
            self.assertEqual(binary.read_bytes(), b"built")
            (root / "Cargo.toml").write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "Frozen source changed"):
                packager.verify_frozen_source(root, sources)
            (root / "Cargo.toml").write_bytes(b"frozen")
            (root / "unexpected").write_bytes(b"new")
            with self.assertRaisesRegex(ValueError, "Frozen source inventory changed"):
                packager.verify_frozen_source(root, sources)

    def test_source_snapshot_omits_private_validation_logs_but_rejects_other_paths(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            contents = {
                "Cargo.lock": b"lock", "LICENSE": b"license",
                "vendor/fasterp/Cargo.toml": b"vendor",
                "validation/checks/run.txt": b"built in /" + b"Users/example/project/",
            }
            for name, data in contents.items():
                path = root / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(data)
            names = ("\0".join(contents) + "\0").encode()
            with mock.patch.object(packager, "ROOT", root), mock.patch.object(
                packager.subprocess, "check_output", return_value=names
            ):
                omitted = []
                sources = packager.source_snapshot(omitted)
                self.assertEqual(omitted, ["validation/checks/run.txt"])
                self.assertNotIn(omitted[0], sources)
                (root / "LICENSE").write_bytes(b"/" + b"home/example/private/")
                with self.assertRaisesRegex(ValueError, "machine-specific path"):
                    packager.source_snapshot()

    def test_every_bundled_source_and_notice_is_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            files = {
                "README.md": b"guide", "docs/other.md": b"more",
                "vendor/fasterp/RNAA_PATCHES.md": b"patches",
                "examples/other.py": b"example", "validation/scientific/reference": b"reference",
                "bin/rnaa": b"binary",
                "RELEASE-MANIFEST.json": b"{}", "DEPENDENCIES.json": b"{}",
                "licenses/lib/LICENSE": b"license",
            }
            for name, data in files.items():
                target = root / name
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(data)
            expected = {name: hashlib.sha256(files[name]).hexdigest() for name in
                        ("README.md", "docs/other.md", "examples/other.py", "validation/scientific/reference",
                         "vendor/fasterp/RNAA_PATCHES.md")}
            dependencies = {"components": [{"license_files": ["licenses/lib/LICENSE"]}]}
            verifier.verify_binary_inventory(root, expected, dependencies)
            (root / "docs/other.md").write_bytes(b"modified")
            with self.assertRaisesRegex(ValueError, "Bundled source file hash mismatch"):
                verifier.verify_binary_inventory(root, expected, dependencies)
            (root / "docs/other.md").write_bytes(files["docs/other.md"])
            (root / "licenses/lib/LICENSE").unlink()
            with self.assertRaisesRegex(ValueError, "Bundled file inventory"):
                verifier.verify_binary_inventory(root, expected, dependencies)

    def test_extract_rejects_duplicates_and_unexpected_roots(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for names, error in [(["expected/a", "expected/./a"], "Duplicate"),
                                 (["other/a"], "Unexpected archive root")]:
                archive = root / "candidate.tar"
                with tarfile.open(archive, "w") as tar:
                    for name in names:
                        info = tarfile.TarInfo(name)
                        info.size = 1
                        tar.addfile(info, io.BytesIO(b"x"))
                with self.assertRaisesRegex(ValueError, error):
                    verifier.extract(archive, root / "out", "expected")


if __name__ == "__main__":
    unittest.main()
