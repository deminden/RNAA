import hashlib
import importlib.util
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

spec = importlib.util.spec_from_file_location(
    "assemble_release", Path(__file__).resolve().parents[1] / "assemble_release.py"
)
release = importlib.util.module_from_spec(spec)
spec.loader.exec_module(release)


class AssembleReleaseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.inputs = self.root / "inputs"
        self.inputs.mkdir()
        for os_name, target in release.PLATFORMS.items():
            candidate = self.inputs / f"rnaa-candidate-{os_name}-rust-stable"
            candidate.mkdir()
            manifest = {"target": target, "git_head": "abc", "dirty": False,
                        "executable_version": "rnaa 1.2.3-rc.1", "source_snapshot_sha256": "same"}
            files = {"RELEASE-MANIFEST.json": json.dumps(manifest).encode(),
                     "DEPENDENCIES.json": b"{}", "source.tar.gz": b"source",
                     f"{target}.tar.gz": b"binary"}
            for name, data in files.items():
                (candidate / name).write_bytes(data)
            (candidate / "SHA256SUMS").write_text("".join(
                f"{hashlib.sha256(data).hexdigest()}  {name}\n" for name, data in files.items()
            ))

    def run_assembly(self, tag="v1.2.3-rc.1", commit="abc"):
        release.assemble(self.inputs, self.root / "out", tag, commit, "1.2.3-rc.1")

    def test_builds_distinct_self_contained_assets(self):
        self.run_assembly()
        output = self.root / "out"
        assets = sorted(output.glob("*.tar.gz"))
        self.assertEqual({asset.name for asset in assets},
                         {f"rnaa-v1.2.3-rc.1-{target}.tar.gz" for target in release.PLATFORMS.values()})
        for asset in assets:
            with tarfile.open(asset) as archive:
                names = archive.getnames()
                self.assertTrue(any(name.endswith("/SHA256SUMS") for name in names))
                self.assertTrue(any(name.endswith("/source.tar.gz") for name in names))
        self.assertEqual(len((output / "SHA256SUMS").read_text().splitlines()), 3)

    def test_accepts_version_tag_without_v_prefix(self):
        self.run_assembly(tag="1.2.3-rc.1")
        self.assertEqual(len(list((self.root / "out").glob("*.tar.gz"))), 3)

    def test_assembles_plain_stable_version(self):
        for candidate in self.inputs.iterdir():
            path = candidate / "RELEASE-MANIFEST.json"
            manifest = json.loads(path.read_text())
            manifest["executable_version"] = "rnaa 0.3.0"
            path.write_text(json.dumps(manifest))
            checksums = candidate / "SHA256SUMS"
            lines = checksums.read_text().splitlines()
            lines = [f"{release.digest(path)}  RELEASE-MANIFEST.json" if line.endswith("  RELEASE-MANIFEST.json") else line for line in lines]
            checksums.write_text("\n".join(lines) + "\n")
        release.assemble(self.inputs, self.root / "out", "0.3.0", "abc", "0.3.0")
        self.assertTrue((self.root / "out" / "rnaa-0.3.0-aarch64-apple-darwin.tar.gz").exists())

    def test_rejects_invalid_tag_or_commit(self):
        with self.assertRaisesRegex(ValueError, "tag must match"):
            self.run_assembly(tag="v9.9.9")
        with self.assertRaisesRegex(ValueError, "identity"):
            self.run_assembly(commit="wrong")

    def test_rejects_changed_or_missing_candidate(self):
        candidate = self.inputs / "rnaa-candidate-macos-15-rust-stable"
        (candidate / "source.tar.gz").write_bytes(b"tampered")
        with self.assertRaisesRegex(ValueError, "checksum mismatch"):
            self.run_assembly()
        (candidate / "source.tar.gz").unlink()
        with self.assertRaisesRegex(ValueError, "inventory"):
            self.run_assembly()

    def test_rejects_differing_snapshot_even_if_checksums_are_valid(self):
        candidate = self.inputs / "rnaa-candidate-macos-15-rust-stable"
        path = candidate / "RELEASE-MANIFEST.json"
        manifest = json.loads(path.read_text())
        manifest["source_snapshot_sha256"] = "different"
        path.write_text(json.dumps(manifest))
        lines = (candidate / "SHA256SUMS").read_text().splitlines()
        lines = [f"{release.digest(path)}  RELEASE-MANIFEST.json" if line.endswith("  RELEASE-MANIFEST.json") else line for line in lines]
        (candidate / "SHA256SUMS").write_text("\n".join(lines) + "\n")
        with self.assertRaisesRegex(ValueError, "same source snapshot"):
            self.run_assembly()


if __name__ == "__main__":
    unittest.main()
