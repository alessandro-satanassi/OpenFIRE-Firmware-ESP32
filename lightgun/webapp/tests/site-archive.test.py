"""Offline test of the WebApp publication step, not of an obsolete build layout.

    python tests/site-archive.test.py        (from lightgun/webapp)

build_site() and build_launcher() produce the publication inputs. The actual
Python step in GLOBAL-build-package-release.yml creates v/<id> and versions.json.
Only disposable directories are changed. No Git, network, PyYAML or hardware.
"""
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest.mock import patch

LIGHTGUN = Path(__file__).resolve().parents[2]
WORKFLOW = LIGHTGUN.parent / ".github/workflows/GLOBAL-build-package-release.yml"
sys.path.insert(0, str(LIGHTGUN / "scripts"))

import webapp_build  # noqa: E402


def publication_script():
    """Extract the real Python heredoc; do not duplicate the publication logic."""
    text = WORKFLOW.read_text(encoding="utf-8-sig")
    blocks = re.findall(
        r"^      - name: Update versioned site, launcher and versions\.json\n"
        r"        run: \|\n"
        r"          python3 - <<'PY'\n(.*?)^          PY[ \t]*$",
        text, re.MULTILINE | re.DOTALL)
    if len(blocks) != 1:
        raise RuntimeError("Expected one WebApp publication step in " + str(WORKFLOW))
    return textwrap.dedent(blocks[0])


def files(folder):
    return {path.relative_to(folder).as_posix(): path.read_bytes()
            for path in folder.rglob("*") if path.is_file()}


class SiteArchiveTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.script = publication_script()
        compile(cls.script, str(WORKFLOW), "exec")

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="of-site-archive-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "webapp_dist"
        self.repo = self.root / "webapp_repo"
        self.repo.mkdir()
        (self.repo / "README.md").write_text("Preserve this documentation.\n", encoding="utf-8")
        (self.repo / "CNAME").write_text("example.invalid\n", encoding="utf-8")
        self.prepare()

    def prepare(self, version_id="6.2", label="6.2.0", kind="stable"):
        self.version = {"id": version_id, "label": label, "type": kind}
        with patch.object(webapp_build, "read_version", return_value=self.version):
            result = webapp_build.build_site(str(LIGHTGUN), str(self.source / "site"))
            webapp_build.build_launcher(str(LIGHTGUN), str(self.source / "launcher"))
        (self.source / "version.json").write_text(
            json.dumps(result["version"]), encoding="utf-8")

    def publish(self, success=True):
        result = subprocess.run(
            [sys.executable, "-c", self.script], cwd=self.root,
            env={**os.environ, "PYTHONUTF8": "1",
                 "GITHUB_ENV": str(self.root / "github-env.txt")},
            capture_output=True, encoding="utf-8", timeout=30)
        self.assertEqual(result.returncode == 0, success, result.stdout + result.stderr)

    def index(self):
        return json.loads((self.repo / "versions.json").read_text(encoding="utf-8"))

    def test_first_publication_matches_build_outputs(self):
        self.publish()
        self.assertEqual(files(self.repo / "v/6.2"), files(self.source / "site"))
        for name in ("index.html", "launcher.js", ".nojekyll"):
            self.assertEqual((self.repo / name).read_bytes(),
                             (self.source / "launcher" / name).read_bytes())
        self.assertFalse((self.repo / "app.js").exists())  # Root is the launcher.
        self.assertFalse((self.repo / "v/versions.json").exists())
        self.assertEqual(self.index(), {"versions": [self.version], "latest": "6.2"})
        self.assertEqual((self.repo / "README.md").read_text(), "Preserve this documentation.\n")
        self.assertEqual((self.repo / "CNAME").read_text(), "example.invalid\n")
        self.assertIn("WEBAPP_VERSION=6.2", (self.root / "github-env.txt").read_text())

    def test_same_version_is_idempotent_and_removes_stale_files(self):
        self.publish()
        before = files(self.repo)
        self.prepare()
        self.publish()
        self.assertEqual(files(self.repo), before)
        (self.repo / "v/6.2/obsolete.js").write_text("// stale", encoding="utf-8")
        self.publish()
        self.assertEqual(files(self.repo), before)

    def test_new_version_preserves_previous_folder(self):
        self.publish()
        (self.repo / "v/6.2/kept.txt").write_text("Keep me.", encoding="utf-8")
        previous = files(self.repo / "v/6.2")
        self.prepare("6.3", "6.3.0", "beta")
        self.publish()
        self.assertEqual(files(self.repo / "v/6.2"), previous)
        self.assertEqual(files(self.repo / "v/6.3"), files(self.source / "site"))
        self.assertEqual(self.index()["latest"], "6.3")
        self.assertEqual([v["id"] for v in self.index()["versions"]], ["6.3", "6.2"])
        self.assertEqual(self.index()["versions"][0], self.version)

    def test_republish_updates_one_entry_and_preserves_extra_metadata(self):
        self.publish()
        self.prepare("6.3", "6.3.0", "beta")
        self.publish()
        previous = files(self.repo / "v/6.3")
        index = self.index()
        index["note"] = "Keep root metadata."
        index["versions"][1]["note"] = "Keep version metadata."
        (self.repo / "versions.json").write_text(json.dumps(index), encoding="utf-8")
        self.prepare("6.2", "6.2.1", "stable")  # Same ID: replace, do not create v/6.2.1.
        self.publish()
        self.assertEqual(files(self.repo / "v/6.3"), previous)
        self.assertEqual(files(self.repo / "v/6.2"), files(self.source / "site"))
        self.assertFalse((self.repo / "v/6.2.1").exists())
        self.assertEqual(self.index()["latest"], "6.2")
        self.assertEqual([v["id"] for v in self.index()["versions"]], ["6.2", "6.3"])
        self.assertEqual(self.index()["versions"][0],
                         {**self.version, "note": "Keep version metadata."})
        self.assertEqual(self.index()["note"], "Keep root metadata.")

    def test_missing_inputs_leave_published_files_untouched(self):
        self.publish()
        before = files(self.repo)
        for name in ("site/app.js", "launcher/.nojekyll", "version.json"):
            with self.subTest(name=name):
                path = self.source / name
                original = path.read_bytes()
                path.unlink()
                self.publish(success=False)
                self.assertEqual(files(self.repo), before)
                path.write_bytes(original)

    def test_invalid_catalog_is_not_overwritten(self):
        self.publish()
        for value in ("{bad", '{"versions":[{"id":"6.2"},{"id":"6.2"}]}'):
            with self.subTest(value=value):
                (self.repo / "versions.json").write_text(value, encoding="utf-8")
                before = files(self.repo)
                self.publish(success=False)
                self.assertEqual(files(self.repo), before)

    def test_invalid_version_metadata_is_rejected_before_copy(self):
        self.publish()
        before = files(self.repo)
        for version in ({"id": "../outside", "label": "bad", "type": "stable"},
                        {"id": "6.2", "label": None, "type": "stable"}):
            with self.subTest(version=version):
                (self.source / "version.json").write_text(json.dumps(version), encoding="utf-8")
                self.publish(success=False)
                self.assertEqual(files(self.repo), before)


if __name__ == "__main__":
    unittest.main(verbosity=2)
