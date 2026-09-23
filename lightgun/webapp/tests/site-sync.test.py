"""Test of the whole WebApp synchronisation: the publication step AND the push step
of GLOBAL-build-package-release.yml, run for real against a local Git repository.

    python tests/site-sync.test.py        (from lightgun/webapp)

It answers one question: what happens to the published repository when something
goes wrong - the copy breaks halfway, or somebody pushes to main while the action
is running. Nothing leaves the temporary directories: the "remote" is a bare
repository on disk, there is no network and no real GitHub.
"""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import unittest

LIGHTGUN = Path(__file__).resolve().parents[2]
WORKFLOW = LIGHTGUN.parent / ".github/workflows/GLOBAL-build-package-release.yml"
sys.path.insert(0, str(LIGHTGUN / "scripts"))

import webapp_build  # noqa: E402


def step_script(name, kind):
    """The real text of a step of the workflow (kind: 'python' or 'shell')."""
    text = WORKFLOW.read_text(encoding="utf-8-sig")
    if kind == "python":
        pattern = (r"^      - name: " + re.escape(name) + r"\n        run: \|\n"
                   r"          python3 - <<'PY'\n(.*?)^          PY[ \t]*$")
    else:
        pattern = (r"^      - name: " + re.escape(name) + r"\n"
                   r"        working-directory: webapp_repo\n        run: \|\n(.*?)(?=^      - name: |\Z)")
    blocks = re.findall(pattern, text, re.MULTILINE | re.DOTALL)
    if len(blocks) != 1:
        raise RuntimeError(f"expected exactly one step {name!r} in {WORKFLOW}")
    return textwrap.dedent(blocks[0])


def git(*args, cwd, check=True):
    return subprocess.run(("git",) + args, cwd=str(cwd), check=check,
                          capture_output=True, text=True)


class SyncTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if shutil.which("git") is None:
            raise unittest.SkipTest("git is not installed: this test needs it for the local repository")
        cls.publish = step_script("Update versioned site, launcher and versions.json", "python")
        cls.push = step_script("Commit and Push to WebApp Repo", "shell")
        compile(cls.publish, str(WORKFLOW), "exec")
        # The app is built once: it is the same for every scenario.
        cls.built = Path(tempfile.mkdtemp(prefix="of-sync-build-"))
        result = webapp_build.build_site(str(LIGHTGUN), str(cls.built / "site"))
        webapp_build.build_launcher(str(LIGHTGUN), str(cls.built / "launcher"))
        (cls.built / "version.json").write_text(
            json.dumps(result["version"], indent=2) + "\n", encoding="utf-8")
        cls.version_id = result["version"]["id"]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.built, ignore_errors=True)

    def setUp(self):
        temporary = tempfile.TemporaryDirectory(prefix="of-sync-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.source = self.root / "webapp_dist"
        shutil.copytree(self.built, self.source)

        # The "remote": a bare repository with a first commit, as the published one.
        self.remote = self.root / "remote.git"
        git("init", "--bare", "-b", "main", str(self.remote), cwd=self.root)
        seed = self.root / "seed"
        git("clone", str(self.remote), str(seed), cwd=self.root)
        git("config", "user.email", "t@t", cwd=seed)
        git("config", "user.name", "t", cwd=seed)
        (seed / "README.md").write_text("La documentazione a mano.\n", encoding="utf-8")
        (seed / "versions.json").write_text(
            json.dumps({"latest": "0.1", "versions": [{"id": "0.1", "label": "0.1.0", "type": "old"}]},
                       indent=2) + "\n", encoding="utf-8")
        git("add", "-A", cwd=seed)
        git("commit", "-m", "inizio", cwd=seed)
        git("push", "origin", "main", cwd=seed)
        shutil.rmtree(seed)

        # The checkout the action works on.
        self.repo = self.root / "webapp_repo"
        git("clone", str(self.remote), str(self.repo), cwd=self.root)
        git("config", "user.email", "a@a", cwd=self.repo)
        git("config", "user.name", "a", cwd=self.repo)

    # ---- the two steps, run as the action runs them -------------------------
    def run_publish(self):
        env = dict(os.environ, GITHUB_ENV=str(self.root / "github_env"))
        (self.root / "github_env").write_text("", encoding="utf-8")
        return subprocess.run([sys.executable, "-c", self.publish],
                              cwd=str(self.root), env=env, capture_output=True, text=True)

    def run_push(self):
        env = dict(os.environ, TAG_NAME="v6.2.1", WEBAPP_VERSION=self.version_id,
                   GIT_AUTHOR_NAME="a", GIT_AUTHOR_EMAIL="a@a",
                   GIT_COMMITTER_NAME="a", GIT_COMMITTER_EMAIL="a@a")
        return subprocess.run(["bash", "-e", "-c", self.push], cwd=str(self.repo),
                              env=env, capture_output=True, text=True)

    def remote_head(self):
        return git("rev-parse", "main", cwd=self.remote).stdout.strip()

    def remote_files(self):
        return sorted(git("ls-tree", "-r", "--name-only", "main", cwd=self.remote).stdout.split())

    # ---- 1. il caso normale -------------------------------------------------
    def test_pubblicazione_normale(self):
        prima = self.remote_head()
        self.assertEqual(self.run_publish().returncode, 0)
        esito = self.run_push()
        self.assertEqual(esito.returncode, 0, esito.stdout + esito.stderr)
        self.assertNotEqual(self.remote_head(), prima, "il commit deve essere arrivato")
        files = self.remote_files()
        for atteso in ("index.html", "launcher.js", ".nojekyll", "versions.json",
                       f"v/{self.version_id}/app.js", f"v/{self.version_id}/boards/pics/rpipico.js"):
            self.assertIn(atteso, files)
        self.assertIn("README.md", files, "la documentazione a mano resta")
        indice = json.loads(git("show", "main:versions.json", cwd=self.remote).stdout)
        self.assertEqual(indice["latest"], self.version_id)
        self.assertEqual([v["id"] for v in indice["versions"]], [self.version_id, "0.1"])

    # ---- 2. la copia si rompe a meta' ---------------------------------------
    def test_copia_interrotta_non_pubblica_niente(self):
        """Il caso peggiore: v/<id> gia' pubblicata, la cartella viene svuotata e la
        copia fallisce. Il danno deve restare nella copia usa-e-getta del runner."""
        self.assertEqual(self.run_publish().returncode, 0)
        self.assertEqual(self.run_push().returncode, 0)
        pubblicato = self.remote_head()
        self.assertIn(f"v/{self.version_id}/app.js", self.remote_files())

        # Il passo vero gira dopo il guasto: copytree fallisce subito dopo il rmtree.
        rotto = ("import shutil\n"
                 "shutil.copytree = lambda *a, **k: (_ for _ in ()).throw(OSError('disco pieno'))\n")
        esito = subprocess.run([sys.executable, "-c", rotto + self.publish], cwd=str(self.root),
                               env=dict(os.environ, GITHUB_ENV=str(self.root / "github_env")),
                               capture_output=True, text=True)
        self.assertNotEqual(esito.returncode, 0, "il passo deve fallire")
        self.assertIn("disco pieno", esito.stderr)

        # La copia di lavoro e' davvero rovinata...
        self.assertFalse((self.repo / "v" / self.version_id).exists(),
                         "la cartella e' stata svuotata: e' il momento scomodo")
        # ...ma il repository pubblicato non se ne accorge, perche' in GitHub Actions
        # un passo fallito ferma il job e il push non parte mai.
        self.assertEqual(self.remote_head(), pubblicato, "il repository pubblicato non e' cambiato")
        self.assertIn(f"v/{self.version_id}/app.js", self.remote_files())

    # ---- 3. qualcuno scrive su main mentre l'action gira --------------------
    def test_modifica_concorrente_non_viene_persa(self):
        altro = self.root / "altro"
        git("clone", str(self.remote), str(altro), cwd=self.root)
        git("config", "user.email", "m@m", cwd=altro)
        git("config", "user.name", "m", cwd=altro)
        (altro / "README.md").write_text("Modifica fatta a mano durante la release.\n", encoding="utf-8")
        git("commit", "-am", "modifica a mano", cwd=altro)
        git("push", "origin", "main", cwd=altro)
        suo = self.remote_head()

        self.assertEqual(self.run_publish().returncode, 0)
        esito = self.run_push()
        self.assertNotEqual(esito.returncode, 0, "il push deve fallire, non forzare")
        self.assertIn("ERRORE FATALE", esito.stdout)
        self.assertEqual(self.remote_head(), suo, "la modifica a mano e' ancora li'")
        self.assertIn("Modifica fatta a mano",
                      git("show", "main:README.md", cwd=self.remote).stdout)
        # E il lavoro non e' perso: sta nel commit locale, basta rilanciare.
        self.assertIn("Auto-sync WebApp", git("log", "-1", "--pretty=%s", cwd=self.repo).stdout)

    # ---- 4. niente da fare --------------------------------------------------
    def test_nessuna_modifica_nessun_commit(self):
        self.assertEqual(self.run_publish().returncode, 0)
        self.assertEqual(self.run_push().returncode, 0)
        dopo = self.remote_head()
        git("clean", "-fdx", cwd=self.repo, check=False)
        self.assertEqual(self.run_publish().returncode, 0)
        esito = self.run_push()
        self.assertEqual(esito.returncode, 0)
        self.assertIn("nessun commit necessario", esito.stdout)
        self.assertEqual(self.remote_head(), dopo)

    # ---- 5. un file che sparisce viene rimosso anche dal repository ---------
    def test_file_sparito_viene_tolto(self):
        self.assertEqual(self.run_publish().returncode, 0)
        self.assertEqual(self.run_push().returncode, 0)
        vittima = f"v/{self.version_id}/boards/pics/rpipico.js"
        self.assertIn(vittima, self.remote_files())
        (self.source / "site" / "boards" / "pics" / "rpipico.js").unlink()
        self.assertEqual(self.run_publish().returncode, 0)
        self.assertEqual(self.run_push().returncode, 0)
        self.assertNotIn(vittima, self.remote_files(), "git add deve registrare anche le rimozioni")


if __name__ == "__main__":
    unittest.main(verbosity=2)
