"""Lo schema della versione: src/OpenFIREversion.h e il job "version" della action.

    python tests/version-scheme.test.py        (from lightgun/webapp)

Prova due cose che devono dire sempre la stessa identica versione:
  - read_version(), che da' il nome alla cartella pubblicata della Web App;
  - il passo della workflow che ne ricava tag, flag di pre-release e titolo.

Il testo del passo viene estratto dal file della workflow, non riscritto qui: se quel
passo cambia nome o forma la prova se ne accorge. Non tocca il repository: lavora su
copie temporanee dell'header.
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
HEADER = LIGHTGUN / "src" / "OpenFIREversion.h"
sys.path.insert(0, str(LIGHTGUN / "scripts"))

import webapp_build  # noqa: E402


def version_step():
    """Il vero heredoc del passo che legge la versione, preso dalla workflow.

    Niente espressioni regolari golose su tutto il file: si cerca il nome del passo,
    poi l'inizio e la fine del suo heredoc."""
    text = WORKFLOW.read_text(encoding="utf-8-sig")
    intestazione = "      - name: Read the version of the firmware\n"
    if text.count(intestazione) != 1:
        raise RuntimeError("expected exactly one version step in " + str(WORKFLOW))
    resto = text[text.index(intestazione) + len(intestazione):]
    apertura = "          python3 - <<'PY'\n"
    if apertura not in resto:
        raise RuntimeError("the version step no longer runs a python heredoc")
    corpo = resto[resto.index(apertura) + len(apertura):]
    chiusura = "          PY\n"
    if chiusura not in corpo:
        raise RuntimeError("the heredoc of the version step is not closed")
    return textwrap.dedent(corpo[:corpo.index(chiusura)])


class SchemeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.step = version_step()
        compile(cls.step, str(WORKFLOW), "exec")
        cls.original = HEADER.read_bytes()

    def setUp(self):
        self.addCleanup(HEADER.write_bytes, self.original)
        temporary = tempfile.TemporaryDirectory(prefix="of-version-")
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    # ---- l'header, con la versione voluta ----------------------------------
    def write_header(self, major=7, minor=0, patch=0, suffix="", short=None):
        text = self.original.decode("utf-8")
        text = re.sub(r"#define OPENFIRE_VERSION_MAJOR \d+", f"#define OPENFIRE_VERSION_MAJOR {major}", text)
        text = re.sub(r"#define OPENFIRE_VERSION_MINOR \d+", f"#define OPENFIRE_VERSION_MINOR {minor}", text)
        text = re.sub(r"#define OPENFIRE_VERSION_PATCH \d+", f"#define OPENFIRE_VERSION_PATCH {patch}", text)
        text = re.sub(r'#define OPENFIRE_VERSION_SUFFIX "[^"]*"', f'#define OPENFIRE_VERSION_SUFFIX "{suffix}"', text)
        text = re.sub(r"#define OPENFIRE_VERSION [\d.]+",
                      f"#define OPENFIRE_VERSION {short if short is not None else f'{major}.{minor}'}", text)
        HEADER.write_bytes(text.encode("utf-8"))

    # ---- il passo della workflow, eseguito davvero -------------------------
    def run_step(self, force_prerelease="false", release_name=""):
        checkout = self.root / "checkout"
        shutil.rmtree(checkout, ignore_errors=True)
        (checkout / "lightgun" / "src").mkdir(parents=True)
        shutil.copytree(LIGHTGUN / "scripts", checkout / "lightgun" / "scripts")
        shutil.copyfile(HEADER, checkout / "lightgun" / "src" / "OpenFIREversion.h")
        output = self.root / "github_output"
        output.write_text("", encoding="utf-8")
        done = subprocess.run(
            [sys.executable, "-c", self.step], cwd=str(checkout),
            env={**os.environ, "PYTHONUTF8": "1", "GITHUB_OUTPUT": str(output),
                 "FORCE_PRERELEASE": force_prerelease, "RELEASE_NAME": release_name},
            capture_output=True, encoding="utf-8", timeout=60)
        values = dict(line.split("=", 1) for line in output.read_text(encoding="utf-8").splitlines() if "=" in line)
        return done, values

    # ---- 1. una versione definitiva ---------------------------------------
    def test_versione_definitiva(self):
        self.write_header(7, 0, 0, "")
        version = webapp_build.read_version(str(LIGHTGUN))
        self.assertEqual(version["id"], "7.0.0")
        self.assertEqual(version["label"], "7.0.0")
        self.assertEqual(version["numbers"], "7.0.0")
        self.assertEqual(version["suffix"], "")
        self.assertEqual(version["type"], "stable")
        self.assertEqual(version["tag"], "v7.0.0")
        self.assertIs(version["prerelease"], False)

        done, values = self.run_step()
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual(values["version"], "7.0.0")
        self.assertEqual(values["tag"], "v7.0.0")
        self.assertEqual(values["prerelease"], "false")
        self.assertEqual(values["release_name"], "v7.0.0 OpenFIRE Firmware ESP32 - GLOBAL Release")

    # ---- 2. con un suffisso ------------------------------------------------
    def test_versioni_col_suffisso(self):
        for suffix, kind in (("beta1", "beta"), ("rc2", "rc"), ("prerelease1", "prerelease"), ("alpha", "alpha")):
            with self.subTest(suffix=suffix):
                self.write_header(7, 1, 3, suffix)
                version = webapp_build.read_version(str(LIGHTGUN))
                self.assertEqual(version["id"], "7.1.3-" + suffix)
                self.assertEqual(version["numbers"], "7.1.3")
                self.assertEqual(version["type"], kind)
                self.assertEqual(version["tag"], "v7.1.3-" + suffix)
                self.assertIs(version["prerelease"], True)

                done, values = self.run_step()
                self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
                self.assertEqual(values["tag"], "v7.1.3-" + suffix)
                self.assertEqual(values["prerelease"], "true", "il suffisso basta da solo")

    # ---- 3. gli scavalchi della maschera -----------------------------------
    def test_scavalchi(self):
        self.write_header(7, 0, 0, "")
        done, values = self.run_step(force_prerelease="true")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual(values["prerelease"], "true", "si puo' forzare anche senza suffisso")

        done, values = self.run_step(release_name="  Il mio titolo  ")
        self.assertEqual(done.returncode, 0, done.stdout + done.stderr)
        self.assertEqual(values["release_name"], "Il mio titolo")
        self.assertEqual(values["tag"], "v7.0.0", "il titolo non tocca il tag")

        # Una versione col suffisso resta pre-release anche senza forzarla.
        self.write_header(7, 0, 0, "rc1")
        done, values = self.run_step(force_prerelease="false")
        self.assertEqual(values["prerelease"], "true")

    # ---- 4. un header che si contraddice si ferma qui ----------------------
    def test_header_incoerente(self):
        casi = [
            ("OPENFIRE_VERSION fermo alla versione di prima", dict(major=7, minor=0, short="6.2"), "6.2"),
            ("OPENFIRE_VERSION con il minore sbagliato", dict(major=7, minor=1, short="7.0"), "7.0"),
            ("MINOR oltre 9", dict(major=7, minor=10, short="7.10"), "MINOR"),
        ]
        for nome, campi, atteso in casi:
            with self.subTest(nome=nome):
                self.write_header(**campi)
                with self.assertRaises(RuntimeError) as errore:
                    webapp_build.read_version(str(LIGHTGUN))
                self.assertIn(atteso, str(errore.exception))
                done, _ = self.run_step()
                self.assertNotEqual(done.returncode, 0, "anche la action deve fermarsi")

    def test_suffisso_malformato(self):
        for suffix in ("1beta", "beta 1", "beta.1", "-beta", "beta/1", "bêta"):
            with self.subTest(suffix=suffix):
                self.write_header(7, 0, 0, suffix)
                with self.assertRaises(RuntimeError) as errore:
                    webapp_build.read_version(str(LIGHTGUN))
                self.assertIn("SUFFIX", str(errore.exception))

    # ---- 5. la versione finisce davvero nella App --------------------------
    def test_la_app_porta_la_versione(self):
        self.write_header(7, 2, 4, "rc1")
        out = self.root / "site"
        result = webapp_build.build_site(str(LIGHTGUN), str(out))
        self.assertEqual(result["version"]["id"], "7.2.4-rc1")
        app = (out / "app.js").read_text(encoding="utf-8")
        build = json.loads(re.search(r"\.BUILD = (\{.*?\});", app).group(1))
        self.assertEqual(build["version"], "7.2.4-rc1", "la App sa di che versione e'")
        self.assertEqual(build["versionLabel"], "7.2.4-rc1")
        # E il firmware manda esattamente la stessa stringa: stessa fonte, stesso formato.
        main = (LIGHTGUN / "src" / "main.cpp").read_text(encoding="utf-8")
        self.assertIn('"%d.%d.%d%s"', main)
        self.assertIn("OPENFIRE_VERSION_TAIL", main)

    # ---- 6. quello che il firmware compila davvero -------------------------
    def test_la_stringa_compilata_dal_firmware(self):
        """Compila l'header con un compilatore vero e guarda cosa ne esce."""
        if shutil.which("g++") is None:
            self.skipTest("g++ non disponibile")
        for major, minor, patch, suffix, atteso in ((7, 0, 0, "", "7.0.0"),
                                                    (7, 0, 0, "beta1", "7.0.0-beta1"),
                                                    (8, 3, 12, "rc2", "8.3.12-rc2")):
            with self.subTest(atteso=atteso):
                self.write_header(major, minor, patch, suffix)
                sorgente = self.root / "prova.cpp"
                sorgente.write_text(
                    '#include <cstdio>\n#include "OpenFIREversion.h"\n'
                    'int main(){ char v[32]; snprintf(v, sizeof v, "%d.%d.%d%s",\n'
                    '  OPENFIRE_VERSION_MAJOR, OPENFIRE_VERSION_MINOR, OPENFIRE_VERSION_PATCH,\n'
                    '  OPENFIRE_VERSION_TAIL); printf("%s", v); return 0; }\n', encoding="utf-8")
                shutil.copyfile(HEADER, self.root / "OpenFIREversion.h")
                binario = self.root / "prova"
                build = subprocess.run(["g++", "-std=c++17", "-Wall", "-Werror", "-o", str(binario), str(sorgente)],
                                       cwd=str(self.root), capture_output=True, text=True, timeout=120)
                self.assertEqual(build.returncode, 0, build.stderr)
                done = subprocess.run([str(binario)], capture_output=True, text=True, timeout=30)
                self.assertEqual(done.stdout, atteso,
                                 "la lightgun deve mandare la stessa stringa che nomina la cartella")
                self.assertEqual(webapp_build.read_version(str(LIGHTGUN))["id"], atteso)


if __name__ == "__main__":
    unittest.main(verbosity=2)
