"""What the build writes for the published site: checks of scripts/webapp_build.py.

    python tests/site-build.test.py        (from lightgun/webapp)

Two outputs, and they are different things:
  dist/site      the app of the firmware being built - one version of it, which is what
                 gets published under v/<version>/ of the site.
  dist/launcher  the home page of the site: the Connect button, which reads the version of
                 the lightgun and opens the app published for it.
"""
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
LIGHTGUN = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(LIGHTGUN, "scripts"))

import webapp_build  # noqa: E402

passed = failed = 0


def ok(condition, message):
    global passed, failed
    if condition:
        passed += 1
        print("PASS " + message)
    else:
        failed += 1
        print("FAIL " + message)


def read(path):
    with open(path, "rb") as f:
        return f.read()


def main():
    version = webapp_build.read_version(LIGHTGUN)
    ok(re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+(-[A-Za-z][A-Za-z0-9]*)?", version["id"]) is not None,
       "the version is three numbers and an optional suffix, and names the folder: " + version["id"])
    ok(version["label"] == version["id"], "shown as it is: " + version["label"])
    ok(version["numbers"].count(".") == 2 and "-" not in version["numbers"],
       "the three numbers alone, for the home page to fall back to: " + version["numbers"])
    ok(version["id"] == version["numbers"] + ("-" + version["suffix"] if version["suffix"] else ""),
       "which is the version without the suffix: " + repr(version["suffix"]))
    ok(version["tag"] == "v" + version["id"], "the Git tag is the version with a v in front: " + version["tag"])
    ok(version["prerelease"] is bool(version["suffix"]),
       "a suffix makes it a pre-release: " + str(version["prerelease"]))
    ok(version["type"] == (version["suffix"].rstrip("0123456789") or "stable"),
       "and the type is the letters of the suffix: " + version["type"])

    out = tempfile.mkdtemp(prefix="of-build-")
    try:
        # ----- the app of this firmware -----
        site = os.path.join(out, "site")
        result = webapp_build.build_site(LIGHTGUN, site)
        names = sorted(os.path.relpath(os.path.join(root, name), site).replace(os.sep, "/")
                       for root, _, files in os.walk(site) for name in files)
        ok("app.js" in names and "index.html" in names and "style.css" in names,
           "the app is there: page, style and script")
        ok(any(n.startswith("boards/pics/") for n in names), "with the board pictures")
        ok(not any(n.startswith("v/") for n in names),
           "and nothing else: no archive inside it, it IS one version")
        ok("versions.json" not in names, "no list of versions: that belongs to the site, not to a version")
        ok(".nojekyll" not in names, "no .nojekyll either: only the root of the site needs it")
        ok(result["version"]["id"] == version["id"], "the build says which version it wrote")

        app = read(os.path.join(site, "app.js")).decode("utf-8")
        ok('"version": "%s"' % version["id"] in app, "the app knows its own version")
        ok('"versionLabel": "%s"' % version["label"] in app, "and the complete number, to show it")
        ok("OF.Version = {" in app, "and carries the check of the version")

        # Ogni rimando al progetto originale, tranne i due dei crediti, manderebbe chi usa
        # la App altrove: a scaricare un firmware che non e' questo, o su una wiki che
        # documenta un altro firmware. Il bundle comprende anche le traduzioni, quindi
        # questo controllo copre sia il codice sia lang/*.json.
        upstream = sorted(set(re.findall(r"https://github\.com/TeamOpenFIRE/[A-Za-z0-9./#_-]*", app)))
        ok(upstream == ["https://github.com/TeamOpenFIRE/OpenFIRE-App",
                        "https://github.com/TeamOpenFIRE/OpenFIRE-Firmware"],
           "only the two credit links still name the original project: " + str(upstream))
        mine = "https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32"
        # Tutte le voci qui sotto sono scritte come REPO + "/qualcosa": senza questa riga
        # cambiare REPO le porterebbe tutte altrove senza che nessun controllo se ne accorga.
        ok("const REPO = '" + mine + "';" in app, "every link is built on this repository")
        for name, url in (("the firmware to download", mine + "/releases/latest"),
                                                    ("the documentation in English",
                           "en: REPO + '/blob/main/lightgun/src/README.md#english-version'"),
                          ("the documentation in Italian",
                           "it: REPO + '/blob/main/lightgun/src/README.md#versione-italiana'"),
                          ("the serial commands page in English",
                           "en: REPO + '/wiki/Serial_Commands_OpenFIRE_EN'"),
                          ("the serial commands page in Italian", "it: REPO + '/wiki/Serial_Commands_OpenFIRE_IT'"),
                          ("where to report issues", mine + "/issues")):
            ok(url in app, "and " + name + " is this project's")
        # E le descrizioni che rimandano ai comandi seriali: l'inglese nelle chiavi, la
        # pagina italiana nelle traduzioni italiane. Nessuna deve piu' citare MAMEHOOKER,
        # che era il nome della pagina sulla wiki del progetto originale.
        ok("MAMEHOOKER" not in app, "nothing points at the page of the original wiki any more")
        ok(app.count("Serial_Commands_OpenFIRE_EN") == 6 and app.count("Serial_Commands_OpenFIRE_IT") == 4,
           "the serial page is named in both languages everywhere it is mentioned: "
           + str((app.count("Serial_Commands_OpenFIRE_EN"), app.count("Serial_Commands_OpenFIRE_IT"))))

        again = webapp_build.build_site(LIGHTGUN, site)
        ok(again["changed"] == [] and again["removed"] == [],
           "building it again writes nothing: a published folder keeps a clean history")

        # Something of ours that no longer belongs to the app is cleaned up.
        stale = os.path.join(site, "boards", "pics", "gone.js")
        with open(stale, "wb") as f:
            f.write(b"// a board that no longer exists")
        webapp_build.build_site(LIGHTGUN, site)
        ok(not os.path.exists(stale), "the folder is kept tidy")

        # ----- the home page of the site -----
        launcher = os.path.join(out, "launcher")
        webapp_build.build_launcher(LIGHTGUN, launcher)
        ok(os.path.isfile(os.path.join(launcher, "index.html")), "the home page is there")
        ok(os.path.isfile(os.path.join(launcher, "launcher.js")), "with its script")
        ok(os.path.isfile(os.path.join(launcher, ".nojekyll")),
           "and .nojekyll, because this is the root of the published site")

        page = read(os.path.join(launcher, "index.html")).decode("utf-8")
        ok('<script src="launcher.js"></script>' in page, "the page loads the joined script")
        ok("OF:SCRIPTS" not in page, "and no longer lists the single ones")
        ok("welcome" not in page and "tab-pins" not in page, "the home page is not the app")

        code = read(os.path.join(launcher, "launcher.js")).decode("utf-8")
        ok("js/core/transport.js" in code and "js/core/protocol.js" in code,
           "it carries the same transport and protocol as the app: one implementation, not two")
        ok("OpenFIREshared.js" in code, "and the shared data it needs to talk to the board")
        ok("getBoardInfo" in code, "it reads only the presentation of the board")
        ok("versions.json" in code, "and looks the version up in the list")
        ok("tab-settings" not in code, "without a line of the app itself")

        empty = webapp_build.build_launcher(LIGHTGUN, launcher)
        ok(empty["changed"] == [], "building the home page again writes nothing")

        # ----- the app embedded in the lightgun -----
        device = webapp_build.bundle(LIGHTGUN, "device", "waveshare-esp32-s3-zero")["app.js"].decode("utf-8")
        ok("OF.Version = {" not in device and "of_version_jump" not in device,
           "the app embedded in the lightgun does not carry the check at all")
        ok('"version"' not in device.split("\n")[1],
           "and no version in its build info: include/web_assets.h does not change because of it")
        ok("js/app/main.js" in device and "js/core/protocol.js" in device,
           "everything else is there, as always")
    finally:
        shutil.rmtree(out, ignore_errors=True)

    print(f"\n{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
