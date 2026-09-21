"""The published site keeps its older versions: checks of scripts/webapp_build.py.

    python tests/site-archive.test.py        (from lightgun/webapp)

Nothing is built for real twice: the version of the firmware is replaced between the two
builds, which is exactly what raising the number in src/OpenFIREversion.h does.
"""
import json
import os
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


def build(out, version_id, label, kind="stable"):
    """One build of the site, as if OpenFIREversion.h said that version."""
    original = webapp_build.read_version
    webapp_build.read_version = lambda project_dir: {"id": version_id, "label": label, "type": kind}
    try:
        return webapp_build.build_site(LIGHTGUN, out)
    finally:
        webapp_build.read_version = original


def main():
    # The version really written in the firmware, before anything is replaced.
    real = webapp_build.read_version(LIGHTGUN)
    ok(real["id"] == "%.1f" % float(real["id"]),
       "the version names the folder as the firmware sends it: " + real["id"])
    ok(real["label"].count(".") == 2, "the complete number is shown in the list: " + real["label"])

    out = tempfile.mkdtemp(prefix="of-site-")
    try:
        # ----- first version -----
        first = build(out, "6.2", "6.2.0")
        root_app = os.path.join(out, "app.js")
        archived_app = os.path.join(out, "v", "6.2", "app.js")
        ok(os.path.isfile(root_app), "the current version is in the root, where everybody goes")
        ok(os.path.isfile(archived_app), "the same app is in v/6.2/")
        ok(read(root_app) == read(archived_app), "root and archive are the same app")
        ok(os.path.isfile(os.path.join(out, ".nojekyll")), "GitHub Pages serves the files as they are")
        ok(not os.path.exists(os.path.join(out, "v", "6.2", ".nojekyll")),
           ".nojekyll only in the root, where GitHub Pages reads it")
        ok(os.path.isfile(os.path.join(out, "v", "6.2", "boards", "pics", "generic.js")),
           "an archived version carries its own board pictures: it does not depend on the others")
        ok(first["version"]["id"] == "6.2", "the build says which version it wrote")

        info = json.loads(read(os.path.join(out, "v", "versions.json")).decode("utf-8"))
        ok(info["latest"] == "6.2", "versions.json: the root is 6.2")
        ok([v["id"] for v in info["versions"]] == ["6.2"], "versions.json: one version")
        ok(info["versions"][0]["label"] == "6.2.0" and info["versions"][0]["type"] == "stable",
           "versions.json: complete number and type")

        page = read(os.path.join(out, "v", "index.html")).decode("utf-8")
        ok('"latest": "6.2"' in page or '"latest":"6.2"' in page, "v/index.html carries the list inside it")
        ok("/*OF:VERSIONS*/" not in page, "v/index.html: the marker was replaced")

        app = read(root_app).decode("utf-8")
        ok('"version": "6.2"' in app, "the app of the site knows its own version")
        ok('"versionLabel": "6.2.0"' in app, "and the complete number, to show it")

        # ----- building again changes nothing -----
        again = build(out, "6.2", "6.2.0")
        ok(again["changed"] == [] and again["removed"] == [],
           "a second identical build writes nothing: the repository keeps a clean history")

        # Something of ours that no longer belongs to the app is cleaned up...
        stale = os.path.join(out, "v", "6.2", "boards", "pics", "gone.js")
        with open(stale, "wb") as f:
            f.write(b"// a board that no longer exists")
        build(out, "6.2", "6.2.0")
        ok(not os.path.exists(stale), "the folder of the version being built is kept tidy")

        # ----- the version number is raised -----
        marker = os.path.join(out, "v", "6.2", "app.js")
        before = read(marker)
        with open(os.path.join(out, "v", "6.2", "LEGGIMI.txt"), "wb") as f:
            f.write(b"kept")
        third = build(out, "6.3", "6.3.0", "beta")

        ok(read(marker) == before, "the previous version stays exactly as it was")
        ok(os.path.exists(os.path.join(out, "v", "6.2", "LEGGIMI.txt")),
           "and nothing of it is removed")
        ok(os.path.isfile(os.path.join(out, "v", "6.3", "app.js")), "the new version has its folder")
        ok(read(os.path.join(out, "app.js")) == read(os.path.join(out, "v", "6.3", "app.js")),
           "the root is now the new version")
        ok('"version": "6.3"' in read(os.path.join(out, "v", "6.3", "app.js")).decode("utf-8"),
           "the archived app carries the version it was built as")
        ok("v/6.3/app.js" in third["changed"], "the report names the folder of the archive")

        info = json.loads(read(os.path.join(out, "v", "versions.json")).decode("utf-8"))
        ok(info["latest"] == "6.3", "versions.json: the root is the new version")
        ok([v["id"] for v in info["versions"]] == ["6.3", "6.2"], "versions.json: newest first")
        ok(info["versions"][1]["label"] == "6.2.0",
           "the old version keeps the number it was published with")
        ok(info["versions"][0]["type"] == "beta", "and the new one says what it is")

        # ----- a version deleted by hand disappears from the list -----
        shutil.rmtree(os.path.join(out, "v", "6.2"))
        build(out, "6.3", "6.3.0", "beta")
        info = json.loads(read(os.path.join(out, "v", "versions.json")).decode("utf-8"))
        ok([v["id"] for v in info["versions"]] == ["6.3"],
           "a version removed from the folder is removed from the list")

        # ----- order of numbers -----
        ok(sorted(["6.9", "6.10", "7.0", "6.2"], key=webapp_build._version_key) == ["7.0", "6.10", "6.9", "6.2"],
           "6.10 is newer than 6.9: the list is ordered by number, not by text")

        # ----- the firmware is not touched -----
        device = webapp_build.bundle(LIGHTGUN, "device", "waveshare-esp32-s3-zero")
        head = device["app.js"].decode("utf-8")[:400]
        ok('"version"' not in head and '"target": "device"' in head,
           "the app embedded in the lightgun has no version added: include/web_assets.h does not change")
    finally:
        shutil.rmtree(out, ignore_errors=True)

    print(f"\n{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
