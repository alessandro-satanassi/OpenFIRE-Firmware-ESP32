"""PlatformIO pre-script: embeds the reduced web app in the firmware (include/web_assets.h)
and, when asked, updates the folder of the published site.

Site (GitHub Pages / Tauri): every build writes it into dist/site. To keep the folder of
the site repository up to date instead, set its path in platformio.ini, under [env] (or in
a single environment):

    custom_webapp_site_dir = ../../OpenFIRE-WebApp

A relative path starts from the lightgun folder; the environment variable
OPENFIRE_WEBAPP_SITE_DIR has priority over platformio.ini, and the value "off" writes no
site at all. The folder receives index.html, style.css, app.js, boards/pics/*.js and
.nojekyll; everything else in it (.git, README, LICENSE, CNAME) is left untouched.

Every build re-reads src/boards (OpenFIREshared.h, boardPics/boards.qrc and the board SVGs) and
webapp/lang/*.json, so changes to boards, pin maps, alternative layouts and pictures are
picked up automatically (the generated files of the webapp folder are refreshed too).

The board is the one the firmware itself selects (OPENFIRE_BOARD #ifdef chain in
src/boards/OpenFIREshared.h) from the -D defines of the environment and board manifest.
Only ESP32 environments serve the web app; the header is rewritten only when it changes,
so an unchanged web app does not trigger a new compile of OpenFIREweb.cpp.
"""
import os
import re
import subprocess
import sys

Import("env")  # noqa: F821 - provided by PlatformIO/SCons

PROJECT_DIR = env.subst("$PROJECT_DIR")  # noqa: F821
sys.path.insert(0, os.path.join(PROJECT_DIR, "scripts"))

import build_board_pics  # noqa: E402
import build_shared_js   # noqa: E402
import webapp_build      # noqa: E402

ENV_NAME = env.subst("$PIOENV")  # noqa: F821
HEADER_PATH = os.path.join(PROJECT_DIR, "include", "web_assets.h")


def site_dir():
    """Folder of the site: dist/site, the configured one, or None ("off")."""
    value = os.environ.get("OPENFIRE_WEBAPP_SITE_DIR", "").strip()
    if not value:
        try:
            value = str(env.GetProjectOption("custom_webapp_site_dir", "") or "").strip()  # noqa: F821
        except Exception:
            value = ""
    if value.lower() in ("off", "no", "none", "0"):
        return None
    if not value:
        return os.path.join(PROJECT_DIR, "dist", "site")
    value = os.path.expanduser(env.subst(value))  # noqa: F821
    return value if os.path.isabs(value) else os.path.abspath(os.path.join(PROJECT_DIR, value))


def update_site():
    """Writes the site into the configured folder (a build must not fail because of it)."""
    target = site_dir()
    if not target:
        return
    try:
        result = webapp_build.build_site(PROJECT_DIR, target)
        changed = len(result["changed"])
        removed = len(result["removed"])
        state = "unchanged" if not changed and not removed else \
            f"{changed} file(s) updated" + (f", {removed} removed" if removed else "")
        print(f"[WebApp] site in {target}: {state}.")
    except Exception as error:
        print(f"[WebApp] WARNING: the site could not be written in {target}: {error}")


def _flatten(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        items = []
        for item in value:
            items += _flatten(item)
        return items
    return [str(value)]


def build_defines():
    """Names defined with -D in build_flags, the board manifest and the environment."""
    texts = []
    try:
        texts += _flatten(env.GetProjectOption("build_flags", ""))  # noqa: F821
    except Exception:
        pass
    texts += _flatten(env.get("BUILD_FLAGS", []))  # noqa: F821
    try:
        texts += _flatten(env.BoardConfig().get("build.extra_flags", ""))  # noqa: F821
    except Exception:
        pass
    defines = set()
    for text in texts:
        try:
            text = env.subst(text)  # noqa: F821
        except Exception:
            pass
        text = re.sub(r';.*', '', text)
        defines.update(re.findall(r'-D\s*([A-Za-z_]\w*)', text))
    for item in _flatten(env.get("CPPDEFINES", [])):  # noqa: F821
        match = re.match(r"\(?'?([A-Za-z_]\w*)", item)
        if match:
            defines.add(match.group(1))
    return defines


def ensure_pillow():
    """Installs Pillow in the PlatformIO Python when a board photo has to be re-compressed."""
    pending = build_board_pics.pictures_needing_pillow(PROJECT_DIR)
    if not pending or build_board_pics.load_pillow():
        return
    python = env.subst("$PYTHONEXE") or sys.executable  # noqa: F821
    print(f"[WebApp] {', '.join(pending)}: installing Pillow in the PlatformIO Python to compress the board photos...")
    try:
        result = subprocess.run([python, "-m", "pip", "install", "--disable-pip-version-check", "-q", "pillow"],
                                timeout=300)
        installed = result.returncode == 0 and build_board_pics.load_pillow()
    except Exception as error:
        print(f"[WebApp] pip failed: {error}")
        installed = False
    if not installed:
        print("[WebApp] WARNING: Pillow could not be installed (no internet?): the new board photos are used "
              "uncompressed for now; the next build tries again.")


def main():
    platform = env.get("PIOPLATFORM", "")  # noqa: F821
    board = build_shared_js.detect_board(PROJECT_DIR, build_defines())
    ensure_pillow()
    webapp_build.refresh_dev_files(PROJECT_DIR)
    update_site()
    if "espressif32" not in platform and "esp32" not in board:
        print(f"[WebApp] {ENV_NAME}: web app not embedded (only ESP32 boards serve it).")
        return
    known = build_shared_js.extract_shared(PROJECT_DIR).get("boardsBoxPositions", {})
    if board not in known or board.startswith("generic"):
        print(f"[WebApp] WARNING: {ENV_NAME}: board '{board}' is not a known board: check the -D ARDUINO_... define "
              f"of the environment and the OPENFIRE_BOARD #ifdef chain in OpenFIREshared.h (generic maps are used).")
    total, changed = webapp_build.write_device_header(PROJECT_DIR, board, HEADER_PATH)
    state = "updated" if changed else "unchanged"
    print(f"[WebApp] {ENV_NAME}: board '{board}', {total / 1024:.1f} KB gzipped in include/web_assets.h ({state}).")


try:
    main()
except Exception as error:  # a broken web app must be visible, not a silent stale header
    print(f"[WebApp] ERROR while packing the web app: {error}")
    env.Exit(1)  # noqa: F821
