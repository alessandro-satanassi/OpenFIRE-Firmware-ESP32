"""PlatformIO pre-script: embeds the reduced web app in the firmware (include/web_assets.h).

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
