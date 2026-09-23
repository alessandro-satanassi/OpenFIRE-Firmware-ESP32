"""Builds the OpenFIRE web app from the webapp folder.

Outputs
  device : the reduced app embedded in the firmware (include/web_assets.h, gzipped):
           only the lightgun's own board, WebSocket connection, its picture inside app.js.
  site   : the complete app (every board, Web Serial) for GitHub Pages and Tauri,
           written to a folder (default: dist/site); one script per board picture in
           boards/pics/ (<picture>.js), so the page also works opened as a local file.
           Only the files of the app are written: the rest of the folder (.git, README,
           LICENSE, CNAME of a site repository) is left untouched, so the site folder
           can be the checkout of the site repository (see scripts/pack_webapp.py).

  launcher : the home page of the published site (default: dist/launcher): the Connect
           button, nothing else. It reads the version of the firmware with the same
           transport and protocol as the app, looks it up in versions.json and opens the
           app published for it (v/<version>/); when that version is not there, it offers
           the ones that are.

What gets published: every firmware has its own app, and dist/site is exactly that app -
one version of it, not a site. Publishing a version means copying dist/site into
v/<version>/ of the site repository and adding it to versions.json; the home page
(dist/launcher) sits at the root and sends everybody to the right version. The version is
the string the firmware itself sends when the app docks, "%.1f" of OPENFIRE_VERSION in
src/OpenFIREversion.h, so 6.2: it is the only version an already installed lightgun can
tell the app about, and it names the folder.

index.html lists the scripts between <!-- OF:SCRIPTS --> and <!-- /OF:SCRIPTS -->: in
the build they are joined, in order, into a single app.js. Two of them are generated:
  boards/OpenFIREshared.js  from src/boards/OpenFIREshared.h (build_shared_js.py)
  lang/translations.js      from webapp/lang/*.json          (build_lang_js.py)
Serving the webapp folder as it is (python -m http.server) runs the unbundled app.

Command line (no PlatformIO needed):
  python scripts/webapp_build.py site   [--out DIR]
  python scripts/webapp_build.py device --board <board> [--header FILE | --out DIR]
  python scripts/webapp_build.py sizes
"""
import argparse
import gzip
import json
import os
import re
import sys

SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

import build_board_pics  # noqa: E402
import build_lang_js     # noqa: E402
import build_shared_js   # noqa: E402

SCRIPTS_BLOCK = re.compile(r'<!--\s*OF:SCRIPTS\b.*?-->(.*?)<!--\s*/OF:SCRIPTS\s*-->', re.DOTALL)
SCRIPT_SRC = re.compile(r'<script\b[^>]*\bsrc="([^"]+)"[^>]*>\s*</script>')
GENERATED_SHARED = "boards/OpenFIREshared.js"
GENERATED_TRANSLATIONS = "lang/translations.js"
# Scripts the app embedded in the lightgun does not need: the gun serves the app of its own
# firmware, so there is no published version to choose. Left out of the device build, which
# is the flash it saves; every call to them in the rest of the app is guarded.
SITE_ONLY_SCRIPTS = ("js/app/version.js",)
DEVICE_ASSETS = (  # (URL path, C variable, source)
    ("/", "web_index_html", "index.html"),
    ("/style.css", "web_style_css", "style.css"),
    ("/app.js", "web_app_js", "app.js"),
)


def project_dir_default():
    return os.path.abspath(os.path.join(SCRIPTS_DIR, ".."))


def _read_text(path):
    with open(path, "r", encoding="utf-8-sig") as f:
        return f.read()


def _define(text, name):
    """Value of a #define, without its trailing comment (None when it is not there)."""
    match = re.search(r'^[ \t]*#define[ \t]+' + name + r'[ \t]+(.+?)[ \t]*(?://.*)?$', text, re.M)
    return match.group(1).strip() if match else None


def read_version(project_dir):
    """Version of the firmware, from src/OpenFIREversion.h.

    id     the string the firmware sends when the app docks ("%.1f" of OPENFIRE_VERSION,
           so "6.2"): it names the folder of the archived site, because it is the only
           version an already installed lightgun can tell the app about.
    label  the complete number, 6.2.0, shown in the list of versions.
    type   stable, beta, rc... as written in the header.
    """
    text = _read_text(os.path.join(project_dir, "src", "OpenFIREversion.h"))
    raw = _define(text, "OPENFIRE_VERSION")
    if raw is None:
        raise RuntimeError("src/OpenFIREversion.h: #define OPENFIRE_VERSION is missing")
    try:
        version_id = "%.1f" % float(raw)  # exactly what the firmware prints when docking
    except ValueError:
        raise RuntimeError(f"src/OpenFIREversion.h: OPENFIRE_VERSION is not a number ({raw!r})")

    parts = []
    for name in ("OPENFIRE_VERSION_MAJOR", "OPENFIRE_VERSION_MINOR", "OPENFIRE_VERSION_PATCH"):
        value = _define(text, name)
        if value is not None and value.isdigit():
            parts.append(value)
    kind = (_define(text, "OPENFIRE_VERSION_TYPE") or "").strip().strip('"').strip()
    return {"id": version_id, "label": ".".join(parts) if len(parts) == 3 else version_id, "type": kind}


def script_list(webapp_dir):
    html = _read_text(os.path.join(webapp_dir, "index.html"))
    block = SCRIPTS_BLOCK.search(html)
    if not block:
        raise RuntimeError("index.html: missing <!-- OF:SCRIPTS --> ... <!-- /OF:SCRIPTS --> block")
    return html, block, SCRIPT_SRC.findall(re.sub(r'<!--.*?-->', '', block.group(1), flags=re.DOTALL))


def bundle(project_dir, target, board=None):
    """Returns {'index.html': bytes, 'style.css': bytes, 'app.js': bytes} for the target."""
    webapp_dir = os.path.join(project_dir, "webapp")
    html, block, scripts = script_list(webapp_dir)

    shared = build_shared_js.extract_shared(project_dir)
    if target == "device":
        shared = build_shared_js.reduce_for_board(shared, board)
    generated = {
        GENERATED_SHARED: build_shared_js.render_shared_js(shared),
        GENERATED_TRANSLATIONS: build_lang_js.render_translations_js(build_lang_js.load_translations(webapp_dir)),
    }

    build_info = {
        "target": target,
        "board": board if target == "device" else None,
    }
    if target == "site":
        # Only the site needs it: the app served by the lightgun is by definition the
        # one of its own firmware, and the device header stays byte for byte the same.
        version = read_version(project_dir)
        build_info["version"] = version["id"]
        build_info["versionLabel"] = version["label"]
    parts = ["/* OpenFIRE Web App - generated by scripts/webapp_build.py, do not edit. */\n"
             "(globalThis.OF = globalThis.OF || {}).BUILD = " + json.dumps(build_info) + ";\n"]
    for src in scripts:
        if target == "device" and src in SITE_ONLY_SCRIPTS:
            continue
        if src in generated:
            code = generated[src]
        else:
            path = os.path.join(webapp_dir, *src.split("/"))
            if not os.path.exists(path):
                raise RuntimeError(f"index.html lists {src}, but webapp/{src} does not exist")
            code = _read_text(path)
        parts.append(f"\n/* ---- {src} ---- */\n{code.rstrip()}\n;\n")

    code = "".join(parts)
    # The published app checks whether it is the one of the connected lightgun; the app of
    # the lightgun must not carry that check at all. Neither is visible by looking at the
    # page, so an old copy of a source file would go unnoticed until somebody connects.
    if target == "site" and "OF.Version = {" not in code:
        raise RuntimeError("site: the joined scripts do not contain js/app/version.js. "
                           "An old webapp/index.html or an old js/app/version.js of another "
                           "source tree leaves the published app without its version check.")
    if target == "device" and "OF.Version = {" in code:
        raise RuntimeError("device: js/app/version.js must not be embedded in the lightgun "
                           "(see SITE_ONLY_SCRIPTS).")

    page = html[:block.start()] + '<script src="app.js"></script>' + html[block.end():]
    return {
        "index.html": page.encode("utf-8"),
        "style.css": _read_text(os.path.join(webapp_dir, "style.css")).encode("utf-8"),
        "app.js": code.encode("utf-8"),
    }


def refresh_dev_files(project_dir):
    """Keeps the generated files of the unbundled webapp folder up to date."""
    webapp_dir = os.path.join(project_dir, "webapp")
    build_shared_js.generate_shared_js(project_dir, webapp_dir)
    build_lang_js.build_lang_js(webapp_dir)
    build_board_pics.all_pictures(project_dir)


def _gzip(data):
    return gzip.compress(data, compresslevel=9, mtime=0)


def device_picture_script(project_dir, board):
    """Script registering the picture of the lightgun's board (generic picture when it has none)."""
    images = build_shared_js.extract_shared(project_dir).get("boardImagesMap", {})
    for key in (board, "generic"):
        if key in images:
            picture = build_board_pics.board_picture(project_dir, images[key])
            if picture is not None:
                return build_board_pics.render_script(images[key], picture)
    return build_board_pics.render_script("generic.svg", b"<svg xmlns='http://www.w3.org/2000/svg'/>")


def device_assets(project_dir, board):
    files = bundle(project_dir, "device", board)
    files["app.js"] += b"\n/* ---- board picture ---- */\n" + device_picture_script(project_dir, board)
    return files


def render_header(files, board, env_name=""):
    lines = [
        "// AUTO-GENERATED by scripts/pack_webapp.py (scripts/webapp_build.py) - do not edit.",
        f"// Web app for board: {board}" + (f"  (environment {env_name})" if env_name else ""),
        "#pragma once",
        "#include <stdint.h>",
        "#include <stddef.h>",
        "",
    ]
    total = 0
    for url, variable, source in DEVICE_ASSETS:
        data = _gzip(files[source])
        total += len(data)
        lines.append(f"// {url}: {len(files[source])} bytes, {len(data)} gzipped")
        lines.append(f"static const uint8_t {variable}_gz[] = {{")
        for offset in range(0, len(data), 24):
            lines.append("    " + ",".join(f"0x{b:02X}" for b in data[offset:offset + 24]) + ",")
        lines.append("};")
        lines.append(f"static const size_t {variable}_gz_len = {len(data)};")
        lines.append("")
    lines.append(f"// Total gzipped web assets: {total} bytes")
    lines.append("")
    return "\n".join(lines), total


def write_device_header(project_dir, board, header_path, env_name=""):
    text, total = render_header(device_assets(project_dir, board), board, env_name)
    changed = build_shared_js.write_if_changed(header_path, text)
    return total, changed


def _write_bytes_if_changed(path, data):
    """Writes only when the content changes (a repository folder keeps its history tidy)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        with open(path, "rb") as f:
            if f.read() == data:
                return False
    with open(path, "wb") as f:
        f.write(data)
    return True


def _write_app(folder, files, pictures, prefix=""):
    """Writes the app (page, style, script, board pictures) into one folder.
    Returns (changed, removed), the names prefixed for the report."""
    changed = []
    for name, data in files.items():
        if _write_bytes_if_changed(os.path.join(folder, name), data):
            changed.append(prefix + name)
    for name, data in pictures.items():
        path = os.path.join(folder, "boards", "pics", *name.split("/"))
        if _write_bytes_if_changed(path, data):
            changed.append(prefix + "boards/pics/" + name)

    # Pictures of boards that no longer exist (and the former .svg copies).
    removed = []
    pics_dir = os.path.join(folder, "boards", "pics")
    if os.path.isdir(pics_dir):
        for name in sorted(os.listdir(pics_dir)):
            if name.lower().endswith((".js", ".svg")) and name not in pictures:
                os.remove(os.path.join(pics_dir, name))
                removed.append(prefix + "boards/pics/" + name)
    return changed, removed


def build_site(project_dir, out_dir):
    """Writes the app of the firmware being built into out_dir: the page, the style, the
    script and the board pictures, and nothing else. It is one version of the app, the one
    that goes with this firmware, and it is what gets published under v/<version>/ of the
    site. Only the files of the app are touched: anything else in the folder (.git, README,
    LICENSE, CNAME of a site repository) is left alone.
    Returns {'files': ..., 'version': ..., 'changed': [names], 'removed': [names]}."""
    files = bundle(project_dir, "site")
    pictures = build_board_pics.all_picture_scripts(project_dir)
    version = read_version(project_dir)

    out_dir = os.path.abspath(out_dir)
    if os.path.abspath(os.path.join(project_dir, "webapp")) == out_dir:
        raise RuntimeError("the site folder cannot be the webapp source folder")

    changed, removed = _write_app(out_dir, files, pictures)
    return {"files": files, "version": version, "changed": changed, "removed": removed}


def launcher_bundle(project_dir):
    """The home page of the published site: {'index.html': bytes, 'launcher.js': bytes}.

    The scripts it lists are joined the same way the app's are, and they are the app's own
    transport and protocol: the page reads the first answer of the lightgun with the same
    code, so there is one implementation and not two."""
    webapp_dir = os.path.join(project_dir, "webapp")
    launcher_dir = os.path.join(webapp_dir, "launcher")
    html = _read_text(os.path.join(launcher_dir, "index.html"))
    block = SCRIPTS_BLOCK.search(html)
    if not block:
        raise RuntimeError("launcher/index.html: missing <!-- OF:SCRIPTS --> ... <!-- /OF:SCRIPTS --> block")
    scripts = SCRIPT_SRC.findall(re.sub(r'<!--.*?-->', '', block.group(1), flags=re.DOTALL))

    # The page lists them as it would need them unbundled, that is relative to itself.
    generated = {"../" + GENERATED_SHARED: build_shared_js.render_shared_js(
        build_shared_js.extract_shared(project_dir))}

    parts = ["/* OpenFIRE Web App home - generated by scripts/webapp_build.py, do not edit. */\n"]
    for src in scripts:
        if src in generated:
            code = generated[src]
        else:
            path = os.path.normpath(os.path.join(launcher_dir, *src.split("/")))
            if not os.path.exists(path):
                raise RuntimeError(f"launcher/index.html lists {src}, but {path} does not exist")
            code = _read_text(path)
        parts.append(f"\n/* ---- {src} ---- */\n{code.rstrip()}\n;\n")

    code = "".join(parts)

    # The home page reads the version with Protocol.getBoardInfo(): built against a
    # protocol.js that does not have it, it would look perfectly fine and never send the
    # dock request at all. Better a build that stops here than a page that fails silently.
    # Each marker must be written in one source only, or the check cannot fail: launcher.js
    # also names WebSerialTransport and the generated board data also names serialCmdTypes_e,
    # so those bare names are always found whatever transport.js or protocol.js is joined.
    for needed in ("_getBoardInfo", "class WebSerialTransport", '"serialCmdTypes_e":'):
        if needed not in code:
            raise RuntimeError(f"launcher: the joined scripts do not contain {needed}. "
                               f"webapp/js/core/protocol.js, transport.js and the data "
                               f"generated from src/boards/OpenFIREshared.h must be the "
                               f"ones of this source tree - an old copy breaks the home page.")

    page = html[:block.start()] + '<script src="launcher.js"></script>' + html[block.end():]
    return {"index.html": page.encode("utf-8"), "launcher.js": code.encode("utf-8")}


def build_launcher(project_dir, out_dir):
    """Writes the home page into out_dir, which is the root of the published site:
    .nojekyll goes here too, because that is where GitHub Pages reads it."""
    files = launcher_bundle(project_dir)
    changed = []
    for name, data in files.items():
        if _write_bytes_if_changed(os.path.join(out_dir, name), data):
            changed.append(name)

    nojekyll = os.path.join(out_dir, ".nojekyll")  # GitHub Pages: serve the files as they are
    if not os.path.exists(nojekyll):
        os.makedirs(out_dir, exist_ok=True)
        with open(nojekyll, "wb"):
            pass
        changed.append(".nojekyll")

    return {"files": files, "changed": changed}


def print_sizes(project_dir):
    shared = build_shared_js.extract_shared(project_dir)
    boards = [b for b in shared["boardsBoxPositions"] if "esp32" in b]
    print(f"{'board':28s} {'html':>6s} {'css':>6s} {'js':>7s} {'picture':>7s} {'total':>8s}  (gzipped bytes; picture inside js)")
    for board in boards:
        files = device_assets(project_dir, board)
        sizes = [len(_gzip(files[source])) for _, _, source in DEVICE_ASSETS]
        picture = len(_gzip(device_picture_script(project_dir, board)))
        print(f"{board:28s} {sizes[0]:6d} {sizes[1]:6d} {sizes[2]:7d} {picture:7d} {sum(sizes):8d}")
    site = bundle(project_dir, "site")
    pictures = sum(len(_gzip(d)) for d in build_board_pics.all_picture_scripts(project_dir).values())
    print(f"{'site (all boards)':28s} {len(_gzip(site['index.html'])):6d} {len(_gzip(site['style.css'])):6d} "
          f"{len(_gzip(site['app.js'])):7d} {pictures:7d}")


def main(argv):
    parser = argparse.ArgumentParser(description="Build the OpenFIRE web app")
    parser.add_argument("target", choices=("site", "launcher", "device", "sizes"))
    parser.add_argument("--project", default=project_dir_default(), help="lightgun folder")
    parser.add_argument("--out", help="output folder (site: <project>/dist/site; "
                                      "launcher: <project>/dist/launcher; device: plain files)")
    parser.add_argument("--board", help="device target board name (e.g. waveshare-esp32-s3-zero)")
    parser.add_argument("--header", help="device header path (default: <project>/include/web_assets.h)")
    args = parser.parse_args(argv)
    project = os.path.abspath(args.project)

    refresh_dev_files(project)
    if args.target == "site":
        out = os.path.abspath(args.out or os.path.join(project, "dist", "site"))
        result = build_site(project, out)
        print(f"App of version {result['version']['label']} written to {out}: "
              f"{len(result['changed'])} file(s) updated"
              + (f", {len(result['removed'])} removed" if result["removed"] else "")
              + f" (publish it as v/{result['version']['id']}/ of the site)")
    elif args.target == "launcher":
        out = os.path.abspath(args.out or os.path.join(project, "dist", "launcher"))
        result = build_launcher(project, out)
        print(f"Home page written to {out}: {len(result['changed'])} file(s) updated")
    elif args.target == "device":
        if not args.board:
            parser.error("device needs --board")
        if args.out:  # plain files, e.g. for tests/sim-server.js
            out = os.path.abspath(args.out)
            os.makedirs(out, exist_ok=True)
            for name, data in device_assets(project, args.board).items():
                with open(os.path.join(out, name), "wb") as f:
                    f.write(data)
            print("Device files written to", out)
        else:
            header = os.path.abspath(args.header or os.path.join(project, "include", "web_assets.h"))
            total, changed = write_device_header(project, args.board, header)
            print(f"{header}: {total} bytes gzipped{'' if changed else ' (unchanged)'}")
    else:
        print_sizes(project)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
