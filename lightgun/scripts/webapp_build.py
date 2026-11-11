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
           the ones that are. Under the button it shows the polling rate test of the app,
           with the app's own script, texts and style (see "The polling rate test").

What gets published: every firmware has its own app, and dist/site is exactly that app -
one version of it, not a site. Publishing a version means copying dist/site into
v/<version>/ of the site repository and adding it to versions.json; the home page
(dist/launcher) sits at the root and sends everybody to the right version. The version is
the full one, MAJOR.MINOR.PATCH plus the optional suffix of src/OpenFIREversion.h (7.0.0,
7.0.0-rc2...), see read_version(): the firmware sends it in its full-version marker when
the app docks, and it names the folder (v/7.0.0/) and the release tag (v7.0.0). The old
"%.1f" field of OPENFIRE_VERSION is only checked for consistency; the launcher keeps a
fallback on it for firmware that does not send the full version.

index.html lists the scripts between <!-- OF:SCRIPTS --> and <!-- /OF:SCRIPTS -->: in
the build they are joined, in order, into a single app.js. Two of them are generated:
  boards/OpenFIREshared.js  from src/boards/OpenFIREshared.h (build_shared_js.py)
  lang/translations.js      from webapp/lang/*.json          (build_lang_js.py)
Serving the webapp folder as it is (python -m http.server) runs the unbundled app.

Command line (no PlatformIO needed):
  python scripts/webapp_build.py site     [--out DIR]
  python scripts/webapp_build.py launcher [--out DIR]
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


VERSION_SUFFIX = re.compile(r"[A-Za-z][A-Za-z0-9]*")
VERSION_KIND = re.compile(r"[A-Za-z]+")   # beta1 -> beta, rc2 -> rc


def _string_define(text, name):
    """Body of a #define holding a C string literal ("" gives "")."""
    raw = _define(text, name)
    if raw is None:
        raise RuntimeError(f"src/OpenFIREversion.h: #define {name} is missing")
    raw = raw.strip()
    if len(raw) < 2 or raw[0] != '"' or raw[-1] != '"':
        raise RuntimeError(f"src/OpenFIREversion.h: {name} must be a string in quotes ({raw!r})")
    return raw[1:-1]


def read_version(project_dir):
    """Version of the firmware, from src/OpenFIREversion.h - and the only place that
    reads it. MAJOR, MINOR, PATCH and SUFFIX are the source; everything else here is
    derived from them, so the folder of the published app, the Git tag and the title of
    the release cannot drift apart from what the firmware actually sends.

    id        the complete version, 7.0.0 or 7.0.0-beta1: it names the folder of the
              archived app (v/<id>/) and it is exactly the string the firmware sends
              when it docks (sVersionFull, see src/main.cpp).
    label     the same, shown in the list of versions.
    numbers   7.0.0, without the suffix: what the home page falls back to.
    suffix    beta1, rc2... empty for a final version.
    type      stable, beta, rc, prerelease...: the letters of the suffix, for the badge.
    tag       v7.0.0-beta1, the Git tag of the release.
    prerelease  True as soon as there is a suffix.

    It also refuses a header that contradicts itself: OPENFIRE_VERSION is still the
    first field the lightgun sends and the only one the desktop App reads, so it has to
    keep saying MAJOR.MINOR - and "%.1f" cannot tell 7.10 from 7.1.
    """
    text = _read_text(os.path.join(project_dir, "src", "OpenFIREversion.h"))

    parts = []
    for name in ("OPENFIRE_VERSION_MAJOR", "OPENFIRE_VERSION_MINOR", "OPENFIRE_VERSION_PATCH"):
        value = _define(text, name)
        if value is None or not value.isdigit():
            raise RuntimeError(f"src/OpenFIREversion.h: {name} is missing or not a whole number ({value!r})")
        parts.append(str(int(value)))
    numbers = ".".join(parts)

    suffix = _string_define(text, "OPENFIRE_VERSION_SUFFIX").strip()
    if suffix and not VERSION_SUFFIX.fullmatch(suffix):
        raise RuntimeError(f"src/OpenFIREversion.h: OPENFIRE_VERSION_SUFFIX must be a letter followed "
                           f"by letters and digits, like beta1 or rc2 ({suffix!r}): it becomes a folder "
                           f"name, a web address and a Git tag.")
    full = numbers + ("-" + suffix if suffix else "")

    raw = _define(text, "OPENFIRE_VERSION")
    if raw is None:
        raise RuntimeError("src/OpenFIREversion.h: #define OPENFIRE_VERSION is missing")
    try:
        short = "%.1f" % float(raw)  # exactly what the firmware prints when docking
    except ValueError:
        raise RuntimeError(f"src/OpenFIREversion.h: OPENFIRE_VERSION is not a number ({raw!r})")
    if int(parts[1]) > 9:
        raise RuntimeError(f"src/OpenFIREversion.h: OPENFIRE_VERSION_MINOR is {parts[1]}: the \"%.1f\" "
                           f"field sent to the desktop App cannot tell {parts[0]}.{parts[1]} from "
                           f"{parts[0]}.{parts[1][0]}.")
    if short != f"{parts[0]}.{parts[1]}":
        raise RuntimeError(f"src/OpenFIREversion.h: OPENFIRE_VERSION is {raw} ({short}) but the version "
                           f"is {full}: it must be {parts[0]}.{parts[1]}, because it is the field the "
                           f"desktop App reads.")

    return {"id": full, "label": full, "numbers": numbers, "suffix": suffix,
            "type": VERSION_KIND.match(suffix).group(0).lower() if suffix else "stable",
            "tag": "v" + full, "prerelease": bool(suffix)}


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
    # Only for the unbundled home page: whatever goes wrong here must not stop the build of
    # the firmware (pack_webapp.py calls this before every build). The launcher build itself
    # stops with the same error, so it does not go unnoticed.
    try:
        _write_bytes_if_changed(os.path.join(webapp_dir, *POLLING_DEV_CSS.split("/")),
                                polling_css(webapp_dir).encode("utf-8"))
    except Exception as error:
        print(f"[WebApp] WARNING: webapp/{POLLING_DEV_CSS} not written: {error}")


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


# ----- The polling rate test of the home page ------------------------------------------------
# "Measure Polling Rate" on the home page shows the polling rate test of the app, before the
# lightgun is docked (it must still be in Run mode, sending its aim as a mouse). It is the
# app's own test and nothing of it is copied by hand into the home page:
#   script  js/app/polling-test.js, listed in launcher/index.html and joined like the others;
#   texts   those of webapp/lang/*.json that the test (and the button) use: polling_translations;
#   style   the rules of webapp/style.css that name "polling", with the colours of the app they
#           use, scoped to the test: polling_css. The page links launcher/polling-test.css
#           unbundled (written by refresh_dev_files) and has the same rules inline once built.
# The test takes el, t and icon from OF.UI (js/app/ui.js), which the home page cannot load:
# launcher.js has the same three, and the build checks that the test needs nothing else.
POLLING_SCRIPT = "js/app/polling-test.js"
POLLING_DEV_CSS = "launcher/polling-test.css"
POLLING_STYLE_BLOCK = re.compile(r'<!--\s*OF:POLLING-STYLE\b.*?-->(.*?)<!--\s*/OF:POLLING-STYLE\s*-->', re.DOTALL)
POLLING_UI = ("el", "t", "icon")
POLLING_SCOPE = ".polling-test"
TEXT_KEY = re.compile(r"""\bt\(\s*(?:'((?:[^'\\\n]|\\.)*)'|"((?:[^"\\\n]|\\.)*)")""")


def _css_items(text):
    """The items of a style sheet (comments already removed): [(prelude, body)], body None
    for a statement such as @import."""
    items = []
    i, n = 0, len(text)
    while i < n:
        j = i
        while j < n and text[j] not in "{;}":
            j += 1
        prelude = " ".join(text[i:j].split())
        if j >= n:
            break
        if text[j] != "{":
            if prelude:
                items.append((prelude, None))
            i = j + 1
            continue
        depth, k = 1, j + 1
        while k < n and depth:
            depth += {"{": 1, "}": -1}.get(text[k], 0)
            k += 1
        items.append((prelude, text[j + 1:k - 1]))
        i = k
    return items


def _split_top(text, separator):
    """Splits at the separators that are not inside brackets or parentheses."""
    parts, depth, start = [], 0, 0
    for index, char in enumerate(text):
        if char in "([":
            depth += 1
        elif char in ")]":
            depth -= 1
        elif char == separator and depth == 0:
            parts.append(text[start:index])
            start = index + 1
    parts.append(text[start:])
    return [" ".join(part.split()) for part in parts if part.strip()]


def _css_rules(text, media=()):
    """Every style rule: (media, selectors, [(property, value)], body), or (media, None, raw,
    raw) for an @keyframes block. Comments must be removed first."""
    for prelude, body in _css_items(text):
        if body is None:
            continue
        lower = prelude.lower()
        if lower.startswith(("@media", "@supports")):
            yield from _css_rules(body, media + (prelude,))
        elif lower.startswith(("@keyframes", "@-webkit-keyframes")):
            raw = prelude + " {" + body + "}"
            yield media, None, raw, raw
        elif not lower.startswith("@"):
            declarations = []
            for part in _split_top(body, ";"):
                name, colon, value = part.partition(":")
                if colon:
                    declarations.append((name.strip(), value.strip()))
            yield media, _split_top(prelude, ","), declarations, body


def _render_css_rule(media, selectors, declarations):
    text = ",\n".join(selectors) + " {\n" + "".join(f"    {name}: {value};\n" for name, value in declarations) + "}\n"
    for prelude in reversed(media):
        text = prelude + " {\n" + "".join("    " + line + "\n" for line in text.splitlines()) + "}\n"
    return text


def _without_comments(text):
    return re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)


def polling_css(webapp_dir):
    """The style of the polling rate test for the home page, from webapp/style.css: every
    rule of it that names "polling" (and only the selectors that do), and before them the
    colours of the app those rules, the test and the home page's own rules for it use,
    defined on the test only, so the rest of the home page keeps its own."""
    style = _without_comments(_read_text(os.path.join(webapp_dir, "style.css")))
    script = _read_text(os.path.join(webapp_dir, *POLLING_SCRIPT.split("/")))
    home = _read_text(os.path.join(webapp_dir, "launcher", "index.html"))

    rules, tokens = [], []
    for media, selectors, declarations, body in _css_rules(style):
        if selectors is None:
            if "polling" in declarations.split("{", 1)[0]:
                rules.append(declarations + "\n")
            continue
        chosen = [selector for selector in selectors if "polling" in selector]
        if chosen:
            # Read from the body, not from the declarations: a nested rule can end up in the
            # name of one (".child { color: red") or in none at all (".child { }").
            if "{" in body:
                raise RuntimeError("launcher: a rule of webapp/style.css for " + ", ".join(chosen) +
                                   " has rules nested inside it, which polling_css() does not read.")
            rules.append(_render_css_rule(media, chosen, declarations))
        elif "{" in body and "polling" in body:
            raise RuntimeError("launcher: webapp/style.css has a rule of the polling rate test nested "
                               "inside " + ", ".join(selectors) + ", which polling_css() does not read.")
        elif all(selector.startswith(":root") for selector in selectors):
            tokens.append((media, selectors, declarations))
    if not any(rule.startswith(POLLING_SCOPE + " {") for rule in rules):
        raise RuntimeError(f"launcher: webapp/style.css has no {POLLING_SCOPE} rule: "
                           f"the polling rate test of the home page would have no style.")

    # The colours used: by the rules taken, by the script (cssColor) and by the rules the
    # home page adds for the test (its buttons), and those they are made of.
    own = ""
    for block in re.findall(r"<style>(.*?)</style>", home, flags=re.DOTALL):
        for media, selectors, declarations, _ in _css_rules(_without_comments(block)):
            if selectors and any("polling" in selector for selector in selectors):
                own += " ".join(value for _, value in declarations) + "\n"
    needed = set(re.findall(r"var\(\s*(--[\w-]+)", "".join(rules) + own))
    needed |= set(re.findall(r"""['"](--[\w-]+)['"]""", script))
    defined = {}
    for _, _, declarations in tokens:
        for name, value in declarations:
            defined.setdefault(name, []).append(value)
    while True:
        more = {name for found in needed for value in defined.get(found, [])
                for name in re.findall(r"var\(\s*(--[\w-]+)", value)} - needed
        if not more:
            break
        needed |= more
    missing = sorted(name for name in needed if name not in defined)
    if missing:
        raise RuntimeError("launcher: the polling rate test uses " + ", ".join(missing) +
                           ", which webapp/style.css does not define in :root.")

    colours = [_render_css_rule(media, [selector + " " + POLLING_SCOPE for selector in selectors],
                                [(name, value) for name, value in declarations if name in needed])
               for media, selectors, declarations in tokens
               if any(name in needed for name, _ in declarations)]
    return ("/* AUTO-GENERATED from webapp/style.css by scripts/webapp_build.py - do not edit.\n"
            "   The style of the polling rate test of the app, for the home page (launcher). */\n"
            + "".join(colours) + "".join(rules))


def polling_translations(webapp_dir, sources):
    """The translations of webapp/lang/*.json for the texts that `sources` (the test and the
    home page) pass to t(): the home page carries a few of them, not the whole app."""
    keys = set()
    for code in sources:
        for single, double in TEXT_KEY.findall(code):
            keys.add(re.sub(r"\\(.)", r"\1", single or double))
    return {code: {key: text for key, text in table.items() if key in keys}
            for code, table in build_lang_js.load_translations(webapp_dir).items()}


def _check_polling_helpers(webapp_dir, test_code, home_code):
    """The test runs here with launcher.js's el, t and icon instead of js/app/ui.js."""
    used = set(re.findall(r"OF\.UI\.(\w+)", test_code))
    for names in re.findall(r"\{([^{}]*)\}\s*=\s*OF\.UI\b", test_code):
        used |= {name.split(":")[0].strip() for name in names.split(",") if name.strip()}
    extra = sorted(used - set(POLLING_UI))
    # A text that is not written out (t(name)) would not be found by polling_translations().
    if len(re.findall(r"\bt\(", test_code)) != len(TEXT_KEY.findall(test_code)):
        raise RuntimeError(f"launcher: {POLLING_SCRIPT} calls t() with a text that is not written "
                           f"in quotes: its translation would not reach the home page.")
    if extra:
        raise RuntimeError(f"launcher: {POLLING_SCRIPT} uses OF.UI." + ", OF.UI.".join(extra) +
                           ", which the home page does not have: add it to the helpers of the "
                           "polling rate test in launcher/launcher.js.")
    icons_block = re.search(r"const ICONS = \{(.*?)\n\s*\};", _read_text(os.path.join(webapp_dir, "js", "app", "ui.js")), re.DOTALL)
    ui_icons = dict(re.findall(r"^\s*(\w+):\s*'([^']*)'", icons_block.group(1) if icons_block else "", re.M))
    for name in sorted(set(re.findall(r"\bicon\(\s*'(\w+)'", test_code))):
        if name not in ui_icons or f"{name}: '{ui_icons[name]}'" not in home_code:
            raise RuntimeError(f"launcher: {POLLING_SCRIPT} draws the icon '{name}': "
                               f"launcher/launcher.js (UI_ICONS) must have it as js/app/ui.js draws it.")


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
    translations_src = "../" + GENERATED_TRANSLATIONS
    generated = {"../" + GENERATED_SHARED: build_shared_js.render_shared_js(
        build_shared_js.extract_shared(project_dir)), translations_src: None}

    sources = {}
    for src in scripts:
        if src in generated:
            continue
        path = os.path.normpath(os.path.join(launcher_dir, *src.split("/")))
        if not os.path.exists(path):
            raise RuntimeError(f"launcher/index.html lists {src}, but {path} does not exist")
        sources[src] = _read_text(path)
    polling_src = "../" + POLLING_SCRIPT
    if polling_src in sources:
        # Only the texts of the polling rate test, not the whole app's.
        generated[translations_src] = build_lang_js.render_translations_js(
            polling_translations(webapp_dir, [sources[polling_src], sources.get("launcher.js", "")]))
        _check_polling_helpers(webapp_dir, sources[polling_src], sources.get("launcher.js", ""))

    parts = ["/* OpenFIRE Web App home - generated by scripts/webapp_build.py, do not edit. */\n"]
    for src in scripts:
        code = generated[src] if src in generated else sources[src]
        if code is None:
            raise RuntimeError(f"launcher/index.html lists {src} without {polling_src}")
        parts.append(f"\n/* ---- {src} ---- */\n{code.rstrip()}\n;\n")

    code = "".join(parts)

    # The home page reads the version with Protocol.getBoardInfo(): built against a
    # protocol.js that does not have it, it would look perfectly fine and never send the
    # dock request at all. Better a build that stops here than a page that fails silently.
    # Each marker must be written in one source only, or the check cannot fail: launcher.js
    # also names WebSerialTransport and the generated board data also names serialCmdTypes_e,
    # so those bare names are always found whatever transport.js or protocol.js is joined.
    for needed in ("_getBoardInfo", "class WebSerialTransport", '"serialCmdTypes_e":', "class PollingTest"):
        if needed not in code:
            raise RuntimeError(f"launcher: the joined scripts do not contain {needed}. "
                               f"webapp/js/core/protocol.js, transport.js and the data "
                               f"generated from src/boards/OpenFIREshared.h must be the "
                               f"ones of this source tree - an old copy breaks the home page.")

    page = html[:block.start()] + '<script src="launcher.js"></script>' + html[block.end():]
    # The style of the polling rate test, inline: the page stays two files.
    style = POLLING_STYLE_BLOCK.search(page)
    if not style:
        raise RuntimeError("launcher/index.html: missing <!-- OF:POLLING-STYLE --> ... <!-- /OF:POLLING-STYLE --> block")
    page = page[:style.start()] + "<style>\n" + polling_css(webapp_dir) + "</style>" + page[style.end():]
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
