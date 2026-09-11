"""Board pictures for the web app: the SVG files listed in src/boards/boardPics/boards.qrc.

The web app inlines the SVG to highlight the pins (elements with id "OF_pin<gpio>"), so
ids and vector content are kept. Two things make the files lighter:
  - editor metadata and indentation are removed;
  - embedded PNG/JPEG pictures are re-encoded as WebP (needs Pillow).

Each picture becomes a script with the same name (rpipico.svg -> boards/pics/rpipico.js):
    OF.BoardPictures["rpipico.svg"] = "<svg ...>";
Scripts load from any page, also opened as a local file (file://), where a fetch() of the
SVG is refused by the browser.

The scripts are cached in webapp/boards/pics/ with the SHA-1 of their source, so a build
machine without Pillow reuses them. The PlatformIO build installs Pillow by itself when a
picture has to be re-compressed (scripts/pack_webapp.py); by hand:
    python scripts/build_board_pics.py
"""
import base64
import hashlib
import io
import json
import os
import re
import sys

SOURCE_RELATIVE = os.path.join("src", "boards", "boardPics")
CACHE_RELATIVE = os.path.join("webapp", "boards", "pics")
WEBP_QUALITY = 85
MAX_RASTER_SIDE = 1024
_STAMP = re.compile(r'^/\* of-source-sha1:([0-9a-f]{40}) raster:(\w+) \*/\r?\n')
_OLD_SVG_STAMP = re.compile(r'^<!-- of-source-sha1:[0-9a-f]{40} raster:\w+ -->\r?\n')  # former .svg cache

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_shared_js  # noqa: E402

try:
    from PIL import Image
except ImportError:  # PlatformIO Python usually has no Pillow: the cache is used instead.
    Image = None


def source_dir(project_dir):
    return os.path.join(project_dir, SOURCE_RELATIVE)


def picture_files(project_dir):
    """Picture files (relative to boardPics, with '/') referenced by boards.qrc that exist."""
    folder = source_dir(project_dir)
    files = set(build_shared_js.read_qrc(project_dir).values())
    return sorted(f for f in files if os.path.isfile(os.path.join(folder, *f.split("/"))))


def script_name(picture):
    """Script of a picture: same path and name, extension .js (rpipico.svg -> rpipico.js)."""
    return os.path.splitext(picture)[0] + ".js"


def render_script(picture, svg_bytes, stamp=""):
    """JavaScript that registers the optimized SVG text of a picture in OF.BoardPictures."""
    text = svg_bytes.decode("utf-8") if isinstance(svg_bytes, bytes) else svg_bytes
    return (stamp + "((globalThis.OF = globalThis.OF || {}).BoardPictures = globalThis.OF.BoardPictures || {})"
            f"[{json.dumps(picture)}] = {json.dumps(text, ensure_ascii=False)};\n").encode("utf-8")


def _script_value(script, start=0):
    """SVG text registered by a script made by render_script, or None when it cannot be read."""
    try:
        begin = script.index("] = ", start) + 4
        end = script.rstrip().rindex(";")
        return json.loads(script[begin:end])
    except ValueError:
        return None


def _cache_path(project_dir, picture):
    return os.path.join(project_dir, CACHE_RELATIVE, *script_name(picture).split("/"))


def load_pillow():
    """Imports Pillow if it became available (e.g. just installed). Returns True when usable."""
    global Image
    if Image is None:
        try:
            import importlib
            importlib.invalidate_caches()
            from PIL import Image as pil_image
            Image = pil_image
        except ImportError:
            return False
    return True


def _source_digest(source):
    # Line endings do not change the optimized picture (git may convert them on checkout).
    return hashlib.sha1(source.replace(b"\r\n", b"\n")).hexdigest()


def pictures_needing_pillow(project_dir):
    """Board SVGs with embedded photos whose compressed copy is missing or out of date."""
    folder = source_dir(project_dir)
    needed = []
    for name in picture_files(project_dir):
        with open(os.path.join(folder, *name.split("/")), "rb") as f:
            source = f.read()
        if re.search(rb'data:image/(?:png|jpeg|jpg);base64,', source) is None:
            continue
        cache_path = _cache_path(project_dir, name)
        stamp = None
        if os.path.exists(cache_path):
            with open(cache_path, "r", encoding="utf-8", errors="replace") as f:
                stamp = _STAMP.match(f.read(200))
        if not stamp or stamp.group(1) != _source_digest(source) or stamp.group(2) != "done":
            needed.append(name)
    return needed


def _reencode_raster(match):
    kind, payload = match.group(1), match.group(2)
    raw = base64.b64decode(re.sub(r'&#10;|&#13;|\s', '', payload))
    image = Image.open(io.BytesIO(raw))
    image.load()
    if max(image.size) > MAX_RASTER_SIDE:
        scale = MAX_RASTER_SIDE / max(image.size)
        image = image.resize((round(image.size[0] * scale), round(image.size[1] * scale)), Image.LANCZOS)
    if image.mode not in ("RGB", "RGBA"):
        image = image.convert("RGBA" if "A" in image.getbands() or "transparency" in image.info else "RGB")
    out = io.BytesIO()
    image.save(out, "WEBP", quality=WEBP_QUALITY, method=6)
    if len(out.getvalue()) >= len(raw):
        return match.group(0)
    return 'data:image/webp;base64,' + base64.b64encode(out.getvalue()).decode('ascii') + '"'


def _minify(svg):
    svg = re.sub(r'<\?xml[^>]*\?>', '', svg)
    svg = re.sub(r'<!--.*?-->', '', svg, flags=re.DOTALL)
    svg = re.sub(r'<metadata\b.*?</metadata>', '', svg, flags=re.DOTALL)
    svg = re.sub(r'<sodipodi:namedview\b[^>]*?(/>|>.*?</sodipodi:namedview>)', '', svg, flags=re.DOTALL)
    svg = re.sub(r'\s(?:inkscape|sodipodi):[\w-]+="[^"]*"', '', svg)
    svg = re.sub(r'[ \t\r\n]+', ' ', svg)
    return svg.strip()


def optimize_svg(source_bytes):
    """Returns (optimized text, raster_done). raster_done is False when Pillow is missing."""
    svg = source_bytes.decode("utf-8")
    has_raster = re.search(r'data:image/(?:png|jpeg|jpg);base64,', svg) is not None
    if has_raster and Image is not None:
        svg = re.sub(r'data:image/(png|jpeg|jpg);base64,([^"]+)"', _reencode_raster, svg)
    return _minify(svg), (not has_raster) or Image is not None


def board_picture(project_dir, filename, update_cache=True):
    """Optimized SVG bytes for a board picture, or None when the source file is missing."""
    source_path = os.path.join(source_dir(project_dir), *filename.split("/"))
    if not os.path.exists(source_path):
        return None
    with open(source_path, "rb") as f:
        source = f.read()
    digest = _source_digest(source)
    cache_path = _cache_path(project_dir, filename)

    if os.path.exists(cache_path):
        with open(cache_path, "r", encoding="utf-8", newline="") as f:
            cached = f.read()
        stamp = _STAMP.match(cached)
        if stamp and stamp.group(1) == digest and (stamp.group(2) == "done" or Image is None):
            text = _script_value(cached, stamp.end())
            if text is not None:
                return text.encode("utf-8")

    text, raster_done = optimize_svg(source)
    if not raster_done:
        print(f"[WebApp] NOTE: {filename}: the embedded photos are not compressed (Pillow is not available).")
    if update_cache:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        stamp = f"/* of-source-sha1:{digest} raster:{'done' if raster_done else 'raw'} */\n"
        with open(cache_path, "wb") as f:
            f.write(render_script(filename, text, stamp))
    return text.encode("utf-8")


def all_pictures(project_dir, update_cache=True):
    """{file: optimized SVG bytes} for every board picture listed in boards.qrc."""
    names = picture_files(project_dir)
    pictures = {}
    for name in names:
        data = board_picture(project_dir, name, update_cache)
        if data is not None:
            pictures[name] = data
    if update_cache:
        _remove_stale_cache(project_dir, names)
    return pictures


def all_picture_scripts(project_dir, update_cache=True):
    """{script file: JavaScript bytes} for every board picture (boards/pics/<name>.js)."""
    return {script_name(name): render_script(name, data) for name, data in all_pictures(project_dir, update_cache).items()}


def _remove_stale_cache(project_dir, source_names):
    """Deletes generated scripts whose source SVG was removed or renamed, and the former .svg cache."""
    cache_dir = os.path.join(project_dir, CACHE_RELATIVE)
    if not os.path.isdir(cache_dir):
        return
    scripts = {script_name(name) for name in source_names}
    for folder, _, files in os.walk(cache_dir):
        for file in files:
            path = os.path.join(folder, file)
            name = os.path.relpath(path, cache_dir).replace(os.sep, "/")
            lower = name.lower()
            if lower.endswith(".js") and name not in scripts:
                stamp = _STAMP
            elif lower.endswith(".svg"):
                stamp = _OLD_SVG_STAMP
            else:
                continue
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                generated = stamp.match(f.read(200)) is not None
            if generated:
                try:
                    os.remove(path)
                except OSError as error:
                    print(f"[WebApp] NOTE: could not remove the old picture cache {name}: {error}")


if __name__ == "__main__":
    import gzip
    project = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
    if Image is None:
        print("Pillow is not installed (pip install pillow): embedded pictures will not be re-encoded.")
    for name, data in all_pictures(project).items():
        with open(os.path.join(source_dir(project), *name.split("/")), "rb") as f:
            original = f.read()
        print(f"{name:32s} {len(gzip.compress(original, 9)) // 1024:5d} KB gz -> {len(gzip.compress(data, 9)) // 1024:4d} KB gz")
