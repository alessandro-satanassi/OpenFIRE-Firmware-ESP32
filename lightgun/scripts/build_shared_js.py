"""Converts src/boards/OpenFIREshared.h (board data shared with the firmware) into JavaScript data.

    extract_shared(project_dir)      -> dict with enums, string maps and board maps
    reduce_for_board(data, board)    -> copy that keeps only one board (firmware web app)
    render_shared_js(data)           -> "const OpenFIREshared = {...};" source text
    generate_shared_js(project, web) -> writes webapp/boards/OpenFIREshared.js (all boards)

Board pictures are not embedded here: see build_board_pics.py.
"""
import json
import os
import re
import sys

HEADER_RELATIVE = os.path.join("src", "boards", "OpenFIREshared.h")
PICTURES_RELATIVE = os.path.join("src", "boards", "boardPics")
QRC_RELATIVE = os.path.join(PICTURES_RELATIVE, "boards.qrc")  # alias (board name) -> picture file

# Maps whose keys are board names (reduced to a single board for the firmware web app).
BOARD_KEYED_MAPS = ("boardsPresetsMap", "boardNames", "boardsBoxPositions", "boardsAltPresets", "boardImagesMap")


def _strip_comments(text):
    text = re.sub(r'/\*.*?\*/', '', text, flags=re.DOTALL)
    return re.sub(r'//.*', '', text)


def _c_string(literal):
    """Body of a C string literal (without quotes) -> Python str."""
    return re.sub(r'\\(.)', lambda m: {'n': '\n', 't': '\t'}.get(m.group(1), m.group(1)), literal)


def _eval(expression, names, problems=None, where=""):
    expression = expression.strip()
    try:
        return eval(expression, {"__builtins__": {}}, names)
    except Exception:
        if problems is not None:
            problems.append(f"OpenFIREshared.h: {where}: '{expression}' not understood")
        return expression


def _warn(message):
    print(f"[WebApp] WARNING: src/boards: {message}")


def extract_shared(project_dir, report=False):
    with open(os.path.join(project_dir, HEADER_RELATIVE), "r", encoding="utf-8") as f:
        content = f.read()

    data = {}
    names = {}
    problems = []

    # --- enum { ... } name; ---
    for match in re.finditer(r'enum\s*\{([^}]+)\}\s*([A-Za-z0-9_]+)\s*;', content):
        values = {}
        current = 0
        for item in (i.strip() for i in _strip_comments(match.group(1)).split(',')):
            if not item:
                continue
            if '=' in item:
                key, expression = (part.strip() for part in item.split('=', 1))
                current = _eval(expression, names)
            else:
                key = item
            values[key] = current
            names[key] = current
            if isinstance(current, int):
                current += 1
        data[match.group(2)] = values

    # --- const char* boardArchs[N] = {...}; (map keys such as boardArchs[boardRP]) ---
    archs = re.search(r'const\s+char\s*\*\s*boardArchs\s*\[\s*\w*\s*\]\s*=\s*\{([^}]*)\}\s*;', content)
    if archs:
        data["boardArchs"] = re.findall(r'"([^"]*)"', _strip_comments(archs.group(1)))
        for index_name, index in data.get("boardArchs_e", {}).items():
            if isinstance(index, int) and 0 <= index < len(data["boardArchs"]):
                names[f"boardArchs[{index_name}]"] = data["boardArchs"][index]

    # --- std::map / std::unordered_map (board maps, *_Strings maps, capabilities) ---
    map_pattern = re.compile(r'const\s+std::(?:unordered_)?map[\s\S]*?\s+([A-Za-z0-9_]+)\s*=\s*\{([\s\S]*?)\};')
    entry_pattern = re.compile(r'\{\s*(?:"([^"]+)"|([A-Za-z0-9_\[\]]+))\s*,\s*([^}]+?)\s*\}')
    for match in map_pattern.finditer(content):
        map_name = match.group(1)
        if (not map_name.startswith('board') and not map_name.endswith('_Strings')
                and map_name not in ('pinCapabilitiesMap', 'mcuCapableMaps')):
            continue
        values = {}
        for entry in entry_pattern.finditer(_strip_comments(match.group(2))):
            key = entry.group(1) or entry.group(2)
            if entry.group(2) and key not in names:
                problems.append(f"OpenFIREshared.h: {map_name}: key '{key}' not understood")
            key = names.get(key, key)
            block = entry.group(3).strip()
            where = f"{map_name}['{key}']"
            if block.startswith('{'):
                values[key] = [_eval(v, names, problems, where) for v in block[1:].split(',') if v.strip()]
            elif block.startswith('"') and block.endswith('"'):
                values[key] = _c_string(block[1:-1])
            else:
                values[key] = _eval(block, names, problems, where)
        data[map_name] = values

    # --- static const <integer type> NAME = value; (e.g. TEMPERATURE_SENSOR_ERROR_VALUE) ---
    scalar_pattern = re.compile(r'static\s+const(?:expr)?\s+(?:unsigned\s+|signed\s+)?(?:int|long|short|char|u?int\d+_t)\s+'
                                r'([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^;]+);')
    for match in scalar_pattern.finditer(_strip_comments(content)):
        value = _eval(re.sub(r'\b(\d+|0[xX][0-9a-fA-F]+)[uUlL]+\b', r'\1', match.group(2).strip()), names, problems, match.group(1))
        if isinstance(value, int):
            data[match.group(1)] = value
            names[match.group(1)] = value

    # --- std::unordered_multimap boardsAltPresets (one board can have several presets) ---
    presets = {}
    alt = re.search(r'unordered_multimap[\s\S]*?\s+boardsAltPresets\s*=\s*\{([\s\S]*?)\n\s*\};', content)
    if alt:
        alt_entry = re.compile(r'\{\s*"([^"]+)"\s*,\s*\{\s*"((?:[^"\\]|\\.)*)"\s*,\s*\{([^}]*)\}\s*\}\s*\}')
        for entry in alt_entry.finditer(_strip_comments(alt.group(1))):
            where = f"boardsAltPresets['{entry.group(1)}'] '{entry.group(2)}'"
            pins = [_eval(v, names, problems, where) for v in entry.group(3).split(',') if v.strip()]
            presets.setdefault(entry.group(1), []).append({"name": _c_string(entry.group(2)), "pin": pins})
    data["boardsAltPresets"] = presets

    # --- Board pictures: src/boards/boardPics/boards.qrc (alias = board name -> file) ---
    data["boardImagesMap"] = read_qrc(project_dir, problems)

    _check(project_dir, data, problems)
    if report:
        for problem in problems:
            _warn(problem)
    return data


def read_qrc(project_dir, problems=None):
    """{alias: file path relative to boardPics, with '/'} from src/boards/boardPics/boards.qrc."""
    qrc_path = os.path.join(project_dir, QRC_RELATIVE)
    if not os.path.exists(qrc_path):
        if problems is not None:
            problems.append("src/boards/boardPics/boards.qrc not found: no board pictures")
        return {}
    with open(qrc_path, "r", encoding="utf-8-sig") as f:
        content = re.sub(r'<!--.*?-->', '', f.read(), flags=re.DOTALL)
    images = {}
    for entry in re.finditer(r'<file\b([^>]*)>\s*([^<]+?)\s*</file>', content):
        alias = re.search(r'\balias\s*=\s*"([^"]+)"', entry.group(1))
        if alias:
            images[alias.group(1)] = entry.group(2).replace("\\", "/")
    return images


def _check(project_dir, data, problems):
    """Board data the web app needs for every board."""
    images = data.get("boardImagesMap", {})
    pictures_dir = os.path.join(project_dir, PICTURES_RELATIVE)
    missing = set()
    for alias, file in list(images.items()):
        if not os.path.exists(os.path.join(pictures_dir, *file.split("/"))):
            problems.append(f"boards.qrc: file '{file}' of '{alias}' does not exist (the generic picture is used)")
            missing.add(alias)
            del images[alias]
    for board, positions in data.get("boardsBoxPositions", {}).items():
        if board not in images and board not in missing:
            problems.append(f"board '{board}' has no picture alias in boards.qrc (the generic picture is used)")
        if board != "generic" and board not in data.get("boardNames", {}):
            problems.append(f"OpenFIREshared.h: board '{board}' is missing from boardNames")
        presets = data.get("boardsPresetsMap", {}).get(board)
        if board != "generic" and presets is None:
            problems.append(f"OpenFIREshared.h: board '{board}' is missing from boardsPresetsMap")
        elif presets is not None and len(presets) != len(positions):
            problems.append(f"OpenFIREshared.h: board '{board}': boardsPresetsMap has {len(presets)} pins, boardsBoxPositions {len(positions)}")
    for board in data.get("boardsPresetsMap", {}):
        if board not in data.get("boardsBoxPositions", {}):
            problems.append(f"OpenFIREshared.h: board '{board}' is missing from boardsBoxPositions")


def detect_board(project_dir, defines):
    """Board name selected by the OPENFIRE_BOARD #ifdef chain of the header for these defines."""
    with open(os.path.join(project_dir, HEADER_RELATIVE), "r", encoding="utf-8") as f:
        content = f.read()
    content = re.sub(r'/\*.*?\*/', '', content, flags=re.DOTALL)
    content = re.sub(r'//[^\n]*', '', content)
    # #ifdef X / #elifdef X / #if defined(X) / #elif defined X, followed by #define OPENFIRE_BOARD "name"
    chain = re.findall(r'#\s*(?:el)?if(?:def\s+|\s+defined\s*\(?\s*)(\w+)\s*\)?\s*\n\s*#\s*define\s+OPENFIRE_BOARD\s+"([^"]+)"', content)
    for define, board in chain:
        if define in defines:
            return board
    fallback = re.search(r'#\s*else\s*\n\s*#\s*define\s+OPENFIRE_BOARD\s+"([^"]+)"', content)
    return fallback.group(1) if fallback else "generic"


def board_arch(data, board):
    """Architecture whose name is part of the board name (e.g. "esp32-s3"), else the first one (RP2040/235X).
    Also works for architectures added later to boardArchs."""
    archs = data.get("boardArchs", ["rp2040_235X", "esp32-s3"])
    matches = [arch for arch in archs[1:] if arch and arch in (board or "")]
    return max(matches, key=len) if matches else archs[0]


def reduce_for_board(data, board):
    reduced = dict(data)
    # An unknown board uses the generic maps.
    keep = {board} if board in data.get("boardsBoxPositions", {}) else {board, "generic"}
    for map_name in BOARD_KEYED_MAPS:
        source = data.get(map_name, {})
        reduced[map_name] = {key: value for key, value in source.items() if key in keep}
    capable = data.get("mcuCapableMaps", {})
    reduced["mcuCapableMaps"] = {key: value for key, value in capable.items()
                                 if key == board or key == board_arch(data, board)}
    return reduced


def render_shared_js(data):
    body = json.dumps(data, ensure_ascii=False, separators=(',', ':'))
    return "// AUTO-GENERATED from src/boards/OpenFIREshared.h - do not edit.\nconst OpenFIREshared = " + body + ";\n"


def write_if_changed(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8", newline="") as f:
            if f.read() == text:
                return False
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write(text)
    return True


def generate_shared_js(project_dir, webapp_dir, report=True):
    path = os.path.join(webapp_dir, "boards", "OpenFIREshared.js")
    write_if_changed(path, render_shared_js(extract_shared(project_dir, report)))
    return path


if __name__ == "__main__":
    project = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
    print("Written", generate_shared_js(project, os.path.join(project, "webapp")))
