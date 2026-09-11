"""Translations of the web app: webapp/lang/<code>.json is the main source.

Each file maps the English text of the app to its translation. Entries identical to
the English text or empty are dropped from the build, the web app shows the English text.

    load_translations(webapp_dir)   -> {code: {source: translation}}
    render_translations_js(data)    -> "OF.TRANSLATIONS = {...};"
    build_lang_js(webapp_dir)       -> writes webapp/lang/translations.js (used when the
                                       unbundled webapp folder is served during development)
"""
import json
import os
import re
import sys

LANG_CODE = re.compile(r'^[a-z]{2,3}(?:-[A-Za-z0-9]{2,8})?$')


def load_translations(webapp_dir):
    lang_dir = os.path.join(webapp_dir, "lang")
    translations = {}
    if not os.path.isdir(lang_dir):
        return translations
    for name in sorted(os.listdir(lang_dir)):
        code, extension = os.path.splitext(name)
        if extension.lower() != ".json" or not LANG_CODE.match(code):
            continue
        try:
            with open(os.path.join(lang_dir, name), "r", encoding="utf-8-sig") as f:
                entries = json.load(f)
        except ValueError as error:
            raise ValueError(f"webapp/lang/{name}: {error}") from None
        translations[code] = {source: text for source, text in entries.items()
                              if isinstance(text, str) and text.strip() and text != source}
    translations.setdefault("en", {})
    return translations


def render_translations_js(translations):
    body = json.dumps(translations, ensure_ascii=False, separators=(',', ':'), sort_keys=True)
    return ("// AUTO-GENERATED from webapp/lang/*.json - do not edit.\n"
            "(globalThis.OF = globalThis.OF || {}).TRANSLATIONS = " + body + ";\n")


def build_lang_js(webapp_dir):
    path = os.path.join(webapp_dir, "lang", "translations.js")
    text = render_translations_js(load_translations(webapp_dir))
    old = None
    if os.path.exists(path):
        with open(path, "r", encoding="utf-8", newline="") as f:
            old = f.read()
    if old != text:
        with open(path, "w", encoding="utf-8", newline="") as f:
            f.write(text)
    return path


if __name__ == "__main__":
    webapp = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "webapp"))
    print("Written", build_lang_js(webapp))
