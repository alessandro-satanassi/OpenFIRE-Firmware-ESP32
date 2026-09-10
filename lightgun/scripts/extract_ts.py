import os
import xml.etree.ElementTree as ET
import json

def extract_ts_to_json(ts_file_path, json_file_path):
    if not os.path.exists(ts_file_path):
        print(f"Skipping {ts_file_path}, file not found.")
        return

    tree = ET.parse(ts_file_path)
    root = tree.getroot()

    translations = {}

    for message in root.findall('.//message'):
        source_el = message.find('source')
        translation_el = message.find('translation')

        if source_el is not None and translation_el is not None:
            source_text = source_el.text
            translation_text = translation_el.text

            # If translation is empty, use the source text (good for English)
            if translation_text is None or translation_text.strip() == "":
                translation_text = source_text

            if source_text:
                translations[source_text] = translation_text

    with open(json_file_path, 'w', encoding='utf-8') as f:
        json.dump(translations, f, ensure_ascii=False, indent=4)
        
    print(f"Extracted {len(translations)} strings to {json_file_path}")

app_dir = "F:/OpenFIREapp/translation"
webapp_lang_dir = "F:/OpenFIREFirmware/lightgun/webapp/lang"

if not os.path.exists(webapp_lang_dir):
    os.makedirs(webapp_lang_dir)

extract_ts_to_json(os.path.join(app_dir, "AppTranslations_it_IT.ts"), os.path.join(webapp_lang_dir, "it.json"))
extract_ts_to_json(os.path.join(app_dir, "AppTranslations_en_US.ts"), os.path.join(webapp_lang_dir, "en.json"))
