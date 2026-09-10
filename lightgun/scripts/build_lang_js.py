import os
import json
import glob

def build_lang_js(webapp_dir):
    lang_dir = os.path.join(webapp_dir, "lang")
    lang_js_path = os.path.join(webapp_dir, "lang.js")
    
    if not os.path.exists(lang_dir):
        return

    translations = {}
    
    # Read all .json files in the lang directory
    for file_path in glob.glob(os.path.join(lang_dir, "*.json")):
        lang_code = os.path.splitext(os.path.basename(file_path))[0]
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                translations[lang_code] = json.load(f)
        except Exception as e:
            print(f"Error loading {file_path}: {e}")

    # Generate the JS file content
    js_content = f"// AUTO-GENERATED from lang/*.json\n"
    js_content += f"const TRANSLATIONS = {json.dumps(translations, ensure_ascii=False, indent=4)};\n\n"
    
    js_content += """
class I18n {
    constructor() {
        const savedLang = localStorage.getItem('of_lang');
        const browserLang = (navigator.language || navigator.userLanguage).substring(0, 2).toLowerCase();
        
        let targetLang = savedLang || browserLang;
        this.currentLang = TRANSLATIONS[targetLang] ? targetLang : "en";
    }

    t(text) {
        if (this.currentLang === "en" || !TRANSLATIONS[this.currentLang]) {
            return text; 
        }
        return TRANSLATIONS[this.currentLang][text] || text;
    }

    translateDOM() {
        const elements = document.querySelectorAll('[data-i18n]');
        elements.forEach(el => {
            if (!el.dataset.originalText) {
                el.dataset.originalText = el.innerText.trim();
            }
            const originalText = el.dataset.originalText;
            el.innerText = this.t(originalText);
        });
    }

    setLanguage(langCode) {
        if (langCode === "en" || TRANSLATIONS[langCode]) {
            this.currentLang = langCode;
            localStorage.setItem('of_lang', langCode);
            this.translateDOM();
            
            const selector = document.getElementById('lang-selector');
            if (selector) selector.value = langCode;
        }
    }
}

const i18n = new I18n();

document.addEventListener("DOMContentLoaded", () => {
    i18n.translateDOM();
    
    const selector = document.getElementById('lang-selector');
    if (selector) {
        selector.value = i18n.currentLang;
        selector.addEventListener('change', (e) => {
            i18n.setLanguage(e.target.value);
        });
    }
});
"""

    with open(lang_js_path, "w", encoding="utf-8") as f:
        f.write(js_content)
        
    print(f"[WebApp Packer] Generated {lang_js_path} with {len(translations)} languages.")

if __name__ == "__main__":
    build_lang_js("F:/OpenFIREFirmware/lightgun/webapp")
