/*  OpenFIRE Web App - translations.

    Keys are the English strings of the app (webapp/lang/<code>.json, bundled as
    OF.TRANSLATIONS). A missing translation shows the English text.

        OF.i18n.t('Save and Send Settings')
        OF.i18n.t('GPIO Pin No. %1.', 12)
        OF.i18n.setLanguage('it')
        OF.i18n.on('change', code => ...)

    Elements with data-i18n are translated by translateDOM(): the original content
    (text, or HTML when the element contains markup) is the key. data-i18n-title and
    data-i18n-placeholder translate those attributes.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};
    const STORAGE_KEY = 'of_lang';

    function storageGet(key) {
        try { return root.localStorage ? root.localStorage.getItem(key) : null; } catch (e) { return null; }
    }

    function storageSet(key, value) {
        try { if (root.localStorage) root.localStorage.setItem(key, value); } catch (e) { /* private mode */ }
    }

    class I18n {
        constructor(translations) {
            this.translations = translations || {};
            if (!this.translations.en) this.translations.en = {};
            this.listeners = new Set();
            this.currentLang = this.pickLanguage();
        }

        pickLanguage() {
            const candidates = [storageGet(STORAGE_KEY)];
            const nav = root.navigator;
            if (nav) candidates.push(...(nav.languages || []), nav.language);
            for (const candidate of candidates) {
                if (!candidate) continue;
                const code = String(candidate).toLowerCase();
                if (this.translations[code]) return code;
                const base = code.split('-')[0];
                if (this.translations[base]) return base;
            }
            return 'en';
        }

        /** Available language codes, English first. */
        languages() {
            return Object.keys(this.translations).sort((a, b) => (a === 'en' ? -1 : b === 'en' ? 1 : a.localeCompare(b)));
        }

        /** Name of a language in that language ("Italiano"). */
        languageName(code) {
            try {
                const name = new Intl.DisplayNames([code], { type: 'language' }).of(code);
                if (name && name !== code) return name.charAt(0).toLocaleUpperCase(code) + name.slice(1);
            } catch (e) { /* old browser */ }
            return code.toUpperCase();
        }

        /** Translation of an English string; %1, %2... are replaced by the extra arguments. */
        t(text, ...args) {
            if (text === undefined || text === null) return '';
            const table = this.translations[this.currentLang];
            let out = (table && table[text]) || String(text);
            if (args.length)
                out = out.replace(/%(\d+)/g, (match, n) => (args[n - 1] !== undefined ? String(args[n - 1]) : match));
            return out;
        }

        setLanguage(code) {
            if (!this.translations[code] || code === this.currentLang) return;
            this.currentLang = code;
            storageSet(STORAGE_KEY, code);
            if (root.document) {
                root.document.documentElement.lang = code;
                this.translateDOM();
            }
            for (const listener of this.listeners) listener(code);
        }

        on(name, listener) {
            if (name === 'change') this.listeners.add(listener);
        }

        off(name, listener) {
            if (name === 'change') this.listeners.delete(listener);
        }

        translateDOM(scope) {
            const base = scope || root.document;
            if (!base) return;
            for (const el of base.querySelectorAll('[data-i18n]')) {
                if (el.dataset.i18nKey === undefined) {
                    const html = el.children.length > 0;
                    el.dataset.i18nKey = (html ? el.innerHTML : el.textContent).trim().replace(/\s+/g, ' ');
                    el.dataset.i18nHtml = html ? '1' : '';
                }
                const value = this.t(el.dataset.i18nKey);
                if (el.dataset.i18nHtml) el.innerHTML = value;
                else el.textContent = value;
            }
            for (const attribute of ['title', 'placeholder']) {
                for (const el of base.querySelectorAll(`[data-i18n-${attribute}]`)) {
                    const keyName = `i18n${attribute.charAt(0).toUpperCase()}${attribute.slice(1)}Key`;
                    if (el.dataset[keyName] === undefined) el.dataset[keyName] = el.getAttribute(attribute) || '';
                    el.setAttribute(attribute, this.t(el.dataset[keyName]));
                }
            }
        }

        /** Fills a <select> with the available languages and keeps it in sync. */
        bindSelector(select) {
            if (!select) return;
            select.textContent = '';
            for (const code of this.languages()) {
                const option = root.document.createElement('option');
                option.value = code;
                option.textContent = this.languageName(code);
                select.appendChild(option);
            }
            select.value = this.currentLang;
            select.addEventListener('change', () => this.setLanguage(select.value));
            this.on('change', (code) => { select.value = code; });
        }
    }

    OF.I18n = I18n;
    OF.i18n = new I18n(OF.TRANSLATIONS);
    if (root.document) root.document.documentElement.lang = OF.i18n.currentLang;
    // Name used by the current app.js.
    root.i18n = OF.i18n;
})(typeof globalThis !== 'undefined' ? globalThis : this);
