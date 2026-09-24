/*  OpenFIRE Web App - the home of the published site.

    One button. It asks the lightgun which firmware it runs - the version is the first
    field of the first answer, so nothing else has to be read - looks that version up in
    versions.json and opens the app published for it, at v/<version>/.

    Nothing is opened by itself when that version is not published: the versions that are
    there are offered instead, so a firmware nobody made an app for can still be tried
    with a neighbouring one (which will say that the versions do not match).

    The lightgun is therefore docked twice: once here, for the version alone, and again by
    the app that opens. The port is handed over in sessionStorage, so the app takes it back
    without asking and the user only sees the page change.

    versions.json sits next to this page and is written when a version is published:
        { "latest": "6.2",
          "versions": [ { "id": "6.2", "label": "6.2.0", "type": "stable" } ] }

    ADDING A LANGUAGE: copy the 'en' block of STRINGS below, change the code and translate.
    .../?lang=it opens the page in that language, and the language is passed on to the app.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const VERSIONS_URL = 'versions.json';
    const SVG_NS = 'http://www.w3.org/2000/svg';
    const THEMES = ['system', 'light', 'dark'];
    const ICONS = {
        system: 'M3.5 5h17a1 1 0 0 1 1 1v9.5a1 1 0 0 1-1 1h-17a1 1 0 0 1-1-1V6a1 1 0 0 1 1-1zM8.5 20h7M12 16.5V20',
        light:  'M12 8a4 4 0 1 0 0 8a4 4 0 1 0 0-8zM12 2.5v2M12 19.5v2M2.5 12h2M19.5 12h2M5.3 5.3l1.4 1.4M17.3 17.3l1.4 1.4M18.7 5.3l-1.4 1.4M6.7 17.3l-1.4 1.4',
        dark:   'M20.8 13.4A8.5 8.5 0 0 1 10.6 3.2a8.5 8.5 0 1 0 10.2 10.2z',
        check:  'M4.5 12.5 9.5 17.5 19.5 6.5',
        caret:  'M6 9l6 6 6-6',
        app:    'M3.5 5h17a1 1 0 0 1 1 1v12a1 1 0 0 1-1 1h-17a1 1 0 0 1-1-1V6a1 1 0 0 1 1-1zM2.5 9h19M5.5 7h0.01M8 7h0.01',
        arrow:  'M5 12h14M13 6l6 6-6 6'
    };
    const JUMP_KEY = 'of_version_jump';     // handed over to the app that is being opened
    const THEME_KEY = 'of_theme';           // the App's own settings, shared with this page
    const LANG_KEY = 'of_lang';

    const STRINGS = {
        en: {
            name: 'English',
            lead: 'Connect the lightgun: this page opens the version of the app that belongs to the firmware it runs.',
            compatBefore: 'Works only with lightguns running firmware 7.x or newer. For the 6.x series use the ',
            toolsLink: 'desktop App',
            compatAfter: '.',
            sixBefore: 'If the lightgun runs firmware of the 6.x series it cannot answer this page: that series needs the ',
            sixAfter: '.',
            connect: 'Connect a Lightgun',
            reading: 'Reading the lightgun...',
            opening: 'Firmware %1: opening the app of that version...',
            noBrowser: 'This browser cannot talk to the lightgun: use Chrome, Edge or Opera on a computer. From here you can still open a version by hand.',
            noAnswer: 'The lightgun did not answer. Check that it is plugged in, and that no other page or program is using it, then try again.',
            busy: 'That port is already in use by another page or another program (the desktop App, a terminal, the Arduino IDE): close it and try again.',
            cannotOpen: 'The port could not be opened. Unplug the lightgun and plug it in again, then try once more.',
            outdated: 'This page is not complete: it was published with files that do not go together. Reload it (Ctrl+F5); if it keeps happening the site has to be published again.',
            theme: 'Theme', theme_system: 'System theme', theme_light: 'Light theme', theme_dark: 'Dark theme',
            language: 'Language',
            noList: 'The list of published versions could not be read. Try again in a moment.',
            notPublished: 'The lightgun runs firmware %1, and no app was published for that version. You can try one of these: it will tell you that the versions do not match.',
            allTitle: 'Published versions',
            allLead: 'Every version of the app stays published. Normally you do not choose: the button above opens the right one.',
            pickOne: 'Choose one to try:',
            showAll: 'See the published versions',
            open: 'Open this version',
            latest: 'latest',
            hub: 'OpenFIRE ESP32: project home page',
            footer: 'OpenFIRE ESP32 - free software, GNU General Public License.'
        },
        it: {
            name: 'Italiano',
            lead: 'Collega la lightgun: questa pagina apre la versione dell’app che corrisponde al firmware che ha dentro.',
            compatBefore: 'Funziona solo con lightgun che hanno il firmware 7.x o superiore. Per la serie 6.x serve l’',
            toolsLink: 'App per computer',
            compatAfter: '.',
            sixBefore: 'Se la lightgun ha un firmware della serie 6.x non può rispondere a questa pagina: per quella serie serve l’',
            sixAfter: '.',
            connect: 'Collega una lightgun',
            reading: 'Leggo la lightgun...',
            opening: 'Firmware %1: apro l’app di quella versione...',
            noBrowser: 'Questo browser non può parlare con la lightgun: usa Chrome, Edge o Opera su un computer. Da qui puoi comunque aprire una versione a mano.',
            noAnswer: 'La lightgun non ha risposto. Controlla che sia collegata e che nessun’altra pagina o programma la stia usando, poi riprova.',
            busy: 'Quella porta è già usata da un’altra pagina o da un altro programma (l’App per computer, un terminale, l’IDE di Arduino): chiudilo e riprova.',
            cannotOpen: 'Non sono riuscito ad aprire la porta. Stacca e riattacca la lightgun, poi riprova.',
            outdated: 'Questa pagina non è completa: è stata pubblicata con file che non stanno insieme. Ricaricala (Ctrl+F5); se continua, il sito va ripubblicato.',
            theme: 'Tema', theme_system: 'Tema di sistema', theme_light: 'Tema chiaro', theme_dark: 'Tema scuro',
            language: 'Lingua',
            noList: 'Non sono riuscito a leggere l’elenco delle versioni pubblicate. Riprova fra un momento.',
            notPublished: 'La lightgun ha il firmware %1, e per quella versione non risulta pubblicata nessuna app. Puoi provarne una di queste: ti dirà che le versioni non coincidono.',
            allTitle: 'Versioni pubblicate',
            allLead: 'Ogni versione dell’app resta pubblicata. Di solito non devi scegliere: il pulsante qui sopra apre quella giusta.',
            pickOne: 'Scegline una da provare:',
            showAll: 'Vedi le versioni pubblicate',
            open: 'Apri questa versione',
            latest: 'ultima',
            hub: 'OpenFIRE ESP32: pagina iniziale del progetto',
            footer: 'OpenFIRE ESP32 - software libero, licenza GNU General Public License.'
        }
    };

    /* Il marchio in alto a sinistra riporta alla pagina iniziale del progetto,
       portandosi dietro la lingua in uso, come fa quella pagina quando manda qui. */
    const HUB_URL = 'https://alessandro-satanassi.github.io/OpenFIRE-ESP32/';
    /* La App per computer, per chi ha ancora un firmware della serie 6.x: il suo
       protocollo e' diverso e questa pagina non riesce nemmeno a leggerlo. */
    const TOOLS_URL = 'https://alessandro-satanassi.github.io/OpenFIRE-ESP32-Tools/';

    let lang = 'en';
    const t = (key, ...args) => {
        let text = (STRINGS[lang] && STRINGS[lang][key]) || STRINGS.en[key] || '';
        args.forEach((value, index) => { text = text.split('%' + (index + 1)).join(String(value)); });
        return text;
    };

    const byId = (id) => root.document.getElementById(id);
    const store = {
        session(key, value) { try { root.sessionStorage.setItem(key, value); } catch (e) { /* private mode */ } },
        local(key) { try { return root.localStorage.getItem(key); } catch (e) { return null; } },
        setLocal(key, value) { try { root.localStorage.setItem(key, value); } catch (e) { /* private mode */ } }
    };

    let theme = 'system';
    let shownList = null;      // the list of versions, while it is on screen

    function icon(name, extra) {
        const svg = root.document.createElementNS(SVG_NS, 'svg');
        svg.setAttribute('class', 'icon' + (extra ? ' ' + extra : ''));
        svg.setAttribute('viewBox', '0 0 24 24');
        svg.setAttribute('aria-hidden', 'true');
        const path = root.document.createElementNS(SVG_NS, 'path');
        path.setAttribute('d', ICONS[name] || '');
        svg.appendChild(path);
        return svg;
    }

    function element(tag, cls, text) {
        const node = root.document.createElement(tag);
        if (cls) node.className = cls;
        if (text !== undefined) node.textContent = text;
        return node;
    }

    function pickLanguage() {
        const codes = Object.keys(STRINGS);
        let asked = null;
        try { asked = new URLSearchParams(root.location.search).get('lang'); } catch (e) { asked = null; }
        if (asked && codes.indexOf(asked) >= 0) return asked;
        const saved = store.local(LANG_KEY);
        if (saved && codes.indexOf(saved) >= 0) return saved;
        const list = root.navigator.languages || [root.navigator.language || ''];
        for (let i = 0; i < list.length; ++i) {
            const short = String(list[i]).toLowerCase().split('-')[0];
            if (codes.indexOf(short) >= 0) return short;
        }
        return 'en';
    }

    /** Same setting as the App (of_theme): choosing it here changes it there too. */
    function applyTheme() {
        if (theme === 'light' || theme === 'dark') root.document.documentElement.dataset.theme = theme;
        else delete root.document.documentElement.dataset.theme;

        const button = byId('theme-button');
        if (button) {
            const path = button.querySelector('.icon path');
            if (path) path.setAttribute('d', ICONS[theme]);
            const label = t('theme') + ': ' + t('theme_' + theme);
            button.setAttribute('aria-label', label);
            button.title = label;
        }
        const items = root.document.querySelectorAll('#theme-menu button');
        for (let i = 0; i < items.length; ++i)
            items[i].setAttribute('aria-checked', String(items[i].dataset.theme === theme));
    }

    /** The theme menu and the language selector, as on the other pages of the project. */
    function buildControls() {
        const bar = root.document.querySelector('.topbar .controls');
        if (!bar) return;
        bar.textContent = '';

        const wrap = element('span', 'theme-wrap');
        const button = element('button', 'control theme-button');
        button.type = 'button';
        button.id = 'theme-button';
        button.setAttribute('aria-haspopup', 'true');
        button.setAttribute('aria-expanded', 'false');
        button.appendChild(icon('system'));
        button.appendChild(icon('caret', 'caret'));

        const menu = element('div', 'theme-menu');
        menu.id = 'theme-menu';
        menu.setAttribute('role', 'menu');
        THEMES.forEach((name) => {
            const entry = root.document.createElement('button');
            entry.type = 'button';
            entry.setAttribute('role', 'menuitemradio');
            entry.dataset.theme = name;
            entry.appendChild(icon(name));
            entry.appendChild(element('span', null, t('theme_' + name)));
            entry.appendChild(icon('check', 'tick'));
            entry.addEventListener('click', () => {
                theme = name;
                store.setLocal(THEME_KEY, name);
                applyTheme();
                menu.dataset.open = 'false';
                button.setAttribute('aria-expanded', 'false');
            });
            menu.appendChild(entry);
        });
        button.addEventListener('click', (event) => {
            event.stopPropagation();
            const open = menu.dataset.open !== 'true';
            menu.dataset.open = String(open);
            button.setAttribute('aria-expanded', String(open));
        });
        wrap.appendChild(button);
        wrap.appendChild(menu);
        bar.appendChild(wrap);

        const select = root.document.createElement('select');
        select.className = 'control';
        select.id = 'lang-select';
        Object.keys(STRINGS).sort((a, b) => STRINGS[a].name.localeCompare(STRINGS[b].name))
            .forEach((code) => {
                const option = root.document.createElement('option');
                option.value = code;
                option.textContent = STRINGS[code].name;
                select.appendChild(option);
            });
        select.value = lang;
        select.setAttribute('aria-label', t('language'));
        select.title = t('language');
        select.addEventListener('change', () => {
            lang = select.value;
            store.setLocal(LANG_KEY, lang);   // chosen here, so the App opens in it too
            render();
        });
        bar.appendChild(select);
    }

    const suffix = () => '?lang=' + encodeURIComponent(lang);

    /** La riga sotto al pulsante. `extra`, se c'e', aggiunge una frase che contiene
        un collegamento: { before, link, href, after }. Niente innerHTML: il testo
        tradotto resta testo, il collegamento e' un nodo costruito qui. */
    function say(text, bad, extra) {
        const node = byId('state');
        node.textContent = text || '';
        node.classList.toggle('bad', !!bad);
        if (!extra) return;
        node.appendChild(root.document.createTextNode(' ' + extra.before));
        const link = element('a', null, extra.link);
        link.href = extra.href;
        node.appendChild(link);
        node.appendChild(root.document.createTextNode(extra.after || ''));
    }

    /** La frase sulla serie 6.x, con il rimando alla pagina degli strumenti. */
    const sixSeries = () => ({
        before: t('sixBefore'), link: t('toolsLink'),
        href: TOOLS_URL + suffix(), after: t('sixAfter')
    });

    /** The versions published beside this page. */
    async function published() {
        try {
            const answer = await fetch(VERSIONS_URL, { cache: 'no-cache' });
            if (!answer.ok) return null;
            const list = await answer.json();
            return list && Array.isArray(list.versions) ? list : null;
        } catch (error) {
            return null;
        }
    }

    function showVersions(list, lead) {
        shownList = { list, lead };
        const box = byId('versions');
        const items = byId('versions-list');
        byId('versions-title').textContent = t('allTitle');
        byId('versions-lead').textContent = lead || t('allLead');
        items.textContent = '';
        list.versions.forEach((entry) => {
            // Un riquadro come quelli della pagina iniziale del progetto: icona,
            // numero di versione, l'eventuale etichetta e la freccia in fondo.
            const link = root.document.createElement('a');
            link.className = 'card';
            link.href = 'v/' + encodeURIComponent(entry.id) + '/' + suffix();

            const top = element('div', 'card-top');
            const mark = element('span', 'card-icon');
            mark.appendChild(icon('app'));
            top.appendChild(mark);
            // 7.0.0-beta1 dice gia' di che versione si tratta: il tipo si aggiunge solo
            // quando l'etichetta non ha un suffisso, cioe' "7.0.0 stable".
            const label = String(entry.label || entry.id);
            top.appendChild(element('h2', null,
                label + (entry.type && label.indexOf('-') < 0 ? ' ' + entry.type : '')));
            if (String(list.latest || '') === String(entry.id))
                top.appendChild(element('span', 'tag', t('latest')));
            link.appendChild(top);

            const go = element('span', 'card-go', t('open'));
            go.appendChild(icon('arrow'));
            link.appendChild(go);

            items.appendChild(link);
        });
        box.classList.add('on');
        byId('show-versions').parentNode.style.display = 'none';
    }

    /** Docks only to read who the board is, then undocks and lets the port go.
        Returns { board } or { error }: which of the ways it can fail matters, because
        each one asks something different of whoever is reading. */
    async function readVersion(port) {
        let probe;
        try {
            probe = new OF.Protocol();
        } catch (error) {
            return { error: 'outdated' };        // the shared data did not come with the page
        }
        // Built against an old protocol.js this page would look perfectly fine and never
        // send the dock request at all: say so instead of blaming the lightgun.
        if (typeof probe.getBoardInfo !== 'function')
            return { error: 'outdated' };

        try {
            await probe.connect(new OF.WebSerialTransport(port));
        } catch (error) {
            return { error: (error && error.code === 'port_busy') ? 'busy' : 'cannotOpen' };
        }
        let board = null;
        try {
            const info = await probe.getBoardInfo();
            board = info.board || null;   // the version is there even when the rest is not
        } catch (error) {
            console.error('[Launcher] reading the board failed:', error);
        }
        try { await probe.disconnect(); } catch (error) { /* the port goes anyway */ }
        return board && String(board.version || '').trim() ? { board } : { error: 'noAnswer' };
    }

    /** The published app for this firmware, from the most precise name to the least:

          7.0.0-beta1   the complete version, the one the app is published under;
          7.0.0         the same without the suffix, so a beta that has no app of its
                        own opens the one of its three numbers;
          7.0           the old short version, which is all a firmware before 7.0 sends.

        Falling back opens a neighbouring app, not a wrong one silently: that app compares
        the versions as soon as the lightgun docks and says they do not match. */
    function entryFor(list, board) {
        const short = String((board && board.version) || '').split('-')[0].trim();
        const full = String((board && board.versionFull) || '').trim();
        const numbers = full.split('-')[0].trim();
        const find = (id) => id && list.versions.find((item) => item && String(item.id) === id);
        return find(full) || find(numbers) || find(short) || null;
    }

    async function connect() {
        const button = byId('connect');
        if (!OF.WebSerialTransport.isSupported()) {
            say(t('noBrowser'), true);
            const list = await published();
            if (list) showVersions(list);
            return;
        }

        let port = null;
        try {
            port = await OF.WebSerialTransport.requestPort();
        } catch (error) {
            port = null;
        }
        if (!port) return;                       // the chooser was closed: nothing to say

        button.disabled = true;
        say(t('reading'));
        try {
            const answer = await readVersion(port);
            if (answer.error) {
                // Chi ha ancora la serie 6.x arriva qui: la lightgun c'e' e il cavo pure,
                // ma parla un protocollo che questa pagina non sa leggere.
                say(t(answer.error), true, answer.error === 'noAnswer' ? sixSeries() : null);
                return;
            }
            const board = answer.board;
            const list = await published();
            if (!list) {
                say(t('noList'), true);
                return;
            }
            const entry = entryFor(list, board);
            const label = String((board.versionFull || '').trim() || board.version);
            if (!entry) {
                say(t('notPublished', label), true);
                showVersions(list, t('pickOne'));
                return;
            }

            // What the app that opens needs: the port to take back, and what to say.
            store.session(JUMP_KEY, JSON.stringify({
                to: String(entry.id),
                latest: String(list.latest || ''),
                latestLabel: String((list.versions.find((item) => String(item.id) === String(list.latest)) || {}).label || list.latest || ''),
                port: OF.WebSerialTransport.describePort(port)
            }));
            say(t('opening', label));
            root.location.href = 'v/' + encodeURIComponent(entry.id) + '/' + suffix();
        } finally {
            button.disabled = false;
        }
    }

    /** Every text of the page, in the language in use: called again when it changes. */
    function render() {
        root.document.documentElement.lang = lang;
        byId('lead').textContent = t('lead');
        const compat = byId('compat');
        if (compat) {
            compat.textContent = t('compatBefore');
            const link = element('a', null, t('toolsLink'));
            link.href = TOOLS_URL + suffix();
            compat.appendChild(link);
            compat.appendChild(root.document.createTextNode(t('compatAfter')));
        }
        byId('connect-label').textContent = t('connect');
        byId('show-versions').textContent = t('showAll');
        byId('footer').textContent = t('footer');

        const brand = byId('brand-link');
        if (brand) {
            brand.href = HUB_URL + suffix();
            brand.title = t('hub');
            brand.setAttribute('aria-label', t('hub'));
        }

        buildControls();
        applyTheme();
        if (shownList) showVersions(shownList.list, shownList.lead);

        // The address says which language is shown: it can be copied and passed on.
        try {
            const url = new URL(root.location.href);
            url.searchParams.set('lang', lang);
            root.history.replaceState(null, '', url);
        } catch (error) { /* file:// or a browser without history */ }
    }

    async function start() {
        theme = store.local(THEME_KEY);
        if (THEMES.indexOf(theme) < 0) theme = 'system';
        lang = pickLanguage();
        render();

        byId('connect').addEventListener('click', () => { connect(); });
        byId('show-versions').addEventListener('click', async () => {
            const list = await published();
            if (list) showVersions(list);
            else say(t('noList'), true);
        });
        // A click anywhere, or Escape, closes the theme menu.
        root.document.addEventListener('click', () => {
            const menu = byId('theme-menu');
            const button = byId('theme-button');
            if (menu) menu.dataset.open = 'false';
            if (button) button.setAttribute('aria-expanded', 'false');
        });
        root.document.addEventListener('keydown', (event) => {
            if (event.key !== 'Escape') return;
            const menu = byId('theme-menu');
            const button = byId('theme-button');
            if (menu) menu.dataset.open = 'false';
            if (button) button.setAttribute('aria-expanded', 'false');
        });
    }

    if (root.document.readyState === 'loading')
        root.document.addEventListener('DOMContentLoaded', start);
    else start();
})(typeof globalThis !== 'undefined' ? globalThis : this);
