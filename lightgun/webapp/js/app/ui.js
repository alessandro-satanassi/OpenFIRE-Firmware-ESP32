/*  OpenFIRE Web App - DOM helpers: elements, dialogs, spin boxes, status bar, files.

        const { el, t } = OF.UI;
        el('button', { class: 'primary', text: t('Save'), on: { click: save } });
        await OF.UI.confirm(t('Commit Confirmation'), t('Are these settings okay?'), t('...'));
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};
    const doc = root.document;

    const t = (text, ...args) => (OF.i18n ? OF.i18n.t(text, ...args) : String(text));

    /** el(tag, { class, text, html, title, attrs, dataset, style, on, hidden, disabled, value }, ...children) */
    function el(tag, props, ...children) {
        const node = tag.includes(':') ?
            doc.createElementNS('http://www.w3.org/2000/svg', tag.split(':')[1]) : doc.createElement(tag);
        if (props) {
            for (const [key, value] of Object.entries(props)) {
                if (value === undefined || value === null) continue;
                switch (key) {
                case 'class': node.setAttribute('class', value); break;
                case 'text': node.textContent = value; break;
                case 'html': node.innerHTML = value; break;
                case 'attrs': for (const [name, v] of Object.entries(value)) if (v !== undefined && v !== null && v !== false) node.setAttribute(name, v === true ? '' : v); break;
                case 'dataset': Object.assign(node.dataset, value); break;
                case 'style':
                    // Custom properties ('--row') need setProperty: a plain assignment is ignored.
                    for (const [name, v] of Object.entries(value)) {
                        if (name.startsWith('--')) node.style.setProperty(name, v);
                        else node.style[name] = v;
                    }
                    break;
                case 'on': for (const [name, handler] of Object.entries(value)) node.addEventListener(name, handler); break;
                default: node[key] = value; break;
                }
            }
        }
        for (const child of children.flat()) {
            if (child === null || child === undefined || child === false) continue;
            node.append(child instanceof root.Node ? child : String(child));
        }
        return node;
    }

    // ----- Icons (24x24 stroke paths) --------------------------------------------

    const ICONS = {
        pins: 'M9 3v3M15 3v3M9 18v3M15 18v3M3 9h3M3 15h3M18 9h3M18 15h3M7 6h10a1 1 0 0 1 1 1v10a1 1 0 0 1-1 1H7a1 1 0 0 1-1-1V7a1 1 0 0 1 1-1zM10 10h4v4h-4z',
        buttons: 'M6 8h12a4 4 0 0 1 4 4v1a4 4 0 0 1-7 2.6l-.6-.6H9.6l-.6.6A4 4 0 0 1 2 13v-1a4 4 0 0 1 4-4zM7 11v3M5.5 12.5h3M15.5 12h.01M18 13.5h.01',
        settings: 'M4 6h9M17 6h3M4 12h3M11 12h9M4 18h11M19 18h1M15 4v4M9 10v4M17 16v4',
        profiles: 'M12 3v4M12 17v4M3 12h4M17 12h4M12 7a5 5 0 1 0 0 10a5 5 0 1 0 0-10zM12 11.5v1',
        tests: 'M9 3h6M10 3v6l-5.5 9.5A1.7 1.7 0 0 0 6 21h12a1.7 1.7 0 0 0 1.5-2.5L14 9V3M7.5 15h9',
        save: 'M5 3h11l3 3v13a2 2 0 0 1-2 2H7a2 2 0 0 1-2-2zM8 3v5h7V3M8 21v-6h8v6',
        edit: 'M4 20h4L19 9l-4-4L4 16zM13.5 6.5l4 4',
        close: 'M6 6l12 12M18 6L6 18',
        menu: 'M4 7h16M4 12h16M4 17h16',
        chevron: 'M6 9l6 6l6-6',
        info: 'M12 3a9 9 0 1 0 0 18a9 9 0 1 0 0-18zM12 11v6M12 7.5v.01',
        warning: 'M12 3.5L2.5 20h19zM12 10v4.5M12 17.5v.01',
        error: 'M12 3a9 9 0 1 0 0 18a9 9 0 1 0 0-18zM9 9l6 6M15 9l-6 6',
        question: 'M12 3a9 9 0 1 0 0 18a9 9 0 1 0 0-18zM9.5 9.5a2.5 2.5 0 1 1 3.5 2.3c-.6.3-1 .9-1 1.6v.3M12 17v.01',
        usb: 'M12 3v14M12 3l-2 3h4zM12 17a2 2 0 1 0 0 4a2 2 0 1 0 0-4zM12 13l-5-2.5V8M12 11l5-2V7M6 6h2v2H6zM16 5.5a1 1 0 1 0 2 0a1 1 0 1 0-2 0',
        wifi: 'M2.5 9a14 14 0 0 1 19 0M5.5 12.5a9.5 9.5 0 0 1 13 0M8.5 16a5 5 0 0 1 7 0M12 19.5v.01',
    };

    function icon(name, extraClass) {
        const svg = el('svg:svg', { class: 'icon' + (extraClass ? ' ' + extraClass : ''), attrs: { viewBox: '0 0 24 24', 'aria-hidden': 'true' } });
        svg.append(el('svg:path', { attrs: { d: ICONS[name] || '' } }));
        return svg;
    }

    // ----- Dialogs ---------------------------------------------------------------

    const dialogRoot = () => doc.getElementById('dialog-root') || doc.body;

    /** Open dialogs: { node, scope, cancel() }. */
    const openDialogs = new Set();

    /**
     * Dialog. Resolves to the value of the pressed button (ESC/close: cancelValue).
     * options: { title, text, html, info, infoHtml, icon, buttons: [{ label, value, primary, danger }],
     *            content: Node, cancelValue, wide, onOpen(dialog), onClose(),
     *            scope: 'session' (closed with cancelValue when the board session ends),
     *            modeless: true (the page stays usable, like a Qt window) }
     */
    function dialog(options) {
        return new Promise((resolve) => {
            const buttons = options.buttons || [{ label: t('OK'), value: true, primary: true }];
            const cancelValue = 'cancelValue' in options ? options.cancelValue : false;
            let result = cancelValue;

            const body = el('div', { class: 'dialog-body' });
            if (options.text) body.append(el('p', { class: 'dialog-text', text: options.text }));
            if (options.html) body.append(el('div', { class: 'dialog-text rich', html: options.html }));
            if (options.info) body.append(el('p', { class: 'dialog-info', text: options.info }));
            if (options.infoHtml) body.append(el('div', { class: 'dialog-info rich', html: options.infoHtml }));
            if (options.content) body.append(options.content);

            const footer = el('div', { class: 'dialog-buttons' });
            const node = el('dialog', { class: 'dialog' + (options.wide ? ' wide' : '') + (options.className ? ' ' + options.className : '') },
                el('div', { class: 'dialog-head' },
                    options.icon ? icon(options.icon, 'dialog-icon ' + options.icon) : null,
                    el('h2', { class: 'dialog-title', text: options.title || '' }),
                    el('button', { class: 'icon-button dialog-close', title: t('Close'), attrs: { 'aria-label': t('Close') },
                        on: { click: () => node.close() } }, icon('close'))),
                body, footer);

            let focusButton = null;
            for (const spec of buttons) {
                const button = el('button', {
                    class: spec.primary ? 'primary' : spec.danger ? 'danger' : '',
                    text: spec.label,
                    on: { click: () => { result = spec.value; node.close(); } },
                });
                if (spec.primary || spec.focus) focusButton = button;
                footer.append(button);
            }
            footer.hidden = !buttons.length;

            const entry = { node, scope: options.scope || null, cancel: () => { result = cancelValue; node.close(); } };
            node.addEventListener('close', () => {
                openDialogs.delete(entry);
                node.remove();
                if (options.onClose) options.onClose();
                resolve(result);
            });
            node.addEventListener('cancel', () => { result = cancelValue; });
            openDialogs.add(entry);
            dialogRoot().append(node);
            if (OF.i18n && options.translate) OF.i18n.translateDOM(node);
            if (options.modeless) {
                node.classList.add('modeless');
                node.addEventListener('keydown', (event) => { if (event.key === 'Escape') { event.preventDefault(); entry.cancel(); } });
            }
            if (options.modeless && typeof node.show === 'function') node.show();
            else if (typeof node.showModal === 'function') node.showModal();
            else node.setAttribute('open', '');
            if (options.onOpen) options.onOpen(node, (value) => { result = value; node.close(); });
            else if (focusButton) focusButton.focus();
        });
    }

    /** Closes the open dialogs of a scope (or all) as if cancelled. */
    function closeDialogs(scope) {
        for (const entry of [...openDialogs])
            if (!scope || entry.scope === scope) entry.cancel();
    }

    /** True while a modal dialog is open. */
    const modalOpen = () => [...openDialogs].some((entry) => !entry.node.classList.contains('modeless'));

    const alert = (title, html, iconName = 'info', options = {}) =>
        dialog({ title, html, icon: iconName, buttons: [{ label: t('OK'), value: true, primary: true }], cancelValue: true, scope: options.scope });

    /** Yes/No question. defaultNo focuses "No" (destructive operations). */
    const confirm = (title, text, info, options = {}) => dialog({
        title, text, info, html: options.html, icon: options.icon || 'question', cancelValue: false, scope: options.scope,
        buttons: [
            { label: t('Yes'), value: true, primary: !options.defaultNo, danger: options.danger },
            { label: t('No'), value: false, primary: !!options.defaultNo },
        ],
    });

    /** Text question. options: { value, maxLength, latin1 (refuse characters above U+00FF, like the Qt name box), scope } */
    function prompt(title, label, options = {}) {
        const input = el('input', { type: 'text', value: options.value || '', attrs: { maxlength: options.maxLength, 'aria-label': label, spellcheck: 'false' } });
        if (options.latin1) {
            let accepted = input.value;
            input.addEventListener('input', () => {
                if ([...input.value].some((c) => c.codePointAt(0) > 255)) {
                    input.value = accepted;
                    input.classList.add('invalid');
                } else {
                    accepted = input.value;
                    input.classList.remove('invalid');
                }
            });
        }
        const content = el('label', { class: 'dialog-field' }, el('span', { text: label }), input);
        return dialog({
            title, content, cancelValue: null, scope: options.scope,
            buttons: [{ label: t('OK'), value: 'ok', primary: true }, { label: t('Cancel'), value: null }],
            onOpen: (node, done) => {
                input.focus();
                input.select();
                input.addEventListener('keydown', (event) => { if (event.key === 'Enter') { event.preventDefault(); done('ok'); } });
            },
        }).then((value) => (value === 'ok' ? input.value : null));
    }

    /** C/Qt "%g" with 6 significant digits (QString::number of a float): 1.23457e+06, 0.0001, nan, inf. */
    function formatG(value) {
        const v = Number(value);
        if (Number.isNaN(v)) return 'nan';
        if (!Number.isFinite(v)) return v < 0 ? '-inf' : 'inf';
        if (v === 0) return Object.is(v, -0) ? '-0' : '0';
        const [mantissa, exponentText] = v.toExponential(5).split('e');
        const exponent = Number(exponentText);
        const trim = (text) => (text.includes('.') ? text.replace(/0+$/, '').replace(/\.$/, '') : text);
        if (exponent < -4 || exponent >= 6)
            return `${trim(mantissa)}e${exponent < 0 ? '-' : '+'}${String(Math.abs(exponent)).padStart(2, '0')}`;
        return trim(v.toFixed(5 - exponent));
    }

    const toHex = (rgb) => '#' + ((rgb >>> 0) & 0xFFFFFF).toString(16).padStart(6, '0');
    const fromHex = (hex) => parseInt(String(hex).replace('#', ''), 16) & 0xFFFFFF;

    /** Colour chooser. Resolves to 0xRRGGBB or null. */
    function pickColor(title, rgb, options = {}) {
        let value = toHex(rgb);
        const swatch = el('div', { class: 'color-swatch large', style: { background: value } });
        const hex = el('input', { type: 'text', value, attrs: { maxlength: 7, spellcheck: 'false', 'aria-label': 'Hex' } });
        const picker = el('input', { type: 'color', value, attrs: { 'aria-label': title } });
        const set = (next, source) => {
            if (!/^#[0-9a-fA-F]{6}$/.test(next)) return;
            value = next.toLowerCase();
            swatch.style.background = value;
            if (source !== hex) hex.value = value;
            if (source !== picker) picker.value = value;
        };
        hex.addEventListener('input', () => set(hex.value.startsWith('#') ? hex.value : '#' + hex.value, hex));
        picker.addEventListener('input', () => set(picker.value, picker));
        const presets = el('div', { class: 'color-presets' },
            ['#ff0000', '#ff8000', '#ffff00', '#00ff00', '#00ffff', '#0000ff', '#8000ff', '#ff00ff', '#ffffff', '#000000']
                .map((color) => el('button', { class: 'color-swatch', title: color, style: { background: color },
                    attrs: { type: 'button', 'aria-label': color }, on: { click: () => set(color) } })));
        const content = el('div', { class: 'color-dialog' }, el('div', { class: 'color-row' }, swatch, picker, hex), presets);
        return dialog({
            title, content, cancelValue: null, scope: options.scope,
            buttons: [{ label: t('OK'), value: 'ok', primary: true }, { label: t('Cancel'), value: null }],
        }).then((result) => (result === 'ok' ? fromHex(value) : null));
    }

    // ----- Inputs ----------------------------------------------------------------

    /**
     * Qt QSpinBox: number input with prefix/suffix and range.
     * spin({ min, max, value, prefix, suffix, hex, onChange }) -> { root, input, get value, set value, set prefix, disabled }
     */
    function spin(options) {
        const min = options.min ?? 0;
        const max = options.max ?? 99;
        const hex = !!options.hex;
        const prefix = el('span', { class: 'affix', text: options.prefix || '' });
        const suffix = el('span', { class: 'affix', text: options.suffix || '' });
        const input = el('input', {
            type: hex ? 'text' : 'number',
            class: 'spin-input' + (hex ? ' hex' : ''),
            attrs: hex ? { maxlength: 4, spellcheck: 'false', autocomplete: 'off', inputmode: 'text' } :
                { min, max, step: 1, inputmode: 'numeric' },
        });
        const wrap = el('span', { class: 'spin' }, prefix, input, suffix);
        let current = options.value ?? min;

        const format = (v) => (hex ? v.toString(16).toUpperCase() : String(v));
        const parse = (text) => {
            const n = hex ? parseInt(text, 16) : Number(text);
            return Number.isFinite(n) ? Math.min(max, Math.max(min, Math.trunc(n))) : current;
        };
        const commit = () => {
            const next = input.value.trim() === '' ? current : parse(input.value);
            input.value = format(next);
            if (next !== current) {
                current = next;
                if (options.onChange) options.onChange(next);
            }
        };
        input.value = format(current);
        input.addEventListener('change', commit);
        input.addEventListener('input', () => {
            if (hex) {
                const clean = input.value.replace(/[^0-9a-fA-F]/g, '');
                if (clean !== input.value) input.value = clean;
            }
            if (input.value.trim() === '') return;
            const n = hex ? parseInt(input.value, 16) : Number(input.value);
            if (Number.isFinite(n) && n >= min && n <= max && Math.trunc(n) === n && n !== current) {
                current = n;
                if (options.onChange) options.onChange(n);
            }
        });
        input.addEventListener('blur', () => { input.value = format(current); });
        input.addEventListener('wheel', (event) => { if (doc.activeElement !== input) event.preventDefault(); }, { passive: false });

        return {
            root: wrap,
            input,
            get value() { return current; },
            set value(v) { current = Math.min(max, Math.max(min, v)); if (doc.activeElement !== input) input.value = format(current); },
            set prefix(text) { prefix.textContent = text; },
            set suffix(text) { suffix.textContent = text; },
            set disabled(on) { input.disabled = !!on; wrap.classList.toggle('disabled', !!on); },
        };
    }

    /** <select> from [{ value, label, disabled }]. */
    function select(items, value, onChange, props = {}) {
        const node = el('select', props);
        for (const item of items)
            node.append(el('option', { value: String(item.value), text: item.label, disabled: !!item.disabled }));
        if (value !== undefined && value !== null) node.value = String(value);
        if (onChange) node.addEventListener('change', () => onChange(node.value));
        // Like the Qt App: the mouse wheel does not change a closed combo box by accident.
        node.addEventListener('wheel', (event) => { if (doc.activeElement !== node) event.preventDefault(); }, { passive: false });
        return node;
    }

    function checkbox(label, checked, onChange, props = {}) {
        const input = el('input', { type: 'checkbox', checked: !!checked });
        if (onChange) input.addEventListener('change', () => onChange(input.checked));
        const node = el('label', Object.assign({ class: 'check' }, props), input, el('span', { class: 'check-label', text: label }));
        node.input = input;
        return node;
    }

    function radio(name, label, checked, onChange, props = {}) {
        const input = el('input', { type: 'radio', name, checked: !!checked });
        if (onChange) input.addEventListener('change', () => { if (input.checked) onChange(); });
        const node = el('label', Object.assign({ class: 'check radio' }, props), input, el('span', { class: 'check-label', text: label }));
        node.input = input;
        return node;
    }

    function group(title, ...children) {
        return el('fieldset', { class: 'group' }, title ? el('legend', { text: title }) : null, ...children);
    }

    /** Enables a control, a checkbox/radio label (node.input) or a whole fieldset. */
    function setEnabled(node, on) {
        if (!node) return;
        if (node.input) node.input.disabled = !on;
        else if ('disabled' in node) node.disabled = !on;
        node.classList.toggle('is-disabled', !on);
    }

    /** Plain text (with \n) as HTML. */
    function escape(text) {
        return String(text).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/\n/g, '<br>');
    }

    // ----- Description box (Qt whatsThis panels) ----------------------------------

    /** Qt whatsThis panel. On phones (no pointer to hover with) it is a bar that opens on a tap. */
    class DescriptionBox {
        constructor(defaultHtml) {
            this.defaultHtml = defaultHtml;
            this.title = el('h3', { class: 'desc-title' });
            this.text = el('div', { class: 'desc-text rich' });
            this.toggleLabel = el('span', { class: 'desc-toggle-label' });
            this.toggle = el('button', {
                class: 'desc-toggle', type: 'button', attrs: { 'aria-expanded': 'false' },
                on: { click: () => this.setOpen(!this.root.classList.contains('open')) },
            }, icon('info'), this.toggleLabel, icon('chevron', 'chevron'));
            this.root = el('section', { class: 'desc-box', attrs: { 'aria-live': 'polite' } }, this.toggle, this.title, this.text);
            this.reset();
        }

        setOpen(open) {
            this.root.classList.toggle('open', open);
            this.toggle.setAttribute('aria-expanded', String(open));
        }

        reset() {
            this.title.textContent = '';
            this.toggleLabel.textContent = t('Description');
            this.text.innerHTML = this.defaultHtml;
        }

        show(title, html) {
            this.title.textContent = title || '';
            this.toggleLabel.textContent = title || t('Description');
            this.text.innerHTML = html || '';
        }

        /** Shows the description while the pointer is over (or focus is in) the element. */
        track(node, title, html) {
            const show = () => this.show(typeof title === 'function' ? title() : title, typeof html === 'function' ? html() : html);
            node.addEventListener('mouseenter', show);
            node.addEventListener('focusin', show);
            return node;
        }
    }

    // ----- Status bar ---------------------------------------------------------------

    class StatusBar {
        constructor(textNode, progressNode) {
            this.textNode = textNode;
            this.progressNode = progressNode;
            this._timer = null;
        }

        /** Qt showMessage: ms = 0 keeps the message until the next one. */
        show(text, ms = 0) {
            clearTimeout(this._timer);
            this._expires = ms > 0 ? Date.now() + ms : 0;
            this.textNode.textContent = text || '';
            this.textNode.title = text || '';
            if (ms > 0) this._timer = setTimeout(() => { this.textNode.textContent = ''; this.textNode.title = ''; }, ms);
        }

        /** Takes over the message, its remaining time and the progress of a status bar being replaced. */
        adopt(previous) {
            if (!previous) return;
            clearTimeout(previous._timer);
            const remaining = previous._expires ? previous._expires - Date.now() : 0;
            if (previous._expires && remaining <= 0) this.show('');
            else this.show(previous.textNode.textContent, remaining);
            this.progressNode.hidden = previous.progressNode.hidden;
            this.progressNode.max = previous.progressNode.max;
            this.progressNode.value = previous.progressNode.value;
        }

        clear() { this.show(''); }

        progressRange(range) {
            if (range) {
                this.progressNode.max = range;
                this.progressNode.value = 0;
                this.progressNode.hidden = false;
            } else {
                this.progressNode.hidden = true;
                this.progressNode.value = 0;
            }
        }

        progress(value) {
            if (!this.progressNode.hidden) this.progressNode.value = value;
        }
    }

    // ----- Files ------------------------------------------------------------------------

    function downloadBytes(name, bytes, type = 'application/octet-stream') {
        const url = URL.createObjectURL(new Blob([bytes], { type }));
        const link = el('a', { href: url, download: name, style: { display: 'none' } });
        doc.body.append(link);
        link.click();
        setTimeout(() => { URL.revokeObjectURL(url); link.remove(); }, 1000);
    }

    /** Opens a file chooser. Resolves to { name, bytes } or null when cancelled. */
    function openFile(accept) {
        return new Promise((resolve) => {
            const input = el('input', { type: 'file', accept, style: { display: 'none' } });
            let settled = false;
            const done = (value) => { if (!settled) { settled = true; input.remove(); resolve(value); } };
            input.addEventListener('change', async () => {
                const file = input.files && input.files[0];
                if (!file) return done(null);
                try { done({ name: file.name, bytes: new Uint8Array(await file.arrayBuffer()) }); } catch (e) { done({ name: file.name, bytes: null }); }
            });
            input.addEventListener('cancel', () => done(null));
            doc.body.append(input);
            input.click();
        });
    }

    // ----- Popup menus ----------------------------------------------------------------

    let openMenu = null;

    function closeMenus() {
        if (openMenu) {
            openMenu.panel.hidden = true;
            openMenu.button.setAttribute('aria-expanded', 'false');
            openMenu = null;
        }
    }

    /** Drop-down menu: items = [{ label, action, checked, disabled, separator, href, hidden }] or a function returning them. */
    function menu(label, items, props = {}) {
        const button = el('button', Object.assign({ class: 'menu-button', attrs: { 'aria-haspopup': 'menu', 'aria-expanded': 'false' } }, props),
            el('span', { text: label }), icon('chevron', 'chevron'));
        const panel = el('div', { class: 'menu-panel', attrs: { role: 'menu' }, hidden: true });
        const wrap = el('div', { class: 'menu' }, button, panel);

        const build = () => {
            panel.textContent = '';
            for (const item of (typeof items === 'function' ? items() : items)) {
                if (!item || item.hidden) continue;
                if (item.separator) { panel.append(el('div', { class: 'menu-separator', attrs: { role: 'separator' } })); continue; }
                const entry = el(item.href ? 'a' : 'button', {
                    class: 'menu-item' + (item.checked !== undefined ? ' checkable' : '') + (item.checked ? ' checked' : ''),
                    attrs: { role: item.checked !== undefined ? 'menuitemcheckbox' : 'menuitem', 'aria-checked': item.checked !== undefined ? String(!!item.checked) : null,
                        href: item.href, target: item.href ? '_blank' : null, rel: item.href ? 'noopener' : null },
                    disabled: !!item.disabled,
                    on: { click: (event) => { closeMenus(); if (item.action) { event.preventDefault(); item.action(); } } },
                }, el('span', { class: 'menu-check', text: item.checked ? '✓' : '' }), el('span', { text: item.label }),
                    item.shortcut ? el('span', { class: 'menu-shortcut', text: item.shortcut }) : null);
                panel.append(entry);
            }
        };

        button.addEventListener('click', (event) => {
            event.stopPropagation();
            const wasOpen = openMenu && openMenu.panel === panel;
            closeMenus();
            if (wasOpen) return;
            build();
            panel.hidden = false;
            button.setAttribute('aria-expanded', 'true');
            openMenu = { panel, button };
            const first = panel.querySelector('.menu-item:not([disabled])');
            if (first && event.detail === 0) first.focus();
        });
        panel.addEventListener('keydown', (event) => {
            const entries = [...panel.querySelectorAll('.menu-item:not([disabled])')];
            const index = entries.indexOf(doc.activeElement);
            if (event.key === 'ArrowDown') { event.preventDefault(); (entries[index + 1] || entries[0]).focus(); }
            else if (event.key === 'ArrowUp') { event.preventDefault(); (entries[index - 1] || entries[entries.length - 1]).focus(); }
            else if (event.key === 'Escape') { closeMenus(); button.focus(); }
        });
        return wrap;
    }

    if (doc) {
        doc.addEventListener('click', (event) => {
            if (openMenu && !openMenu.panel.contains(event.target)) closeMenus();
        });
        // Links of the Qt texts (description boxes, dialogs, previewer) open outside the App,
        // like QDesktopServices: the page and its board session stay open.
        doc.addEventListener('click', (event) => {
            if (event.defaultPrevented || event.button !== 0) return;
            const link = event.target.closest && event.target.closest('a[href]');
            if (!link || link.hasAttribute('download') || link.target === '_blank') return;
            let url;
            try { url = new URL(link.getAttribute('href'), root.location.href); } catch (e) { return; }
            if (!/^https?:$/.test(url.protocol)) return;
            event.preventDefault();
            root.open(url.href, '_blank', 'noopener');
        });
        doc.addEventListener('keydown', (event) => { if (event.key === 'Escape') closeMenus(); });
    }

    OF.UI = {
        t, el, icon, dialog, closeDialogs, modalOpen, alert, confirm, prompt, pickColor, toHex, fromHex, formatG,
        spin, select, checkbox, radio, group, setEnabled, escape,
        DescriptionBox, StatusBar, downloadBytes, openFile, menu, closeMenus,
    };
})(typeof globalThis !== 'undefined' ? globalThis : this);
