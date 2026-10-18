/*  OpenFIRE Web App - "Board Layout" tab (Qt: pinsTab).

    One box per GPIO placed around the board picture as described by boardsBoxPositions
    (left: [box][label], right: [label][box], middle: label above box). Pointing at a box
    shows the circle drawn over that pin in the board picture.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const GPIO_TOOLTIP = 'GPIO Pin No. %1.\n\n' +
        'Blue pin numbers are members of I2C0.\n' +
        'Orange are members of I2C1.\n' +
        'Purple pin numbers can automatically select any I2C channel in software.\n' +
        'Gray cannot use I2C devices.';

    const INVERT_TOOLTIP = 'Inverted logic: tick it for a normally closed (NC) contact, read as released while closed.';

    function build(app) {
        const { el, t, select, checkbox, menu } = OF.UI;
        const S = app.state.S;
        const P = S.boardBoxPositions_e;
        const state = app.state;

        const left = el('div', { class: 'pins-col pins-left' });
        const right = el('div', { class: 'pins-col pins-right' });
        const middle = el('div', { class: 'pins-middle' });
        const picture = el('div', { class: 'board-picture', attrs: { role: 'img', 'aria-label': t('Board Layout') } });
        const layout = el('div', { class: 'pins-layout' },
            left, el('div', { class: 'pins-center' }, picture, middle), right);

        const custom = checkbox(t('Use Custom Pins'), false, (on) => {
            report(state.setCustomPins(on));
            app.refresh();
        }, { title: t('Enable to use custom pins mapping, or disable to use the default mappings for the current board.') });

        // Unlocks the NC checkboxes of the pins, against changes made by mistake. Never saved:
        // it starts off at every load and goes off when custom pins or Show Unsafe Settings are turned off.
        let invertEditing = false;
        const invertEdit = checkbox(t('Edit NC Inputs'), false, (on) => {
            invertEditing = on;
            update();
        }, { title: t('Enable to change the NC (normally closed) checkboxes of the pins. This box is not saved: it only protects the inversions from changes made by mistake.') });

        const layoutsMenu = menu(t('User Layouts'), () => [
            { label: t('Import Custom Layout...'), action: importLayout },
            { label: t('Export Custom Layout...'), action: exportLayout, disabled: !state.customPins },
        ], { class: 'menu-button tool-button' });
        layoutsMenu.classList.add('up');

        const presets = el('select', { class: 'presets-box', attrs: { 'aria-label': t('Select Layout Preset') } });
        presets.addEventListener('change', () => {
            const index = Number(presets.value);
            if (index >= 0) {
                report(state.applyAltPreset(index));
                app.refresh();
            }
        });

        const scroll = el('div', { class: 'tab-scroll' }, layout);
        const rootNode = el('div', { class: 'tab-body pins-tab' },
            scroll,
            // Qt PinsBottomHalf: checkbox on the left, User Layouts and the presets box on the right.
            el('div', { class: 'pins-bottom' }, el('div', { class: 'pins-bottom-left' }, custom, invertEdit),
                el('div', { class: 'pins-bottom-right' }, layoutsMenu, presets)));

        let boxes = [];
        let inverts = [];
        let loadedBoard = null;

        // The board layout is shown whole: it shrinks instead of being scrolled
        // (on narrow screens the single column is scrolled as usual).
        const wide = root.matchMedia ? root.matchMedia('(min-width: 761px)') : null;
        const fit = () => OF.UI.fitToHeight(scroll, layout, 0.6, !wide || wide.matches);
        OF.UI.fitOnResize(scroll, fit);
        if (wide && wide.addEventListener) wide.addEventListener('change', fit);

        function report(messages) {
            if (messages && messages.length) app.status.show(t(messages[messages.length - 1]), 10000);
        }

        function highlight(gpio, on) {
            OF.Boards.highlightPin(picture, gpio, on);
        }

        /** (Re)creates the boxes for the loaded board. */
        function populate() {
            left.textContent = '';
            right.textContent = '';
            middle.textContent = '';
            boxes = [];
            inverts = [];
            loadedBoard = state.board ? state.board.type : null;
            if (!state.loaded) return;

            const names = OF.Maps.functionNames(S);
            const positions = state.boxPositions;
            const placed = [];
            let maxLeft = 0;
            let maxRight = 0;

            for (let gpio = 0; gpio < positions.length; ++gpio) {
                const position = positions[gpio];
                const side = position & P.posCheck;
                if (side === P.posNothing) { boxes.push(null); inverts.push(null); continue; }
                const slot = position ^ side;

                const items = names.map((name, index) => ({
                    value: index,
                    label: t(name),
                    disabled: index > 0 && state.functionDisabled(gpio, index - 1),
                }));
                const box = select(items, 0, (value) => {
                    report(state.changePin(gpio, Number(value)));
                    app.refresh();
                }, { class: 'pin-box', attrs: { 'aria-label': `GPIO${gpio}` } });
                box.addEventListener('mouseenter', () => highlight(gpio, true));
                box.addEventListener('mouseleave', () => highlight(gpio, false));
                box.addEventListener('focus', () => highlight(gpio, true));
                box.addEventListener('blur', () => highlight(gpio, false));

                const label = el('span', {
                    class: 'gpio-label ' + state.pinI2CClass(gpio),
                    text: `«GPIO${gpio}»`,
                    title: t(GPIO_TOOLTIP, gpio),
                });
                // Shown only for the inputs that can be inverted (see update()).
                const invert = checkbox(t('NC'), false, (on) => {
                    const fn = state.functionAt(gpio);
                    if (state.invertible(fn)) state.setInverted(fn, on);
                    app.refresh();
                }, { class: 'check pin-invert', title: t(INVERT_TOOLTIP) });
                boxes.push(box);
                inverts.push(invert);

                if (side === P.posLeft) {
                    placed.push([left, slot, el('div', { class: 'pin-item', style: { '--row': slot } }, invert, box, label)]);
                    maxLeft = Math.max(maxLeft, slot);
                } else if (side === P.posRight) {
                    placed.push([right, slot, el('div', { class: 'pin-item', style: { '--row': slot } }, label, box, invert)]);
                    maxRight = Math.max(maxRight, slot);
                } else if (side === P.posMiddle) {
                    placed.push([middle, slot, el('div', { class: 'pin-item vertical', style: { order: slot } }, label, box, invert)]);
                }
            }
            // Row of boardsBoxPositions (Qt grid row, empty rows keep their space); page order follows it too.
            placed.sort((a, b) => a[1] - b[1]).forEach(([column, , item]) => column.append(item));
            left.style.setProperty('--rows', Math.max(1, maxLeft));
            right.style.setProperty('--rows', Math.max(1, maxRight));
            middle.hidden = !middle.children.length;

            OF.Boards.showPicture(picture, state.board.type);

            presets.textContent = '';
            presets.append(el('option', { value: '-1', text: t('Select Layout Preset'), disabled: true, hidden: true }));
            state.altPresets.forEach((preset, index) => presets.append(el('option', { value: String(index), text: preset.name })));
            presets.hidden = !state.altPresets.length;
        }

        function update() {
            if (!state.loaded) return;
            if (state.board.type !== loadedBoard || !boxes.length) { populate(); fit(); }
            const customOn = state.customPins;
            if (!customOn || !app.showUnsafe) invertEditing = false;
            boxes.forEach((box, gpio) => {
                if (!box) return;
                const fn = state.functionAt(gpio);
                box.value = String(fn + 1);
                box.disabled = !customOn;
                const invert = inverts[gpio];
                const allowed = state.invertible(fn);
                const inverted = allowed && state.inverted(fn);
                // Locked: only the inverted inputs show it, as an orange "NC" tag over the box (style.css).
                invert.classList.toggle('off', !allowed || (!invertEditing && !inverted));
                invert.classList.toggle('badge', !invertEditing && inverted);
                invert.input.disabled = !allowed || !invertEditing;
                invert.input.checked = inverted;
            });
            custom.input.checked = customOn;
            invertEdit.hidden = !customOn || !app.showUnsafe;
            invertEdit.input.checked = invertEditing;
            // Qt: applying a preset changes the pin boxes, which reset the presets box to its placeholder.
            presets.selectedIndex = 0;
        }

        async function importLayout() {
            const session = app.session;
            const file = await OF.UI.openFile('.ofl');
            // The board may have been disconnected (or another one docked) while choosing the file.
            if (session !== app.session || !state.loaded || app.busy) return;
            if (!file) {
                app.status.show(t('Canceled custom layout load operation.'), 5000);
                return;
            }
            if (!file.bytes) {
                await OF.UI.alert(t('File Read Error'), OF.UI.escape(t('Custom layout file could not be read.')), 'warning', { scope: 'session' });
                return;
            }
            let result;
            try {
                result = state.importLayout(file.bytes);
            } catch (error) {
                result = 'corrupt';   // not a layout file (for example too big to read)
            }
            if (result === 'corrupt') {
                await OF.UI.alert(t('File Read Error'), OF.UI.escape(t('Custom layout file could not be read.')), 'warning', { scope: 'session' });
                return;
            }
            if (result !== 'ok') {
                await OF.UI.alert(t("Board Doesn't Match"), OF.UI.escape(t('Custom layout file is not compatible with this board.')), 'warning', { scope: 'session' });
                return;
            }
            app.refresh();
            app.status.show(t('Successfully imported custom layout!'), 5000);
        }

        function exportLayout() {
            try {
                OF.UI.downloadBytes(`${state.board.type}.ofl`, state.exportLayout());
                app.status.show(t('Custom layout export successful!'), 5000);
            } catch (error) {
                OF.UI.alert(t('File Write Error'), OF.UI.escape(t('Custom layout file could not be written.')), 'warning');
            }
        }

        return {
            root: rootNode,
            update,
            reload() { loadedBoard = null; invertEditing = false; update(); },
        };
    }

    OF.Tabs = OF.Tabs || {};
    OF.Tabs.pins = { id: 'pins', label: 'Board Layout', icon: 'pins', build };
})(typeof globalThis !== 'undefined' ? globalThis : this);
