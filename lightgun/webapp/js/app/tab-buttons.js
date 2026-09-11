/*  OpenFIRE Web App - "Button Mapping" tab (Qt: buttonFuncTab).

    For the 13 digital inputs: output type + output for on-screen, off-screen and gamepad mode
    (the gamepad slot only has gamepad outputs). A row is enabled when its input has a pin.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const SERIAL_LINK = "<a href='https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/wiki/MAMEHOOKER-Documentation#m---mode-commands'><span style=' text-decoration: underline; color:#8ab4f8;'>Serial command</span></a> ";

    const TEXTS = [
        {
            typeName: 'Onscreen Button Output Type for %1',
            typeHelp: '<p>Select the type of button that <i>%1</i> will function as <b>when the gun is pointing at the screen.</b></p>' +
                '<p>Each available input defined in the current <i>Board Layout</i> can be defined as a button press for one of the ' +
                'three available device outputs that the gun presents to the connected device.</p>',
            valueName: 'Onscreen Button Mapping for %1',
            valueHelp: '<p>Select the output that <i>%1</i> will send to the connected device <b>when the gun is pointing at the screen.</b></p>' +
                '<p>Each available input defined in the current <i>Board Layout</i> can be defined as a button press for one of the ' +
                'three available device outputs that the gun presents to the connected device.</p>',
        },
        {
            typeName: 'Offscreen Button Output Type for %1',
            typeHelp: '<p>Select the type of button that <i>%1</i> will function as <b>when the gun is pointing outside of the screen.</b></p>' +
                '<p>Each available input defined in the current <i>Board Layout</i> can be defined as a button press for one of the ' +
                'three available device outputs that the gun presents to the connected device.</p>',
            valueName: 'Offscreen Button Mapping for %1',
            valueHelp: '<p>Select the output that <i>%1</i> will send to the connected device <b>when the gun is pointing outside of the screen.</b></p>' +
                '<p>Each available input defined in the current <i>Board Layout</i> can be defined as a button press for one of the ' +
                'three available device outputs that the gun presents to the connected device.</p>',
        },
        {
            typeName: 'Gamepad Mode Output Type for %1',
            typeHelp: '<p>Select the type of button that <i>%1</i> will function as <b>when the gun set to Gamepad Output Mode.</b></p>' +
                '<p>Only Gamepad-type outputs are available for Gamepad Output Mode, which can be set via ' + SERIAL_LINK +
                '<tt>M0x1</tt>.</p>',
            valueName: 'Gamepad Mode Button Mapping for %1',
            valueHelp: '<p>Select the output that <i>%1</i> will send to the connected device <b>when the gun is set to Gamepad Output Mode.</b></p>' +
                '<p>Only Gamepad buttons are available to be mapped for Gamepad Output Mode, which can be set via ' + SERIAL_LINK +
                '<tt>M0x1</tt>.</p>' +
                '<p><b>NOTE:</b> When connected to the MiSTer FPGA device, these mappings will NOT be reflected, ' +
                'as OpenFIRE has a hard-coded button layout specifically optimized for the MiSTer ecosystem.</p>',
        },
    ];

    const TAB_HELP = '<html><head/><body><p>This tab shows the currently loaded Button Mapping settings.</p><p>Any enabled button input in the current <span style=" font-style:italic;">Board Layout</span> can be assigned to output any button from any of the three Input Devices that OpenFIRE presents to any connected device (a Mouse, Keyboard, and Xbox-like Gamepad), depending on the condition that the input is pressed.</p><p>Hover over an option to view detailed info about it here.</p></body></html>';
    const ASTICK_HELP = '<html><head/><body><p>If an analog stick is enabled in the current <span style=" font-style:italic;">Board Layout,</span> this determines which type of Button Output it uses.</p><p><span style=" font-weight:700;">Gamepad Analog Stick</span> will use either the Left or Right analog stick of the Gamepad device, depending on the circumstances and/or whether <span style=" font-style:italic;">Gamepad Output Mode</span> is set via Serial. The other two settings will translate the stick\'s analog movements into either digital <span style=" font-weight:700;">Gamepad D-Pad</span> or <span style=" font-weight:700;">Keyboard Arrow Key</span> presses.</p><p><span style=" font-weight:700; font-style:italic;">Default:</span><span style=" font-style:italic;"> Gamepad Analog Stick</span></p></body></html>';

    function build(app) {
        const { el, t, select } = OF.UI;
        const state = app.state;
        const S = state.S;
        const E = state.E;
        const M = OF.Maps;
        const names = M.functionNames(S);
        const desc = new OF.UI.DescriptionBox(t(TAB_HELP));

        const slotTitles = ['On-Screen Function', 'Off-Screen Function', 'Gamepad Mode'];
        const header = el('div', { class: 'btn-row btn-header', attrs: { 'aria-hidden': 'true' } },
            el('span'), ...slotTitles.map((title) => el('span', { class: 'btn-slot-title', text: t(title) })));

        const rows = [];
        const list = el('div', { class: 'btn-list' }, header);

        const typeItems = M.inputTypes.map((name, value) => ({ value, label: t(name) }));
        const outputItems = (type) => (M.outputs[type] || []).map((entry) => ({ value: entry.value, label: t(entry.name) }));

        for (let button = 0; button < M.BUTTON_COUNT - 1; ++button) {
            const name = t(names[button + 1]);
            const row = el('fieldset', { class: 'btn-row' }, el('span', { class: 'btn-name', text: name }));
            const slots = [];
            for (let slot = 0; slot < 3; ++slot) {
                const texts = TEXTS[slot];
                const typeBox = select(typeItems, null, (value) => {
                    state.setButtonType(button, slot, Number(value));
                    fillOutputs(slots[slot]);
                    app.refresh();
                }, { class: 'btn-type', disabled: slot >= 2, attrs: { 'aria-label': t(texts.typeName, name) } });
                const valueBox = select([], null, (value) => {
                    state.setButtonValue(button, slot, Number(value));
                    app.refresh();
                }, { class: 'btn-value', attrs: { 'aria-label': t(texts.valueName, name) } });
                desc.track(typeBox, t(texts.typeName, name), t(texts.typeHelp, name));
                desc.track(valueBox, t(texts.valueName, name), t(texts.valueHelp, name));
                const entry = { button, slot, typeBox, valueBox, type: -1 };
                slots.push(entry);
                row.append(el('div', { class: 'btn-slot', dataset: { caption: t(slotTitles[slot]) } }, typeBox, valueBox));
            }
            rows.push({ row, slots });
            list.append(row);
        }

        function fillOutputs(entry) {
            const type = state.buttonType(entry.button, entry.slot);
            if (entry.type !== type) {
                entry.valueBox.textContent = '';
                for (const item of outputItems(type))
                    entry.valueBox.append(el('option', { value: String(item.value), text: item.label }));
                entry.type = type;
            }
            entry.typeBox.value = String(type);
            entry.valueBox.value = String(state.buttonValue(entry.button, entry.slot));
            if (entry.valueBox.selectedIndex < 0 && entry.valueBox.options.length) entry.valueBox.selectedIndex = 0;
        }

        const aStickMode = select(M.analogModes.map((label, value) => ({ value, label: t(label) })), 0, (value) => {
            state.setSetting(state.T.analogMode, Number(value));
            app.refresh();
        }, { attrs: { 'aria-label': t('Analog Stick Output Mode') } });
        desc.track(aStickMode, t('Analog Stick Output Mode'), t(ASTICK_HELP));
        const aStickBox = el('fieldset', { class: 'group' }, el('legend', { text: t('Analog Stick') }),
            el('div', { class: 'row center' }, el('label', { text: t('Send Analog Stick As:') }), aStickMode));

        const rootNode = el('div', { class: 'tab-body' },
            el('div', { class: 'tab-scroll' },
                el('div', { class: 'tab-content buttons' },
                    el('fieldset', { class: 'group' }, el('legend', { text: t('Digital Inputs') }), list),
                    aStickBox)),
            desc.root);

        function update() {
            if (!state.loaded) return;
            for (const { row, slots } of rows) {
                const button = slots[0].button;
                row.disabled = !state.pinMapped(button);
                slots.forEach(fillOutputs);
            }
            aStickBox.disabled = !(state.pinMapped(E.analogX) && state.pinMapped(E.analogY));
            aStickMode.value = String(state.setting(state.T.analogMode));
        }

        return { root: rootNode, update, onShow() { desc.reset(); } };
    }

    OF.Tabs = OF.Tabs || {};
    OF.Tabs.buttons = { id: 'buttons', label: 'Button Mapping', icon: 'buttons', build };
})(typeof globalThis !== 'undefined' ? globalThis : this);
