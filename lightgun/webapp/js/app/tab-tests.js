/*  OpenFIRE Web App - "Gun Tests" tab (Qt: testsTab).

    Live button/switch/temperature/analog stick readings sent by the docked board, feedback tests,
    the IR camera tester, the aim test and the board actions (bootloader, clear save memory).
    It describes the board as it was loaded or last saved (state.testPins).
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    // Last switch update from the board: { read, on } bit masks, null until one arrives (older
    // firmware never sends it). Outside build(), so a rebuilt page (language) keeps it.
    let switchState = null;

    function build(app) {
        const { el, t } = OF.UI;
        const state = app.state;
        const S = state.S;
        const E = state.E;
        const T = state.T;
        const C = S.serialCmdTypes_e;
        const names = OF.Maps.functionNames(S);

        // ----- Inputs test ---------------------------------------------------------------
        const buttonLabels = [];
        const buttonGrid = el('div', { class: 'test-buttons' });
        for (let i = 0; i < OF.Maps.BUTTON_COUNT; ++i) {
            const label = el('div', { class: 'test-button' });
            buttonLabels.push(label);
            buttonGrid.append(label);
        }
        const buttonsGroup = el('fieldset', { class: 'group buttons-group' }, el('legend', { text: t('Buttons') }), buttonGrid);

        // Rumble, solenoid and autofire switches: on, off, or dashed without a pin. The board
        // reads them with the normally closed setting applied, as in Run mode (sSwitchStateUpd).
        const SWITCHES = [
            { fn: E.rumbleSwitch, label: 'Rumble' },
            { fn: E.solenoidSwitch, label: 'Solenoid' },
            { fn: E.autofireSwitch, label: 'Autofire' },
        ];
        const switchItems = SWITCHES.map((sw) => {
            const name = el('span', { class: 'switch-name' });
            const item = el('div', { class: 'test-switch', attrs: { role: 'img' } }, name, el('span', { class: 'switch-pill' }));
            return { ...sw, item, name };
        });
        const switchesGroup = el('fieldset', { class: 'group switches-group' }, el('legend', { text: t('Switch') }),
            el('div', { class: 'test-switches' }, ...switchItems.map((sw) => sw.item)));

        const temperature = el('div', { class: 'temperature' });
        const temperatureBox = el('div', { class: 'temperature-box' }, temperature);

        const positionLabel = el('div', { class: 'analog-label mono' });
        const dot = el('svg:circle', { attrs: { cx: 128, cy: 128, r: 8 } });
        const analogView = el('svg:svg', { class: 'analog-view', attrs: { viewBox: '0 0 256 256', role: 'img', 'aria-label': t('Analog Stick') } },
            el('svg:line', { class: 'axis', attrs: { x1: 128, y1: 0, x2: 128, y2: 256 } }),
            el('svg:line', { class: 'axis', attrs: { x1: 0, y1: 128, x2: 256, y2: 128 } }),
            dot);
        // The stick is as tall as the buttons, its square view shrinks to fit.
        const analogViewBox = el('div', { class: 'analog-view-box' }, analogView);
        const analogGroup = el('fieldset', { class: 'group analog-group' }, el('legend', { text: t('Analog Stick') }), positionLabel,
            analogViewBox);

        // Two rows: buttons and analog stick, then the switches with the temperature beside them.
        const inputsGrid = el('div', { class: 'inputs-test' }, buttonsGroup, analogGroup, switchesGroup, temperatureBox);
        const inputsTest = el('fieldset', { class: 'group' }, el('legend', { text: t('Inputs Test') }), inputsGrid);

        // ----- Feedback tests -----------------------------------------------------------------
        const testButton = (label, command, message) => el('button', { text: t(label), on: { click: async () => {
            if (await app.sendCommand(command)) app.status.show(t(message), 2500);
        } } });
        const rumbleTest = testButton('Test Rumble Motor', C.sTestRumble, 'Sent a rumble test pulse.');
        const solenoidTest = testButton('Test Solenoid', C.sTestSolenoid, 'Sent a solenoid test pulse.');
        const redTest = testButton('Test Red LED', C.sTestLEDR, 'Set LED to Red.');
        const greenTest = testButton('Test Green LED', C.sTestLEDG, 'Set LED to Green.');
        const blueTest = testButton('Test Blue LED', C.sTestLEDB, 'Set LED to Blue.');
        const feedbackTests = el('fieldset', { class: 'group' }, el('legend', { text: t('Feedback Tests') }),
            el('div', { class: 'row fill feedback-row' }, rumbleTest, solenoidTest),
            el('div', { class: 'row fill feedback-row' }, redTest, greenTest, blueTest));

        // The IR camera tester and the aim test share one row (same test mode of the gun).
        const irTest = el('button', { text: t('Open IR Camera Tester...'), on: { click: () => app.openIRTest() } });
        const aimTest = el('button', { text: t('Open Aim Test...'), on: { click: () => app.openIRTest(true) } });
        const testWindows = el('div', { class: 'row fill test-windows' }, irTest, aimTest);

        // ----- Board actions --------------------------------------------------------------------
        const reboot = el('button', { on: { click: () => app.rebootToBootloader() } });
        const clear = el('button', { class: 'danger-text', text: t('Clear Save Memory [!]'), on: { click: () => app.clearSaveMemory() } });
        const boardActions = el('fieldset', { class: 'group' }, el('legend', { text: t('Board Actions') }),
            el('div', { class: 'row fill board-actions' }, reboot, clear));

        const content = el('div', { class: 'tab-content narrow' }, inputsTest, feedbackTests, testWindows, el('hr'), boardActions);
        const scroll = el('div', { class: 'tab-scroll' }, content);
        const rootNode = el('div', { class: 'tab-body' }, scroll);

        // Like the board layout: the whole tab is shown, it shrinks instead of being scrolled
        // (on narrow screens the single column is scrolled as usual).
        const wide = root.matchMedia ? root.matchMedia('(min-width: 761px)') : null;
        /** The stick column as wide as the square view is tall (the height the buttons leave it),
            so the view stays square and the buttons take the rest. Not in the single column. */
        function squareAnalog() {
            if (!analogViewBox.clientHeight) return; // tab not shown
            if (wide && !wide.matches) {
                inputsGrid.style.removeProperty('--analog-width');
            } else {
                for (let pass = 0; pass < 3; ++pass) {
                    const sides = analogGroup.offsetWidth - analogViewBox.clientWidth;
                    const width = Math.ceil(analogViewBox.clientHeight + sides);
                    if (Math.abs(width - analogGroup.offsetWidth) <= 1) break;
                    inputsGrid.style.setProperty('--analog-width', width + 'px');
                }
            }
            // Square whatever is left (the legend may make the column wider than the view is tall).
            const side = Math.min(analogViewBox.clientWidth, analogViewBox.clientHeight);
            analogView.style.width = analogView.style.height = side + 'px';
        }
        const fit = () => {
            squareAnalog();
            OF.UI.fitToHeight(scroll, content, 0.6, !wide || wide.matches);
        };
        OF.UI.fitOnResize(scroll, fit);
        if (wide && wide.addEventListener) wide.addEventListener('change', fit);

        const mapped = (fn) => (state.testPins[fn] ?? -1) >= 0;

        /**
         * Name of an input box: the plain name when it has a pin; '(N/C)' when it has none.
         * The two pedals also work without a pin, over the air: when the board is in that mode
         * they say so, and add that nothing answered when no pedal is paired.
         */
        function inputLabel(index) {
            // Short name where the box is small: in Italian "Pedale alternativo (wireless, disconnesso)"
            // took four lines and made the tab shrink to fit the window.
            const name = names[index + 1] === 'Alt Pedal' ? t('Alt. Pedal') : t(names[index + 1]);
            if (mapped(index))
                return { text: name, connected: true };
            if (state.isPedal(index) && state.pedalWirelessMode)
                return state.pedalWireless ?
                    { text: t('%1 (wireless)', name), connected: true } :
                    { text: t('%1 (wireless, disconnected)', name), connected: false };
            return { text: `${name} (N/C)`, connected: false };
        }

        /** The three switches, from the pins at the last load/save and the board's last update. */
        function updateSwitches() {
            switchItems.forEach((sw, i) => {
                const connected = mapped(sw.fn);
                const known = connected && switchState !== null;
                // Without that bit the board does not read the switch (a build without it).
                const read = known && !!((switchState.read >> i) & 1);
                const on = read && !!((switchState.on >> i) & 1);
                const off = !connected || (known && !read);
                sw.name.textContent = t(sw.label);
                sw.item.classList.toggle('on', on);
                sw.item.classList.toggle('not-connected', off);
                sw.item.classList.toggle('unknown', connected && !known);
                const state = off ? t('Not Connected') : !known ? t('Switch state not available: update the firmware.') :
                    on ? t('Enabled') : t('Disabled');
                sw.item.title = state;
                sw.item.setAttribute('aria-label', `${t(sw.label)}: ${state}`);
            });
        }

        /** Qt LabelsUpdate: on load and after a save. */
        function resetReadings() {
            buttonLabels.forEach((label, i) => {
                label.classList.remove('pressed');
                const { text, connected } = inputLabel(i);
                label.textContent = text;
                label.classList.toggle('not-connected', !connected);
            });
            temperature.className = 'temperature';
            if (mapped(E.tempPin)) {
                temperature.textContent = t('Temperature Read...');
                temperature.classList.remove('not-connected');
            } else {
                temperature.textContent = t('Temperature Sensor (N/C)');
                temperature.classList.add('not-connected');
            }
            const analog = mapped(E.analogX) && mapped(E.analogY);
            positionLabel.textContent = analog ? '' : t('Not Connected');
            dot.setAttribute('cx', 128);
            dot.setAttribute('cy', 128);
            // The last switch states stay: the board repeats them twice a second, with the
            // pins and the normally closed setting it has after a save.
            updateSwitches();
        }

        /** Buttons only (Qt CaliWindowExiting): the board does not send the releases made
         *  during calibration, so no button stays shown as pressed when it ends. */
        function clearButtons() {
            buttonLabels.forEach((label) => label.classList.remove('pressed'));
        }

        /** A new board session: nothing known about its switches yet. */
        function onLoad() {
            switchState = null;
            resetReadings();
        }

        function update() {
            if (!state.loaded) return;
            const recovery = app.commitNeedsRetry;
            const testing = app.irTestActive;
            inputsTest.disabled = recovery || testing;
            feedbackTests.disabled = recovery || testing;
            irTest.disabled = recovery;
            aimTest.disabled = recovery;
            boardActions.disabled = testing;
            analogGroup.disabled = !(mapped(E.analogX) && mapped(E.analogY));
            rumbleTest.disabled = !state.toggle(state.B.rumble);
            solenoidTest.disabled = !state.toggle(state.B.solenoid);
            redTest.disabled = !mapped(E.ledR);
            greenTest.disabled = !mapped(E.ledG);
            blueTest.disabled = !mapped(E.ledB);
            reboot.textContent = t(state.isRP ? 'Reboot to Bootloader' : 'Restart Microcontroller in Firmware Update Mode');
            fit();
        }

        /** Board events for this tab. Returns true when handled. */
        function onEvent(command, payload) {
            switch (command) {
            case C.sBtnPressed:
            case C.sBtnReleased:
                if (payload.length === 1 && payload[0] < buttonLabels.length)
                    buttonLabels[payload[0]].classList.toggle('pressed', command === C.sBtnPressed);
                return true;
            case C.sSwitchStateUpd:
                if (payload.length >= 2) {
                    switchState = { read: payload[0], on: payload[1] };
                    updateSwitches();
                }
                return true;
            case C.sTemperatureUpd: {
                if (payload.length !== 1) return true;
                const temp = payload[0];
                temperature.classList.remove('not-connected', 'hot', 'warm', 'ok', 'fault');
                if (temp === S.TEMPERATURE_SENSOR_ERROR_VALUE) {
                    temperature.textContent = t('Temperature: FAULT');
                    temperature.classList.add('fault');
                    return true;
                }
                temperature.textContent = t('Temperature: %1°C', temp);
                if (temp > state.orig.settings[T.tempShutdown]) temperature.classList.add('hot');
                else if (temp > state.orig.settings[T.tempWarning]) temperature.classList.add('warm');
                else temperature.classList.add('ok');
                return true;
            }
            case C.sAnalogPosUpd: {
                if (payload.length !== 4) return true;
                const view = new DataView(payload.buffer, payload.byteOffset, 4);
                const x = view.getUint16(0, true);
                const y = view.getUint16(2, true);
                // Qt: uint8_t pos = ~(value / 16)
                dot.setAttribute('cx', (~Math.floor(x / 16)) & 0xFF);
                dot.setAttribute('cy', (~Math.floor(y / 16)) & 0xFF);
                positionLabel.textContent = `${x} , ${y}`;
                return true;
            }
            default:
                return false;
            }
        }

        return { root: rootNode, update, onEvent, resetReadings, clearButtons, onLoad };
    }

    OF.Tabs = OF.Tabs || {};
    OF.Tabs.tests = { id: 'tests', label: 'Gun Tests', icon: 'tests', build };
})(typeof globalThis !== 'undefined' ? globalThis : this);
