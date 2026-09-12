/*  OpenFIRE Web App - "Gun Tests" tab (Qt: testsTab).

    Live button/temperature/analog stick readings sent by the docked board, feedback tests,
    the IR camera tester and the board actions (bootloader, clear save memory).
    It describes the board as it was loaded or last saved (state.testPins).
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

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
        const temperature = el('div', { class: 'temperature' });
        const buttonsGroup = el('fieldset', { class: 'group' }, el('legend', { text: t('Buttons') }), buttonGrid, temperature);

        const positionLabel = el('div', { class: 'analog-label mono' });
        const dot = el('svg:circle', { attrs: { cx: 128, cy: 128, r: 8 } });
        const analogView = el('svg:svg', { class: 'analog-view', attrs: { viewBox: '0 0 256 256', role: 'img', 'aria-label': t('Analog Stick') } },
            el('svg:line', { class: 'axis', attrs: { x1: 128, y1: 0, x2: 128, y2: 256 } }),
            el('svg:line', { class: 'axis', attrs: { x1: 0, y1: 128, x2: 256, y2: 128 } }),
            dot);
        const analogGroup = el('fieldset', { class: 'group analog-group' }, el('legend', { text: t('Analog Stick') }), positionLabel, analogView);

        const inputsTest = el('fieldset', { class: 'group' }, el('legend', { text: t('Inputs Test') }),
            el('div', { class: 'inputs-test' }, buttonsGroup, analogGroup));

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
            el('div', { class: 'row fill' }, rumbleTest, solenoidTest),
            el('div', { class: 'row fill' }, redTest, greenTest, blueTest));

        const irTest = el('button', { class: 'wide-button', text: t('Open IR Camera Tester...'), on: { click: () => app.openIRTest() } });

        // ----- Board actions --------------------------------------------------------------------
        const reboot = el('button', { on: { click: () => app.rebootToBootloader() } });
        const clear = el('button', { class: 'danger-text', text: t('Clear Save Memory [!]'), on: { click: () => app.clearSaveMemory() } });
        const boardActions = el('fieldset', { class: 'group' }, el('legend', { text: t('Board Actions') }),
            el('div', { class: 'row fill board-actions' }, reboot, clear));

        const content = el('div', { class: 'tab-content narrow' }, inputsTest, feedbackTests, irTest, el('hr'), boardActions);
        const scroll = el('div', { class: 'tab-scroll' }, content);
        const rootNode = el('div', { class: 'tab-body' }, scroll);

        // Like the board layout: the whole tab is shown, it shrinks instead of being scrolled
        // (on narrow screens the single column is scrolled as usual).
        const wide = root.matchMedia ? root.matchMedia('(min-width: 761px)') : null;
        const fit = () => OF.UI.fitToHeight(scroll, content, 0.6, !wide || wide.matches);
        OF.UI.fitOnResize(scroll, fit);
        if (wide && wide.addEventListener) wide.addEventListener('change', fit);

        const mapped = (fn) => (state.testPins[fn] ?? -1) >= 0;

        /** Qt LabelsUpdate: on load and after a save. */
        function resetReadings() {
            buttonLabels.forEach((label, i) => {
                label.classList.remove('pressed');
                label.textContent = mapped(i) ? t(names[i + 1]) : `${t(names[i + 1])} (N/C)`;
                label.classList.toggle('not-connected', !mapped(i));
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
        }

        function update() {
            if (!state.loaded) return;
            const recovery = app.commitNeedsRetry;
            const testing = app.irTestActive;
            inputsTest.disabled = recovery || testing;
            feedbackTests.disabled = recovery || testing;
            irTest.disabled = recovery;
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

        return { root: rootNode, update, onEvent, resetReadings, onLoad: resetReadings };
    }

    OF.Tabs = OF.Tabs || {};
    OF.Tabs.tests = { id: 'tests', label: 'Gun Tests', icon: 'tests', build };
})(typeof globalThis !== 'undefined' ? globalThis : this);
