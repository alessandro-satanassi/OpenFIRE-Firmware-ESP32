/*  OpenFIRE Web App - "Calibration Profiles" tab (Qt: profilesTab).

    One row per profile: rename, active profile (selected immediately on the board),
    calibration values, IR sensitivity, run mode, layout, aspect ratio and colour.
    "Calibrate Profile N" opens the fullscreen calibration window.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const TAB_HELP = '<html><head/><body><p>This tab shows information and settings to tweak about your current calibration profiles.</p><p>The table above represents the amount of profiles the current board can store, including currently selected profile, names shown for each profile when using a compatible <span style=" font-style:italic;">I2C Display,</span> and colors used for lighting devices when switching to them from <span style=" font-style:italic;">Pause Mode.</span></p><p>Hover over an option to view detailed info about it here.</p></body></html>';
    const SERIAL_LINK = "<a href='https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/wiki/MAMEHOOKER-Documentation#m---mode-commands'><span style=' text-decoration: underline; color:#8ab4f8;'>Serial command</span></a> ";

    const HELP = {
        rename: ['Rename Profile %1',
            '<p>Click to rename this Calibration Profile.</p>' +
            '<p>Aside from differentiating between different profiles for different displays, ' +
            'Cali Profile names are displayed in Pause Mode when using a compatible <i>I2C Display.</i></p>'],
        irSensitivity: ['Camera Sensitivity for Profile %1',
            '<p>This setting determines the sensitivity of the IR Camera for this Calibration Profile.</p>' +
            '<p>If the camera seems to have trouble picking up IR emitters (and is causing coarse cursor movement), ' +
            'adjusting this setting higher might fix issues with tracking.<br>' +
            'Conversely, setting sensitivity too high may cause indirect IR sources ' +
            '(such as sunlight or IR bouncing off of reflective surfaces) ' +
            'to be picked up instead, causing the cursor to jitter or erratically jump across the screen.</p>'],
        runMode: ['Camera Position Averaging Mode for Profile %1',
            '<p>This setting determines the cursor Averaging Mode for this Calibration Profile.</p>' +
            '<p>The movement of the aiming cursor can be smoothed out by averaging a select number of frames, ' +
            'at the cost of a small increase in latency; conversely, disabling this position averaging can ' +
            'reduce latency, at the cost of some added jitter in mouse movement.</p>' +
            '<p>The default is <b>1-Frame Avg</b>, which should be the preferred balance for most people.</p>'],
        layoutType: ['IR Emitter Layout for Profile %1',
            '<p>This setting determines the IR Layout to be used with this Calibration Profile.</p>' +
            '<p>Each Cali Profile can be set to use either the <i>Square Layout,</i> ' +
            'which uses two pairs of emitters on the top and bottom, and <i>Diamond Layout,</i> ' +
            'which uses one emitter at the center of each side of the display.</p>' +
            '<p><i>Square Layout</i> generally has much higher accuracy at any angle and allows for ' +
            'playing closer to the screen or using external Fish Eye lenses without viewport distortion, ' +
            'while <i>Diamond Layout</i> is for screen compatibility with certain legacy lightgun systems ' +
            '(allowing OpenFIRE guns to play with such other lightgun systems on the same display).</p>' +
            '<p>If unsure, use <b>Square Layout</b> ' +
            '(unless you also use a different brand of lightgun that needs a diamond IR layout to function).</p>'],
        aspectRatio: ['Aspect Ratio Correction for Profile %1',
            '<p>This setting determines the type of Aspect Ratio Correction used for this Profile when <i>4:3 Mode</i> is set.</p>' +
            '<p>When ' + SERIAL_LINK +
            '<tt>M3x1</tt> is received, the firmware stretches the effective range for fullscreen applications in Windows ' +
            'that runs in resolutions <b>narrower</b> than the full display width/height; ' +
            'this setting determines the stretch factor that will be used for these 4:3 applications.</p>' +
            '<p>Do note that this restriction exclusively applies to the <b>Windows Operating System ONLY ' +
            "for legacy applications that DON'T support the monitor's full resolution;</b> <i>Linux-based systems</i> and games run via <i>Wine/Proton</i> " +
            '<b>do not need this workaround,</b> except for certain applications like <i>CXBX-Reloaded</i> that don\'t scale down the effective range correctly for 4:3 content.</p>' +
            '<p>If unsure, set to <b>the aspect ratio of your display.</b> Setting to <i>4:3</i> effectively disables range stretching when toggled.</p>'],
        color: ['Profile Menu Color for Cali Profile %1',
            '<p>Open a window to select the color used to represent this profile in <i>Pause Mode.</i></p>' +
            '<p>Each profile can be assigned a color used to identify them when switching profiles on the lightgun itself, ' +
            'which is emitted by a 4-pin RGB LED and/or an active NeoPixel strand.</p>'],
        calibrate: ['Open Calibration Window for Cali Profile %1', 'Click to start the calibration process for this profile.'],
    };

    const VALUE_COLUMNS = [
        ['topOffset', 'Top'], ['bottomOffset', 'Bottom'], ['leftOffset', 'Left'],
        ['rightOffset', 'Right'], ['TLled', 'TLled'], ['TRled', 'TRled'],
    ];
    const OPTION_COLUMNS = [
        ['irSensitivity', 'Sensitivity', 'irSensitivity'], ['runMode', 'Run Mode', 'runModes'],
        ['layoutType', 'Layout', 'layouts'], ['aspectRatio', 'Display Ratio', 'aspectRatios'],
    ];

    /** Qt QString::number: offsets are integers, TLled/TRled floats ("%g", 6 significant digits). */
    function formatNumber(value, isFloat) {
        return isFloat ? OF.UI.formatG(value) : String(Math.trunc(Number(value)));
    }

    function build(app) {
        const { el, t, icon, select } = OF.UI;
        const state = app.state;
        const desc = new OF.UI.DescriptionBox(t(TAB_HELP));
        const help = (key, n) => [t(HELP[key][0], n), t(HELP[key][1], n)];

        const table = el('table', { class: 'profiles-table' });
        const tableScroll = el('div', { class: 'table-scroll' }, table);
        const buttons = el('div', { class: 'cali-buttons' });

        // The whole table is shown: it shrinks instead of being scrolled sideways
        // (on narrow screens it is scrolled as usual).
        const wide = root.matchMedia ? root.matchMedia('(min-width: 761px)') : null;
        const fit = () => OF.UI.fitToWidth(tableScroll, table, 0.6, !wide || wide.matches);
        OF.UI.fitOnResize(tableScroll, fit);
        if (wide && wide.addEventListener) wide.addEventListener('change', fit);
        const rootNode = el('div', { class: 'tab-body' },
            el('div', { class: 'tab-scroll' },
                el('div', { class: 'tab-content' },
                    el('fieldset', { class: 'group' }, el('legend', { text: t('Calibration Profiles') }),
                        tableScroll, buttons))),
            desc.root);

        let rows = [];
        let builtCount = -1;

        function populate() {
            table.textContent = '';
            buttons.textContent = '';
            rows = [];
            builtCount = state.profileCount;

            table.append(el('thead', null, el('tr', null,
                el('th', { attrs: { colspan: 2 } }),
                ...VALUE_COLUMNS.map(([, title]) => el('th', { text: t(title) })),
                ...OPTION_COLUMNS.map(([, title]) => el('th', { text: t(title) })),
                el('th', { text: t('Color') }))));

            const body = el('tbody');
            for (let i = 0; i < builtCount; ++i) {
                const n = i + 1;
                const rename = el('button', { class: 'icon-button', title: t(HELP.rename[0], n), attrs: { 'aria-label': t(HELP.rename[0], n) },
                    on: { click: () => renameProfile(i) } }, icon('edit'));
                desc.track(rename, ...help('rename', n));

                const radioInput = el('input', { type: 'radio', name: 'of-profile' });
                const radioText = el('span', { class: 'mono' });
                const radio = el('label', { class: 'check radio profile-radio' }, radioInput, radioText);
                radio.input = radioInput;
                radioInput.addEventListener('change', () => { if (radioInput.checked) selectProfile(i); });

                const values = VALUE_COLUMNS.map(() => el('td', { class: 'mono value' }));
                const options = OPTION_COLUMNS.map(([key, , list]) => {
                    const box = select(OF.Maps[list].map((label, value) => ({ value, label: t(label) })), 0, (value) => {
                        state.setProfileField(i, key, Number(value));
                        app.refresh();
                    }, { attrs: { 'aria-label': t(HELP[key][0], n) } });
                    desc.track(box, ...help(key, n));
                    return box;
                });

                const swatch = el('span', { class: 'color-swatch' });
                const color = el('button', { class: 'color-button only', title: t(HELP.color[0], n), attrs: { 'aria-label': t(HELP.color[0], n) },
                    on: { click: () => chooseColor(i) } }, swatch);
                desc.track(color, ...help('color', n));

                body.append(el('tr', null,
                    el('td', null, rename), el('td', null, radio),
                    ...values,
                    ...options.map((box) => el('td', null, box)),
                    el('td', null, color)));

                const calibrate = el('button', { text: t('Calibrate Profile %1', n), on: { click: () => app.calibrate(i) } });
                desc.track(calibrate, ...help('calibrate', n));
                buttons.append(calibrate);

                rows.push({ radio, radioText, values, options, swatch, calibrate });
            }
            table.append(body);
        }

        async function renameProfile(i) {
            const session = app.session;
            const name = await OF.UI.prompt(t('Input Name'), t('Set name for Calibration Profile %1', i + 1),
                { value: state.cur.profiles[i].name, maxLength: OF.AppState.NAME_MAX, latin1: true, scope: 'session' });
            if (session !== app.session || !state.loaded || app.busy) return;
            if (name) state.renameProfile(i, name);
            app.refresh();
        }

        async function chooseColor(i) {
            const session = app.session;
            const color = await OF.UI.pickColor(t(HELP.color[0], i + 1), state.cur.profiles[i].color, { scope: 'session' });
            if (session !== app.session || !state.loaded || app.busy) return;
            if (color !== null) {
                state.setProfileField(i, 'color', color);
                app.refresh();
            }
        }

        function selectProfile(i) {
            if (i !== state.cur.selectedProfile) app.selectProfile(i);
        }

        function update() {
            if (!state.loaded) return;
            if (builtCount !== state.profileCount) populate();
            fit();
            const locked = app.commitNeedsRetry;
            rows.forEach((row, i) => {
                const profile = state.cur.profiles[i];
                row.radioText.textContent = `${i + 1}. ${profile.name}`;
                row.radio.input.checked = i === state.cur.selectedProfile;
                row.radio.input.disabled = locked;
                VALUE_COLUMNS.forEach(([key], column) => { row.values[column].textContent = formatNumber(profile[key], key === 'TLled' || key === 'TRled'); });
                OPTION_COLUMNS.forEach(([key], column) => { row.options[column].value = String(profile[key]); });
                row.swatch.style.background = OF.UI.toHex(profile.color);
                row.calibrate.disabled = locked;
            });
        }

        return { root: rootNode, update, onShow() { desc.reset(); }, reload() { builtCount = -1; update(); } };
    }

    OF.Tabs = OF.Tabs || {};
    OF.Tabs.profiles = { id: 'profiles', label: 'Calibration Profiles', icon: 'profiles', build, formatNumber };
})(typeof globalThis !== 'undefined' ? globalThis : this);
