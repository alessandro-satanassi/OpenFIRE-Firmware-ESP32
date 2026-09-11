/*  OpenFIRE Web App: editing rules of the configuration (js/app/state.js), ported from the Qt App.

    Run from lightgun/webapp:   node --test --test-concurrency=1 tests/*.test.js
*/
'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const WEBAPP = path.join(__dirname, '..');
const load = (relative) => vm.runInThisContext(fs.readFileSync(path.join(WEBAPP, relative), 'utf8'), { filename: relative });

load('boards/OpenFIREshared.js');
vm.runInThisContext('globalThis.OpenFIREshared = OpenFIREshared;');
load('js/core/boards.js');
load('js/core/protocol.js');
load('js/app/maps.js');
load('js/app/state.js');
load('js/app/ui.js');
load('js/app/testfont.js');
load('js/app/fullscreen.js');

const OF = globalThis.OF;
const S = globalThis.OpenFIREshared;
const E = S.boardInputs_e;
const B = S.boolTypes_e;
const T = S.settingsTypes_e;

function makeConfig({ customPins = false, pins = null } = {}) {
    const codec = new OF.ProtocolCodec(S);
    const config = codec.createConfig();
    config.toggles[B.customPins] = customPins;
    config.toggles[B.rumble] = true;
    config.toggles[B.solenoid] = true;
    if (pins) config.pins = pins.slice();
    config.settings[T.customLEDcount] = 3;
    config.settings[T.customLEDstatic] = 2;
    config.settings[T.tempWarning] = 42;
    config.settings[T.tempShutdown] = 50;
    config.buttons[0] = [0, 1, 0, 1, 2, 9];
    for (let i = 0; i < 4; ++i) {
        const profile = codec.createProfile();
        profile.name = `Profile ${i + 1}`;
        config.profiles.push(profile);
    }
    config.selectedProfile = 0;
    config.tinyUSB = { id: 1, name: 'FIRECon P1' };
    return config;
}

function loaded(type = 'rpipico', options) {
    const state = new OF.AppState(S);
    state.load({ version: '6.2', type, arch: OF.Boards.arch(type), cameraError: false }, makeConfig(options));
    return state;
}

test('load: custom pins off shows the board default layout, nothing to save', () => {
    const state = loaded();
    const preset = S.boardsPresetsMap.rpipico;
    preset.forEach((fn, gpio) => { if (fn > E.btnUnmapped) assert.equal(state.pinOf(fn), gpio); });
    assert.deepEqual(state.testPins, state.cur.pins);
    assert.equal(state.isDirty(), false);
    assert.equal(state.customPins, false);
    assert.equal(state.prettyName(state.cur.tinyUSB.name), 'FIRECon P1 | Raspberry Pi Pico (RP2040)');
    assert.equal(state.prettyName(''), 'Unnamed Device | Raspberry Pi Pico (RP2040)');
});

test('load: custom pins from the board, duplicates resolved like the Qt pin boxes', () => {
    const pins = new Array(E.boardInputsCount).fill(-1);
    pins[E.btnTrigger] = 5;
    pins[E.btnGunA] = 5; // same GPIO twice: the later function wins
    pins[E.rumblePin] = 16;
    const state = loaded('rpipico', { customPins: true, pins });
    assert.equal(state.functionAt(5), E.btnGunA);
    assert.equal(state.pinOf(E.btnTrigger), -1);
    assert.equal(state.pinOf(E.rumblePin), 16);
});

test('pins: a function moves to its new GPIO and the previous box is unmapped', () => {
    const state = loaded();
    state.setCustomPins(true);
    assert.equal(state.isDirty(), true);
    const triggerPin = state.pinOf(E.btnTrigger);
    state.changePin(0, E.btnTrigger + 1);
    assert.equal(state.pinOf(E.btnTrigger), 0);
    assert.equal(state.functionAt(triggerPin), -1);
    state.changePin(0, 0);
    assert.equal(state.pinOf(E.btnTrigger), -1);
});

test('pins: I2C channel rules of the Qt App (pair on another channel, camera/peripheral clash)', () => {
    const state = loaded();
    state.setCustomPins(true);
    assert.equal(state.pinOf(E.camSDA), 20);    // I2C0
    assert.equal(state.pinOf(E.periphSDA), 18); // I2C1
    assert.equal(state.pinOf(E.periphSCL), 19); // I2C1
    const messages = state.changePin(4, E.periphSDA + 1); // GPIO4 is I2C0 SDA
    assert.deepEqual(messages, [
        "Peripheral pins are not on the same I2C channel! Please check peripheral pins' mappings.",
        'Camera and Peripheral Data pins clashed! Please remap Camera SDA.',
    ]);
    assert.equal(state.pinOf(E.periphSDA), 4);
    assert.equal(state.pinOf(E.periphSCL), -1);
    assert.equal(state.pinOf(E.camSDA), -1);
});

test('pins: capabilities disable analog and I2C functions (RP channels, ESP32 any pin)', () => {
    const rp = loaded();
    assert.equal(rp.functionDisabled(0, E.analogX), true);      // no ADC
    assert.equal(rp.functionDisabled(26, E.analogX), false);    // ADC0
    assert.equal(rp.functionDisabled(4, E.camSDA), false);      // I2C0 SDA
    assert.equal(rp.functionDisabled(4, E.camSCL), true);
    assert.equal(rp.functionDisabled(5, E.periphSDA), true);    // I2C0 SCL
    assert.equal(rp.pinI2CClass(4), 'i2c0');
    assert.equal(rp.pinI2CClass(18), 'i2c1');
    const esp = loaded('waveshare-esp32-s3-zero');
    assert.equal(esp.pinI2CClass(4), 'any');
    assert.equal(esp.functionDisabled(4, E.camSCL), false);
    assert.equal(esp.pinI2CClass(45), 'none');
    assert.equal(esp.functionDisabled(45, E.camSDA), true);
});

test('pins: unmapping the rumble or solenoid pin turns their feedback off', () => {
    const state = loaded();
    state.setCustomPins(true);
    state.setToggle(B.solenoid, false);
    state.setToggle(B.rumbleFF, true);
    state.setToggle(B.autofire, true);
    assert.equal(state.toggle(B.autofire), true);
    state.changePin(state.pinOf(E.rumblePin), 0);
    assert.equal(state.toggle(B.rumble), false);
    assert.equal(state.toggle(B.rumbleFF), false);
    assert.equal(state.toggle(B.autofire), false);
    state.setToggle(B.solenoid, true);
    state.changePin(state.pinOf(E.solenoidPin), 0);
    assert.equal(state.toggle(B.solenoid), false);
});

test('toggles: solenoid and rumble force feedback exclude each other, autofire needs one of them', () => {
    const state = loaded();
    assert.equal(state.autofireAvailable, true);
    state.setToggle(B.autofire, true);
    state.setToggle(B.rumbleFF, true);
    assert.equal(state.toggle(B.solenoid), false);
    assert.equal(state.toggle(B.autofire), true);   // rumble + rumbleFF
    state.setToggle(B.solenoid, true);
    assert.equal(state.toggle(B.rumbleFF), false);
    state.setToggle(B.solenoid, false);
    assert.equal(state.toggle(B.autofire), false);
    assert.equal(state.autofireAvailable, false);
});

test('custom pins off restores the default layout; alternative presets turn custom pins on', () => {
    const state = loaded();
    state.setCustomPins(true);
    state.changePin(0, E.btnPedal + 1);
    state.setCustomPins(false);
    assert.equal(state.functionAt(0), S.boardsPresetsMap.rpipico[0]);
    assert.equal(state.toggle(B.rumble), true); // feedback pins are mapped again

    const presets = state.altPresets;
    assert.ok(presets.length > 0);
    state.applyAltPreset(0);
    assert.equal(state.customPins, true);
    assert.equal(state.altPresetIndex, 0);
    presets[0].pin.forEach((fn, gpio) => { if (fn > E.btnUnmapped && gpio < state.pinCount) assert.equal(state.functionAt(gpio), fn); });
    state.changePin(1, 0);
    assert.equal(state.altPresetIndex, -1);
});

test('custom layout files (.ofl): export and import, other boards refused', () => {
    const state = loaded();
    state.setCustomPins(true);
    state.changePin(0, E.btnPedal + 1);
    const bytes = state.exportLayout();
    assert.equal(String.fromCharCode(...bytes.subarray(0, 8)), 'rpipico\n');

    const other = loaded();
    assert.equal(other.importLayout(bytes), 'ok');
    assert.equal(other.customPins, true);
    assert.deepEqual(other.cur.pins, state.cur.pins);

    const esp = loaded('waveshare-esp32-s3-zero');
    assert.equal(esp.importLayout(bytes), 'mismatch');
});

test('save: commitDone makes the edits the saved state (pins cleared when custom pins are off)', () => {
    const state = loaded();
    state.setSetting(T.rumbleStrength, 128);
    state.renameProfile(1, 'A very long profile name');
    assert.equal(state.cur.profiles[1].name, 'A very long pro');
    assert.equal(state.isDirty(), true);
    state.commitDone();
    assert.equal(state.isDirty(), false);
    assert.ok(state.orig.pins.every((pin) => pin === -1));
    state.setSelectedProfile(2);
    assert.equal(state.isDirty(), true);
});

test('button mapping: changing the output type picks its first output, the saved type gets back its output', () => {
    const state = loaded();
    state.setButtonType(0, 0, OF.Maps.INPUT_KEYBOARD);
    assert.equal(state.buttonValue(0, 0), 0xFF);
    state.setButtonValue(0, 0, 'a'.charCodeAt(0));
    assert.equal(state.isDirty(), true);
    state.setButtonType(0, 0, OF.Maps.INPUT_MOUSE);
    assert.equal(state.buttonValue(0, 0), 1);
    assert.equal(state.isDirty(), false);
});

test('settings: static NeoPixels never exceed the strand length', () => {
    const state = loaded();
    state.setSetting(T.customLEDstatic, 3);
    assert.equal(state.setting(T.customLEDstatic), 3);
    state.setSetting(T.customLEDcount, 1);
    assert.equal(state.setting(T.customLEDstatic), 1);
    state.setSetting(T.customLEDstatic, 3);
    assert.equal(state.setting(T.customLEDstatic), 1);
    assert.equal(state.pixelsChanged(), true);
});

test('TinyUSB identifier: player presets, Latin-1 names, id with an empty name', () => {
    const state = loaded();
    state.setUsbPreset(3);
    assert.deepEqual(state.cur.tinyUSB, { id: 3, name: 'FIRECon P3' });
    assert.equal(state.usbPreset, 3);
    assert.equal(state.setUsbName('Pistola ñ'), true);
    assert.equal(state.setUsbName('Ω gun'), false);
    assert.equal(state.cur.tinyUSB.name, 'Pistola ñ');
    state.setUsbName('');
    state.setUsbId(2);
    assert.deepEqual(state.cur.tinyUSB, { id: 2, name: 'FIRECon P2' });
    state.setUsbId(0x1234);
    assert.equal(state.usbPreset, 0);
});

test('calibration results: cancelled, successful and malformed values', () => {
    const state = loaded();
    const none = { topOffset: -1, bottomOffset: -1, leftOffset: -1, rightOffset: -1, TLled: -1, TRled: -1 };
    assert.equal(state.applyCalibration(1, none), 'cancelled');
    assert.equal(state.isDirty(), false);
    assert.equal(state.applyCalibration(1, { topOffset: 10, bottomOffset: 20, leftOffset: 30, rightOffset: 40, TLled: 1023.5, TRled: 6656.25 }), 'ok');
    assert.equal(state.cur.profiles[1].TRled, 6656.25);
    assert.equal(state.applyCalibration(1, { topOffset: 99999, bottomOffset: 20, leftOffset: 30, rightOffset: 40, TLled: 1, TRled: 2 }), 'malformed');
    state.setProfileField(0, 'layoutType', 1);
    assert.equal(state.layoutChangedUnsaved(), true);
});

test('maps: output lists of the Qt App', () => {
    const M = OF.Maps;
    assert.equal(M.outputs[M.INPUT_MOUSE].length, 5);
    assert.equal(M.outputs[M.INPUT_KEYBOARD].length, 54);
    assert.equal(M.outputs[M.INPUT_GAMEPAD].length, 16);
    for (const list of M.outputs) assert.equal(new Set(list.map((entry) => entry.value)).size, list.length);
    assert.equal(M.outputs[M.INPUT_KEYBOARD].find((entry) => entry.name === 'F12').value, 0xCD);
    const names = M.functionNames(S);
    assert.equal(names[0], 'Unmapped');
    assert.equal(names[E.btnTrigger + 1], 'Trigger');
    assert.equal(names.length, E.boardInputsCount + 1);
});

test('bitmap typeface: glyphs, composed accents, missing characters', () => {
    const glyph = OF.FullscreenWindow.glyphPixels;
    const e = glyph('e');
    const eGrave = glyph('è');
    assert.equal(e.length, 64);
    assert.ok(e.subarray(0, 16).every((v) => v === 0));
    assert.ok(eGrave.subarray(0, 16).some((v) => v === 1));
    assert.deepEqual(Array.from(eGrave.subarray(16)), Array.from(e.subarray(16)));
    assert.deepEqual(Array.from(glyph('È')), Array.from(glyph('E'))); // no room above capitals
    assert.ok(glyph(' ').every((v) => v === 0));
    assert.equal(glyph('日'), null);
});

test('load: values the Qt widgets cannot show are corrected and can be saved', () => {
    const config = makeConfig();
    config.settings[T.tempWarning] = 5;          // Qt spin box 20..99
    config.settings[T.rumbleStrength] = 1000;    // 0..255
    config.settings[T.customLEDcount] = 0;       // 1..100
    config.settings[T.customLEDstatic] = 3;      // not above the strand length
    config.toggles[B.rumbleFF] = true;           // with the solenoid on
    config.toggles[B.autofire] = true;
    const state = new OF.AppState(S);
    state.load({ version: '6.2', type: 'rpipico', arch: OF.Boards.arch('rpipico') }, config);
    assert.equal(state.setting(T.tempWarning), 20);
    assert.equal(state.setting(T.rumbleStrength), 255);
    assert.equal(state.setting(T.customLEDcount), 1);
    assert.equal(state.setting(T.customLEDstatic), 1);
    assert.equal(state.toggle(B.solenoid), false);
    assert.equal(state.toggle(B.rumbleFF), true);
    assert.equal(state.toggle(B.autofire), true);      // rumble used as force feedback
    assert.equal(state.orig.settings[T.tempWarning], 5); // the board keeps its value until Save
    assert.equal(state.isDirty(), true);
    for (const [name, [min, max]] of Object.entries(OF.AppState.SETTING_RANGES))
        assert.ok(T[name] !== undefined && min <= max, name);
});

test('unknown board: generic boxes without a default layout (Qt boardsPresetsMap.count)', () => {
    const state = loaded('not-a-known-board');
    assert.equal(state.pinCount, S.boardsBoxPositions.generic.length);
    assert.ok(state.cur.pins.every((pin) => pin === -1));
    assert.deepEqual(state.presetPins, []);
    assert.ok(loaded().presetPins.length > 0);
});

test('lost session: unsaved edits and RAM-only calibrations can be restored', () => {
    const before = loaded();
    assert.equal(before.snapshotDiffers(before.snapshot()), false);
    const values = { topOffset: 10, bottomOffset: 20, leftOffset: 30, rightOffset: 40, TLled: 1.5, TRled: 2.5 };
    assert.equal(before.applyCalibration(1, values), 'ok');
    const snapshot = before.snapshot();

    // The board docks again with the calibration in its RAM: nothing differs from the board now.
    const config = makeConfig();
    Object.assign(config.profiles[1], values);
    const after = new OF.AppState(S);
    after.load({ version: '6.2', type: 'rpipico', arch: OF.Boards.arch('rpipico') }, config);
    assert.equal(after.isDirty(), false);
    assert.equal(after.snapshotDiffers(snapshot), true);
    assert.equal(after.restoreSnapshot(snapshot), true);
    assert.equal(after.isDirty(), true);           // Save available to write it to the flash
    assert.equal(after.cur.profiles[1].TLled, 1.5);
    after.commitDone();
    assert.equal(after.isDirty(), false);

    const other = loaded('waveshare-esp32-s3-zero');
    assert.equal(other.snapshotDiffers(snapshot), false);
    assert.equal(other.restoreSnapshot(snapshot), false);
});

test('profile names: Latin-1 only, like the product name', () => {
    const state = loaded();
    assert.equal(state.renameProfile(0, 'Schermo 1'), true);
    assert.equal(state.renameProfile(0, 'Écran café'), true);
    assert.equal(state.cur.profiles[0].name, 'Écran café');
    assert.equal(state.renameProfile(0, '画面'), false);
    assert.equal(state.renameProfile(0, ''), false);
    assert.equal(state.cur.profiles[0].name, 'Écran café');
});

test('layout import: rumble and solenoid get back their saved toggles (Qt import)', () => {
    const source = loaded();
    source.setCustomPins(true); // the default layout becomes the custom one
    assert.ok(source.pinMapped(E.rumblePin) && source.pinMapped(E.solenoidPin));
    const bytes = source.exportLayout();

    const state = loaded();
    state.setToggle(B.rumble, false);
    state.setToggle(B.solenoid, false);
    assert.equal(state.importLayout(bytes), 'ok');
    assert.equal(state.toggle(B.solenoid), true);
    assert.equal(state.toggle(B.rumble), true);
    assert.equal(state.toggle(B.rumbleFF), false);
});

test('numbers: Qt QString::number "%g" with 6 significant digits', () => {
    const g = OF.UI.formatG;
    assert.equal(g(0), '0');
    assert.equal(g(1.5), '1.5');
    assert.equal(g(-12.25), '-12.25');
    assert.equal(g(123456), '123456');
    assert.equal(g(1234567), '1.23457e+06');
    assert.equal(g(0.0001), '0.0001');
    assert.equal(g(0.00001234), '1.234e-05');
    assert.equal(g(2.5000001), '2.5');
    assert.equal(g(NaN), 'nan');
    assert.equal(g(-Infinity), '-inf');
    assert.equal(g(Math.fround(0.1)), '0.1');
});
