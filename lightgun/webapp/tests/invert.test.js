/*  OpenFIRE Web App: inverted inputs (normally closed contacts), settings[inputInvertMask].

    Run from lightgun/webapp:   node --test --test-concurrency=1 tests/*.test.js
*/
'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const { MockFirmware, MemoryLink, createMemoryTransport } = require('./mock-firmware.js');

const WEBAPP = path.join(__dirname, '..');
const load = (relative) => vm.runInThisContext(fs.readFileSync(path.join(WEBAPP, relative), 'utf8'), { filename: relative });

load('boards/OpenFIREshared.js');
vm.runInThisContext('globalThis.OpenFIREshared = OpenFIREshared;');
load('js/core/boards.js');
load('js/core/transport.js');
load('js/core/protocol.js');
load('js/app/maps.js');
load('js/app/state.js');

const OF = globalThis.OF;
const S = globalThis.OpenFIREshared;
const E = S.boardInputs_e;
const B = S.boolTypes_e;
const T = S.settingsTypes_e;

if (!process.env.OF_TEST_VERBOSE) {
    console.warn = () => {};
    console.debug = () => {};
}

function loaded({ customPins = true, mask = 0 } = {}) {
    const codec = new OF.ProtocolCodec(S);
    const config = codec.createConfig();
    config.toggles[B.customPins] = customPins;
    config.pins = S.boardsPresetsMap.rpipico.reduce((pins, fn, gpio) => {
        if (fn > E.btnUnmapped) pins[fn] = gpio;
        return pins;
    }, new Array(E.boardInputsCount).fill(-1));
    config.settings[T.customLEDcount] = 1;
    config.settings[T.tempWarning] = 42;
    config.settings[T.tempShutdown] = 50;
    config.settings[T.inputInvertMask] = mask;
    for (let i = 0; i < 4; ++i) config.profiles.push(Object.assign(codec.createProfile(), { name: `P${i + 1}` }));
    config.tinyUSB = { id: 1, name: 'FIRECon P1' };
    const state = new OF.AppState(S);
    state.load({ version: '7.0', type: 'rpipico', arch: OF.Boards.arch('rpipico'), cameraError: false }, config);
    return state;
}

test('shared data: invertibleInputs and the InputInvertMask setting come from OpenFIREshared.h', () => {
    assert.deepEqual(S.invertibleInputs, [
        E.btnTrigger, E.btnGunA, E.btnGunB, E.btnStart, E.btnSelect, E.btnGunC,
        E.btnGunUp, E.btnGunDown, E.btnGunLeft, E.btnGunRight,
        E.btnPedal, E.btnPedal2, E.btnPump, E.btnHome,
        E.rumbleSwitch, E.solenoidSwitch, E.autofireSwitch]);
    assert.equal(S.settingsTypes_Strings.InputInvertMask, T.inputInvertMask);
    assert.ok(S.invertibleInputs.every((fn) => fn >= 0 && fn < 32), 'every bit fits the 32-bit setting');
});

test('invertible: only the digital inputs, only with custom pins', () => {
    const state = loaded();
    for (let fn = -1; fn < E.boardInputsCount; ++fn)
        assert.equal(state.invertible(fn), S.invertibleInputs.includes(fn), `function ${fn}`);
    state.setCustomPins(false);
    assert.ok(S.invertibleInputs.every((fn) => !state.invertible(fn)), 'no checkbox with the default layout');
});

test('nothing inverted: mask 0 loads clean, nothing to save', () => {
    const state = loaded();
    assert.equal(state.isDirty(), false);
    assert.ok(S.invertibleInputs.every((fn) => !state.inverted(fn)));
});

test('ticking and unticking: one bit per function, Save follows the changes', () => {
    const state = loaded();
    state.setInverted(E.btnTrigger, true);
    assert.equal(state.cur.settings[T.inputInvertMask], 1 << E.btnTrigger);
    assert.equal(state.isDirty(), true, 'Save enabled');
    state.setInverted(E.autofireSwitch, true);
    assert.equal(state.cur.settings[T.inputInvertMask], (1 << E.btnTrigger) | (1 << E.autofireSwitch));
    assert.equal(state.inverted(E.autofireSwitch), true);
    assert.equal(state.inverted(E.btnGunA), false);
    state.setInverted(E.btnTrigger, false);
    state.setInverted(E.autofireSwitch, false);
    assert.equal(state.cur.settings[T.inputInvertMask], 0);
    assert.equal(state.isDirty(), false, 'back to the saved value: nothing to save');
});

test('saved inversions are read back; the inversion follows the function to another GPIO', () => {
    const state = loaded({ mask: 1 << E.btnPedal });
    assert.equal(state.inverted(E.btnPedal), true);
    assert.equal(state.isDirty(), false);
    const free = [...Array(state.pinCount).keys()].find((g) => state.functionAt(g) === -1 && !state.functionDisabled(g, E.btnPedal));
    state.changePin(free, E.btnPedal + 1);
    assert.equal(state.functionAt(free), E.btnPedal);
    assert.equal(state.inverted(E.btnPedal), true);
});

test('commitDone: the edited mask becomes the saved one', () => {
    const state = loaded();
    state.setInverted(E.btnHome, true);
    state.commitDone();
    assert.equal(state.isDirty(), false);
    assert.equal(state.orig.settings[T.inputInvertMask], 1 << E.btnHome);
});

test('protocol: the mask reaches the firmware and comes back unchanged', async () => {
    const link = new MemoryLink();
    const firmware = new MockFirmware(S, link);
    firmware.start();
    const protocol = new OF.Protocol({ shared: S });
    let transport = createMemoryTransport(OF, link);
    try {
        await protocol.connect(transport);
        const first = await protocol.getSettings();
        assert.equal(first.ok, true);
        assert.equal(first.config.settings[T.inputInvertMask], 0, 'default: nothing inverted');
        const current = protocol.codec.cloneConfig(first.config);
        const mask = ((1 << E.btnTrigger) | (1 << E.btnPedal) | (1 << E.solenoidSwitch)) >>> 0;
        current.settings[T.inputInvertMask] = mask;
        const saved = await protocol.commitSettings(current, first.config);
        assert.equal(saved.ok, true);
        assert.equal(firmware.getSetting(T.inputInvertMask), mask);

        await protocol.disconnect();
        transport = createMemoryTransport(OF, link);
        await protocol.connect(transport);
        const again = await protocol.getSettings();
        assert.equal(again.ok, true);
        assert.equal(again.config.settings[T.inputInvertMask], mask);
    } finally {
        link.connected = false;
        await firmware.stop();
        if (transport.isOpen) await transport.close();
    }
});

test('.ofl layouts: the inversions travel with the file and replace the ones already set', () => {
    const source = loaded();
    const free = S.boardsPresetsMap.rpipico.indexOf(E.btnUnmapped);
    source.changePin(free, E.rumbleSwitch + 1);
    source.setInverted(E.btnTrigger, true);
    source.setInverted(E.rumbleSwitch, true);
    source.setInverted(E.btnPedal2, true);   // not on any GPIO: not written in the file
    const bytes = source.exportLayout();

    const target = loaded({ mask: (1 << E.btnGunA) >>> 0 });
    assert.equal(target.importLayout(bytes), 'ok');
    assert.deepEqual(target.cur.pins, source.cur.pins);
    assert.equal(target.inverted(E.btnTrigger), true);
    assert.equal(target.inverted(E.rumbleSwitch), true);
    assert.equal(target.inverted(E.btnGunA), false, 'the target inversion not in the file is cleared');
    assert.equal(target.inverted(E.btnPedal2), false);

    // A file without NC records (older files, nothing inverted) clears every inversion.
    const plain = loaded();
    const again = loaded({ mask: (1 << E.btnTrigger) >>> 0 });
    assert.equal(again.importLayout(plain.exportLayout()), 'ok');
    assert.equal(again.cur.settings[T.inputInvertMask], 0);
});

test('.ofl layouts: older readers (Qt app, previous WebApp) still read the pins of a file with NC records', () => {
    const source = loaded();
    S.invertibleInputs.forEach((fn) => source.setInverted(fn, true));
    const bytes = source.exportLayout();
    const expected = new Array(source.pinCount).fill(E.btnUnmapped);
    source.cur.pins.forEach((gpio, fn) => { if (gpio > E.btnUnmapped) expected[gpio] = fn; });

    // Qt on_actionImport_Custom_Layout_triggered: read up to and including NUL; known name -> read 1 byte.
    const qt = new Array(source.pinCount).fill(E.btnUnmapped);
    let at = bytes.indexOf(0x0A) + 1, rounds = 0;
    while (at < bytes.length) {
        assert.ok(++rounds < 1000, 'the Qt loop ends');
        const nul = bytes.indexOf(0, at);
        const name = String.fromCharCode(...bytes.subarray(at, nul < 0 ? at : nul));
        at = nul < 0 ? at : nul + 1;
        const fn = S.boardInputs_Strings[name];
        if (fn !== undefined) qt[bytes[at++]] = fn;
    }
    assert.deepEqual(qt, expected);

    // Previous WebApp importLayout: (name, NUL, byte) records, unknown names skipped.
    const web = new Array(source.pinCount).fill(E.btnUnmapped);
    at = bytes.indexOf(0x0A) + 1;
    while (at < bytes.length) {
        const nul = bytes.indexOf(0, at);
        if (nul < 0 || nul + 1 >= bytes.length) break;
        const fn = S.boardInputs_Strings[String.fromCharCode(...bytes.subarray(at, nul))];
        if (fn !== undefined) web[bytes[nul + 1]] = fn;
        at = nul + 2;
    }
    assert.deepEqual(web, expected);
});

test('.ofl layouts: a truncated or damaged file is refused and changes nothing', () => {
    const source = loaded();
    source.changePin(S.boardsPresetsMap.rpipico.indexOf(E.btnUnmapped), E.rumbleSwitch + 1);
    source.setInverted(E.btnTrigger, true);
    source.setInverted(E.rumbleSwitch, true);
    const bytes = source.exportLayout();
    const header = bytes.indexOf(0x0A) + 1;
    const text = (from, to) => String.fromCharCode(...bytes.subarray(from, to));
    assert.equal(text(bytes.length - 15, bytes.length - 6), 'Checksum\0', 'checksum record before the end');
    assert.equal(text(bytes.length - 5, bytes.length), 'End\0\0', 'the file ends with the end record');

    const refused = (data, why) => {
        const target = loaded({ mask: (1 << E.btnGunA) >>> 0 });
        const before = JSON.stringify(target.cur);
        assert.equal(target.importLayout(data), 'corrupt', why);
        assert.equal(JSON.stringify(target.cur), before, `${why}: nothing changed`);
        assert.equal(target.isDirty(), false, why);
    };
    // Every cut, between records too.
    for (let length = header; length < bytes.length; ++length) refused(bytes.slice(0, length), `cut at ${length}`);
    // Every single byte changed after the board line.
    for (let at = header; at < bytes.length; ++at)
        for (const flip of [0x01, 0x80, 0xFF]) {
            const damaged = bytes.slice();
            damaged[at] ^= flip;
            refused(damaged, `byte ${at} ^ ${flip}`);
        }
    refused(Uint8Array.from([...bytes, 0]), 'bytes after the end record');
    // A file of the previous format (no checksum and end records).
    refused(bytes.slice(0, bytes.length - 15), 'previous format');

    // The whole file still imports, the other board is still told apart.
    assert.equal(loaded().importLayout(bytes), 'ok');
    assert.equal(loaded().importLayout(Uint8Array.from([...'other\n'].map((c) => c.charCodeAt(0)))), 'mismatch');
});

test('.ofl layouts: records a later version may add are skipped', () => {
    const source = loaded();
    source.setInverted(E.btnTrigger, true);
    const bytes = Array.from(source.exportLayout());
    const body = bytes.slice(0, bytes.length - 15);   // without the checksum and end records
    for (const c of 'Future Record') body.push(c.charCodeAt(0));
    body.push(0, 42);
    // Same CRC-8 as the WebApp (polynomial 0x07).
    let crc = 0;
    for (const byte of body) { crc ^= byte; for (let bit = 0; bit < 8; ++bit) crc = (crc & 0x80) ? ((crc << 1) ^ 0x07) & 0xFF : (crc << 1) & 0xFF; }
    const file = Uint8Array.from([...body, ...[...'Checksum'].map((c) => c.charCodeAt(0)), 0, crc, ...[...'End'].map((c) => c.charCodeAt(0)), 0, 0]);
    const target = loaded();
    assert.equal(target.importLayout(file), 'ok');
    assert.deepEqual(target.cur.pins, source.cur.pins);
    assert.equal(target.inverted(E.btnTrigger), true);
});

test('taking a function off its pin clears its inversion; moving it keeps it', () => {
    const state = loaded({ mask: ((1 << E.btnTrigger) | (1 << E.btnGunA) | (1 << E.btnGunB)) >>> 0 });
    // Unmapped: cleared, and ticking it again is needed.
    state.changePin(state.pinOf(E.btnTrigger), 0);
    assert.equal(state.inverted(E.btnTrigger), false);
    // Replaced by another function: the one taken off is cleared, the new one is not inverted.
    const gpioA = state.pinOf(E.btnGunA);
    state.changePin(gpioA, E.btnHome + 1);
    assert.equal(state.inverted(E.btnGunA), false);
    assert.equal(state.inverted(E.btnHome), false);
    // Moved to another pin: kept.
    const free = [...Array(state.pinCount).keys()].find((g) => state.functionAt(g) === -1 && !state.functionDisabled(g, E.btnGunB));
    state.changePin(free, E.btnGunB + 1);
    assert.equal(state.inverted(E.btnGunB), true);
    // Mapping the function again does not bring the old inversion back.
    state.changePin(gpioA, E.btnGunA + 1);
    assert.equal(state.inverted(E.btnGunA), false);
});

test('presets and Use Custom Pins: the pins change, so the inversions are cleared', () => {
    const mask = ((1 << E.btnTrigger) | (1 << E.btnPedal)) >>> 0;
    const state = loaded({ mask });
    state.setCustomPins(false);
    assert.equal(state.cur.settings[T.inputInvertMask], 0, 'default layout: nothing inverted');
    state.setCustomPins(true);
    assert.equal(state.cur.settings[T.inputInvertMask], mask, 'back to the saved custom layout: saved inversions');
    assert.equal(state.isDirty(), false);

    const plain = loaded({ customPins: false });
    plain.setCustomPins(true);
    assert.equal(plain.cur.settings[T.inputInvertMask], 0);

    if (state.altPresets.length) {
        state.applyAltPreset(0);
        assert.ok(S.invertibleInputs.every((fn) => !state.inverted(fn)), 'preset: nothing inverted');
    }
});

test('.ofl layouts: records in function order, then the NC records in the same order', () => {
    const state = loaded();
    S.invertibleInputs.forEach((fn) => state.setInverted(fn, true));
    const bytes = state.exportLayout();
    const names = [];
    let at = bytes.indexOf(0x0A) + 1;
    while (at < bytes.length) {
        const nul = bytes.indexOf(0, at);
        names.push(String.fromCharCode(...bytes.subarray(at, nul)));
        at = nul + 2;
    }
    const byValue = [];
    for (const [name, value] of Object.entries(S.boardInputs_Strings)) byValue[value] = name;
    const mapped = [...Array(E.boardInputsCount).keys()].filter((fn) => state.pinMapped(fn));
    assert.deepEqual(names, [
        ...mapped.map((fn) => byValue[fn]),
        ...mapped.filter((fn) => S.invertibleInputs.includes(fn)).map((fn) => 'NC:' + byValue[fn]),
        'Checksum', 'End']);
});

test('.ofl layouts: a big file that is not a layout is refused without errors', () => {
    const big = new Uint8Array(4 * 1024 * 1024).fill(0x41);   // no newline
    assert.equal(loaded().importLayout(big), 'mismatch');
    const longLine = big.slice(); longLine[longLine.length - 1] = 0x0A;
    assert.equal(loaded().importLayout(longLine), 'mismatch');
    const header = [...'rpipico\n'].map((c) => c.charCodeAt(0));
    const longName = new Uint8Array(header.length + big.length + 2);
    longName.set(header); longName.set(big, header.length);
    assert.equal(loaded().importLayout(longName), 'corrupt');
});

test('.ofl layouts: an NC record of a function with no pin is ignored', () => {
    const source = loaded();
    source.changePin(source.pinOf(E.btnHome), 0);
    const body = Array.from(source.exportLayout()).slice(0, -15);
    for (const c of 'NC:Home Button') body.push(c.charCodeAt(0));
    body.push(0, 0);
    let crc = 0;
    for (const byte of body) { crc ^= byte; for (let bit = 0; bit < 8; ++bit) crc = (crc & 0x80) ? ((crc << 1) ^ 0x07) & 0xFF : (crc << 1) & 0xFF; }
    const file = Uint8Array.from([...body, ...[...'Checksum'].map((c) => c.charCodeAt(0)), 0, crc, ...[...'End'].map((c) => c.charCodeAt(0)), 0, 0]);
    const target = loaded();
    assert.equal(target.importLayout(file), 'ok');
    assert.equal(target.pinMapped(E.btnHome), false);
    assert.equal(target.inverted(E.btnHome), false);
});
