/*  OpenFIRE Web App: translations and board helpers.

    Run from lightgun/webapp:   node --test --test-concurrency=1 tests/*.test.js
    Requires boards/OpenFIREshared.js and lang/translations.js (python scripts/webapp_build.py sizes).
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
load('lang/translations.js');
load('js/core/i18n.js');
load('js/core/boards.js');
load('js/core/transport.js');
load('js/core/protocol.js');

const OF = globalThis.OF;

test('i18n: translation, English fallback and %n arguments', () => {
    const i18n = new OF.I18n({ en: {}, it: { 'About': 'Informazioni', 'GPIO Pin No. %1.': 'Pin GPIO n. %1.' } });
    assert.equal(i18n.t('About'), 'About');
    i18n.currentLang = 'it';
    assert.equal(i18n.t('About'), 'Informazioni');
    assert.equal(i18n.t('Not translated'), 'Not translated');
    assert.equal(i18n.t('GPIO Pin No. %1.', 12), 'Pin GPIO n. 12.');
    assert.equal(i18n.t('%1 of %2', 1), '1 of %2');
    assert.deepEqual(new OF.I18n({ it: {}, en: {}, de: {} }).languages(), ['en', 'de', 'it']);
});

test('i18n: bundled translations contain Italian and drop untranslated entries', () => {
    const it = OF.TRANSLATIONS.it;
    assert.ok(Object.keys(it).length > 100);
    for (const [source, text] of Object.entries(it))
        assert.notEqual(source, text);
});

test('boards: site build lists the compatible boards and resolves unknown ones', () => {
    const B = OF.Boards;
    assert.equal(B.isDevice, false);
    const boards = B.list();
    assert.ok(boards.includes('rpipico') && boards.includes('waveshare-esp32-s3-zero'));
    assert.ok(!boards.some((board) => board.includes('generic')));
    assert.equal(B.resolve('some-future-board'), 'generic');
    assert.equal(B.displayName('rpipico2'), 'Raspberry Pi Pico 2 (RP2350)');
    assert.equal(B.pictureName('rpipico'), 'rpipico.svg');
    assert.equal(B.pictureUrl('rpipico'), 'boards/pics/rpipico.js');
    assert.equal(B.pictureUrl('some-future-board'), 'boards/pics/generic.js');
});

test('boards: every picture is a script with the same name that registers the SVG (works from file://)', () => {
    const B = OF.Boards;
    const S = OpenFIREshared;
    for (const board of B.list()) {
        const file = path.join(WEBAPP, B.pictureUrl(board));
        assert.ok(fs.existsSync(file), `${board}: ${B.pictureUrl(board)} (run scripts/webapp_build.py)`);
        vm.runInThisContext(fs.readFileSync(file, 'utf8'), { filename: B.pictureUrl(board) });
        const svg = OF.BoardPictures[B.pictureName(board)];
        assert.ok(svg.startsWith('<svg') && svg.includes('OF_pin'), board);
    }
    assert.ok(Object.keys(S.boardImagesMap).length >= B.list().length);
});

test('boards: architecture and pin capabilities', () => {
    const B = OF.Boards;
    const S = OpenFIREshared;
    assert.equal(B.arch('esp32-s3-devkitc-1'), 'esp32-s3');
    assert.equal(B.arch('rpipico'), 'rp2040_235X');
    assert.equal(B.isEsp32('waveshare-esp32-s3-pico'), true);
    assert.equal(B.capabilities('rpipico'), S.mcuCapableMaps.rpipico);
    assert.equal(B.capabilities('waveshareZero'), S.mcuCapableMaps.rp2040_235X);
    assert.equal(B.capabilities('esp32-s3-devkitc-1'), S.mcuCapableMaps['esp32-s3']);
});

test('boards: architecture from the board name, also for architectures added to boardArchs', () => {
    const B = OF.Boards;
    const S = OpenFIREshared;
    const saved = S.boardArchs;
    S.boardArchs = [...saved, 'esp32-p4', 'esp32-s3-lp'];
    try {
        assert.equal(B.arch('devkit-esp32-p4'), 'esp32-p4');
        assert.equal(B.arch('mini-esp32-s3-lp'), 'esp32-s3-lp');
        assert.equal(B.arch('waveshare-esp32-s3-zero'), 'esp32-s3');
        assert.equal(B.arch('adafruitKB2040'), 'rp2040_235X');
        assert.equal(OF.ProtocolUtils ? OF.ProtocolUtils.boardArch(S, 'devkit-esp32-p4') : 'esp32-p4', 'esp32-p4');
    } finally {
        S.boardArchs = saved;
    }
});

test('boards: device build uses the picture inside app.js', async () => {
    const saved = OF.BoardPictures;
    OF.BUILD = { target: 'device', board: 'waveshare-esp32-s3-zero' };
    OF.BoardPictures = { 'waveshare-esp32-s3-zero.svg': '<svg id="own"/>' };
    try {
        assert.equal(OF.Boards.isDevice, true);
        assert.equal(await OF.Boards.loadPicture('waveshare-esp32-s3-zero'), '<svg id="own"/>');
    } finally {
        delete OF.BUILD;
        OF.BoardPictures = saved;
    }
});
