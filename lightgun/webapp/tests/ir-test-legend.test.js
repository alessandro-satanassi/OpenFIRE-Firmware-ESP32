/* IR-test graphics regressions. Run with node --test; no lightgun is required. */
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const source = fs.readFileSync(path.join(__dirname, '../js/app/fullscreen.js'), 'utf8');

function setup(language = 'en', mode = 'irtest') {
    const translations = JSON.parse(fs.readFileSync(path.join(__dirname, '../lang/', language + '.json'), 'utf8'));
    const root = { OF: { i18n: { currentLang: language, t: key => translations[key] || key } },
        document: {}, location: { search: '' }, Image: class {}, Date, setTimeout, clearTimeout };
    vm.runInNewContext(source, root);
    const win = new root.OF.FullscreenWindow(mode);
    win.width = 1920; win.height = 1080;
    const calls = [];
    const ctx = new Proxy({}, { get: (obj, key) => key in obj ? obj[key] : (...args) => calls.push([key, ...args]),
        set: (obj, key, value) => { calls.push(['property', key, value]); obj[key] = value; return true; } });
    return { root, win, calls, ctx };
}

test('IR aim retains a 50-pixel circle and the six calibration crosshair strokes, without animation', () => {
    const { win, ctx, calls } = setup();
    win._drawIRTestAim(ctx, 715, -42);
    assert.deepEqual(calls.filter(c => c[0] === 'translate'), [['translate', 715, -42]]);
    assert.deepEqual(calls.filter(c => c[0] === 'scale'), [['scale', 25 / 24.42, 25 / 24.42]]);
    assert.deepEqual(calls.filter(c => c[0] === 'ellipse'), [['ellipse', 0, 0, 24.42, 24.42, 0, 0, Math.PI * 2]]);
    assert.equal(calls.filter(c => c[0] === 'moveTo').length, 6);
    assert.equal(calls.filter(c => c[0] === 'lineTo').length, 6);
    assert.equal(calls.filter(c => c[0] === 'stroke').length, 2);
    assert.equal(calls.filter(c => c[0] === 'save').length, 1);
    assert.equal(calls.filter(c => c[0] === 'restore').length, 1);
    assert(!calls.some(c => ['rotate', 'setLineDash', 'fill'].includes(c[0])));
});

test('aim and emitter coordinates keep their original mapping on 16:9, 4:3 and ultrawide screens', () => {
    for (const [width, height] of [[1920, 1080], [1024, 768], [2560, 1080]]) {
        const { win, ctx, calls } = setup(); win.width = width; win.height = height;
        win.coords = [1200, 300, 2641, 300, -399, 780, 2640, 780, -15, 1100, 955, 545];
        const original = Array.from(win.coords);
        win._background = win._centered = win._updateIRTestLegend = () => {};
        let aim, polygon;
        win._drawIRTestAim = (context, x, y) => { aim = [x, y]; };
        win._polyline = (context, points) => { polygon = JSON.parse(JSON.stringify(points)); };
        win._drawIRTest(ctx);
        assert.deepEqual(aim, [-15, 1100]); assert.deepEqual(Array.from(win.coords), original);
        assert(calls.some(c => c[0] === 'scale' && c[1] === width / 1920 && c[2] === height / 1080));
        const scale = Math.min(width / 1920, height / 1080);
        const ox = (width - 1920 * scale) / 2, oy = (height - 1080 * scale) / 2;
        assert.deepEqual(polygon, [[600, 300], [1320, 300], [1320, 780], [-200, 780]].map(([x, y]) => [ox + x * scale, oy + y * scale]));
        const rings = calls.filter(c => c[0] === 'ellipse');
        assert.equal(rings.length, 5); // four classic emitters plus unchanged red centre
        assert.deepEqual(rings.at(-1).slice(1, 5), [955, 545, 25, 25]);
    }
});

test('legend is still drawn without coordinates; no artificial samples are required', () => {
    const { win, ctx } = setup(); let updated = 0;
    win._background = win._centered = () => {};
    win._updateIRTestLegend = () => ++updated;
    win._drawIRTest(ctx); assert.equal(updated, 1); assert.equal(win.coords, null);
});

test('IR test clears its transparent live layer; the opaque-background helper is unchanged', () => {
    const { win, ctx, calls } = setup();
    win._centered = win._updateIRTestLegend = () => {};
    win._drawIRTest(ctx);
    assert.deepEqual(calls.filter(c => c[0] === 'clearRect'), [['clearRect', 0, 0, 1920, 1080]]);
    assert(!calls.some(c => c[0] === 'fillRect'));
    calls.length = 0;
    win._background(ctx, 'dimgray');
    assert.deepEqual(calls.filter(c => c[0] === 'fillRect'), [['fillRect', 0, 0, 1920, 1080]]);
    const css = fs.readFileSync(path.join(__dirname, '../style.css'), 'utf8');
    assert.match(css, /\.fullscreen-window\.mode-irtest\s*\{\s*background:\s*midnightblue;/);
    assert.match(css, /\.mode-irtest\s*>\s*\.ir-test-legend\s*\{\s*z-index:\s*0;/);
    assert.match(css, /\.mode-irtest\s*>\s*\.fullscreen-canvas\s*\{\s*position:\s*relative;\s*z-index:\s*1;/);
    assert.match(css, /\.mode-irtest\s*>\s*\.fullscreen-retry\s*\{\s*z-index:\s*2;/);
});

test('IR layout setter cannot change calibration, alignment or closed windows', () => {
    for (const mode of ['calibrate', 'alignment', 'irtest']) {
        const { win } = setup('en', mode); let updated = 0;
        win._updateIRTestLegend = () => ++updated;
        win.setIRTestLayout(true);
        assert.equal(updated, mode === 'irtest' ? 1 : 0);
        win.closed = true; win.setIRTestLayout(false);
        assert.equal(updated, mode === 'irtest' ? 1 : 0);
    }
});

test('legend uses the applied profile layout and adds no protocol commands', () => {
    const main = fs.readFileSync(path.join(__dirname, '../js/app/main.js'), 'utf8');
    const start = main.indexOf('async openIRTest()');
    const body = main.slice(start, main.indexOf('/** Qt CaliWindowExiting.', start));
    assert(body.includes('diamond: this.state.orig.profiles[this.state.cur.selectedProfile].layoutType'));
    assert.equal((body.match(/this\.sendCommand\(/g) || []).length, 1);
    assert(body.includes('this.sendCommand(this.C.sIRTest, [1])'));
});

test('every static legend phrase is explicitly translated in English and Italian', () => {
    const keys = Array.from(source.matchAll(/data-ir-key="([^"]+)"/g), m => m[1]);
    assert.equal(new Set(keys).size, 13);
    for (const language of ['en', 'it']) {
        const json = JSON.parse(fs.readFileSync(path.join(__dirname, '../lang/', language + '.json'), 'utf8'));
        keys.forEach(key => { assert.equal(typeof json[key], 'string', key); assert(json[key].length); });
        assert(!json['    your aim, and the gray crosshair should be  '].includes('circle'));
    }
});

test('vertical legend keeps the lower-left anchor, sample spacing and translations on resize', () => {
    const css = fs.readFileSync(path.join(__dirname, '../style.css'), 'utf8');
    assert.match(css, /\.ir-test-legend\s*\{[^}]*width:\s*250px;\s*height:\s*540px;/);
    assert.match(css, /\.ir-test-legend \.ir-measures\s*\{[^}]*grid-template-columns:\s*1fr;/);
    assert.match(css, /\.ir-test-legend \.ir-signal-samples, \.ir-test-legend \.ir-size-samples\s*\{[^}]*width:\s*164\.5px;/);
    const keys = Array.from(source.matchAll(/data-ir-key="([^"]+)"/g), m => m[1]);
    for (const language of ['en', 'it']) {
        const { win, root } = setup(language);
        const translated = keys.map(key => ({ dataset: { irKey: key }, textContent: '' }));
        const attrs = {}, legend = { style: {}, setAttribute: (key, value) => { attrs[key] = value; },
            querySelectorAll: selector => selector === '[data-ir-key]' ? translated : [] };
        win.irTestLegend = legend;
        for (const [width, height] of [[1920, 1080], [1280, 720], [1024, 768], [800, 600], [640, 480], [2560, 1080]]) {
            win.width = width; win.height = height;
            win._updateIRTestLegend();
            const scale = Number(legend.style.transform.match(/scale\(([^)]+)\)/)[1]);
            const left = parseFloat(legend.style.left), bottom = parseFloat(legend.style.bottom);
            assert(scale > 0 && scale <= 1);
            assert.equal(left, Math.min(24, width / 80));
            assert(left + 250 * scale <= width && bottom + 540 * scale <= height);
            if (width === 1920) { assert.equal(scale, 1); assert.equal(bottom, 24); }
            assert.equal(legend.lang, language);
            assert.equal(attrs['aria-label'], root.OF.i18n.t('IR Camera Test Legend'));
            translated.forEach(el => assert.equal(el.textContent, root.OF.i18n.t(el.dataset.irKey)));
        }
    }
});
