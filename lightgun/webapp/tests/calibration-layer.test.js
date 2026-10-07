/* Calibration panel stacking regressions. Run with node --test. */
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const source = fs.readFileSync(path.join(__dirname, '../js/app/fullscreen.js'), 'utf8');

function setup(language, stage) {
    const translations = JSON.parse(fs.readFileSync(path.join(__dirname, '../lang/', language + '.json'), 'utf8'));
    const root = { OF: { i18n: { currentLang: language, t: key => translations[key] || key } },
        document: {}, location: { search: '' }, Date,
        Image: class { constructor() { this.complete = true; this.naturalWidth = 61; } },
        matchMedia: () => ({ matches: false }), setTimeout: () => 1, clearTimeout: () => {} };
    vm.runInNewContext(source, root);
    const win = new root.OF.FullscreenWindow('calibrate');
    win.stage = stage; win.width = 1920; win.height = 1080;
    win._centered = () => ({ width: 0, height: 0 });
    const calls = [];
    win._drawCaliIrPanel = () => calls.push(['ir-panel']);
    const progress = win._drawCalibrationProgress;
    win._drawCalibrationProgress = function(...args) { calls.push(['progress']); return progress.apply(this, args); };
    win._updateCaliLegend = () => calls.push(['legend-position']);
    const ctx = new Proxy({}, { get: (obj, key) => key in obj ? obj[key] : (...args) => calls.push([key, ...args]),
        set: (obj, key, value) => { obj[key] = value; return true; } });
    return { win, calls, ctx };
}

test('IR panel is painted before the crosshair, in all six shots and verification, in EN/IT', () => {
    for (const language of ['en', 'it']) for (let stage = 0; stage <= 6; ++stage) {
        const { win, calls, ctx } = setup(language, stage);
        const [x, y] = Array.from(win._crosshairPosition());
        const size = 61.44 * win.textScale('crosshair');
        const values = JSON.stringify(win.values);
        win._drawCalibration(ctx);
        assert.deepEqual(calls.filter(c => c[0] === 'clearRect'), [['clearRect', 0, 0, 1920, 1080]]);
        assert(!calls.some(c => c[0] === 'fillRect'));
        const panel = calls.findIndex(c => c[0] === 'ir-panel');
        const crosshair = calls.findIndex(c => c[0] === 'drawImage');
        assert(panel >= 0 && crosshair > panel);
        assert(calls.findIndex(c => c[0] === 'progress') < crosshair);
        assert.deepEqual(calls[crosshair].slice(2), [x - size/2, y - size/2, size, size]);
        assert.equal(win.stage, stage); assert.equal(JSON.stringify(win.values), values);
    }
});

test('verification coordinates remain unchanged over either panel and near screen edges', () => {
    for (const point of [[100,900], [1800,900], [0,540], [1920,540], [960,0], [960,1080]]) {
        const { win, ctx, calls } = setup('en', 6);
        win.mouse = { x: point[0], y: point[1] };
        win._drawCalibration(ctx);
        const image = calls.find(c => c[0] === 'drawImage'), size = 61.44 * win.textScale('crosshair');
        assert.deepEqual(image.slice(2), [point[0] - size/2, point[1] - size/2, size, size]);
        assert.deepEqual(Array.from(win._crosshairPosition()), point);
    }
});

test('calibration stacking is mode-specific and preserves the fullscreen retry above the canvas', () => {
    const css = fs.readFileSync(path.join(__dirname, '../style.css'), 'utf8');
    assert.match(css, /\.fullscreen-window\.mode-calibrate\s*\{\s*background:\s*dimgray;/);
    assert.match(css, /\.mode-calibrate\s*>\s*\.calibration-legend\s*\{\s*z-index:\s*0;/);
    assert.match(css, /\.mode-calibrate\s*>\s*\.fullscreen-canvas\s*\{\s*position:\s*relative;\s*z-index:\s*1;/);
    assert.match(css, /\.mode-calibrate\s*>\s*\.fullscreen-retry\s*\{\s*z-index:\s*2;/);
});
