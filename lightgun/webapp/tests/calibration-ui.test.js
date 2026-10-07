/* Calibration UI regressions: no device is needed; run with node --test. */
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const source = fs.readFileSync(path.join(__dirname, '../js/app/fullscreen.js'), 'utf8');

function setup(language = 'en') {
    let now = 1000000, reduce = false, id = 0;
    const timers = new Map(), frames = new Map();
    const translations = JSON.parse(fs.readFileSync(path.join(__dirname, '../lang/', language + '.json'), 'utf8'));
    const noop = () => {};
    const root = {
        OF: { i18n: { currentLang: language, t: (key) => translations[key] || key } },
        document: { removeEventListener: noop, getElementById: () => null },
        location: { search: '' }, Date: { now: () => now },
        matchMedia: () => ({ matches: reduce }), removeEventListener: noop,
        Image: class { constructor() { this.complete = false; this.naturalWidth = 0; } },
        setTimeout(fn, delay) { const token = ++id; timers.set(token, { fn, delay }); return token; },
        clearTimeout(token) { timers.delete(token); },
        requestAnimationFrame(fn) { const token = ++id; frames.set(token, fn); return token; },
        cancelAnimationFrame(token) { frames.delete(token); },
    };
    vm.runInNewContext(source, root);
    const win = new root.OF.FullscreenWindow('calibrate');
    win.width = 1280; win.height = 800;
    function draw() {
        const arcs = [];
        const ctx = new Proxy({
            arc(x, y, radius) { arcs.push({ x, y, radius }); },
            fill() { arcs.at(-1).fill = this.fillStyle; },
            stroke() { arcs.at(-1).stroke = this.strokeStyle; },
        }, { get: (obj, key) => key in obj ? obj[key] : noop });
        win._drawCalibrationProgress(ctx, 188);
        return arcs;
    }
    function drawRing(color) {
        const calls = { arcs: [], translations: [], rotations: [], dashes: [], strokes: [], fills: 0, saved: 0, restored: 0 };
        const ctx = new Proxy({
            arc(...args) { calls.arcs.push(args); },
            translate(...args) { calls.translations.push(args); },
            rotate(angle) { calls.rotations.push(angle); },
            setLineDash(dash) { calls.dashes.push(Array.from(dash)); },
            stroke() { calls.strokes.push({ color: this.strokeStyle, width: this.lineWidth, cap: this.lineCap }); },
            fill() { ++calls.fills; },
            save() { ++calls.saved; }, restore() { ++calls.restored; },
            createRadialGradient() { throw new Error('a target ring must not draw a halo'); },
        }, { get: (obj, key) => key in obj ? obj[key] : noop });
        win._drawTargetRing(ctx, color);
        return calls;
    }
    return { win, timers, frames, draw, drawRing, time: (value) => { now = value; }, reduced: (value) => { reduce = value; } };
}

test('six dots track every target and verification, independently of IR quality', () => {
    for (const [width, height, scale] of [[800, 600, 1], [1280, 800, 1], [1920, 1080, 1.5], [2560, 1440, 2]]) {
        const s = setup(); s.win.width = width; s.win.height = height;
        for (let stage = 0; stage <= 6; ++stage) {
            s.win.setStage(stage);
            const arcs = s.draw();
            assert.equal(arcs.length, 6); // progress never draws a shot-confirmation ring
            arcs.slice(0, 6).forEach((dot, i) => {
                assert.equal(dot.x, width / 2 + (i - 2.5) * 22 * scale);
                assert.equal(dot.y, 188 - 52 * scale);
                assert.equal(dot.radius, (i === stage ? 8 : 6) * scale);
                assert.equal(dot.fill, i === stage ? '#00dc00' : i < stage ? 'rgb(160,160,164)' : undefined);
                assert.equal(dot.stroke, i > stage ? 'rgb(160,160,164)' : undefined);
            });
        }
    }
});

test('only consecutive forward stages acknowledge a shot; no duplicate, reset, skipped or refused acknowledgment', () => {
    const s = setup(), w = s.win;
    w.setStage(1); const shot = w.acceptedShot;
    assert.equal(shot.stage, 0);
    s.time(1000100); w.setStage(1);
    assert.equal(w.acceptedShot, shot); // no restart
    w.showIrWarning(1);
    assert.equal(w.stage, 1); assert.equal(w.acceptedShot, shot);
    w.setStage(0); assert.equal(w.acceptedShot, null);
    w.setStage(3); assert.equal(w.acceptedShot, null);
    w.setStage(4); assert.equal(w.acceptedShot.stage, 3);
    w.setStage(5); assert.equal(w.acceptedShot.stage, 4); // only the most recent ring
    w.setStage(-1); w.setStage(99); assert.equal(w.stage, 5);
    w.setStage(0); assert.equal(w.acceptedShot, null);
    assert.deepEqual(JSON.parse(JSON.stringify(w.values)), { topOffset: -1, bottomOffset: -1, leftOffset: -1, rightOffset: -1, TLled: -1, TRled: -1 });
});

test('confirmation expires during verification, including reduced motion', () => {
    for (const reduce of [false, true]) {
        const s = setup(); s.reduced(reduce); s.win.setStage(5); s.win.setStage(6);
        assert.equal(s.draw().length, 6);
        assert.equal(s.drawRing().strokes.at(-1).color, 'rgba(255,255,255,1.000)');
        s.time(1000080);
        assert.equal(s.drawRing().strokes.at(-1).color, reduce ? 'rgba(255,255,255,1.000)' : 'rgba(255,255,255,0.500)');
        assert.equal(s.timers.size, 1); // repeated redraw does not add timers
        s.time(1000120); assert.equal(s.drawRing().arcs.length, 0); assert.equal(s.win.acceptedShot, null);
        assert.equal(s.timers.size, 0); assert.equal(s.win._shotTimer, 0);
    }
});

test('all six calibration values retain their types and values; final info does not duplicate the shot', () => {
    for (const language of ['en', 'it']) {
        const s = setup(language), w = s.win;
        w.setStage(5);
        const values = [-25, 36, -47, 58, 1234.5, 2345.75];
        values.forEach((value, index) => {
            const data = new Uint8Array(5); data[0] = index + 1;
            const view = new DataView(data.buffer);
            if (index < 4) view.setInt32(1, value, true); else view.setFloat32(1, value, true);
            w.setInfo(data);
        });
        assert.deepEqual(Object.values(w.values), values);
        assert.equal(w.stage, 6); assert.equal(w.acceptedShot.stage, 5);
        const shot = w.acceptedShot; s.time(1000100); w.setStage(6); assert.equal(w.acceptedShot, shot);
        assert.equal(w.infoPrefixes[0].trim(), language === 'it' ? 'Offset superiore:' : 'Top Offset:');
        assert.equal(new Set(w.infoPrefixes.map((text) => text.length)).size, 1);
        w.mouse = { x: -20, y: 900 };
        assert.deepEqual(Array.from(w._crosshairPosition()), [640, 400]); // held until the final flash ends
        s.time(1000120); s.drawRing();
        assert.deepEqual(Array.from(w._crosshairPosition()), [-20, 900]); // latest outside-screen aiming retained
    }
});

test('panel adds a translated central label for Square/Diamond, preserving emitter geometry and missing marks', () => {
    for (const language of ['en', 'it']) for (const diamond of [false, true]) {
        const s = setup(language), w = s.win;
        w.options.diamond = diamond; w.coordsTime = 1000000;
        w.coords = [0, 0, 0, 0, 1, 0, 0, 0]; // BL / bottom missing
        const emitters = [], curves = [], labels = [];
        w._drawEmitter = (...args) => emitters.push(args);
        const ctx = new Proxy({ quadraticCurveTo(...args) { curves.push(args); },
            fillText(text, x, y) { labels.push({ text, x, y, color: this.fillStyle, align: this.textAlign, font: this.font }); }
        }, { get: (o, k) => k in o ? o[k] : () => {} });
        w._drawCaliIrPanel(ctx, 960);
        const side = 1280 * 0.14, margin = 24, x = 1280 - margin - side, y = 800 - margin - side;
        const d = side * 0.25;
        const places = diamond ? [[side / 2, d], [d, side / 2], [side / 2, side - d], [side - d, side / 2]] :
            [[d, d], [side - d, d], [d, side - d], [side - d, side - d]];
        assert.equal(curves.length, 4); assert.equal(emitters.length, 4);
        assert.deepEqual(labels.map(label => label.text), language === 'it' ? ['LED', 'IR'] : ['IR', 'LEDs']);
        labels.forEach((label, i) => {
            assert.equal(label.x, x + side / 2);
            assert.equal(label.y, y + side / 2 + (i * 16 - 8) * side / 268.8);
            assert.equal(label.align, 'center'); assert.equal(label.color, '#dedede');
            assert.ok(label.font.includes('Segoe UI'));
        });
        emitters.forEach((args, i) => {
            assert.equal(args[1], x + places[i][0]); assert.equal(args[2], y + places[i][1]);
            assert.equal(args[5].seen, i !== 2);
            if (i === 2) { assert.equal(args[4], '#ff3030'); assert.equal(args[3], side * 0.09 / 25); }
        });
    }
});

test('raised progress dots remain below the upper target at the supported desktop scales', () => {
    for (const [width, height] of [[800, 600], [1280, 720], [1920, 1080], [2560, 1440]]) {
        const s = setup(), w = s.win;
        w.width = width; w.height = height; w.stage = 1;
        const stageTop = height * .25 - 8 * w.textScale('heading') / 2;
        const arcs = [], noop = () => {};
        const ctx = new Proxy({ arc(x, y, r) { arcs.push({ x, y, r }); } }, { get: (o, k) => k in o ? o[k] : noop });
        w._drawCalibrationProgress(ctx, stageTop);
        const radius = 24.42 * w.textScale('crosshair') * 48 / 39;
        const stroke = 2.1 * w.textScale('crosshair') / 2;
        assert.equal(arcs.length, 6);
        assert(arcs.every(dot => dot.y - dot.r > radius + stroke + 10));
    }
});

test('shutdown cancels confirmation and queued rendering; late callbacks cannot rearm animation', () => {
    const s = setup(), w = s.win;
    w.canvas = {}; w.overlay = { remove() {} };
    w._draw = () => { throw new Error('drawing after shutdown'); };
    w.setStage(1); s.drawRing();
    assert.equal(s.frames.size, 1); assert.equal(s.timers.size, 1);
    const queuedFrame = [...s.frames.values()][0], queuedTimer = [...s.timers.values()][0].fn;
    w.shutdown();
    assert.equal(s.frames.size, 0); assert.equal(s.timers.size, 0); assert.equal(w.acceptedShot, null);
    queuedFrame(); queuedTimer();
    assert.equal(s.frames.size, 0); assert.equal(s.timers.size, 0);
});

test('outer ring leaves every target centre and original crosshair geometry unchanged at every scale', () => {
    for (const [width, height, scale] of [[800, 600, 2], [1280, 800, 3], [1920, 1080, 3], [2560, 1440, 4]]) {
        const s = setup(); s.win.width = width; s.win.height = height;
        const positions = [[width / 2, height / 2], [width / 2, 0], [width / 2, height], [0, height / 2], [width, height / 2], [width / 2, height / 2]];
        positions.forEach((position, stage) => {
            s.win.setStage(stage);
            s.time(1000000 + stage * 1000 + 120);
            s.drawRing(); // expire the preceding flash before inspecting the next normal target
            const ring = s.drawRing([0, 224, 0]);
            assert.deepEqual(ring.translations[0], position);
            assert.equal(ring.arcs.length, 1); assert.equal(ring.strokes.length, 1);
            assert.ok(Math.abs(ring.arcs[0][2] - 24.42 * scale * 48 / 39) < 1e-10);
            assert.equal(ring.strokes[0].color, '#00e000');
            assert.equal(ring.strokes[0].width, 2.1 * scale); assert.equal(ring.strokes[0].cap, 'round');
            const circumference = 2 * Math.PI * ring.arcs[0][2];
            assert.ok(Math.abs(ring.dashes[0][0] / circumference * 360 - 21) < 1e-10);
            assert.ok(Math.abs(ring.dashes[0][1] / circumference * 360 - 15) < 1e-10);
            assert.equal(ring.fills, 0); assert.equal(ring.saved, ring.restored);
        });
    }
});

test('rotation is constant, clockwise and time-based with an eight-second full turn', () => {
    const s = setup();
    for (const elapsed of [0, 250, 1000, 2000, 4000, 7999, 8000, 9000, 86400000]) {
        s.time(1000000 + elapsed);
        const ring = s.drawRing();
        assert.ok(Math.abs(ring.rotations[0] - (elapsed % 8000) * 2 * Math.PI / 8000) < 1e-12);
        assert.equal(ring.strokes[0].color, '#ff8758'); // old firmware / no IR data
        assert.equal(s.timers.size, 1); // no growing timer queue
    }
    const originalStart = s.win.targetTime;
    s.win.setStage(0); assert.equal(s.win.targetTime, originalStart); // duplicate does not restart
    s.win.setStage(1); s.time(s.win.acceptedShot.time + 120); s.drawRing(); assert.equal(s.drawRing().rotations[0], 0);
});

test('ring retains the supplied crosshair colour without pulsing, even with missing or weak emitters', () => {
    const s = setup();
    for (const color of [[255, 24, 24], [255, 192, 88], [168, 255, 144], [0, 224, 0]]) {
        const stroke = s.drawRing(color).strokes[0];
        s.time(1005000);
        assert.deepEqual(s.drawRing(color).strokes[0], stroke);
        assert.ok(!stroke.color.includes('rgba'));
    }
    const before = Object.values(s.win.values); s.win.showIrWarning(1);
    assert.equal(s.win.stage, 0); assert.deepEqual(Object.values(s.win.values), before);
});

test('reduced motion stops the ring timer, but still draws an unfilled static ring', () => {
    const s = setup(); s.drawRing(); assert.equal(s.timers.size, 1);
    s.reduced(true); s.time(1002000);
    const ring = s.drawRing();
    assert.equal(s.timers.size, 0); assert.equal(s.win._targetRingTimer, 0);
    assert.equal(ring.rotations[0], 0); assert.equal(ring.strokes.length, 1); assert.equal(ring.fills, 0);
    s.reduced(false); s.drawRing(); assert.equal(s.timers.size, 1);
});

test('verification stops the coloured ring but briefly confirms the last centre, including final-info fallback', () => {
    for (const fallback of [false, true]) {
        const s = setup(); s.win.setStage(5); s.drawRing(); assert.equal(s.timers.size, 1);
        if (fallback) {
            const data = new Uint8Array(5); data[0] = 6; new DataView(data.buffer).setFloat32(1, 1234.5, true);
            s.win.setInfo(data);
        } else s.win.setStage(6);
        assert.equal(s.win._targetRingTimer, 0); assert.equal(s.timers.size, 0);
        const confirmation = s.drawRing();
        assert.equal(confirmation.arcs.length, 1);
        assert.equal(confirmation.strokes[0].color, 'rgba(255,255,255,1.000)');
        assert.deepEqual(confirmation.translations, [[640, 400]]);
        s.win.mouse = { x: -40, y: 1200 }; assert.deepEqual(Array.from(s.win._crosshairPosition()), [640, 400]);
        assert.deepEqual(s.drawRing().translations, [[640, 400]]); // confirmation does not follow the cursor
        s.time(1000120); assert.equal(s.drawRing().arcs.length, 0); assert.equal(s.timers.size, 0);
        assert.deepEqual(Array.from(s.win._crosshairPosition()), [-40, 1200]);
        s.win.setStage(0); assert.equal(s.drawRing().strokes.length, 1);
    }
});

test('late target-ring callbacks cannot draw or rearm after shutdown', () => {
    const s = setup(); s.drawRing(); const callback = [...s.timers.values()][0].fn;
    s.win.overlay = { remove() {} };
    s.win.shutdown(); callback();
    assert.equal(s.drawRing().arcs.length, 0); assert.equal(s.timers.size, 0); assert.equal(s.frames.size, 0);
});

test('solid stationary flash holds the original crosshair until 120 ms, then shows the next target', () => {
    const s = setup(), positions = [[640,400], [640,0], [640,800], [0,400], [1280,400], [640,400]];
    for (let stage = 1; stage <= 6; ++stage) {
        const start = 1000000 + (stage - 1) * 1000;
        s.time(start + 400); s.win.setStage(stage);
        const ring = s.drawRing([0,224,0]);
        assert.equal(s.win.stage, stage);
        assert.deepEqual(ring.translations.at(-1), positions[stage - 1]);
        assert.deepEqual(Array.from(s.win._crosshairPosition()), positions[stage - 1]);
        assert.equal(ring.strokes.at(-1).color, 'rgba(255,255,255,1.000)');
        assert.equal(ring.rotations.length, 0); assert.deepEqual(ring.dashes, [[]]);
        assert.equal(ring.strokes.length, 1);
        s.time(start + 480);
        const faded = s.drawRing([0,224,0]);
        assert.equal(faded.strokes.at(-1).color, 'rgba(255,255,255,0.500)');
        assert.deepEqual(faded.translations.at(-1), positions[stage - 1]);
        assert.deepEqual(Array.from(s.win._crosshairPosition()), positions[stage - 1]);
        s.time(start + 520);
        const next = s.drawRing([0,224,0]);
        assert.equal(s.win.acceptedShot, null);
        assert.deepEqual(Array.from(s.win._crosshairPosition()), positions[Math.min(stage,5)]);
        assert.equal(next.strokes.length, stage <= 5 ? 1 : 0);
        if(stage <= 5) assert.equal(next.strokes[0].color, '#00e000');
    }
});

test('confirmation scales and clips at every edge like the coloured ring; resizing preserves target anchors', () => {
    const s = setup(); s.win.setStage(1); s.win.setStage(2);
    for (const [width, height] of [[800,600], [1280,800], [1920,1080], [2560,1440]]) {
        s.win.width = width; s.win.height = height;
        const ring = s.drawRing();
        assert.deepEqual(ring.translations, [[width / 2,height / 2]]); // overlapping shots keep the visible target still
        assert.deepEqual(ring.dashes, [[]]);
        assert.equal(ring.strokes[0].width, 2.1 * s.win.textScale('crosshair'));
        assert.deepEqual(Array.from(s.win._crosshairPosition()), [width/2,height/2]);
    }
});

test('fast consecutive shots replace rather than queue confirmations, and reset cancels only the old confirmation', () => {
    const s = setup(); s.win.setStage(1); s.drawRing();
    const old = s.win._shotTimer;
    s.time(1000010); s.win.setStage(2); s.drawRing();
    assert.ok(!s.timers.has(old)); assert.equal(s.timers.size, 1);
    assert.equal(s.win.acceptedShot.stage, 1);
    assert.equal(s.win.acceptedShot.positionStage, 0);
    s.win.setStage(0); assert.equal(s.win.acceptedShot, null); assert.equal(s.win._shotTimer, 0);
    assert.equal(s.drawRing().strokes.length, 1); assert.equal(s.timers.size, 1);
});

test('late white-confirmation callbacks cannot draw or rearm after shutdown', () => {
    const s = setup(); s.win.setStage(1); s.drawRing();
    const callbacks = [...s.timers.values()].map(timer => timer.fn);
    s.win.overlay = { remove() {} }; s.win.shutdown();
    callbacks.forEach(callback => callback());
    assert.equal(s.drawRing().arcs.length, 0); assert.equal(s.timers.size, 0); assert.equal(s.frames.size, 0);
    assert.equal(s.win._shotTimer, 0); assert.equal(s.win.acceptedShot, null);
});

test('flash is fully white for 40 ms and fades for 80 ms, with no dash or rotation at any point', () => {
    const s = setup(); s.win.setStage(1);
    for(const [elapsed, alpha] of [[0,1], [39,1], [40,1], [60,0.75], [80,0.5], [100,0.25], [119,0.0125]]) {
        s.time(1000000 + elapsed);
        const flash = s.drawRing();
        assert.deepEqual(flash.translations, [[640,400]]);
        assert.deepEqual(flash.dashes, [[]]); assert.equal(flash.rotations.length, 0);
        assert.equal(flash.strokes[0].color, `rgba(255,255,255,${alpha.toFixed(3)})`);
        assert.deepEqual(Array.from(s.win._crosshairPosition()), [640,400]);
    }
    s.time(1000120); s.drawRing();
    assert.deepEqual(Array.from(s.win._crosshairPosition()), [640,0]);
    assert.equal(s.win._shotTimer, 0);
});

test('an expired flash cannot hold an older target when a new stage arrives before its redraw', () => {
    const s = setup(); s.win.setStage(1);
    s.time(1000120); s.win.setStage(2);
    assert.equal(s.win.acceptedShot.positionStage, 1);
    assert.deepEqual(Array.from(s.win._crosshairPosition()), [640,0]);
    s.time(1000240); s.drawRing();
    assert.deepEqual(Array.from(s.win._crosshairPosition()), [640,800]);
});

test('skipped or backwards stages cancel the held position immediately without changing acquired values', () => {
    const s = setup(); s.win.setStage(1); s.drawRing();
    const payload = new Uint8Array(5); payload[0] = 1; new DataView(payload.buffer).setInt32(1, -25, true);
    s.win.setInfo(payload); s.win.setStage(4);
    assert.equal(s.win.acceptedShot, null); assert.equal(s.win._shotTimer, 0);
    assert.deepEqual(Array.from(s.win._crosshairPosition()), [1280,400]); assert.equal(s.win.values.topOffset, -25);
    s.win.setStage(5); s.drawRing(); s.win.setStage(3);
    assert.equal(s.win.acceptedShot, null); assert.deepEqual(Array.from(s.win._crosshairPosition()), [0,400]);
});

test('a late redraw finishes the flash immediately; final expiry timer never extends the deadline', () => {
    const s = setup(); s.win.setStage(1); s.time(1000119); s.drawRing();
    assert.equal([...s.timers.values()][0].delay, 1);
    s.time(1000900); s.drawRing();
    assert.equal(s.win.acceptedShot, null); assert.equal(s.win._shotTimer, 0);
    assert.deepEqual(Array.from(s.win._crosshairPosition()), [640,0]);
    assert.equal(s.timers.size, 1); // just the new target's rotation
});

function drawCalibration(s) {
    const texts = [], boxes = [], panels = [];
    s.win._centered = (ctx, lines, y, scale, tint) => {
        const width = Math.max(0, ...lines.map(line => [...line].length)) * 8 * scale;
        const height = lines.length * 8 * scale;
        texts.push({ lines: Array.from(lines), x: (s.win.width - width) / 2, y, width, height, scale, tint });
        return { width, height };
    };
    s.win._centeredSegments = () => {};
    s.win._drawCaliIrPanel = (ctx, textRight) => panels.push(textRight);
    const noop = () => {};
    const ctx = new Proxy({ fillRect(...rect) { boxes.push(rect); } },
        { get: (obj, key) => key in obj ? obj[key] : noop });
    s.win._drawCalibration(ctx);
    const reminder = texts.find(text => text.lines[1] === (s.language === 'it' ? 'senza ruotare la lightgun.' : 'without rotating the lightgun.'));
    return { texts, boxes, panels, reminder };
}

test('two-line posture reminder stays above unchanged bottom instructions in every stage, in English and Italian', () => {
    for (const language of ['en', 'it']) {
        for (const [width, height] of [[1280, 720], [1920, 1080], [2560, 1440], [3840, 2160]]) {
            const s = setup(language); s.language = language; s.win.width = width; s.win.height = height;
            s.win.values = { topOffset: 0, bottomOffset: 0, leftOffset: 0, rightOffset: 0, TLled: 100, TRled: 3000 };
            for (let stage = 0; stage <= 6; ++stage) {
                s.win.stage = stage;
                const { texts, reminder } = drawCalibration(s);
                assert.deepEqual(reminder.lines, language === 'it' ?
                    ['Posizionati centralmente di fronte allo schermo,', 'senza ruotare la lightgun.'] :
                    ['Stand in front of the centre of the screen,', 'without rotating the lightgun.']);
                const tutorial = texts.find(text => text !== reminder && Math.abs(text.y + text.height / 2 - height * 0.8) < 1e-7);
                assert.ok(tutorial, `bottom instructions remain centred at 80%: stage ${stage}`);
                assert.equal(reminder.y + reminder.height + 8 * reminder.scale, tutorial.y);
                assert.equal(reminder.tint, undefined, 'posture reminder remains white');
                assert.ok(reminder.x >= 0 && reminder.x + reminder.width <= width);
            }
        }
    }
});

test('posture reminder survives refused shots and leaves a full row above the warning box', () => {
    for (const language of ['en', 'it']) {
        const s = setup(language); s.language = language; s.win.width = 1920; s.win.height = 1080;
        s.win.showIrWarning(1);
        const { boxes, reminder } = drawCalibration(s);
        assert.ok(reminder);
        const box = boxes.find(rect => rect[1] > s.win.height * 0.6);
        assert.ok(box);
        assert.ok(reminder.y + reminder.height + 8 * reminder.scale <= box[1]);
        assert.equal(s.win.stage, 0); assert.equal(s.win.acceptedShot, null);
    }
});

test('posture reminder remains white in both incomplete and malformed verification branches', () => {
    for (const values of [{ TLled: -1, TRled: -1 }, { TLled: 40000, TRled: 40000 }]) {
        const s = setup(); s.win.stage = 6; s.win.values = values;
        const { reminder } = drawCalibration(s);
        assert.ok(reminder); assert.equal(reminder.tint, undefined);
        assert.equal(s.win.values, values);
    }
});

test('static legend follows the existing panel size without requiring IR samples or adding timers', () => {
    for (const language of ['en', 'it']) {
        const s = setup(language), w = s.win;
        const nodes = ['Legend', 'IR signal intensity', 'Weak', 'Strong', 'IR point size', 'Small', 'Large',
            'IR not detected', 'The crosshair colour indicates the emitter with the weakest signal.',
            'All 4 emitters must be detected to proceed.'].map(key => ({ dataset: { legendKey: key }, textContent: '' }));
        const attributes = {};
        w.caliLegend = { style: {}, setAttribute: (key, value) => { attributes[key] = value; }, querySelectorAll: () => nodes };
        for (const [width, height] of [[1280, 720], [1920, 1080], [2560, 1440]]) {
            w.width = width; w.height = height;
            for (let stage = 0; stage <= 6; ++stage) {
                w.stage = stage;
                const values = Object.assign({}, w.values), timers = s.timers.size;
                w._updateCaliLegend(width * .8);
                const margin = 12 * w.textScale('small');
                const side = Math.max(140, Math.min(width * .14, 300, width * .2 - 2 * margin));
                assert.ok(Math.abs(Number(w.caliLegend.style.transform.slice(6, -1)) * 268.8 - side) < 1e-7);
                assert.equal(w.caliLegend.style.left, margin + 'px');
                assert.equal(w.caliLegend.style.bottom, margin + 'px');
                assert.equal(attributes.lang, language);
                assert.equal(nodes[0].textContent, language === 'it' ? 'Legenda' : 'Legend');
                assert.equal(nodes[1].textContent, language === 'it' ? 'Intensità segnale IR' : 'IR signal intensity');
                assert.equal(nodes[7].textContent, language === 'it' ? 'IR non rilevato' : 'IR not detected');
                assert.equal(s.timers.size, timers);
                assert.deepEqual(JSON.parse(JSON.stringify(w.values)), values);
                assert.equal(w.coords, null);
            }
        }
    }
});

test('the legend is not updated in alignment or IR-test mode', () => {
    for (const mode of ['alignment', 'irtest']) {
        const s = setup(), w = s.win;
        w.mode = mode;
        w.overlay = { getBoundingClientRect: () => ({ width: 1920, height: 1080 }) };
        w.canvas = { width: 1920, height: 1080, getContext: () => ({ setTransform() {} }) };
        let drawn = false;
        w._drawAlignment = w._drawIRTest = () => { drawn = true; };
        w._updateCaliLegend = () => { throw Error('legend must remain calibration-only'); };
        w._draw();
        assert.equal(drawn, true);
    }
});
