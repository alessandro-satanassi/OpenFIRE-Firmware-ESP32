/* Lens refinement solver (js/core/lensfit.js) on simulated cameras. Run with node --test. */
'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');

const root = {};
vm.runInNewContext(fs.readFileSync(path.join(__dirname, '../js/core/lensfit.js'), 'utf8'), root);
const L = root.OF.LensFit;

// Camera records as in OpenFIRECameraProfile.h (internal camera units).
const CAMS = {
    R3: { resX: 4096, resY: 4096, record: { k1: -0.128946, k2: 0.155087, max: 1.5, min: 0.8, focal: 0.378 / 0.539 } },
    R2: { resX: 4096, resY: 4096, record: { k1: 0, k2: 0, max: 1.2, min: 0.8, focal: 1 / 0.34727 } },
};

function rng(seed) {
    let s = seed >>> 0;
    const next = () => { s = (s * 1664525 + 1013904223) >>> 0; return s / 4294967296; };
    const normal = () => Math.sqrt(-2 * Math.log(next() + 1e-12)) * Math.cos(2 * Math.PI * next());
    return { next, normal, uniform: (a, b) => a + (b - a) * next() };
}

function gain(k1, k2, q) { return 1 + k1 * q + k2 * q * q; }

/** Observed point whose firmware correction (true lens) is the ideal point u. */
function distort(cam, k, u) {
    const cx = cam.resX / 2, cy = cam.resY / 2, dx = u[0] - cx, dy = u[1] - cy;
    let lo = 0, hi = 2;
    for (let i = 0; i < 50; ++i) {
        const m = (lo + hi) / 2, q = (dx * m / cx) ** 2 + (dy * m / cy) ** 2;
        if (m * gain(k.k1, k.k2, q) < 1) lo = m; else hi = m;
    }
    const m = (lo + hi) / 2;
    return [cx + dx * m, cy + dy * m];
}

function look(fwd, roll) {
    const n = Math.hypot(...fwd), f = fwd.map((v) => v / n), up = [0, -1, 0];
    let x = [up[1] * f[2] - up[2] * f[1], up[2] * f[0] - up[0] * f[2], up[0] * f[1] - up[1] * f[0]];
    const nx = Math.hypot(...x); x = x.map((v) => v / nx);
    const y = [f[1] * x[2] - f[2] * x[1], f[2] * x[0] - f[0] * x[2], f[0] * x[1] - f[1] * x[0]];
    const c = Math.cos(roll), s = Math.sin(roll);
    return [[c * x[0] - s * y[0], c * x[1] - s * y[1], c * x[2] - s * y[2]],
        [s * x[0] + c * y[0], s * x[1] + c * y[1], s * x[2] + c * y[2]], f];
}

/** Square tracker samples (TL, TR, BL, BR) of a 55" screen, true lens/focal given. */
function simulate(cam, truth, n, noise, seed, opts = {}) {
    const r = rng(seed);
    const W = 1.21, H = 0.68;
    const layout = opts.diamond ? [[W / 2, 0], [W, H / 2], [W / 2, H], [0, H / 2]] : [[0, 0], [W, 0], [0, H], [W, H]];
    const leds = layout.map(([x, y]) => [x + r.normal() * 0.003, y + r.normal() * 0.003]);
    const f = cam.resX / 2 * truth.focal, cx = cam.resX / 2, cy = cam.resY / 2;
    const samples = [];
    let guard = 0;
    while (samples.length < n && guard++ < n * 200) {
        const pos = [W / 2 + r.uniform(-0.6, 0.6), H / 2 + r.uniform(-0.3, 0.3), -r.uniform(opts.dmin || 0.8, opts.dmax || 3)];
        const tgt = [r.uniform(-0.05, 1.05) * W, r.uniform(-0.05, 1.05) * H, 0];
        const roll = opts.roll !== undefined ? opts.roll + r.normal() * 0.01 : r.normal() * 0.25;
        const R = look([tgt[0] - pos[0], tgt[1] - pos[1], tgt[2] - pos[2]], roll);
        const pts = [];
        let ok = true;
        for (const [X, Y] of leds) {
            const d = [X - pos[0], Y - pos[1], -pos[2]];
            const pc = R.map((row) => row[0] * d[0] + row[1] * d[1] + row[2] * d[2]);
            if (pc[2] <= 0) { ok = false; break; }
            const p = distort(cam, truth, [cx + f * pc[0] / pc[2], cy + f * pc[1] / pc[2]]);
            p[0] += r.normal() * noise; p[1] += r.normal() * noise;
            // The scene above is seen upside down: turn the image half round about its centre
            // (the lens is radial, nothing else changes), so that a level gun reads level.
            p[0] = 2 * cx - p[0]; p[1] = 2 * cy - p[1];
            if (opts.mirror) p[0] = cam.resX - p[0];   // tracker X reversed (camera mounted mirrored)
            if (p[0] < 0 || p[1] < 0 || p[0] > cam.resX - 1 || p[1] > cam.resY - 1) { ok = false; break; }
            pts.push(p);
        }
        if (ok) samples.push({ diamond: !!opts.diamond, complete: true, points: pts });
    }
    return samples;
}

/** One Square view from pos (metres, z < 0 in front of the screen) aiming at a screen fraction. */
function view(cam, truth, pos, aim, roll, r, noise) {
    const W = 1.21, H = 0.68;
    const leds = [[0, 0], [W, 0], [0, H], [W, H]];
    const f = cam.resX / 2 * truth.focal, cx = cam.resX / 2, cy = cam.resY / 2;
    const R = look([aim[0] * W - pos[0], aim[1] * H - pos[1], -pos[2]], roll);
    const pts = [];
    for (const [X, Y] of leds) {
        const d = [X - pos[0], Y - pos[1], -pos[2]];
        const pc = R.map((row) => row[0] * d[0] + row[1] * d[1] + row[2] * d[2]);
        if (pc[2] <= 0) return null;
        const p = distort(cam, truth, [cx + f * pc[0] / pc[2], cy + f * pc[1] / pc[2]]);
        p[0] = 2 * cx - p[0] - r.normal() * noise; p[1] = 2 * cy - p[1] - r.normal() * noise;
        if (p[0] < 0 || p[1] < 0 || p[0] > cam.resX - 1 || p[1] > cam.resY - 1) return null;
        pts.push(p);
    }
    return { diamond: false, complete: true, points: pts };
}

function info(cam, active) { return { resX: cam.resX, resY: cam.resY, record: cam.record, active: active || cam.record }; }

test('settings bits: float round trip and the unset value', () => {
    assert.equal(L.toBits(null), 0xFFFFFFFF);
    assert.equal(L.fromBits(0xFFFFFFFF), null);
    assert.equal(L.toBits(0), 0);
    for (const v of [0, -0.128946, 0.155087, 0.7012987])
        assert(Math.abs(L.fromBits(L.toBits(v)) - v) < 1e-7);
});

test('physical check of the curve, as in the firmware', () => {
    assert(L.curveValid(-0.128946, 0.155087).ok);
    assert(L.curveValid(0, 0).ok);
    assert(!L.curveValid(-0.6, 0).ok);    // factor below zero at the corners
    assert(!L.curveValid(0.1, -0.2).ok);  // radius no longer growing
    const c = L.curveValid(-0.128946, 0.155087);
    assert(c.min > 0.97 && c.min < 1 && c.max > 1.36 && c.max < 1.37);
    // Radius growth that turns back between two of the old 41 sampled points (q = 0.775).
    assert(!L.curveValid(-0.8606451613, 0.3331529657).ok);
    // Exact extremes: the vertex of the parabola, not the nearest sample.
    const v = L.curveValid(-0.3, 0.2);   // vertex at q = 0.75: 1 - 0.1125
    assert(Math.abs(v.min - (1 - 0.3 * 0.75 + 0.2 * 0.75 * 0.75)) < 1e-6);
});

test('sLensInfo and sLensSample payloads', () => {
    const p = new Uint8Array(49), v = new DataView(p.buffer);
    p[0] = 1; v.setInt32(1, 4096, true); v.setInt32(5, 3072, true);
    [0.006, 0, 1.2, 0.8, 3.37, 0.01, 0.02, 1.2, 0.8, 3.3].forEach((x, i) => v.setFloat32(9 + 4 * i, x, true));
    const i = L.parseInfo(p);
    assert.equal(i.resX, 4096); assert.equal(i.resY, 3072);
    assert(Math.abs(i.record.k1 - 0.006) < 1e-7 && Math.abs(i.active.focal - 3.3) < 1e-6);
    assert.equal(L.parseInfo(p.subarray(0, 40)), null);
    const s = new Uint8Array(34), w = new DataView(s.buffer);
    s[0] = 1; s[1] = 3;
    for (let k = 0; k < 8; ++k) w.setInt32(2 + 4 * k, 100 * k - 50, true);
    const sample = L.parseSample(s);
    assert(sample.diamond && sample.complete);
    assert.equal(JSON.stringify(sample.points), JSON.stringify([[-50, 50], [150, 250], [350, 450], [550, 650]]));
});

test('collector: paints the screen, needs the whole sweep, ignores incomplete and repeated frames', () => {
    const cam = CAMS.R3;
    const col = new L.Collector(info(cam));
    const samples = simulate(cam, cam.record, 1500, 2, 7);
    assert.equal(col.add({ diamond: false, complete: false, points: samples[0].points }), false);
    // An incomplete view is not kept, but it moves the aim on screen (the cursor follows the gun).
    assert.equal(col.frames.length, 0);
    assert(col.aim && Number.isFinite(col.aim[0]) && Number.isFinite(col.aim[1]));
    const aim0 = col.aim.slice();
    assert.equal(col.add({ diamond: false, complete: false, points: samples[1].points }), false);
    assert.equal(col.frames.length, 0);
    assert.notDeepEqual(col.aim, aim0);
    const aim1 = col.aim.slice();
    assert.equal(col.add({ diamond: false, complete: false, points: [[-1, -1], [-1, -1], [-1, -1], [-1, -1]] }), false);
    assert.deepEqual(col.aim, aim1);   // nothing usable: the aim stays where it was
    assert.equal(col.add(samples[0]), true);
    assert.equal(col.add(samples[0]), false);   // the gun did not move
    for (const s of samples) col.add(s);
    const p = col.progress();
    assert(p.count >= 250 && p.count <= 900, String(p.count));
    assert(p.screen > 0.9, String(p.screen));
    assert(p.total > 0.8, JSON.stringify(p));
});

test('collector never gets stuck: a long straight sweep, then tilted right, then left', () => {
    const cam = CAMS.R3;
    const col = new L.Collector(info(cam));
    // Straight: the top edge of the LED frame within 6 degrees (no tilt counted at all).
    const level = (s) => Math.abs(Math.atan2(s.points[1][1] - s.points[0][1], s.points[1][0] - s.points[0][0])) < 0.1;
    for (const s of simulate(cam, cam.record, 12000, 2, 41, { roll: 0 }).filter(level)) col.add(s);
    const straight = col.progress();
    // Every screen cell and size band already full before tilting (the old collector stopped here).
    assert(straight.count >= 600 && straight.roll === 0, JSON.stringify(straight));
    for (const s of simulate(cam, cam.record, 800, 2, 42, { roll: 0.3 })) col.add(s);
    for (const s of simulate(cam, cam.record, 800, 2, 43, { roll: -0.3 })) col.add(s);
    const p = col.progress();
    assert(p.count <= 900, String(p.count));
    assert.equal(p.roll, 1, JSON.stringify(p));
    assert(p.done, JSON.stringify(p));
    // The counts stay those of the frames actually kept.
    const screen = new Array(col.screen.length).fill(0);
    for (const m of col.meta) if (m.cell >= 0) screen[m.cell]++;
    assert.equal(JSON.stringify(screen), JSON.stringify(Array.from(col.screen)));
    assert.equal(col.meta.length, col.frames.length);
});

test('guide: a user who only follows the steps on screen finishes, and all green means done', () => {
    for (const cam of [CAMS.R3, CAMS.R2]) {
        const r = rng(51);
        const col = new L.Collector(info(cam));
        let hand = 1, dist = 2.6, steps = new Set(), closer = 0, farther = 0, target = 0, lost = 0, n = 0;
        for (; n < 20000; ++n) {
            const g = col.guide();
            steps.add(g.step);
            if (g.step === 4) break;
            let aim, roll = r.normal() * 0.02;
            if (g.step === 1) {
                aim = [r.uniform(0, 1), r.uniform(0, 1)];                           // aims everywhere
                if (lost > 12) { dist = Math.min(5, dist * 1.05); lost = 0; }      // red warning often: steps back
            } else if (g.step === 2) {
                if (g.closer) { dist = Math.max(0.7, dist * 0.97); ++closer; aim = [0.5, 0.5]; }
                else if (g.farther) { dist = Math.min(5, dist * 1.03); ++farther; aim = [0.5, 0.5]; }
                else { ++target; aim = [g.target[0] + r.normal() * 0.02, g.target[1] + r.normal() * 0.02]; }
            } else {
                aim = [0.5 + r.normal() * 0.1, 0.5 + r.normal() * 0.1];
                // Tilts towards the side still missing, watching the needle: turns the other way
                // when the needle goes the wrong way.
                const want = g.left < 1 ? -1 : 1;
                if (Math.abs(g.roll) > g.limit && Math.sign(g.roll) !== want) hand = -hand;
                roll = want * hand * 0.3;
            }
            const pos = [0.605 + r.normal() * 0.05, 0.34 + r.normal() * 0.03, -dist];
            const s = view(cam, cam.record, pos, aim, roll, r, 2);
            if (s) { col.add(s); lost = Math.max(0, lost - 0.5); } else { ++lost; col.add({ diamond: false, complete: false, points: [[0, 0], [0, 0], [0, 0], [0, 0]] }); }
        }
        const p = col.progress();
        if (process.env.GUIDE_LOG) console.log(cam === CAMS.R3 ? 'R3' : 'R2', { n, closer, farther, target, dist: dist.toFixed(2), field: Array.from(col.field), count: p.count });
        assert(p.done, `${cam === CAMS.R3 ? 'R3' : 'R2'} ${n} ${JSON.stringify(p)} steps ${[...steps]} closer ${closer} farther ${farther} target ${target} dist ${dist.toFixed(2)} field ${JSON.stringify(Array.from(col.field))}`);
        assert(col.guide().step === 4);
        assert(steps.has(1));
        // Screen all green implies the frames the fit needs.
        assert(p.frames === 1 && p.screen === 1);
    }
});

test('R3: K1/K2 different from the record are recovered, error lower than with the record values', () => {
    const cam = CAMS.R3;
    const truth = { k1: cam.record.k1 + 0.05, k2: cam.record.k2 - 0.04, focal: cam.record.focal };
    const col = new L.Collector(info(cam));
    for (const s of simulate(cam, truth, 3000, 4, 11)) col.add(s);
    const res = L.solve(info(cam), col.frames);
    assert(res.ok, res.reason);
    assert(res.tier === 'k1k2' || res.tier === 'k1k2focal', res.tier);
    assert(Math.abs(res.k1 - truth.k1) < 0.02, `K1 ${res.k1}`);
    assert(Math.abs(res.k2 - truth.k2) < 0.02, `K2 ${res.k2}`);
    assert(res.rms < res.rmsRecord * 0.97, `${res.rms} vs ${res.rmsRecord}`);
    assert(L.curveValid(res.k1, res.k2).ok);
    assert.equal(res.validated, true);   // also better on the views left out of the fit
    assert.equal(res.activeCustom, false);
});

test('values already refined: kept when the new fit is not clearly better, replaced when worse', async () => {
    const cam = CAMS.R3;
    const truth = { k1: cam.record.k1 + 0.05, k2: cam.record.k2 - 0.04, focal: cam.record.focal };
    const col = new L.Collector(info(cam));
    for (const s of simulate(cam, truth, 3000, 4, 14)) col.add(s);
    const good = L.solve(info(cam, { ...cam.record, k1: truth.k1, k2: truth.k2 }), col.frames);
    assert(good.ok && good.activeCustom && good.keepActive, JSON.stringify({ tier: good.tier, rmsActive: good.rmsActive, rms: good.rms }));
    const bad = L.solve(info(cam, { ...cam.record, k1: cam.record.k1 - 0.05 }), col.frames);
    assert(bad.ok && bad.activeCustom && !bad.keepActive && bad.rmsActive > bad.rms);
});

test('good values in use and noisy views: kept, never replaced by the record values by chance', () => {
    const cam = CAMS.R3;
    const truth = { k1: cam.record.k1 + 0.05, k2: cam.record.k2 - 0.04, focal: cam.record.focal };
    for (const seed of [43, 45]) {   // with these views a zero margin went back to the record values
        const col = new L.Collector(info(cam));
        for (const s of simulate(cam, truth, 600, 32, seed)) col.add(s);
        const res = L.solve(info(cam, { ...cam.record, k1: truth.k1, k2: truth.k2 }), col.frames);
        assert(res.ok && res.activeCustom, res.reason);
        assert(!(res.tier === 'record' && !res.keepActive), `seed ${seed}: back to the record values, ${JSON.stringify({ rmsActive: res.rmsActive, rms: res.rms })}`);
    }
});

test('solveAsync gives the page back during the fit and can be cancelled', async () => {
    const cam = CAMS.R3;
    const col = new L.Collector(info(cam));
    for (const s of simulate(cam, cam.record, 1500, 4, 15)) col.add(s);
    root.setTimeout = setTimeout;
    let calls = 0;
    const cancelled = await root.OF.LensFit.solveAsync(info(cam), col.frames, { onProgress: () => ++calls, isCancelled: () => calls >= 2 });
    assert.equal(cancelled, null);
    const progress = [];
    const res = await root.OF.LensFit.solveAsync(info(cam), col.frames, { onProgress: (p) => progress.push(p) });
    assert(res.ok);
    assert(progress.length > 3 && progress.every((p, i) => p >= 0 && p <= 1 && (i === 0 || p >= progress[i - 1])));
    assert.equal(JSON.stringify(res.k1), JSON.stringify(L.solve(info(cam), col.frames).k1));
});

test('R3: a wrong nominal focal ratio is found when the data allow it', () => {
    const cam = CAMS.R3;
    const truth = { k1: cam.record.k1 + 0.05, k2: cam.record.k2 - 0.04, focal: cam.record.focal * 1.05 };
    const col = new L.Collector(info(cam));
    for (const s of simulate(cam, truth, 3000, 4, 12)) col.add(s);
    const res = L.solve(info(cam), col.frames);
    assert(res.ok, res.reason);
    assert.equal(res.tier, 'k1k2focal');
    assert(Math.abs(res.focal / truth.focal - 1) < 0.015, `focal ${res.focal} vs ${truth.focal}`);
});

test('record values already right: nothing (or nothing meaningful) changes', () => {
    const cam = CAMS.R3;
    const col = new L.Collector(info(cam));
    for (const s of simulate(cam, cam.record, 3000, 4, 13)) col.add(s);
    const res = L.solve(info(cam), col.frames);
    assert(res.ok, res.reason);
    assert(res.rms <= res.rmsRecord + 1e-9);
    if (res.k1 !== null) assert(Math.abs(res.k1 - cam.record.k1) < 0.02, `K1 ${res.k1}`);
    if (res.k2 !== null) assert(Math.abs(res.k2 - cam.record.k2) < 0.02, `K2 ${res.k2}`);
});

test('R2 (narrow field): only K1 is refined, K2 stays the record value', () => {
    const cam = CAMS.R2;
    const truth = { k1: 0.03, k2: 0, focal: cam.record.focal };
    const col = new L.Collector(info(cam));
    for (const s of simulate(cam, truth, 4000, 4, 21, { dmin: 1.9, dmax: 3.5 })) col.add(s);
    const res = L.solve(info(cam), col.frames);
    assert(res.ok, res.reason);
    assert.equal(res.tier, 'k1');
    assert.equal(res.k2, null);
    assert.equal(res.focal, null);
    assert(Math.abs(res.k1 - 0.03) < 0.005, `K1 ${res.k1}`);
    assert.equal(res.near, null);   // only when the record values stay
});

test('record values right (R2): the record is kept', () => {
    const cam = CAMS.R2;
    const col = new L.Collector(info(cam));
    for (const s of simulate(cam, cam.record, 4000, 4, 22, { dmin: 1.9, dmax: 3.5 })) col.add(s);
    const res = L.solve(info(cam), col.frames);
    assert(res.ok, res.reason);
    assert.equal(res.tier, 'record');
    assert.equal(res.k1, null);
    // The closest calculated values are still given, very close to the record ones.
    assert(res.near, 'near values given');
    assert(res.near.k1 === null || Math.abs(res.near.k1 - cam.record.k1) < 0.01, `near K1 ${res.near.k1}`);
    assert(res.near.k2 === null || Math.abs(res.near.k2 - cam.record.k2) < 0.02, `near K2 ${res.near.k2}`);
});

test('Diamond layout seen through a mirrored tracker X: same lens values', () => {
    const cam = CAMS.R3;
    const truth = { k1: cam.record.k1 + 0.05, k2: cam.record.k2 - 0.04, focal: cam.record.focal };
    const col = new L.Collector(info(cam));
    for (const s of simulate(cam, truth, 4000, 4, 31, { diamond: true, mirror: true })) col.add(s);
    assert(col.progress().screen > 0.8, JSON.stringify(col.progress()));
    const res = L.solve(info(cam), col.frames);
    assert(res.ok, res.reason);
    assert(res.k1 !== null && res.k2 !== null, res.tier);
    assert(Math.abs(res.k1 - truth.k1) < 0.02, `K1 ${res.k1}`);
    assert(Math.abs(res.k2 - truth.k2) < 0.02, `K2 ${res.k2}`);
});
