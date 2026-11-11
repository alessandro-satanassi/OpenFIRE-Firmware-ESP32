/*  OpenFIRE Web App: statistics of the polling rate test (welcome screen).

    Run from lightgun/webapp:   node --test tests/polling-test.test.js
*/
'use strict';

const test = require('node:test');
const assert = require('node:assert/strict');
const fs = require('fs');
const path = require('path');
const vm = require('vm');

const WEBAPP = path.join(__dirname, '..');
vm.runInThisContext(fs.readFileSync(path.join(WEBAPP, 'js/app/polling-test.js'), 'utf8'), { filename: 'js/app/polling-test.js' });
const OF = globalThis.OF;

const PERIOD = 1000 / 209;

/** Reports every `period` ms from `start`, `count` of them, sampled like the display does
    (every 16 ms); returns the time of the last one. */
function feed(stats, start, count, period = PERIOD, skip = () => false) {
    let time = start;
    let nextFrame = start;
    for (let i = 0; i < count; i++) {
        time = start + i * period;
        if (!skip(i)) stats.addReport(time);
        if (time >= nextFrame) { stats.sample(time); nextFrame = time + 16; }
    }
    return time;
}

const near = (actual, expected, tolerance, label) =>
    assert.ok(Math.abs(actual - expected) <= tolerance, `${label || ''} ${actual} != ${expected}`);

test('polling stats: steady camera rate', () => {
    const stats = new OF.PollingStats();
    const end = feed(stats, 1000, 1000);
    const s = stats.snapshot(end);
    near(s.average, 209, 0.01, 'average');
    near(s.current, 209, 0.5, 'current');
    near(s.min, 209, 0.5, 'min');
    near(s.max, 209, 0.5, 'max');
    assert.equal(s.reports, 1000);
    assert.equal(s.moving, true);
});

test('polling stats: the average covers only the last 3 s of movement', () => {
    const stats = new OF.PollingStats();
    let end = feed(stats, 0, 2000, 1000 / 100);          // 20 s at 100 Hz
    near(stats.snapshot(end).average, 100, 0.01);
    end = feed(stats, end + PERIOD, 209 * 12);           // then 12 s at 209 Hz
    const s = stats.snapshot(end);
    near(s.average, 209, 0.01, 'average');                // the 100 Hz part left the window
    near(s.max, 209, 0.5, 'max');
    near(s.min, 209, 0.5, 'min');                         // min and max: the graph (last 10 s)
});

test('polling stats: a pause does not push the movement out of the average', () => {
    const stats = new OF.PollingStats();
    let end = feed(stats, 0, 300, 1000 / 150);            // 2 s at 150 Hz
    end = feed(stats, end + 30000, 100);                  // 30 s still, then 0.5 s at 209 Hz
    const s = stats.snapshot(end);
    // Window: 2 s at 150 Hz and 0.5 s at 209 Hz, the 30 s pause is not part of it.
    const expected = (299 + 99) * 1000 / (299 * 1000 / 150 + 99 * PERIOD);
    near(s.average, expected, 0.01, 'average');
});

test('polling stats: the graph stops while the gun is still', () => {
    const stats = new OF.PollingStats();
    let end = feed(stats, 0, 300);                        // ~1.4 s of movement
    stats.sample(end + 100);                              // the pause is not seen yet (120 ms)
    const before = stats.samples.length;
    assert.ok(before >= 12 && before <= 15, before);
    for (let t = end; t < end + 3000; t += 16) stats.sample(t);   // 3 s still
    assert.equal(stats.samples.length, before);
    assert.equal(stats.moving(end + 3000), false);
    assert.equal(stats.current(end + 3000), null);
    end = feed(stats, end + 3000, 300);                   // moves again: the graph goes on
    assert.ok(stats.samples.length > before + 10);
    for (const sample of stats.samples) near(sample.hz, 209, 1, 'sample');
    assert.equal(stats.samples.filter((sample) => sample.resumed).length, 1);
});

test('polling stats: the graph keeps the last 10 s of movement', () => {
    const stats = new OF.PollingStats();
    feed(stats, 0, 209 * 30);
    assert.equal(stats.samples.length, 100);
});

test('polling stats: missing reports lower the rate', () => {
    const stats = new OF.PollingStats();
    const end = feed(stats, 0, 1000, PERIOD, (i) => i % 10 === 5);
    const s = stats.snapshot(end);
    near(s.average, 209 * 0.9, 0.5, 'average');
});

test('polling stats: repeated or backward times are ignored, reset clears everything', () => {
    const stats = new OF.PollingStats();
    stats.addReport(10);
    stats.addReport(10);
    stats.addReport(9);
    stats.addReport(NaN);
    stats.addReport(15);
    assert.equal(stats.snapshot(15).reports, 2);
    assert.equal(stats.count, 1);
    stats.reset();
    const s = stats.snapshot(20);
    assert.equal(s.reports, 0);
    assert.equal(s.average, null);
    assert.equal(s.min, null);
    assert.equal(s.moving, false);
    assert.equal(stats.samples.length, 0);
});

test('polling stats: a 1000 Hz mouse is measured as such', () => {
    const stats = new OF.PollingStats();
    const end = feed(stats, 0, 2000, 1);
    near(stats.snapshot(end).average, 1000, 0.01);
});

test('polling stats: min and max follow the graph, a mouse used before goes away', () => {
    const stats = new OF.PollingStats();
    let end = feed(stats, 0, 5000, 1);                    // 5 s of a 1000 Hz mouse
    near(stats.snapshot(end).max, 1000, 1, 'max with the mouse');
    end = feed(stats, end + 1000, 209 * 5);               // 5 s of lightgun: the mouse still shows
    near(stats.snapshot(end).max, 1000, 1, 'max, mouse still in the graph');
    near(stats.snapshot(end).min, 209, 0.5, 'min');
    end = feed(stats, end + 1000, 209 * 11);              // 11 more seconds: only the lightgun
    const s = stats.snapshot(end);
    near(s.max, 209, 0.5, 'max');
    near(s.min, 209, 0.5, 'min');
});

test('polling test: colour of the average (green from 180 Hz, yellow from 60 Hz, red below)', () => {
    const level = OF.PollingTest.averageLevel;
    assert.equal(level(null), null);
    assert.equal(level(209), 'ok');
    assert.equal(level(180), 'ok');
    assert.equal(level(179), 'warn');
    assert.equal(level(150), 'warn');
    assert.equal(level(60), 'warn');
    assert.equal(level(59), 'bad');
});

test('polling test: the graph changes colour exactly where the line crosses 180 and 60 Hz', () => {
    const split = OF.PollingTest.splitAtThresholds;
    const level = OF.PollingTest.averageLevel;
    // From a mouse at 1000 Hz down to 30 Hz in one step: green, orange and red pieces.
    const pieces = split(0, 1000, 1, 30);
    assert.equal(pieces.length, 3);
    assert.deepEqual(pieces.map((p) => level((p[1] + p[3]) / 2)), ['ok', 'warn', 'bad']);
    near(pieces[0][3], 180, 1e-9, 'first cut');
    near(pieces[1][3], 60, 1e-9, 'second cut');
    near(pieces[0][2], (1000 - 180) / (1000 - 30), 1e-9, 'where');
    // Within one band: a single piece.
    assert.equal(split(3, 200, 4, 950).length, 1);
    // Upwards too.
    assert.deepEqual(split(0, 50, 1, 209).map((p) => level((p[1] + p[3]) / 2)), ['bad', 'warn', 'ok']);
});

