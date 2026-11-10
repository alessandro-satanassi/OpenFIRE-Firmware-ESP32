/*  OpenFIRE Web App - browser checks of the aim test (Gun Tests) against the simulated lightgun.

    The button beside the IR camera tester, the gun in IR test mode while the window is open,
    the trigger firing the solenoid (or the rumble when the solenoid is off), the other buttons,
    Esc leaving the test mode.

    Needs Playwright (npm install -g playwright, NODE_PATH set to the global modules).
    Run from lightgun/webapp:   node tests/aim-test-e2e.js
    OF_E2E_SCREENSHOTS=<folder> also saves screenshots.
*/
'use strict';

const path = require('path');
let chromium;
try {
    ({ chromium } = require('playwright'));
} catch (error) {
    console.log('Playwright is not installed: npm install -g playwright (and set NODE_PATH to the global modules).');
    process.exit(0);
}
const { startServer } = require('./sim-server.js');
const SHOTS = process.env.OF_E2E_SCREENSHOTS || null;
const shot = async (page, name) => { if (SHOTS) await page.screenshot({ path: path.join(SHOTS, name + '.png') }); };
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
async function waitFor(fn, ms = 6000) { const t = Date.now(); while (!(await fn())) { if (Date.now() - t > ms) return false; await sleep(30); } return true; }
let failures = 0;
const ok = (c, m) => { console.log((c ? 'PASS ' : 'FAIL ') + m); if (!c) failures++; };

(async () => {
    const sim = await startServer({ port: 8132, root: path.join(__dirname, '..'), board: 'rpipico' });
    const fw = sim.firmware;
    const S = globalThis.OpenFIREshared, E = S.boardInputs_e, B = S.boolTypes_e, C = S.serialCmdTypes_e;
    const browser = await chromium.launch();
    for (const lang of ['it', 'en']) {
        // Italian with the solenoid on, English with the solenoid off and the rumble on.
        fw.toggles[B.solenoid] = lang === 'it' ? 1 : 0;
        fw.toggles[B.rumble] = 1;
        const page = await browser.newPage({ viewport: { width: 1280, height: 820 } });
        const errors = [];
        page.on('pageerror', (e) => errors.push(e.message));
        page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
        await page.goto('http://localhost:8132/?ws&lang=' + lang);
        ok(await waitFor(() => page.evaluate(() => !!(window.OF && OF.app && OF.app.state.loaded))), lang + ': docked');
        await page.click('.tab-button[aria-controls="tab-tests"]');
        await sleep(200);
        const row = await page.$$eval('#tab-tests .test-windows > button', (b) => b.map((e) => Math.round(e.getBoundingClientRect().top)));
        ok(row.length === 2 && row[0] === row[1], lang + ': IR tester and aim test on one row');
        await shot(page, 'aim-tab-' + lang);

        await page.click('#tab-tests .test-windows > button:nth-child(2)');
        ok(await waitFor(async () => await page.locator('.fullscreen-window.mode-aimtest').count() === 1), lang + ': aim test open');
        ok(await waitFor(() => fw.runMode === 'processing'), lang + ': gun in test mode');
        ok(await waitFor(() => page.evaluate(() => !!(OF.app.window && OF.app.window.aim))), lang + ': aim point received');

        const solenoid = fw.toggles[B.solenoid];
        fw.lastTest = 0;
        fw.queueEvent(C.sBtnPressed, Uint8Array.of(E.btnTrigger));
        fw.queueEvent(C.sBtnReleased, Uint8Array.of(E.btnTrigger));
        fw.queueEvent(C.sBtnPressed, Uint8Array.of(E.btnGunA));
        await sleep(60);
        await shot(page, 'aim-shot-' + lang);
        ok(await waitFor(() => fw.lastTest === (solenoid ? C.sTestSolenoid : fw.toggles[B.rumble] ? C.sTestRumble : 0), 2000),
            lang + ': trigger -> ' + (solenoid ? 'solenoid' : 'rumble') + ' (' + fw.lastTest + ')');
        ok(await page.evaluate(() => OF.app.window.shots.length === 1 && !!OF.app.window.button), lang + ': hit mark and button name');
        fw.queueEvent(C.sBtnReleased, Uint8Array.of(E.btnGunA));
        // F: no recoil on the trigger, then on again.
        await page.keyboard.press('f');
        fw.lastTest = 0;
        fw.queueEvent(C.sBtnPressed, Uint8Array.of(E.btnTrigger));
        fw.queueEvent(C.sBtnReleased, Uint8Array.of(E.btnTrigger));
        await sleep(300);
        ok(fw.lastTest === 0 && await page.evaluate(() => OF.app.window.shots.length === 2), lang + ': F off -> shot without recoil');
        await page.keyboard.press('f');
        ok(await page.evaluate(() => OF.app.window.feedback === true), lang + ': F on again');
        if (SHOTS) { await sleep(600); await shot(page, 'aim-rest-' + lang); }

        // Emitters (sTestBlobs before sTestCoords): 3 seen = estimated aim, 1 seen = lost, no hit mark.
        ok(await page.evaluate(() => OF.app.window.aim.seen === 3 && OF.app.window._aimFresh(Date.now())), lang + ': 3 emitters -> aim valid');
        const lost = await page.evaluate(() => {
            const win = OF.app.window;
            const coords = new Uint8Array(48);
            new DataView(coords.buffer).setInt32(32, 300, true);
            const blobs = new Uint8Array(21);
            blobs[0] = 1; blobs[1] = 1; // only the first emitter seen
            const shots = win.shots.length;
            win.setTestBlobs(blobs);
            win.drawTest(coords);
            win.onButton(0, true);
            return !win._aimFresh(Date.now()) && win.shots.length === shots;
        });
        ok(lost, lang + ': 1 emitter -> aim lost, shot not marked');
        // Tracked flag (bit 2): 2 emitters but the tracker kept the old point (after a dark moment,
        // also right after opening) -> lost; 2 emitters and updated -> valid.
        const tracked = await page.evaluate(() => {
            const win = OF.app.window;
            const coords = new Uint8Array(48);
            const send = (flags) => {
                const blobs = new Uint8Array(21);
                blobs[0] = 1;
                flags.forEach((f, i) => { blobs[1 + i * 5] = f; });
                win.setTestBlobs(blobs);
                win.drawTest(coords);
                return win._aimFresh(Date.now());
            };
            const stale = send([1, 1, 0, 0]);               // 2 seen, not updated
            const updated = send([5, 5, 4, 4]);             // 2 seen, updated
            const again = send([1, 1, 1, 0]);               // 3 seen, not updated yet
            // Emitters seen but not tracked: "Aim unavailable", not "IR emitters not seen".
            const texts = []; const proto = CanvasRenderingContext2D.prototype; const fill = proto.fillText;
            proto.fillText = function (text, ...rest) { texts.push(String(text)); return fill.call(this, text, ...rest); };
            win._draw();
            proto.fillText = fill;
            const message = texts.includes(OF.i18n.t('Aim unavailable')) && !texts.includes(OF.i18n.t('IR emitters not seen'));
            return !stale && updated && !again && message;
        });
        ok(tracked, lang + ': the tracked flag decides');
        // Redraw timer: every 40 ms only while something fades, otherwise once when the aim gets old.
        const delays = await page.evaluate(() => {
            const win = OF.app.window;
            const seen = []; const set = window.setTimeout;
            window.setTimeout = (fn, ms) => { seen.push(ms); return set(fn, ms); };
            win.shots = []; win.button = null; win.flash = 0;
            win.aim = { x: 0.5, y: 0.5, time: Date.now(), seen: 4, tracked: true };
            win._draw();
            win.flash = Date.now();
            win._draw();
            window.setTimeout = set;
            return seen;
        });
        ok(delays.length === 2 && delays[0] > 400 && delays[1] === 40, lang + ': redraw timer ' + JSON.stringify(delays));
        await page.keyboard.press('Escape');
        ok(await waitFor(async () => await page.locator('.fullscreen-window').count() === 0), lang + ': closed with Esc');
        ok(await waitFor(() => fw.runMode !== 'processing'), lang + ': test mode left');
        ok(await waitFor(() => page.evaluate(() => !OF.app.irTestActive && !document.getElementById('app').inert)), lang + ': page usable again');
        ok(errors.length === 0, lang + ': no errors ' + errors.join(' | '));
        await page.close();
    }
    await browser.close();
    if (sim.close) await sim.close();
    console.log(failures ? failures + ' FAILED' : 'ALL PASSED');
    process.exit(failures ? 1 : 0);
})();
