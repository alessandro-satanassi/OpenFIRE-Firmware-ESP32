/*  OpenFIRE Web App - browser checks of the switches in Gun Tests against the simulated lightgun.

    Rumble, solenoid and autofire switches (sSwitchStateUpd): on/off with normally open and normally
    closed contacts, without a pin, after a save that changes the normally closed setting, with a
    firmware that does not send them, and the layout (analog stick as tall as the buttons, temperature
    beside the switches).

    Needs Playwright (npm install -g playwright, NODE_PATH set to the global modules).
    Run from lightgun/webapp:   node tests/gun-tests-switches-e2e.js
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
const WEBAPP = path.join(__dirname, '..');
const SHOTS = process.env.OF_E2E_SCREENSHOTS || null;
const shot = async (page, name) => { if (SHOTS) await page.screenshot({ path: path.join(SHOTS, name + '.png') }); };
const { startServer } = require('./sim-server.js');
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
async function waitFor(fn, ms = 6000) { const t = Date.now(); while (!(await fn())) { if (Date.now() - t > ms) return false; await sleep(30); } return true; }
let failures = 0;
const ok = (c, m) => { console.log((c ? 'PASS ' : 'FAIL ') + m); if (!c) failures++; };
const loaded = (p) => p.evaluate(() => !!(window.OF && OF.app && OF.app.state.loaded));

(async () => {
    const sim = await startServer({ port: 8131, root: WEBAPP, board: 'rpipico' });
    const fw = sim.firmware;
    const S = globalThis.OpenFIREshared, E = S.boardInputs_e, B = S.boolTypes_e, T = S.settingsTypes_e;
    const setSetting = (i, v) => new DataView(fw.settings.buffer).setUint32(i * 4, v >>> 0, true);
    // Custom pins: rumble switch on GPIO22, solenoid switch on GPIO26 normally closed, autofire without pin.
    fw.toggles[B.customPins] = 1;
    fw.pins[E.rumbleSwitch] = 22;
    fw.pins[E.solenoidSwitch] = 26;
    setSetting(T.inputInvertMask, 1 << E.solenoidSwitch);

    const browser = await chromium.launch();
    const page = await browser.newPage({ viewport: { width: 1280, height: 820 } });
    const errors = [];
    page.on('pageerror', (e) => errors.push(e.message));
    page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
    page.on('dialog', (d) => d.accept().catch(() => {}));
    await page.goto('http://localhost:8131/?ws');
    ok(await waitFor(() => loaded(page)), 'docked');
    await page.click('.tab-button[aria-controls="tab-tests"]');
    const cls = () => page.$$eval('#tab-tests .test-switch', (n) => n.map((e) => e.className.replace('test-switch', '').trim() || 'off'));
    const names = await page.$$eval('#tab-tests .test-switch .switch-name', (n) => n.map((e) => e.textContent));
    ok(JSON.stringify(names) === '["Rumble","Solenoid","Autofire"]', 'names ' + names);
    ok(await waitFor(async () => JSON.stringify(await cls()) === '["off","on","not-connected"]'),
        'contacts open: rumble NO off, solenoid NC on, autofire N/C -> ' + JSON.stringify(await cls()));
    await shot(page, 'switches-light');

    fw.switchClosed = [true, true, false];
    ok(await waitFor(async () => JSON.stringify(await cls()) === '["on","off","not-connected"]'),
        'contacts closed: rumble on, solenoid off -> ' + JSON.stringify(await cls()));
    fw.switchClosed = [false, false, true];
    ok(await waitFor(async () => JSON.stringify(await cls()) === '["off","on","not-connected"]'), 'autofire contact without pin stays N/C');

    // Layout: edges lined up.
    const r = await page.evaluate(() => {
        const b = (s) => document.querySelector('#tab-tests ' + s).getBoundingClientRect();
        const scroll = document.querySelector('#tab-tests .tab-scroll');
        return { buttons: b('.buttons-group'), analog: b('.analog-group'), sw: b('.switches-group'), temp: b('.temperature'),
                 view: b('.analog-view'), legend: b('.switches-group > legend'), fits: scroll.scrollHeight <= scroll.clientHeight + 1 };
    });
    const near = (a, b) => Math.abs(a - b) <= 1.5;
    ok(near(r.buttons.top, r.analog.top) && near(r.buttons.bottom, r.analog.bottom), 'analog stick as tall as the buttons');
    ok(near(r.sw.bottom, r.temp.bottom), 'temperature bottom = switch box bottom');
    ok(near(r.temp.top, r.legend.top + r.legend.height / 2) && near(r.sw.top - r.buttons.bottom + r.legend.height / 2, 12),
        `temperature top on the switch box border (${r.temp.top} vs ${r.legend.top + r.legend.height / 2})`);
    ok(near(r.buttons.left, r.sw.left) && near(r.buttons.right, r.sw.right) && near(r.analog.left, r.temp.left) && near(r.analog.right, r.temp.right), 'columns aligned');
    ok(near(r.view.width, r.view.height) && r.view.height > 60 && r.view.bottom <= r.analog.bottom, `square stick view inside its box ${r.view.width}x${r.view.height}`);
    ok(r.fits, 'tab fits the window');

    // Save with a new normally-closed setting: the board applies it, the view follows.
    await page.evaluate((fn) => { OF.app.state.setInverted(fn, true); OF.app.refresh(); }, E.rumbleSwitch);
    await page.click('.save-button');
    await page.click('dialog button.primary');
    ok(await waitFor(async () => (await page.locator('.status-text').innerText()).includes('Sent settings successfully!')), 'save completed');
    ok((fw.getSetting(T.inputInvertMask) >>> E.rumbleSwitch) & 1, 'firmware got the new mask');
    ok(await waitFor(async () => JSON.stringify(await cls()) === '["on","on","not-connected"]'),
        'rumble now normally closed: open contact = on -> ' + JSON.stringify(await cls()));

    // Dark theme + narrow.
    await page.emulateMedia({ colorScheme: 'dark' });
    await shot(page, 'switches-dark');
    await page.setViewportSize({ width: 420, height: 900 });
    await sleep(300);
    await shot(page, 'switches-narrow');
    const overflow = await page.evaluate(() => document.documentElement.scrollWidth > window.innerWidth + 1);
    ok(!overflow, 'no horizontal overflow at phone width');

    // Firmware without sSwitchStateUpd: mapped switches show an unknown state.
    await page.setViewportSize({ width: 1280, height: 820 });
    await page.emulateMedia({ colorScheme: 'light' });
    fw.switchReports = false;
    await page.reload();
    ok(await waitFor(() => loaded(page)), 'docked again');
    await page.click('.tab-button[aria-controls="tab-tests"]');
    await sleep(800);
    ok(JSON.stringify(await cls()) === '["unknown","unknown","not-connected"]', 'old firmware: unknown -> ' + JSON.stringify(await cls()));
    const title = await page.getAttribute('#tab-tests .test-switch', 'title');
    ok(/firmware/.test(title), 'tooltip: ' + title);
    await shot(page, 'switches-old-firmware');

    ok(errors.length === 0, 'no page errors ' + JSON.stringify(errors));
    await browser.close();
    await sim.close();
    console.log(failures ? `${failures} FAILED` : 'ALL PASSED');
    process.exit(failures ? 1 : 0);
})().catch((e) => { console.error(e); process.exit(2); });
