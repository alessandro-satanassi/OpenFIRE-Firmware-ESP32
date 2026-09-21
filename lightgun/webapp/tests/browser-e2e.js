/*  OpenFIRE Web App - browser checks of the built app against the simulated lightgun.

    Needs Playwright (npm install -g playwright, NODE_PATH set to the global modules) and the built folders:
      python scripts/webapp_build.py device --board waveshare-esp32-s3-zero --out dist/device
      python scripts/webapp_build.py site
    Run from lightgun/webapp:   node tests/browser-e2e.js
    OF_E2E_SCREENSHOTS=<folder> also saves screenshots.

    The site is tested with a simulated Web Serial port (navigator.serial) connected to the
    serial link of the simulated firmware; the lightgun page uses its WebSocket.
*/
'use strict';

const fs = require('fs');
const path = require('path');
let chromium;
try {
    ({ chromium } = require('playwright'));
} catch (error) {
    console.log('Playwright is not installed: npm install -g playwright (and set NODE_PATH to the global modules).');
    process.exit(0);
}
const LIGHTGUN = path.join(__dirname, '..', '..');
const { startServer } = require('./sim-server.js');
const SHOTS = process.env.OF_E2E_SCREENSHOTS || null;
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
async function waitFor(fn, ms = 8000) { const t = Date.now(); while (!(await fn())) { if (Date.now() - t > ms) return false; await sleep(50); } return true; }
const ok = (cond, msg) => { console.log((cond ? 'PASS ' : 'FAIL ') + msg); if (!cond) process.exitCode = 1; };
const shot = async (page, name) => { if (SHOTS) await page.screenshot({ path: path.join(SHOTS, name + '.png') }); };

function watchErrors(page) {
    const errors = [];
    page.on('pageerror', (e) => errors.push(e.message));
    page.on('console', (m) => { if (m.type() === 'error') errors.push(m.text()); });
    return errors;
}

const loaded = (page) => page.evaluate(() => !!(window.OF && OF.app && OF.app.state.loaded));
const statusText = (page) => page.locator('.status-text').innerText();
const saveLabel = (page) => page.locator('.save-label').innerText();
const tab = (page, id) => page.click(`.tab-button[aria-controls="tab-${id}"]`);

/** Simulated Web Serial port bridged to the serial link of the simulated firmware. */
async function installFakeSerial(context, sim) {
    const state = { page: null, opens: [], busy: false, granted: false, closes: [] };
    await context.exposeBinding('__ofSerialWrite', (source, bytes) => sim.serialLink.appWrite(Uint8Array.from(bytes)));
    // Order in which the page tears the port down: a real port only guarantees that the
    // bytes left once the writable stream is closed, so 'flush' must come before 'port'.
    await context.exposeBinding('__ofSerialClose', (source, what) => { state.closes.push(what); });
    await context.exposeBinding('__ofSerialOpen', (source, baud) => {
        state.page = source.page;
        state.opens.push(baud);
        if (baud === 1200) sim.firmware.rebootedToBootloader = true;
        return state.busy;
    });
    sim.serialLink.toApp = (bytes) => {
        if (state.page) state.page.evaluate((data) => window.__ofSerialPush && window.__ofSerialPush(data), Array.from(bytes)).catch(() => {});
    };
    await context.addInitScript(() => {
        const listeners = new Set();
        class FakePort {
            getInfo() { return { usbVendorId: 0xF143, usbProductId: 0x0001 }; }
            async open(options) {
                if (await window.__ofSerialOpen(options.baudRate)) { const e = new Error('busy'); e.name = 'NetworkError'; throw e; }
                this.readable = new ReadableStream({ start: (controller) => { window.__ofSerialPush = (data) => { try { controller.enqueue(new Uint8Array(data)); } catch (e) { /* closed */ } }; this._controller = controller; } });
                this.writable = new WritableStream({
                    write: (chunk) => window.__ofSerialWrite(Array.from(chunk)),
                    close: () => window.__ofSerialClose('flush'),
                });
            }
            async setSignals() {}
            async close() {
                window.__ofSerialClose('port');
                window.__ofSerialPush = null;
                try { this._controller.close(); } catch (e) { /* already closed */ }
                this.readable = null;
                this.writable = null;
            }
        }
        const port = new FakePort();
        // Like a real browser, the permission belongs to the site and not to the page:
        // after a jump to another page of the same site the port is still there.
        let granted = false;
        try { granted = localStorage.getItem('__of_fake_serial_granted') === '1'; } catch (e) { /* private mode */ }
        const serial = {
            requestPort: async () => {
                granted = true;
                try { localStorage.setItem('__of_fake_serial_granted', '1'); } catch (e) { /* private mode */ }
                return port;
            },
            getPorts: async () => (granted ? [port] : []),
            addEventListener: (name, fn) => listeners.add(fn),
            removeEventListener: (name, fn) => listeners.delete(fn),
        };
        Object.defineProperty(Navigator.prototype, 'serial', { get: () => serial, configurable: true });
    });
    return state;
}

(async () => {
    const browser = await chromium.launch();

    // ===================== lightgun page (WebSocket) =====================
    let sim = await startServer({ port: 8123, root: LIGHTGUN + '/dist/device', board: 'waveshare-esp32-s3-zero' });
    const ctx = await browser.newContext({ locale: 'en-US', viewport: { width: 1280, height: 860 } });
    const page = await ctx.newPage();
    const errors = watchErrors(page);
    await page.goto('http://localhost:8123/');
    ok(await waitFor(() => loaded(page)), 'device page docks by itself');
    ok((await page.locator('.board-title').innerText()) === 'FIRECon P1 | Waveshare ESP32-S3-Zero', 'board title: ' + await page.locator('.board-title').innerText());
    ok(await waitFor(() => page.evaluate(() => !!document.querySelector('.board-picture').shadowRoot?.querySelector('svg'))), 'board picture shown as inline SVG (inside app.js)');
    await page.hover('.pin-box[aria-label="GPIO4"]');
    ok(await page.evaluate(() => getComputedStyle(document.querySelector('.board-picture').shadowRoot.getElementById('OF_pin4')).opacity === '1'), 'pointing at a pin box lights its circle');
    const rowsOk = (scope) => page.evaluate((scope) => {
        const P = OF.Boards.shared.boardBoxPositions_e;
        const positions = OF.app.state.boxPositions;
        return [...document.querySelectorAll(`${scope} .pins-col`)].every((column) => {
            const items = [...column.children];
            return items.every((item, k) => {
                const gpio = Number(item.querySelector('.gpio-label').textContent.match(/\d+/)[0]);
                const slot = positions[gpio] & ~P.posCheck;
                const placed = getComputedStyle(item).gridRowStart === String(slot);
                return placed && (k === 0 || item.getBoundingClientRect().top > items[k - 1].getBoundingClientRect().top);
            });
        });
    }, scope);
    ok(await rowsOk('#tab-pins'), 'pin boxes on the rows of boardsBoxPositions (empty rows kept, top to bottom)');
    const fits = async (size) => {
        if (size) await page.setViewportSize(size);
        await sleep(400);
        return page.evaluate(() => {
            const scroll = document.querySelector('#tab-pins .tab-scroll');
            return scroll.scrollHeight <= scroll.clientHeight + 1;
        });
    };
    ok(await fits(), 'board layout fits the window without scrolling');
    ok(await fits({ width: 1280, height: 720 }), 'board layout still fits a shorter window');
    await page.setViewportSize({ width: 1280, height: 860 });
    await tab(page, 'profiles');
    await sleep(400);
    ok(await page.evaluate(() => {
        const scroll = document.querySelector('#tab-profiles .table-scroll');
        return scroll.scrollWidth <= scroll.clientWidth + 1;
    }), 'profiles table fits the window without scrolling sideways');
    await tab(page, 'pins');
    // The User Layouts button sits on the right of the bottom bar: its menu must stay in the window.
    await page.click('#tab-pins .pins-bottom .menu-button');
    ok(await page.evaluate(() => {
        const box = document.querySelector('#tab-pins .menu-panel').getBoundingClientRect();
        return box.left >= 0 && box.top >= 0 && box.right <= innerWidth + 1 && box.bottom <= innerHeight + 1;
    }), 'the User Layouts menu opens inside the window');
    await page.click('#tab-pins .pins-bottom .menu-button');
    ok(await page.locator('.device-bar').count() === 0, 'no device selector on the lightgun page');
    ok(await page.locator('.link-state.on').count() === 1, 'the lightgun page shows the plug as connected too');
    ok(await page.locator('.disconnect-button').count() === 0, 'no Disconnect button on the lightgun page');
    ok(await page.locator('.menu-button', { hasText: 'Board Previews' }).count() === 0, 'board previews hidden on the lightgun page');
    ok(sim.firmware.gunMode === 'docked' && sim.firmware.sessionActive, 'firmware docked');
    await shot(page, 'device-docked');
    const s1 = sim.firmware.sessionCounter;

    sim.dropClient();
    ok(await waitFor(async () => (await statusText(page)).includes('reconnecting')), 'lost -> reconnecting');
    ok(await waitFor(() => loaded(page)), 'docked again after WebSocket drop');
    ok(sim.firmware.sessionCounter > s1, 'new firmware session');

    const s2 = sim.firmware.sessionCounter;
    await page.reload();
    ok(await waitFor(() => loaded(page)), 'docked again after page reload');
    ok(sim.firmware.sessionCounter > s2, 'new firmware session after reload');

    // Two pages of the lightgun open: they must not take the gun from each other forever.
    const page2 = await ctx.newPage();
    const errors2 = watchErrors(page2);
    await page2.goto('http://localhost:8123/');
    const stopped = async (p) => p.evaluate(() => OF.app.connection.autoStopped);
    ok(await waitFor(async () => (await stopped(page)) || (await stopped(page2)), 30000), 'two lightgun pages: one of them stops retrying');
    const [winner, loser] = (await stopped(page)) ? [page2, page] : [page, page2];
    ok(await waitFor(() => loaded(winner)), 'two lightgun pages: the other one stays docked');
    const sessions = sim.firmware.sessionCounter;
    await sleep(4000);
    ok(sim.firmware.sessionCounter === sessions && await loaded(winner), 'two lightgun pages: no more takeovers');
    ok(await loser.locator('.welcome button', { hasText: 'Reconnect' }).isVisible(), 'stopped page offers Reconnect');
    ok((await statusText(loser)).includes('taken over by another page'), 'stopped page explains why');
    await winner.close();
    await loser.click('.welcome button:has-text("Reconnect")');
    ok(await waitFor(() => loaded(loser)), 'Reconnect docks the page again');
    const survivor = loser;
    ok(errors.length === 0 && errors2.length === 0, 'no page errors with two pages ' + JSON.stringify(errors.concat(errors2)));

    await survivor.selectOption('.lang-selector', 'it');
    ok((await survivor.locator('.tab-button[aria-controls="tab-pins"] .tab-label').innerText()) === 'Layout Scheda', 'language switch translates the tabs');
    ok((await survivor.locator('.lang-selector option').allInnerTexts()).join(',') === 'English,Italiano', 'language names');
    ok(await loaded(survivor), 'still docked after the language change');
    await survivor.selectOption('.lang-selector', 'en');

    await survivor.close();
    ok(await waitFor(() => !sim.firmware.sessionActive && sim.firmware.gunMode === 'run'), 'closing the page returns the gun to Run mode');
    ok(errors.length === 0 && errors2.length === 0, 'no page errors ' + JSON.stringify(errors.concat(errors2)));
    await sim.close();

    // ===================== unbundled webapp folder: editing flows =====================
    sim = await startServer({ port: 8125, root: LIGHTGUN + '/webapp', board: 'rpipico' });
    const fw = sim.firmware;
    const E = await (async () => { const p = await ctx.newPage(); await p.goto('http://localhost:8125/?ws'); await waitFor(() => loaded(p)); const e = await p.evaluate(() => OF.Boards.shared.boardInputs_e); await p.close(); return e; })();
    await waitFor(() => !fw.sessionActive);
    const dev = await ctx.newPage();
    const devErrors = watchErrors(dev);
    // Reloading with unsaved edits: the browser asks; the checks leave the page.
    dev.on('dialog', (dialog) => dialog.accept().catch(() => {}));
    await dev.goto('http://localhost:8125/?ws');
    ok(await waitFor(() => loaded(dev)), 'unbundled webapp docks with ?ws');
    ok((await saveLabel(dev)) === '[Nothing To Save]', 'nothing to save after loading');

    // Board Layout
    ok(await dev.locator('.pin-box[aria-label="GPIO0"]').isDisabled(), 'pin boxes locked while custom pins are off');
    await dev.click('text=Use Custom Pins');
    await dev.selectOption('.pin-box[aria-label="GPIO0"]', String(E.btnPedal + 1));
    ok((await saveLabel(dev)) === 'Save and Send Settings', 'edits enable Save');
    await dev.selectOption('.pin-box[aria-label="GPIO4"]', String(E.periphSDA + 1));
    ok((await statusText(dev)).includes('Camera and Peripheral Data pins clashed'), 'I2C clash reported in the status bar');
    await shot(dev, 'dev-pins');

    // Button Mapping
    await tab(dev, 'buttons');
    const typeBox = dev.locator('.btn-row:not(.btn-header) >> nth=3').locator('.btn-type').first();
    await typeBox.hover();
    await typeBox.selectOption('0');
    ok(await dev.locator('#tab-buttons .desc-title').innerText().then((t) => t.length > 0), 'description box follows the pointer');

    // Gun Settings
    await tab(dev, 'settings');
    const intensity = dev.locator('.spin', { hasText: '/ 255' }).locator('input');
    await intensity.fill('128');
    await intensity.press('Tab');
    ok(await dev.locator('text=Temperature Warning Threshold:').isHidden(), 'unsafe settings hidden by default');
    await dev.getByRole('button', { name: 'View', exact: true }).click();
    await dev.click('.menu-item:has-text("Show Unsafe Settings")');
    ok(await dev.locator('text=Temperature Warning Threshold:').isVisible(), 'View > Show Unsafe Settings shows the temperature thresholds');
    await dev.click('text=Advanced View');
    await dev.fill('.text-input', 'Gun Ω');
    ok((await dev.inputValue('.text-input')) !== 'Gun Ω', 'product name refuses characters outside Latin-1');
    await dev.fill('.text-input', 'My Gun');
    await shot(dev, 'dev-settings');

    // Calibration Profiles: rename and select
    await tab(dev, 'profiles');
    await dev.click('.profiles-table .icon-button >> nth=1');
    await dev.fill('dialog input[type="text"]', 'Living Room');
    await dev.click('dialog button.primary');
    ok(await waitFor(async () => (await dev.locator('.profile-radio >> nth=1').innerText()).includes('2. Living Room')), 'profile renamed');
    await dev.click('.profile-radio >> nth=2');
    ok(await waitFor(() => fw.currentProfile === 2), 'selecting a profile selects it on the board');
    await dev.click('.profiles-table .icon-button >> nth=0');
    await dev.fill('dialog input[type="text"]', 'Salón');
    await dev.fill('dialog input[type="text"]', 'Salón 日');
    ok((await dev.inputValue('dialog input[type="text"]')) === 'Salón' && await dev.locator('dialog input.invalid').count() === 1,
        'profile name refuses characters outside Latin-1');
    await dev.click('dialog button:has-text("Cancel")');

    // Links of the Qt texts open outside the page: the session stays.
    await dev.locator('.profiles-table select >> nth=3').hover();
    const popup = dev.context().waitForEvent('page', { timeout: 5000 }).catch(() => null);
    await dev.locator('#tab-profiles .desc-text a').first().click();
    const opened = await popup;
    ok(!!opened && await loaded(dev) && dev.url().startsWith('http://localhost:8125/'), 'description link opens in a new tab, the page keeps its board');
    if (opened) await opened.close();

    // Save
    await dev.click('.save-button');
    await dev.click('dialog button.primary');
    ok(await waitFor(async () => (await statusText(dev)).includes('Sent settings successfully!'), 8000), 'save completed');
    ok((await saveLabel(dev)) === '[Nothing To Save]', 'nothing to save after saving');
    ok(fw.toggles[0] === 1 && fw.pins[E.btnPedal] === 0, 'firmware received the custom pins');
    ok(fw.getSetting(0) === 128, 'firmware received rumble intensity');
    ok(fw.getProfileName(1) === 'Living Room', 'firmware received the profile name');
    ok(fw.getUSB().name === 'My Gun', 'firmware received the product name');
    ok((await dev.locator('.board-title').innerText()).startsWith('My Gun |'), 'title shows the saved name');

    // Calibration
    await dev.click('text=Calibrate Profile 2');
    ok(await waitFor(() => dev.locator('.fullscreen-window.mode-calibrate').count().then((n) => n === 1)), 'calibration window opens');
    await sleep(300);
    await shot(dev, 'dev-cali-init');
    for (let i = 0; i < 6; ++i) { fw.trigger = true; await sleep(250); }
    await dev.mouse.move(600, 300);
    await shot(dev, 'dev-cali-verify');
    fw.trigger = true;
    ok(await waitFor(() => dev.locator('.fullscreen-window').count().then((n) => n === 0)), 'calibration window closes at the end');
    ok((await statusText(dev)).includes('Calibration for Profile 2 successful'), 'calibration result in the status bar');
    ok((await dev.locator('.profiles-table tbody tr >> nth=1').locator('td.value').first().innerText()) === '90', 'calibrated offsets in the table');

    await dev.click('text=Calibrate Profile 1');
    await sleep(300);
    await dev.keyboard.press('Escape');
    ok(await waitFor(() => dev.locator('.fullscreen-window').count().then((n) => n === 0)), 'ESC cancels the calibration through the board');
    ok((await statusText(dev)).includes('Cancelled Calibration for Profile 1'), 'cancelled calibration reported');

    // Link lost with the Save question open and a calibration not saved yet.
    await dev.click('.save-button');
    ok(await dev.locator('dialog', { hasText: 'Are these settings okay?' }).count() === 1, 'save question open');
    sim.dropClient();
    ok(await waitFor(() => dev.locator('dialog', { hasText: 'Are these settings okay?' }).count().then((n) => n === 0)), 'session dialogs close when the board is lost');
    ok(await waitFor(() => loaded(dev)), 'docked again after the loss');
    ok(await waitFor(() => dev.locator('dialog', { hasText: 'The connection was lost before your changes were saved.' }).count().then((n) => n === 1)), 'unsaved edits offered back');
    await dev.click('dialog button:has-text("Yes")');
    ok(await waitFor(async () => (await saveLabel(dev)) === 'Save and Send Settings'), 'restored calibration can be saved');
    ok((await dev.locator('.profiles-table tbody tr >> nth=1').locator('td.value').first().innerText()) === '90', 'restored calibration values');
    ok(!(await dev.evaluate(() => OF.app.busy)), 'no save left running');

    // Gun Tests
    await tab(dev, 'tests');
    fw.queueEvent(fw.C.sBtnPressed, Uint8Array.of(E.btnTrigger));
    fw.queueEvent(fw.C.sTemperatureUpd, Uint8Array.of(40));
    ok(await waitFor(() => dev.locator('.test-button.pressed').count().then((n) => n === 1)), 'pressed button highlighted');
    ok(await waitFor(() => dev.locator('.temperature.warm').count().then((n) => n === 1)), 'temperature above the warning threshold');
    fw.queueEvent(fw.C.sBtnReleased, Uint8Array.of(E.btnTrigger));
    ok(await waitFor(() => dev.locator('.test-button.pressed').count().then((n) => n === 0)), 'released button');
    await dev.click('text=Test Rumble Motor');
    ok(await waitFor(() => fw.lastTest === fw.C.sTestRumble), 'rumble test command');
    await sleep(400);
    ok(await dev.evaluate(() => {
        const scroll = document.querySelector('#tab-tests .tab-scroll');
        const buttons = document.querySelector('#tab-tests .test-buttons').getBoundingClientRect();
        const temp = document.querySelector('#tab-tests .temperature').getBoundingClientRect();
        const box = document.querySelector('#tab-tests .test-buttons').closest('.group').getBoundingClientRect();
        return scroll.scrollHeight <= scroll.clientHeight + 1 && buttons.height > 120 && temp.bottom <= box.bottom;
    }), 'gun tests fit the window, the buttons fill the box with the temperature at the bottom');

    await dev.dblclick('text=Open IR Camera Tester...');
    ok(await waitFor(() => fw.runMode === 'processing'), 'IR test mode on the board');
    await sleep(300);
    ok(await dev.locator('.fullscreen-window').count() === 1, 'double click opens a single IR test window');
    ok(await waitFor(async () => (await saveLabel(dev)) === '[Disabled while in Test Mode]'), 'Save disabled during the IR test');
    await sleep(400);
    await shot(dev, 'dev-irtest');
    await dev.keyboard.press('Escape');
    ok(await waitFor(() => fw.runMode === 'normal'), 'IR test mode ended with ESC');
    ok(await waitFor(async () => (await dev.locator('.fullscreen-window').count()) === 0 && !(await dev.evaluate(() => OF.app.irTestActive))), 'IR test window and lock gone');
    ok(await dev.evaluate(() => !document.getElementById('app').inert), 'page usable again after the fullscreen window');

    await dev.getByRole('button', { name: 'Emitter Alignment', exact: true }).click();
    ok(await dev.locator('.fullscreen-window.mode-alignment').count() === 1, 'alignment assistant opens');
    await dev.keyboard.press('Escape');
    ok(await waitFor(() => dev.locator('.fullscreen-window').count().then((n) => n === 0)), 'alignment assistant closes');

    await dev.click('text=Clear Save Memory [!]');
    await dev.click('dialog button:has-text("No")');
    ok(await waitFor(async () => (await statusText(dev)).includes('Clear operation canceled.')), 'clear save memory can be cancelled');

    // Recoverable save error, then a successful retry
    await tab(dev, 'settings');
    await intensity.fill('100');
    await intensity.press('Tab');
    fw.failNextSave = true;
    await dev.click('.save-button');
    await dev.click('dialog button.primary');
    ok(await waitFor(() => dev.locator('dialog', { hasText: 'Save not confirmed' }).count().then((n) => n === 1), 10000), 'save error dialog');
    await dev.click('dialog button.primary');
    ok((await saveLabel(dev)) === 'Save and Send Settings', 'Save stays available after a failed save');
    // Reconnecting: the board reports that the save is still pending (Qt DiffUpdate on reconnect).
    await dev.reload();
    ok(await waitFor(() => loaded(dev)), 'docked again with a pending save');
    ok(await waitFor(() => dev.locator('dialog', { hasText: 'Save not confirmed' }).count().then((n) => n === 1)), 'pending save reported after reconnecting');
    await dev.click('dialog button.primary');
    ok(await waitFor(async () => (await saveLabel(dev)) === 'Save and Send Settings'), 'Save available for the pending save');
    await tab(dev, 'settings');
    await dev.locator('.spin', { hasText: '/ 255' }).locator('input').fill('100');
    await dev.locator('.spin', { hasText: '/ 255' }).locator('input').press('Tab');
    await tab(dev, 'tests');
    ok(await dev.locator('text=Open IR Camera Tester...').isDisabled(), 'hardware tests disabled until a save is confirmed');
    await dev.click('.save-button');
    await dev.click('dialog button.primary');
    ok(await waitFor(async () => (await statusText(dev)).includes('Sent settings successfully!'), 10000), 'retry saves');
    ok(fw.getSetting(0) === 100 && await dev.locator('text=Open IR Camera Tester...').isEnabled(), 'retry reached the board, tests enabled again');

    // Camera error reported by the board
    fw.camNotAvailable = true;
    await dev.click('text=Open IR Camera Tester...');
    ok(await waitFor(() => dev.locator('dialog', { hasText: 'Device Error: Camera not available!' }).count().then((n) => n === 1)), 'camera error dialog');
    await dev.click('dialog button.primary');
    ok(await waitFor(() => dev.locator('.fullscreen-window').count().then((n) => n === 0)), 'no test window left open after the camera error');
    ok(await waitFor(async () => (await saveLabel(dev)) === '[Nothing To Save]' && fw.runMode === 'normal'), 'IR test mode ended after the camera error');
    fw.camNotAvailable = false;
    await dev.keyboard.press('Escape');

    // Theme: its own button beside the language selector, not the View menu any more
    await dev.click('.theme-menu .menu-button');
    await dev.click('.menu-item:has-text("Dark Theme")');
    ok(await dev.evaluate(() => document.documentElement.dataset.theme === 'dark'), 'dark theme selected');
    ok(await dev.evaluate(() => document.querySelector('.theme-button').title === 'Dark Theme'), 'the theme button names the theme in use');
    await dev.getByRole('button', { name: 'View', exact: true }).click();
    ok(await dev.locator('.menu-panel:not([hidden]) .menu-item:has-text("Theme")').count() === 0, 'no themes left in the View menu');
    await dev.keyboard.press('Escape');
    await dev.reload();
    ok(await dev.evaluate(() => document.documentElement.dataset.theme === 'dark'), 'theme remembered');
    ok(await waitFor(() => loaded(dev)), 'docked after reload');

    // Phone width: nothing wider than the screen
    await dev.setViewportSize({ width: 390, height: 844 });
    for (const id of ['pins', 'buttons', 'settings', 'profiles', 'tests']) {
        await tab(dev, id);
        const overflow = await dev.evaluate(() => document.documentElement.scrollWidth - window.innerWidth);
        ok(overflow <= 0, `phone width: ${id} fits the screen (${overflow})`);
    }
    ok(await dev.locator('#tab-tests .board-actions button').evaluateAll((nodes) => nodes.every((n) => n.getBoundingClientRect().right <= innerWidth)),
        'phone width: board action buttons inside the screen');
    await tab(dev, 'settings');
    ok(await dev.locator('#tab-settings .desc-toggle').isVisible() && await dev.locator('#tab-settings .desc-text').isHidden(),
        'phone width: the description is a closed bar');
    await dev.click('#tab-settings .desc-toggle');
    ok(await dev.locator('#tab-settings .desc-text').isVisible(), 'phone width: tapping the bar opens the description');
    await shot(dev, 'dev-phone-tests');
    ok(devErrors.length === 0, 'no dev errors ' + JSON.stringify(devErrors));
    await dev.close();
    await sim.close();

    // ===================== site (Web Serial) =====================
    sim = await startServer({ port: 8124, root: LIGHTGUN + '/dist/site', board: 'esp32-s3-devkitc-1' });
    const siteContext = await browser.newContext({ locale: 'en-US', viewport: { width: 1280, height: 860 } });
    const serial = await installFakeSerial(siteContext, sim);
    const site = await siteContext.newPage();
    const siteErrors = watchErrors(site);
    await site.goto('http://localhost:8124/');
    ok(await site.locator('.welcome .big-button').isVisible(), 'site shows the Connect button');

    await site.click('.menu-button:has-text("Board Previews")');
    const options = await site.locator('.preview-select option:not([disabled])').count();
    ok(options === 11, 'compatible boards list: ' + options);
    const order = await site.locator('.preview-select option:not([disabled])').evaluateAll((nodes) => nodes.map((n) => n.value));
    ok(order.join() === order.slice().sort().join(), 'boards listed by board name like the Qt std::map');
    ok(await site.evaluate(() => !document.querySelector('.preview-dialog').matches(':modal')), 'Boards Previewer is a window, not modal');
    for (const board of await site.evaluate(() => OF.Boards.list())) {
        await site.selectOption('.preview-select', board);
        await waitFor(() => site.evaluate((b) => { const h = document.querySelector('.previewer .board-picture'); return h && h.dataset.board === b && !!h.shadowRoot?.querySelector('svg'); }, board));
        const r = await site.evaluate(() => {
            const host = document.querySelector('.previewer .board-picture');
            const items = [...document.querySelectorAll('.previewer .pin-item')];
            const dark = [];
            for (const item of items) {
                const gpio = Number(item.querySelector('.gpio-label').textContent.match(/\d+/)[0]);
                item.dispatchEvent(new MouseEvent('mouseenter'));
                const pin = host.shadowRoot.getElementById('OF_pin' + gpio);
                const lit = pin && getComputedStyle(pin).opacity === '1';
                item.dispatchEvent(new MouseEvent('mouseleave'));
                if (!lit || getComputedStyle(pin).opacity !== '0') dark.push(gpio);
            }
            return { items: items.length, dark };
        });
        ok(r.dark.length < r.items, `${board}: ${r.items - r.dark.length}/${r.items} pin circles light up` + (r.dark.length ? ` (no OF_pin circle for GPIO ${r.dark.join(',')})` : ''));
    }
    await site.selectOption('.preview-select', 'esp32-s3-devkitc-1');
    await waitFor(() => site.evaluate(() => !!document.querySelector('.previewer .board-picture').shadowRoot?.querySelector('svg')));
    const previewRects = () => site.evaluate(() => [...document.querySelectorAll('.previewer .pin-item, .previewer .gpio-label')]
        .map((node) => { const r = node.getBoundingClientRect(); return [r.left, r.top, r.width].map(Math.round).join(','); }).join(';'));
    const previewBefore = await previewRects();
    let previewMoves = 0;
    for (const item of await site.$$('.previewer .pins-col .pin-item')) {
        await item.hover();
        if ((await previewRects()) !== previewBefore) previewMoves++;
    }
    ok(previewMoves === 0, 'previewer: pointing at a pin does not move the columns (' + previewMoves + ')');
    ok(await site.evaluate(() => {
        const S = OF.Boards.shared;
        const P = S.boardBoxPositions_e;
        const positions = S.boardsBoxPositions['esp32-s3-devkitc-1'];
        return [...document.querySelectorAll('.previewer .pins-col .pin-item')].every((item) => {
            const gpio = Number(item.querySelector('.gpio-label').textContent.match(/\d+/)[0]);
            return getComputedStyle(item).gridRowStart === String(positions[gpio] & ~P.posCheck);
        });
    }), 'previewer: pins on the rows of boardsBoxPositions');
    await site.selectOption('.preview-select', 'waveshare-esp32-s3-zero');
    await shot(site, 'site-preview');
    await site.keyboard.press('Escape');

    // A wireless pedal answered at start-up and the user asked for one: it has no pin,
    // and must be neither described as not connected nor left unconfigurable.
    sim.firmware.pedalWireless = true;
    sim.firmware.toggles[sim.firmware.S.boolTypes_e.pedalWireless] = 1;

    await site.click('.welcome .big-button');
    ok(await waitFor(() => loaded(site)), 'site docks over Web Serial');
    ok(sim.firmware.sessionLink === 'serial', 'session on the serial link');
    ok(serial.opens[0] === 9600, 'port opened at 9600 baud');
    ok(await site.locator('.device-bar').count() === 0, 'no device selector row on the site either');
    ok(await site.locator('.disconnect-button').isVisible(), 'Disconnect appears in the status bar while docked');
    ok(await site.locator('.link-state.on').count() === 1, 'the plug shows the link as connected');
    await shot(site, 'site-docked');
    await site.click('.disconnect-button');
    ok(await waitFor(async () => !(await loaded(site))), 'Disconnect undocks');
    ok(await waitFor(() => !sim.firmware.sessionActive && sim.firmware.gunMode === 'run'), 'gun back to Run mode');
    ok(await site.locator('.disconnect-button').isHidden(), 'Disconnect hidden when no board is docked');
    ok(await site.locator('.link-state').isVisible() && await site.locator('.link-state.on').count() === 0,
        'the plug stays in the status bar and shows the link as gone');

    serial.busy = true;
    await site.click('.welcome .big-button');
    ok(await waitFor(() => site.locator('dialog', { hasText: 'Serial port is already in use!' }).count().then((n) => n === 1)), 'busy port dialog');
    await site.click('dialog button.primary');
    serial.busy = false;

    await site.click('.welcome .big-button');
    ok(await waitFor(() => loaded(site)), 'Connect docks again after the busy port');
    await tab(site, 'tests');
    const pedalText = async () => (await site.locator('.test-button', { hasText: 'Pedal' }).first().innerText()).replace(/\s+/g, ' ');
    ok((await pedalText()) === 'Pedal (wireless)', 'a paired wireless pedal says so: ' + await pedalText());
    // The same pedal with nothing answering: it says it is not connected, not '(N/C)'.
    await site.evaluate(() => { OF.app.state.board.pedalWireless = false; OF.app.tabs.tests.resetReadings(); });
    ok((await pedalText()) === 'Pedal (wireless, disconnected)', 'an absent wireless pedal says so: ' + await pedalText());
    await site.evaluate(() => { OF.app.state.board.pedalWireless = true; OF.app.tabs.tests.resetReadings(); });
    await tab(site, 'buttons');
    ok(await site.evaluate(() => {
        const rows = [...document.querySelectorAll('#tab-buttons .btn-row')];
        const row = rows.find((r) => (r.querySelector('.btn-name') || {}).textContent === 'Pedal');
        return !!row && !row.disabled;
    }), 'the wireless pedal can be configured although it has no pin');

    // Output the board starts from. Value 0 is the absolute mouse, i.e. what every
    // build did before this setting existed, so that is where it has to start.
    await tab(site, 'settings');
    const bootBox = 'select[aria-label="Startup Mode"]';
    const bootIndex = await site.evaluate(() => OF.Boards.shared.settingsTypes_e.bootOutputMode);
    ok(await site.locator(bootBox).count() === 1, 'the startup output box is in the Input group');
    ok((await site.locator(bootBox).inputValue()) === '0', 'it starts on Absolute Mouse');
    await site.selectOption(bootBox, '2');
    ok(await waitFor(async () => (await saveLabel(site)) === 'Save and Send Settings'), 'changing the startup output enables Save');
    await site.click('.save-button');
    await site.click('dialog button.primary');
    ok(await waitFor(async () => (await statusText(site)).includes('Sent settings successfully!'), 8000), 'startup output saved');
    ok(sim.firmware.getSetting(bootIndex) === 2,
        'firmware received the startup output mode: ' + sim.firmware.getSetting(bootIndex));
    await site.selectOption(bootBox, '0');
    await site.click('.save-button');
    await site.click('dialog button.primary');
    ok(await waitFor(() => sim.firmware.getSetting(bootIndex) === 0, 8000), 'and back to Absolute Mouse');
    await tab(site, 'tests');
    serial.closes.length = 0;
    await site.click('text=Restart Microcontroller in Firmware Update Mode');
    ok(await waitFor(() => sim.firmware.rebootedToBootloader), 'ESP32 restart command reaches the board');
    // The port used to be closed straight after the write, which on a real board threw the
    // restart command away every now and then.
    ok(await waitFor(() => serial.closes.indexOf('flush') >= 0 &&
        (serial.closes.indexOf('port') < 0 || serial.closes.indexOf('flush') < serial.closes.indexOf('port'))),
        'the written bytes are flushed before the port is closed (' + serial.closes.join(' -> ') + ')');
    ok((await site.locator('.status-text').innerText()).length > 0, 'the restart is reported in the status bar');
    ok(await waitFor(async () => !(await loaded(site))), 'page undocked after the restart');
    ok(siteErrors.length === 0, 'no site errors ' + JSON.stringify(siteErrors));
    await site.close();
    await sim.close();

    // ===================== the published site: home and versions =====================
    // What is published is one folder per version of the App (v/<version>/, a copy of what
    // the build writes in dist/site) and, at the root, the home page with the Connect
    // button, which reads the version of the firmware and opens the App published for it.
    // The site is put together here out of the two builds, as publishing does.
    const SITE = path.join(LIGHTGUN, 'dist', 'test-site');
    const copyInto = (from, to) => {
        fs.mkdirSync(to, { recursive: true });
        for (const entry of fs.readdirSync(from, { withFileTypes: true })) {
            const source = path.join(from, entry.name);
            const target = path.join(to, entry.name);
            if (entry.isDirectory()) copyInto(source, target);
            else fs.copyFileSync(source, target);
        }
    };
    fs.rmSync(SITE, { recursive: true, force: true });
    copyInto(path.join(LIGHTGUN, 'dist', 'launcher'), SITE);
    copyInto(path.join(LIGHTGUN, 'dist', 'site'), path.join(SITE, 'v', '6.2'));
    // An older version, as it would have been published in its day.
    copyInto(path.join(LIGHTGUN, 'dist', 'site'), path.join(SITE, 'v', '6.1'));
    const oldApp = path.join(SITE, 'v', '6.1', 'app.js');
    fs.writeFileSync(oldApp, fs.readFileSync(oldApp, 'utf8')
        .replace('"version": "6.2"', '"version": "6.1"')
        .replace('"versionLabel": "6.2.0"', '"versionLabel": "6.1.0"'));
    fs.writeFileSync(path.join(SITE, 'versions.json'), JSON.stringify({
        latest: '6.2',
        versions: [{ id: '6.2', label: '6.2.0', type: 'stable' },
                   { id: '6.1', label: '6.1.0', type: 'stable' }]
    }));

    sim = await startServer({ port: 8126, root: SITE, board: 'esp32-s3-devkitc-1' });
    const verContext = await browser.newContext({ locale: 'en-US', viewport: { width: 1280, height: 860 } });
    await installFakeSerial(verContext, sim);
    const mismatch = (page) => page.locator('dialog', { hasText: 'Versions do not match' });
    const CMD = globalThis.OpenFIREshared.serialCmdTypes_e;
    const settingsAsked = () => sim.firmware.log.filter((e) => e.dir === 'app->fw' &&
        [CMD.sGetToggles, CMD.sGetSettings, CMD.sGetProfile, CMD.sGetBtns].indexOf(e.command) >= 0).length;

    // ----- the home page opens the App of the firmware -----
    sim.firmware.version = '6.2-abcdef0';
    sim.firmware.versionFull = '6.2.0-stable';
    let home = await verContext.newPage();
    const homeErrors = watchErrors(home);
    await home.goto('http://localhost:8126/');
    ok(!(await home.locator('.welcome .big-button').count()), 'the home of the site is not the App');
    ok(await home.locator('#connect').isVisible(), 'it is a page with the Connect button');
    ok(await home.locator('#theme-button').isVisible() && await home.locator('#lang-select').isVisible(),
        'with the theme and the language, like the other pages of the project');
    await home.selectOption('#lang-select', 'it');
    ok(await waitFor(async () => (await home.locator('#connect-label').innerText()).includes('Collega')),
        'the language changes the page there and then');
    ok(home.url().includes('lang=it'), 'and travels in the address: ' + home.url());
    await home.click('#theme-button');
    await home.click('#theme-menu button[data-theme="dark"]');
    ok(await waitFor(async () => await home.evaluate(() => getComputedStyle(document.body).backgroundColor) === 'rgb(20, 22, 27)'),
        'the theme too, and it is the App\'s own setting');
    ok(await home.evaluate(() => localStorage.getItem('of_theme')) === 'dark', 'remembered in of_theme');
    await home.click('#theme-button');
    await home.click('#theme-menu button[data-theme="system"]');
    await home.selectOption('#lang-select', 'en');
    sim.firmware.log.length = 0;
    await home.click('#connect');
    ok(await home.waitForURL(/\/v\/6\.2\//, { timeout: 15000 }).then(() => true).catch(() => false),
        'Connect opens the App of the firmware: ' + home.url());
    ok(await waitFor(() => loaded(home), 15000), 'and that App docks by itself, without another click');
    ok(await home.evaluate(() => OF.BUILD.version) === '6.2', 'it is the App of 6.2');
    ok(await mismatch(home).count() === 0, 'with nothing to warn about');
    ok(homeErrors.length === 0, 'no errors on the home ' + JSON.stringify(homeErrors));
    await home.close();

    // ----- an older firmware gets its own App -----
    sim.firmware.version = '6.1-abcdef0';
    sim.firmware.versionFull = '6.1.0-stable';
    home = await verContext.newPage();
    await home.goto('http://localhost:8126/');
    sim.firmware.log.length = 0;
    await home.click('#connect');
    ok(await home.waitForURL(/\/v\/6\.1\//, { timeout: 15000 }).then(() => true).catch(() => false),
        'an older firmware opens the App published for it: ' + home.url());
    ok(/[?&]lang=/.test(home.url()), 'and the language travels with it: ' + home.url());
    ok(await waitFor(() => loaded(home), 15000), 'which docks by itself');
    ok(await home.evaluate(() => OF.BUILD.version) === '6.1', 'it really is the App of 6.1');
    ok(await mismatch(home).count() === 0, 'and says nothing, because now they match');
    const notice = await statusText(home);
    ok(notice.includes('6.1.0') && notice.includes('6.2.0'),
        'the status bar says which App this is and that a newer firmware exists: ' + JSON.stringify(notice));
    await shot(home, 'site-home-opened-6.1');
    await home.close();

    // ----- a firmware nobody published an App for -----
    sim.firmware.version = '6.9-abcdef0';
    sim.firmware.versionFull = '6.9.0-beta';
    home = await verContext.newPage();
    await home.goto('http://localhost:8126/');
    await home.click('#connect');
    ok(await waitFor(async () => (await home.locator('#state').innerText()).includes('6.9.0-beta'), 15000),
        'a firmware without a published App is reported: ' + JSON.stringify((await home.locator('#state').innerText()).slice(0, 120)));
    await sleep(700);
    ok(!/\/v\//.test(home.url()), 'nothing is opened by itself: ' + home.url());
    ok(await home.locator('#versions .item').count() === 2, 'the published versions are offered instead');
    await shot(home, 'site-home-not-published');
    // and one of them can be tried by hand: it is the one that then says the versions differ
    await home.click('#versions .item >> nth=0');
    ok(await home.waitForURL(/\/v\/6\.2\//, { timeout: 15000 }).then(() => true).catch(() => false),
        'choosing one by hand opens it: ' + home.url());
    await home.click('.welcome .big-button');
    ok(await waitFor(() => mismatch(home).count().then((n) => n === 1), 15000),
        'and that App says the versions do not match');
    await home.close();

    // ----- the question an App of another version asks -----
    sim.firmware.version = '6.1-abcdef0';
    sim.firmware.versionFull = '6.1.0-stable';
    let app = await verContext.newPage();
    const appErrors = watchErrors(app);
    await app.goto('http://localhost:8126/v/6.2/');          // opened by its own address
    ok(await waitFor(() => app.evaluate(() => !!(window.OF && OF.app))), 'a published version opens on its own');
    await app.click('.welcome .big-button');
    ok(await waitFor(() => mismatch(app).count().then((n) => n === 1), 15000),
        'and asks what to do when the firmware is of another version');
    const question = (await mismatch(app).innerText()).replace(/\s+/g, ' ');
    ok(question.includes('6.2.0') && question.includes('6.1.0-stable'),
        'naming both versions: ' + JSON.stringify(question.slice(0, 130)));
    ok(await mismatch(app).locator('button', { hasText: 'Carry on' }).count() === 1 &&
       await mismatch(app).locator('button', { hasText: 'Go back' }).count() === 1,
        'with Carry on and Go back');
    await shot(app, 'site-version-question');
    await app.click('dialog button:has-text("Carry on")');
    ok(await loaded(app), 'Carry on keeps the App and the connection');
    await sleep(500);
    ok(/\/v\/6\.2\//.test(app.url()), 'and stays where it was: ' + app.url());
    // it does not ask again at every reconnection
    await app.click('.disconnect-button');
    ok(await waitFor(async () => !(await loaded(app))), 'undocked');
    await app.click('.welcome .big-button');
    ok(await waitFor(() => loaded(app)), 'docked again');
    await sleep(700);
    ok(await mismatch(app).count() === 0, 'and it is asked once, not at every connection');
    await app.close();

    // Go back: the lightgun is undocked and the home page is opened again.
    app = await verContext.newPage();
    await app.goto('http://localhost:8126/v/6.2/');
    await app.click('.welcome .big-button');
    ok(await waitFor(() => mismatch(app).count().then((n) => n === 1), 15000), 'asked again on a new page');
    await app.click('dialog button:has-text("Go back")');
    ok(await app.waitForURL(/8126\/($|\?)/, { timeout: 15000 }).then(() => true).catch(() => false),
        'Go back returns to the home page: ' + app.url());
    ok(await app.locator('#connect').isVisible(), 'where the Connect button is');
    await app.close();

    // Closing the window is going back too.
    app = await verContext.newPage();
    await app.goto('http://localhost:8126/v/6.2/');
    await app.click('.welcome .big-button');
    ok(await waitFor(() => mismatch(app).count().then((n) => n === 1), 15000), 'asked again');
    await app.click("dialog .dialog-close");
    ok(await app.waitForURL(/8126\/($|\?)/, { timeout: 15000 }).then(() => true).catch(() => false),
        'the X of the window goes back as well: ' + app.url());
    ok(appErrors.length === 0, 'no errors in the published App ' + JSON.stringify(appErrors));
    await app.close();

    // ----- when the versions match, nothing is asked -----
    sim.firmware.version = '6.2-abcdef0';
    sim.firmware.versionFull = '6.2.0-stable';
    app = await verContext.newPage();
    await app.goto('http://localhost:8126/v/6.2/');
    await app.click('.welcome .big-button');
    ok(await waitFor(() => loaded(app)), 'the App of the same version docks');
    await sleep(700);
    ok(await mismatch(app).count() === 0, 'and nothing is asked');
    ok(settingsAsked() > 0, 'the settings are read, as always');
    await app.close();
    await verContext.close();

    // The language chosen on the other pages of the project travels in the address.
    for (const [query, label, expected] of [['?lang=it', 'italian', 'Collega una lightgun'],
                                            ['', 'the browser language', 'Connect a Lightgun']]) {
        const langContext = await browser.newContext({ locale: 'en-US', viewport: { width: 1000, height: 800 } });
        const langPage = await langContext.newPage();
        await langPage.goto('http://localhost:8126/' + query);
        ok(await waitFor(async () => (await langPage.locator('#connect').innerText()).includes(expected)),
            `?lang: the home opens in ${label} (${JSON.stringify(query)})`);
        await langContext.close();
    }
    await sim.close();

    // ===================== opened as local files (file://) =====================
    for (const [label, folder] of [['site', LIGHTGUN + '/dist/site'], ['webapp folder', LIGHTGUN + '/webapp']]) {
        const local = await siteContext.newPage();
        const localErrors = watchErrors(local);
        await local.goto('file://' + folder + '/index.html');
        ok(await waitFor(() => local.evaluate(() => !!(window.OF && OF.app))), `${label} from file://: the app starts`);
        await local.click('.menu-button:has-text("Board Previews")');
        await local.selectOption('.preview-select', 'rpipico');
        ok(await waitFor(() => local.evaluate(() => !!document.querySelector('.previewer .board-picture').shadowRoot?.querySelector('svg'))),
            `${label} from file://: board picture loaded from boards/pics/rpipico.js`);
        const lit = await local.evaluate(() => {
            const item = [...document.querySelectorAll('.previewer .pin-item')].find((node) => /«GPIO0»/.test(node.textContent));
            item.dispatchEvent(new MouseEvent('mouseenter'));
            const pin = document.querySelector('.previewer .board-picture').shadowRoot.getElementById('OF_pin0');
            return !!pin && getComputedStyle(pin).opacity === '1';
        });
        ok(lit, `${label} from file://: pointing at a pin lights its circle`);
        ok(localErrors.length === 0, `${label} from file://: no errors ` + JSON.stringify(localErrors));
        await local.close();
    }

    await browser.close();
    process.exit(process.exitCode || 0);
})().catch((e) => { console.error(e); process.exit(1); });
