/*  OpenFIRE Web App - polling rate test ("Measure Polling Rate", welcome screen).

    Measures how often the lightgun really sends its aim to the computer, before the App
    connects: the gun must be in Run mode (a docked gun stops reporting the aim) and its output
    must be the mouse (absolute pointer). Nothing goes to the gun: the page only listens to the
    pointer events of the browser.

    Every new aim is an HID report. The browser joins the reports of one display frame into one
    pointermove event; getCoalescedEvents() gives them back one by one with their own time, so
    the rate measured is the one of the reports (about the frame rate of the IR camera, 209 fps,
    or what reaches the computer through the dongle), not the refresh rate of the display.

    A gun held still sends no new aim, and an aim outside the window or the test area sends none
    to the page: a gap longer than PAUSE_MS ends a "burst" of movement and is not counted, so
    the averages stay those of the movement.

        const stats = new OF.PollingStats();      // pure, also used by tests/polling-test.test.js
        stats.addReport(timeMs);                  // each HID report (coalesced event)
        stats.sample(nowMs);                      // every display frame: graph, min and max
        stats.snapshot(nowMs);                    // { current, average (3 s), min, max }

        const view = new OF.PollingTest({ onExit });
        container.append(view.root);
        view.start();  ...  view.stop();
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};
    const doc = root.document;

    const CAMERA_FPS = 209;                       // OpenFIRECameraProfile.h: DFRobot and PAJ7025
    const PAUSE_MS = 60;                          // longer gap: the gun stopped or left the area
    const WINDOW_MS = 500;                        // "current" rate: reports of the last half second
    const AVERAGE_MS = 3000;                      // average: the last 3 s of movement
    const SETTLE_MS = 150;                        // a burst needs this much movement before it is shown
    const SAMPLE_MS = 100;                        // one point of the graph / of min and max
    const HISTORY_S = 10;                         // seconds of movement shown by the graph
    const HISTORY_SAMPLES = HISTORY_S * 1000 / SAMPLE_MS;

    /** Statistics of the report times (ms, performance.now() timeline). No DOM. */
    class PollingStats {
        constructor(options = {}) {
            this.pauseMs = options.pauseMs || PAUSE_MS;
            this.reset();
        }

        reset() {
            this.reports = 0;               // all the reports, also the first of each burst
            this.count = 0;                 // intervals ever counted (within the bursts)
            this.last = null;               // time of the last report
            this.recent = [];               // report times of the last WINDOW_MS, for "current"
            this.window = [];               // intervals of the last AVERAGE_MS of movement
            this.windowSum = 0;
            this.samples = [];              // { hz, resumed } every SAMPLE_MS of movement: the graph
            this.lastSample = null;
            this.resume = false;
        }

        /** One HID report at `time` (ms). Times that go back or repeat are ignored. */
        addReport(time) {
            if (!Number.isFinite(time)) return;
            if (this.last !== null && time <= this.last) return;
            this.reports++;
            if (this.last === null || time - this.last > this.pauseMs) {
                // First report of a burst: the gap before it is a pause, not an interval.
                this.recent = [];
            } else {
                this._addInterval(time - this.last);
            }
            this.last = time;
            this.recent.push(time);
            let drop = 0;
            while (drop < this.recent.length && this.recent[drop] < time - WINDOW_MS) drop++;
            if (drop) this.recent.splice(0, drop);
        }

        _addInterval(interval) {
            this.count++;
            this.window.push(interval);
            this.windowSum += interval;
            // Only time spent moving fills the window: a pause does not push old reports out.
            let drop = 0;
            while (this.windowSum - (this.window[drop] || 0) >= AVERAGE_MS && drop < this.window.length - 1) {
                this.windowSum -= this.window[drop];
                drop++;
            }
            if (drop) this.window.splice(0, drop);
        }

        /** Whether the gun is sending reports now. */
        moving(now) {
            return this.last !== null && now - this.last <= this.pauseMs * 2;
        }

        /** Rate of the reports of the last WINDOW_MS of the current burst (Hz), null if too few. */
        current(now) {
            if (!this.moving(now) || this.recent.length < 3) return null;
            // The first reports of a burst are too few for a rate (it would jump up and down).
            const span = this.recent[this.recent.length - 1] - this.recent[0];
            if (span < SETTLE_MS) return null;
            return (this.recent.length - 1) * 1000 / span;
        }

        /** Called by the view at every display frame: one point every SAMPLE_MS while moving.
            While the gun is still nothing is added, so the graph stops and goes on afterwards. */
        sample(now) {
            const hz = this.current(now);
            if (hz === null) {
                if (this.lastSample !== null) this.resume = true;   // next point: after a pause
                this.lastSample = null;
                return;
            }
            if (this.lastSample !== null && now - this.lastSample < SAMPLE_MS) return;
            this.lastSample = now;
            // Only windows that cover enough of the movement count for min and max.
            const full = this.recent[this.recent.length - 1] - this.recent[0] >= WINDOW_MS * 0.6;
            this.samples.push({ hz, full, resumed: !!this.resume && this.samples.length > 0 });
            this.resume = false;
            if (this.samples.length > HISTORY_SAMPLES) this.samples.splice(0, this.samples.length - HISTORY_SAMPLES);
        }

        /** Lowest and highest rate of the graph (last 10 s of movement): a peak or a drop goes
            away once it leaves the graph, e.g. a mouse used before the lightgun. */
        range() {
            let min = null;
            let max = null;
            for (const sample of this.samples) {
                if (!sample.full) continue;
                min = min === null ? sample.hz : Math.min(min, sample.hz);
                max = max === null ? sample.hz : Math.max(max, sample.hz);
            }
            return { min, max };
        }

        snapshot(now) {
            const n = this.window.length;
            const mean = n ? this.windowSum / n : null;
            const range = this.range();
            return {
                moving: this.moving(now),
                current: this.current(now),
                average: mean ? 1000 / mean : null,      // last 3 s of movement
                min: range.min,               // last 10 s of movement, like the graph
                max: range.max,
                reports: this.reports,
            };
        }
    }

    // ----- View -----------------------------------------------------------------------

    const CURRENT_MS = 150;                       // the big number: lively but readable
    const NUMBERS_MS = 250;                       // the boxes: 4 times a second, not at every frame

    // Average (Hz): green from 180 Hz (about 15% under the 209 fps of the camera), orange from 60 Hz (the refresh of an ordinary display),
    // red below. Indicative thresholds.
    const AVERAGE_OK_HZ = 180;
    const AVERAGE_WARN_HZ = 60;
    const averageLevel = (hz) => hz === null || !Number.isFinite(hz) ? null :
        hz >= AVERAGE_OK_HZ ? 'ok' : hz >= AVERAGE_WARN_HZ ? 'warn' : 'bad';

    /** Splits the stretch (i0, hz0)-(i1, hz1) where it crosses the colour thresholds:
        [[fromIndex, fromHz, toIndex, toHz], ...], each piece of one colour only. */
    function splitAtThresholds(i0, hz0, i1, hz1) {
        const cuts = [AVERAGE_OK_HZ, AVERAGE_WARN_HZ]
            .filter((limit) => (hz0 - limit) * (hz1 - limit) < 0)
            .map((limit) => (limit - hz0) / (hz1 - hz0))
            .sort((a, b) => a - b);
        const points = [[i0, hz0], ...cuts.map((f) => [i0 + (i1 - i0) * f, hz0 + (hz1 - hz0) * f]), [i1, hz1]];
        const pieces = [];
        for (let k = 1; k < points.length; k++) pieces.push([points[k - 1][0], points[k - 1][1], points[k][0], points[k][1]]);
        return pieces;
    }

    const fixed = (value, digits) => value === null || value === undefined || !Number.isFinite(value) ? '—' : value.toFixed(digits);

    function cssColor(node, name, fallback) {
        const value = root.getComputedStyle ? root.getComputedStyle(node).getPropertyValue(name).trim() : '';
        return value || fallback;
    }

    class PollingTest {
        constructor(options = {}) {
            this.onExit = options.onExit || (() => {});
            this.stats = options.stats || new PollingStats();   // kept when the page is rebuilt (language)
            this.running = false;
            this.inside = false;
            this.mouseSeen = false;
            this.otherPointer = false;
            this._frame = null;
            this._onMove = (event) => this.onPointerMove(event);
            this._onLeave = () => { this.inside = false; };
            this._onEnter = () => { this.inside = true; };
            this._onKey = (event) => this.onKey(event);
            // The buttons of the gun must not trigger the browser while testing: no context menu
            // (right button, e.g. A or the off-screen shot), no auto-scroll (middle button), no
            // page back or forward (buttons 4 and 5). The left button stays as it is.
            this._onContext = (event) => event.preventDefault();
            this._onButton = (event) => { if (event.button !== 0) event.preventDefault(); };
            this.build();
        }

        build() {
            const { el, t, icon } = OF.UI;
            const tile = (key, label, unit) => {
                const value = el('span', { class: 'polling-tile-value', text: '—' });
                this.tiles[key] = value;
                return el('div', { class: 'polling-tile' },
                    el('span', { class: 'polling-tile-label', text: label }),
                    el('span', { class: 'polling-tile-row' }, value, unit ? el('span', { class: 'polling-tile-unit', text: unit }) : null));
            };
            this.tiles = {};

            const header = el('div', { class: 'polling-header' },
                el('h2', { class: 'polling-title' }, icon('gauge'), el('span', { text: t('Polling Rate Test') })),
                el('div', { class: 'polling-tools' },
                    el('button', { class: 'polling-tool', on: { click: () => this.reset() } }, el('span', { text: t('Reset') })),
                    el('button', { class: 'polling-tool primary', on: { click: () => this.onExit() } }, el('span', { text: t('Back') }))));

            this.hint = el('p', { class: 'polling-hint' },
                t('Aim the lightgun inside this area and keep moving it. The lightgun must be in Run mode (not connected to the App) with mouse output.'),
                ' ', el('strong', { text: t('Expected value with OpenFIRE: about %1 Hz, the frame rate of the camera (with the USB cable and with the dongle).', CAMERA_FPS) }));

            this.state = el('span', { class: 'polling-state' });
            this.currentValue = el('span', { class: 'polling-current-value', text: '—' });
            const hero = el('div', { class: 'polling-hero' },
                el('div', { class: 'polling-current' }, this.currentValue, el('span', { class: 'polling-current-unit', text: 'Hz' })),
                this.state);

            const tiles = el('div', { class: 'polling-tiles' },
                tile('average', t('Average (last 3 s)'), 'Hz'),
                tile('max', t('Maximum (last 10 s)'), 'Hz'),
                tile('min', t('Minimum (last 10 s)'), 'Hz'));
            this.tiles.average.parentNode.parentNode.title = t('Rate of the last 3 s of movement: green from 180 Hz, orange from 60 Hz, red below. Indicative thresholds.');

            this.rateCanvas = el('canvas', { class: 'polling-chart' });
            const charts = el('div', { class: 'polling-charts' },
                el('figure', { class: 'polling-figure' },
                    el('figcaption', { text: t('Rate over the last %1 s of movement (line: %2 fps of the camera)', HISTORY_S, CAMERA_FPS) }),
                    el('div', { class: 'polling-chart-box' }, this.rateCanvas)));

            // Decided once, before the test starts, so the layout never moves during the test.
            const coalesced = !!(root.PointerEvent && root.PointerEvent.prototype &&
                typeof root.PointerEvent.prototype.getCoalescedEvents === 'function');
            this.note = el('p', { class: 'polling-note', hidden: coalesced,
                text: t('This browser does not give the single reports (coalesced events): the rate measured is limited by the refresh of the display. Use Chrome or Edge.') });

            this.area = el('div', { class: 'polling-area', attrs: { tabindex: '-1' } }, el('div', { class: 'polling-top' }, hero, tiles), charts);
            this.root = el('div', { class: 'polling-test' }, header, this.hint, this.note, this.area);
        }

        start() {
            if (this.running) return;
            this.running = true;
            this.render(root.performance.now());
            // The whole test (header included) listens: the aim does not have to stay in the box.
            this.root.addEventListener('pointermove', this._onMove);
            this.root.addEventListener('pointerleave', this._onLeave);
            this.root.addEventListener('pointerenter', this._onEnter);
            // Capture phase: the test looks at the menus and dialogs before their own Escape closes them.
            doc.addEventListener('keydown', this._onKey, true);
            this.root.addEventListener('contextmenu', this._onContext);
            this.root.addEventListener('mousedown', this._onButton);
            this.root.addEventListener('mouseup', this._onButton);
            this.root.addEventListener('auxclick', this._onButton);
            this._loop();
        }

        stop() {
            if (!this.running) return;
            this.running = false;
            this.root.removeEventListener('pointermove', this._onMove);
            this.root.removeEventListener('pointerleave', this._onLeave);
            this.root.removeEventListener('pointerenter', this._onEnter);
            doc.removeEventListener('keydown', this._onKey, true);
            this.root.removeEventListener('contextmenu', this._onContext);
            this.root.removeEventListener('mousedown', this._onButton);
            this.root.removeEventListener('mouseup', this._onButton);
            this.root.removeEventListener('auxclick', this._onButton);
            if (this._frame !== null) root.cancelAnimationFrame(this._frame);
            this._frame = null;
        }

        reset() {
            this.stats.reset();
            this.mouseSeen = false;
            this.otherPointer = false;
            this.render(root.performance.now(), true);
        }

        /** Esc and R are the test's only when nothing is over it: a dialog (a modal one, or the one
            the key is pressed in), an open menu, a fullscreen window (it makes the page behind
            inert) or a text field keep their keys. */
        ownsKey(event) {
            if (event.defaultPrevented || this.root.closest('[inert]')) return false;
            if (doc.querySelector('dialog[open]:not(.modeless)')) return false;   // a modal dialog is open
            const target = event.target;
            if (target && typeof target.closest === 'function' &&
                target.closest('dialog, [role="menu"], input, textarea, select, [contenteditable]:not([contenteditable="false"])')) return false;
            return !doc.querySelector('[aria-haspopup][aria-expanded="true"]');
        }

        onKey(event) {
            if (!this.ownsKey(event)) return;
            if (event.key === 'Escape') { event.preventDefault(); this.onExit(); }
            else if ((event.key === 'r' || event.key === 'R') && !event.ctrlKey && !event.metaKey && !event.altKey) this.reset();
        }

        onPointerMove(event) {
            this.inside = true;
            if (event.pointerType && event.pointerType !== 'mouse') { this.otherPointer = true; return; }
            this.mouseSeen = true;
            if (doc.hidden) return;
            const reports = typeof event.getCoalescedEvents === 'function' ? event.getCoalescedEvents() : [];
            if (reports.length) for (const report of reports) this.stats.addReport(report.timeStamp);
            else this.stats.addReport(event.timeStamp);
        }

        _loop() {
            if (!this.running) return;
            const now = root.performance.now();
            this.stats.sample(now);
            this.render(now);
            this._frame = root.requestAnimationFrame(() => this._loop());
        }

        render(now, force) {
            const t = OF.UI.t;
            const s = this.stats.snapshot(now);
            let state;
            let kind;
            if (s.moving) { state = t('Measuring'); kind = 'ok'; }
            else if (!s.reports) {
                state = this.otherPointer ? t('Touch and pen are not measured: use the lightgun') : t('Waiting for the lightgun: move it inside the area');
                kind = 'idle';
            } else if (!this.inside) { state = t('Aim outside the area'); kind = 'warn'; }
            else { state = t('Paused: keep moving the lightgun'); kind = 'warn'; }
            this.state.textContent = state;
            this.state.dataset.kind = kind;

            // Numbers: a few times a second, or they flicker too fast to read; at once when the
            // measure starts or stops (the graph is still drawn at every frame).
            const measuring = s.current !== null;
            const changed = force || measuring !== this._shownMeasuring;
            this._shownMeasuring = measuring;
            if (changed || this._currentAt === undefined || now - this._currentAt >= CURRENT_MS) {
                this._currentAt = now;
                this.currentValue.textContent = s.current === null ? '—' : String(Math.round(s.current));
            }
            if (changed || this._shownAt === undefined || now - this._shownAt >= NUMBERS_MS) {
                this._shownAt = now;
                this.renderNumbers(s);
            }
            this.drawRate();
        }

        renderNumbers(s) {
            const t = OF.UI.t;
            const tiles = this.tiles;
            tiles.average.textContent = fixed(s.average, 0);
            const averageTile = tiles.average.parentNode.parentNode;
            const average = averageLevel(s.average);
            if (average) averageTile.dataset.level = average;
            else delete averageTile.dataset.level;
            tiles.max.textContent = fixed(s.max, 0);
            tiles.min.textContent = fixed(s.min, 0);
        }

        _prepare(canvas) {
            const ratio = root.devicePixelRatio || 1;
            const width = Math.max(1, Math.round(canvas.clientWidth));
            const height = Math.max(1, Math.round(canvas.clientHeight));
            if (canvas.width !== Math.round(width * ratio) || canvas.height !== Math.round(height * ratio)) {
                canvas.width = Math.round(width * ratio);
                canvas.height = Math.round(height * ratio);
            }
            const ctx = canvas.getContext('2d');
            if (!ctx) return null;
            ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
            ctx.clearRect(0, 0, width, height);
            return { ctx, width, height };
        }

        drawRate() {
            const area = this._prepare(this.rateCanvas);
            if (!area) return;
            const { ctx, width, height } = area;
            const accent = cssColor(this.root, '--accent', '#2fa6e4');
            const faint = cssColor(this.root, '--faint', '#767f8f');
            const border = cssColor(this.root, '--border', '#363b47');
            const samples = this.stats.samples;
            const top = Math.max(300, ...samples.map((sample) => sample.hz)) * 1.1;
            const pad = { left: 34, right: 6, top: 6, bottom: 14 };
            const w = width - pad.left - pad.right;
            const h = height - pad.top - pad.bottom;
            const y = (hz) => pad.top + h - hz / top * h;
            // Newest point on the right edge, one step per SAMPLE_MS of movement.
            const x = (index) => pad.left + w - (samples.length - 1 - index) / (HISTORY_SAMPLES - 1) * w;

            ctx.font = '10px ' + cssColor(this.root, '--font', 'sans-serif');
            ctx.fillStyle = faint;
            ctx.strokeStyle = border;
            ctx.lineWidth = 1;
            ctx.textAlign = 'right';
            ctx.textBaseline = 'middle';
            const step = top > 1200 ? 500 : top > 600 ? 200 : 100;
            for (let hz = 0; hz <= top; hz += step) {
                ctx.beginPath();
                ctx.moveTo(pad.left, Math.round(y(hz)) + 0.5);
                ctx.lineTo(width - pad.right, Math.round(y(hz)) + 0.5);
                ctx.stroke();
                ctx.fillText(String(hz), pad.left - 4, y(hz));
            }
            // Reference: frame rate of the camera.
            ctx.save();
            ctx.setLineDash([4, 4]);
            ctx.strokeStyle = faint;
            ctx.beginPath();
            ctx.moveTo(pad.left, Math.round(y(CAMERA_FPS)) + 0.5);
            ctx.lineTo(width - pad.right, Math.round(y(CAMERA_FPS)) + 0.5);
            ctx.stroke();
            ctx.restore();

            ctx.strokeStyle = accent;
            ctx.lineWidth = 2;
            ctx.lineJoin = 'round';
            // Where the gun was still: a faint mark, the line goes on without the pause.
            ctx.save();
            ctx.strokeStyle = border;
            ctx.setLineDash([2, 3]);
            for (let index = 0; index < samples.length; index++) {
                if (!samples[index].resumed) continue;
                const px = Math.round(x(index - 0.5)) + 0.5;
                ctx.beginPath();
                ctx.moveTo(px, pad.top);
                ctx.lineTo(px, pad.top + h);
                ctx.stroke();
            }
            ctx.restore();
            // The line takes the colours of the average (green, orange, red) by its height: a
            // stretch that crosses a threshold (180 or 60 Hz) changes colour right where it
            // crosses it, whatever the scale of the graph.
            const colors = {
                ok: cssColor(this.root, '--ok', '#4cc96f'),
                warn: cssColor(this.root, '--mid', '#fb923c'),
                bad: cssColor(this.root, '--danger', '#ff5a6a'),
            };
            let current = null;
            const stroke = (level, fromIndex, fromHz, toIndex, toHz) => {
                if (level !== current) {
                    if (current) ctx.stroke();
                    current = level;
                    ctx.strokeStyle = colors[level] || accent;
                    ctx.beginPath();
                    ctx.moveTo(x(fromIndex), y(fromHz));
                }
                ctx.lineTo(x(toIndex), y(toHz));
            };
            for (let index = 1; index < samples.length; index++) {
                for (const piece of splitAtThresholds(index - 1, samples[index - 1].hz, index, samples[index].hz))
                    stroke(averageLevel((piece[1] + piece[3]) / 2) || 'ok', ...piece);
            }
            if (current) ctx.stroke();
        }
    }

    OF.PollingStats = PollingStats;
    OF.PollingTest = PollingTest;
    OF.PollingTest.CAMERA_FPS = CAMERA_FPS;
    OF.PollingTest.averageLevel = averageLevel;
    OF.PollingTest.splitAtThresholds = splitAtThresholds;
})(typeof window !== 'undefined' ? window : globalThis);
