/*  OpenFIRE Web App - fullscreen window of the aim test ("Open Aim Test...", Gun Tests).

    Checks the aim against a grid with nine reference targets (centre, edge midpoints, corners),
    the buttons and the recoil, in one screen. The gun is in the same test mode as the IR camera
    tester (sIRTest): main.js hands over the aim point (sTestCoords, fields 8-9, in the 1920x1080
    test space) and the buttons (sBtnPressed / sBtnReleased). The trigger flashes the crosshair,
    leaves a fading hit mark and asks main.js for the recoil; any other button shows its name in a
    short semi-transparent circle around the crosshair. Drawn like the GunTester of IAmOrion
    (OpenFIRE-Tools, GunTester/src/Program.cs).

        const win = new OF.AimTestWindow({ onExit, onTrigger, buttonName });
        win.open();                       // inside the click handler (fullscreen needs a user gesture)
        win.drawTest(payload);            // sTestCoords
        win.onButton(index, pressed);     // sBtnPressed / sBtnReleased
        win.shutdown();

    Keys: Esc closes, H hides the key hints, F turns the recoil (solenoid / rumble) off and on.
    The grid always stays.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};
    const doc = root.document;

    const MODE_AIMTEST = 'aimtest';
    const TRIGGER = 0;              // boardInputs_e.btnTrigger

    const FONT = '"Segoe UI", system-ui, -apple-system, Roboto, sans-serif';
    const COLOR_BACKGROUND = 'rgb(14, 14, 18)';
    const COLOR_GRID_MINOR = 'rgb(38, 38, 48)';     // 1 px
    const COLOR_GRID_MAJOR = 'rgb(70, 70, 86)';     // 1.5 px, every tenth of the screen
    const COLOR_BORDER = 'rgb(110, 110, 130)';      // 2 px
    const COLOR_TARGET = 'rgb(150, 150, 170)';      // 2 px
    const COLOR_HINTS = 'rgba(200, 200, 200, 0.47)';
    const COLOR_MESSAGE = 'gainsboro';
    const COLOR_AIM = [255, 70, 70];
    const COLOR_AIM_PARTIAL = [255, 165, 0]; // 2-3 emitters seen: the aim is estimated
    const COLOR_EDGE = 'orange';

    const GRID_LINES = 20;          // lines every twentieth of the screen, major every other one
    const TARGET_INSET = 0.05;      // targets 5% in from the edges
    const TEST_W = 1920;            // sTestCoords space of the aim point
    const TEST_H = 1080;
    const COORDS_MAX_AGE = 500;     // ms without sTestCoords: the camera does not see the emitters
    const BLOBS_MAX_AGE = 250;      // sTestBlobs come right before sTestCoords (fullscreen.js IRTEST_BLOBS_*)
    const BLOBS_VERSION = 1;
    const BLOBS_LENGTH = 21;
    const BLOB_TRACKED = 0x04;      // sTestBlobs flag (all entries): the tracker updated the aim this frame
    const FLASH_TIME = 120;         // crosshair flash on a shot
    const SHOT_TIME = 2000;         // hit marks fade out in this time
    const SHOT_MARKS = 12;
    const BUTTON_TIME = 400;        // circle with the name of a button

    const rgba = (c, a) => `rgba(${c[0]}, ${c[1]}, ${c[2]}, ${a})`;

    class AimTestWindow {
        constructor(options = {}) {
            this.mode = MODE_AIMTEST;
            this.options = options;
            this.closed = false;
            this.showHints = true;
            this.feedback = true;   // F: recoil on the trigger
            this.aim = null;        // { x, y, time, seen, tracked }: last aim point, 0..1 of the screen
            this.blobs = null;      // { seen, tracked, time }: the last sTestBlobs
            this.shots = [];        // { x, y, time }
            this.flash = 0;         // time of the last shot
            this.button = null;     // { name, time }
            this._frame = 0;
            this._timer = 0;
            this._grid = null;      // the grid drawn once per window size (offscreen canvas)
        }

        t(text) { return OF.i18n ? OF.i18n.t(text) : text; }

        open() {
            const canvas = doc.createElement('canvas');
            canvas.className = 'fullscreen-canvas';
            const overlay = doc.createElement('div');
            overlay.className = 'fullscreen-window mode-aimtest';
            overlay.tabIndex = -1;
            overlay.setAttribute('role', 'dialog');
            overlay.setAttribute('aria-label', this.t('Aim Test'));
            const retry = doc.createElement('button');
            retry.className = 'fullscreen-retry';
            retry.type = 'button';
            retry.hidden = true;
            retry.textContent = this.t('Click here to show this window in fullscreen');
            retry.addEventListener('click', () => this._requestFullscreen());
            overlay.append(canvas, retry);
            (doc.getElementById('overlay-root') || doc.body).append(overlay);
            this.overlay = overlay;
            this.canvas = canvas;
            this.retryButton = retry;

            const app = doc.getElementById('app');
            this._appWasInert = app ? app.inert : false;
            if (app) app.inert = true;

            this._onKey = (event) => {
                if (OF.UI && OF.UI.modalOpen && OF.UI.modalOpen()) return; // the message box first
                if (event.ctrlKey || event.altKey || event.metaKey) return;
                if (event.key === 'Escape') this.shutdown();
                else if (event.key === 'h' || event.key === 'H') { this.showHints = !this.showHints; this.render(); }
                // The hints come back to show the new state.
                else if (event.key === 'f' || event.key === 'F') { this.feedback = !this.feedback; this.showHints = true; this.render(); }
                else return;
                event.preventDefault();
                event.stopPropagation();
            };
            this._onResize = () => this.render();
            this._onFullscreen = () => {
                if (!doc.fullscreenElement && this._wasFullscreen && !this.closed) this.shutdown();
                this._wasFullscreen = !!doc.fullscreenElement;
                this.render();
            };
            doc.addEventListener('keydown', this._onKey, true);
            root.addEventListener('resize', this._onResize);
            doc.addEventListener('fullscreenchange', this._onFullscreen);
            overlay.addEventListener('contextmenu', (event) => event.preventDefault());
            this._wasFullscreen = false;
            this._requestFullscreen();
            overlay.focus();
            this.render();
        }

        _requestFullscreen() {
            const overlay = this.overlay;
            if (this.closed || !overlay.requestFullscreen) return;
            this.retryButton.hidden = true;
            overlay.requestFullscreen({ navigationUI: 'hide' })
                .then(() => { this._wasFullscreen = true; overlay.focus(); })
                .catch(() => { if (!this.closed) this.retryButton.hidden = false; });
        }

        shutdown() {
            if (this.closed) return;
            this.closed = true;
            root.cancelAnimationFrame(this._frame);
            root.clearTimeout(this._timer);
            this._grid = null;
            doc.removeEventListener('keydown', this._onKey, true);
            root.removeEventListener('resize', this._onResize);
            doc.removeEventListener('fullscreenchange', this._onFullscreen);
            if (doc.fullscreenElement === this.overlay && doc.exitFullscreen) doc.exitFullscreen().catch(() => {});
            this.overlay.remove();
            const app = doc.getElementById('app');
            if (app) app.inert = this._appWasInert;
            if (this.options.onExit) this.options.onExit(this.mode, null);
        }

        escape() { this.shutdown(); }

        // ----- Board events ------------------------------------------------------------------

        /**
         * sTestCoords: 12 int32, the aim point in fields 8 and 9. The firmware sends them also
         * without emitters (the old point stays). From the sTestBlobs just before: seen = count of
         * the emitters seen, tracked = the tracker updated the point (both null without sTestBlobs).
         */
        drawTest(payload) {
            if (this.closed || payload.length !== 48) return;
            const now = Date.now();
            const view = new DataView(payload.buffer, payload.byteOffset, 48);
            const blobs = this.blobs && now - this.blobs.time <= BLOBS_MAX_AGE ? this.blobs : null;
            this.aim = {
                x: view.getInt32(32, true) / TEST_W, y: view.getInt32(36, true) / TEST_H, time: now,
                seen: blobs ? blobs.seen : null,
                tracked: blobs ? blobs.tracked : null,
            };
            this.render();
        }

        /** sTestBlobs: version, then for each emitter flags (bit 0 seen), average, max, area (uint16). */
        setTestBlobs(payload) {
            if (this.closed || payload.length < BLOBS_LENGTH || payload[0] !== BLOBS_VERSION) return;
            let seen = 0;
            for (let i = 0; i < 4; ++i) seen += payload[1 + i * 5] & 1;
            this.blobs = { seen, tracked: (payload[1] & BLOB_TRACKED) !== 0, time: Date.now() };
        }

        onButton(index, pressed) {
            if (this.closed || !pressed) return;
            const now = Date.now();
            if (index === TRIGGER) {
                this.flash = now;
                if (this._aimFresh(now)) {
                    this.shots.push({ x: this.aim.x, y: this.aim.y, time: now });
                    if (this.shots.length > SHOT_MARKS) this.shots.shift();
                }
                if (this.feedback && this.options.onTrigger) this.options.onTrigger();
            } else {
                this.button = { name: this.options.buttonName ? this.options.buttonName(index) : String(index), time: now };
            }
            this.render();
        }

        /** Aim point recent and updated by the tracker (no sTestBlobs: the coordinates are trusted).
            A shot is marked only then. */
        _aimFresh(now) {
            const aim = this.aim;
            return !!aim && now - aim.time <= COORDS_MAX_AGE && aim.tracked !== false;
        }

        // ----- Drawing ---------------------------------------------------------------------

        render() {
            if (this.closed || this._frame) return;
            this._frame = root.requestAnimationFrame(() => {
                this._frame = 0;
                if (!this.closed) this._draw();
            });
        }

        _draw() {
            const rect = this.overlay.getBoundingClientRect();
            const dpr = root.devicePixelRatio || 1;
            const w = rect.width;
            const h = rect.height;
            const canvas = this.canvas;
            const pw = Math.round(w * dpr);
            const ph = Math.round(h * dpr);
            if (canvas.width !== pw || canvas.height !== ph) {
                canvas.width = pw;
                canvas.height = ph;
            }
            const ctx = canvas.getContext('2d');
            const unit = Math.max(8, h / 90); // the text and the targets follow the screen height
            const now = Date.now();
            const fresh = this._aimFresh(now);

            // The grid never changes: drawn once for this size, then copied on every frame.
            if (!this._grid || this._grid.width !== pw || this._grid.height !== ph) {
                this._grid = doc.createElement('canvas');
                this._grid.width = pw;
                this._grid.height = ph;
                const gridCtx = this._grid.getContext('2d');
                gridCtx.setTransform(dpr, 0, 0, dpr, 0, 0);
                this._drawGrid(gridCtx, w, h, unit);
            }
            ctx.setTransform(1, 0, 0, 1, 0, 0);
            ctx.drawImage(this._grid, 0, 0);
            ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            this.shots = this.shots.filter((shot) => now - shot.time < SHOT_TIME);
            for (const shot of this.shots) this._drawShot(ctx, shot.x * w, shot.y * h, unit, now - shot.time);

            ctx.textAlign = 'center';
            ctx.textBaseline = 'top';
            if (this.showHints) {
                ctx.font = `${unit * 1.15}px ${FONT}`;
                ctx.fillStyle = COLOR_HINTS;
                ctx.fillText(`Esc ${this.t('quit')}    H ${this.t('hide this line')}    F ${this.t(this.feedback ? 'feedback on' : 'feedback off')}`,
                    w / 2, unit * 0.6);
            }
            if (!fresh) {
                // Emitters seen but the tracker not updated (after a dark moment): the aim is not usable yet.
                const seen = this.aim && now - this.aim.time <= COORDS_MAX_AGE && this.aim.seen > 0;
                ctx.font = `bold ${unit * 1.5}px ${FONT}`;
                ctx.fillStyle = COLOR_MESSAGE;
                ctx.fillText(this.t(seen ? 'Aim unavailable' : 'IR emitters not seen'), w / 2, h * 0.40);
            }
            if (this.aim) this._drawAim(ctx, w, h, unit, now, fresh);

            // Redraw every 40 ms only while something fades (hit marks, flash, button circle);
            // otherwise once, when the aim gets old. New coordinates redraw anyway (drawTest).
            const fading = this.shots.map((shot) => shot.time + SHOT_TIME);
            if (this.button) fading.push(this.button.time + BUTTON_TIME);
            fading.push(this.flash + FLASH_TIME);
            let delay = fading.some((end) => end > now) ? 40 : Infinity;
            if (fresh) delay = Math.min(delay, this.aim.time + COORDS_MAX_AGE + 1 - now);
            root.clearTimeout(this._timer);
            if (delay !== Infinity) this._timer = root.setTimeout(() => this.render(), delay);
        }

        /** Background, grid, border and the nine targets. */
        _drawGrid(ctx, w, h, unit) {
            ctx.fillStyle = COLOR_BACKGROUND;
            ctx.fillRect(0, 0, w, h);
            const line = (x0, y0, x1, y1) => {
                ctx.beginPath();
                ctx.moveTo(x0, y0);
                ctx.lineTo(x1, y1);
                ctx.stroke();
            };
            for (let i = 1; i < GRID_LINES; ++i) {
                const major = i % 2 === 0;
                ctx.strokeStyle = major ? COLOR_GRID_MAJOR : COLOR_GRID_MINOR;
                ctx.lineWidth = major ? 1.5 : 1;
                line(w * i / GRID_LINES, 0, w * i / GRID_LINES, h);
                line(0, h * i / GRID_LINES, w, h * i / GRID_LINES);
            }
            ctx.strokeStyle = COLOR_BORDER;
            ctx.lineWidth = 2;
            ctx.strokeRect(1, 1, w - 3, h - 3);

            const r = unit * 2.2;
            ctx.strokeStyle = COLOR_TARGET;
            for (const fx of [TARGET_INSET, 0.5, 1 - TARGET_INSET]) {
                for (const fy of [TARGET_INSET, 0.5, 1 - TARGET_INSET]) {
                    const cx = w * fx;
                    const cy = h * fy;
                    ctx.beginPath();
                    ctx.arc(cx, cy, r, 0, Math.PI * 2);
                    ctx.moveTo(cx + r / 3, cy);
                    ctx.arc(cx, cy, r / 3, 0, Math.PI * 2);
                    ctx.stroke();
                    line(cx - r * 1.5, cy, cx + r * 1.5, cy);
                    line(cx, cy - r * 1.5, cx, cy + r * 1.5);
                }
            }
        }

        /** A hit mark: a ring that opens on the shot, then a dot that fades out. */
        _drawShot(ctx, x, y, unit, age) {
            const k = age / SHOT_TIME;
            if (age < 250) {
                ctx.beginPath();
                ctx.arc(x, y, unit * (1 + 4 * age / 250), 0, Math.PI * 2);
                ctx.strokeStyle = `rgba(255, 255, 255, ${0.8 * (1 - age / 250)})`;
                ctx.lineWidth = 2;
                ctx.stroke();
            }
            ctx.beginPath();
            ctx.arc(x, y, unit * 0.45, 0, Math.PI * 2);
            ctx.fillStyle = rgba(COLOR_AIM, 0.9 * (1 - k));
            ctx.fill();
        }

        /** Crosshair (white flash on a shot), the name of the last button, OFF-SCREEN at the edge. */
        _drawAim(ctx, w, h, unit, now, fresh) {
            const outside = this.aim.x < 0 || this.aim.x > 1 || this.aim.y < 0 || this.aim.y > 1;
            const cx = Math.min(Math.max(this.aim.x, 0), 1) * (w - 1);
            const cy = Math.min(Math.max(this.aim.y, 0), 1) * (h - 1);
            const r = unit * 2.8;

            if (this.button && now - this.button.time < BUTTON_TIME) {
                const fade = 1 - (now - this.button.time) / BUTTON_TIME;
                ctx.beginPath();
                ctx.arc(cx, cy, r * 2.6, 0, Math.PI * 2);
                ctx.fillStyle = `rgba(80, 150, 255, ${0.35 * fade})`;
                ctx.fill();
                ctx.font = `bold ${unit * 1.4}px ${FONT}`;
                ctx.textAlign = 'center';
                ctx.textBaseline = 'top';
                ctx.fillStyle = `rgba(255, 255, 255, ${fade})`;
                ctx.fillText(this.button.name, cx, Math.min(cy + r * 2.0, h - unit * 2));
            }

            const flash = now - this.flash < FLASH_TIME;
            ctx.globalAlpha = fresh ? 1 : 0.35; // aim lost: the last point stays, dimmed
            // Shape of the IR camera tester crosshair (fullscreen.js _drawIRTestAim: ring 24.42,
            // ticks 18.83-30, centre cross 3.8) scaled to the same footprint as before (1.8 r).
            const k = r * 1.8 / 30;
            const ring = 24.42 * k;
            const tick = [18.83 * k, 30 * k];
            const dot = 3.8 * k;
            const partial = this.aim.seen !== null && this.aim.seen < 4;
            const color = flash ? 'white' : rgba(partial ? COLOR_AIM_PARTIAL : COLOR_AIM, 1);
            for (const [stroke, ringWidth, lineWidth] of [['rgba(0, 0, 0, 0.63)', 6, 5], [color, 3, 2]]) {
                ctx.strokeStyle = stroke;
                ctx.lineWidth = ringWidth;
                ctx.beginPath();
                ctx.arc(cx, cy, ring, 0, Math.PI * 2);
                ctx.stroke();
                ctx.lineWidth = lineWidth;
                ctx.beginPath();
                ctx.moveTo(cx - tick[1], cy); ctx.lineTo(cx - tick[0], cy);
                ctx.moveTo(cx + tick[0], cy); ctx.lineTo(cx + tick[1], cy);
                ctx.moveTo(cx, cy - tick[1]); ctx.lineTo(cx, cy - tick[0]);
                ctx.moveTo(cx, cy + tick[0]); ctx.lineTo(cx, cy + tick[1]);
                ctx.moveTo(cx - dot, cy); ctx.lineTo(cx + dot, cy);
                ctx.moveTo(cx, cy - dot); ctx.lineTo(cx, cy + dot);
                ctx.stroke();
            }
            if (flash) {
                ctx.beginPath();
                ctx.arc(cx, cy, r * 0.6, 0, Math.PI * 2);
                ctx.fillStyle = rgba(COLOR_AIM, 0.6);
                ctx.fill();
            }
            if (outside && fresh) {
                ctx.font = `bold ${unit * 1.3}px ${FONT}`;
                ctx.textAlign = 'center';
                ctx.textBaseline = 'top';
                const text = this.t('OFF-SCREEN');
                const half = ctx.measureText(text).width / 2;
                const tx = Math.min(Math.max(cx, half + 4), w - half - 4);
                const ty = Math.min(Math.max(cy - r * 2.2 - unit * 1.6, 4), h - unit * 1.6 - 4);
                ctx.fillStyle = COLOR_EDGE;
                ctx.fillText(text, tx, ty);
            }
            ctx.globalAlpha = 1;
        }
    }

    AimTestWindow.MODE_AIMTEST = MODE_AIMTEST;

    OF.AimTestWindow = AimTestWindow;
})(typeof globalThis !== 'undefined' ? globalThis : this);
