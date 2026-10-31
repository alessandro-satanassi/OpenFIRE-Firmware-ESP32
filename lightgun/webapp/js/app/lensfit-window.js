/*  OpenFIRE Web App - fullscreen window of the lens refinement ("Refine lens", Gun Settings).

    The user sweeps the gun over the whole screen. The screen is a grid that gets "painted"
    where the gun points (grey = empty, yellow = some views, green = complete); an inset shows
    the camera view, painted where the LEDs landed: its edges matter most for the lens. When
    everything is green the lens values are fitted (OF.LensFit) and shown with Apply / Retry /
    Cancel. Nothing is sent to the gun here: main.js puts the accepted values in the settings,
    saved with the normal Save.

        const win = new OF.LensFitWindow({ onExit: (mode, values) => ... }); // values null = cancelled
        win.open();                       // inside the click handler (fullscreen needs a user gesture)
        win.setLensInfo(payload);         // sLensInfo
        win.addLensSample(payload);       // sLensSample
        win.shutdown();
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};
    const doc = root.document;

    const MODE_LENSFIT = 'lensfit';
    const INFO_TIMEOUT = 3000;      // firmware without the lens refinement never sends sLensInfo
    const LOST_TIME = 400;          // ms without a complete view before the red warning

    const COLOR_EMPTY = [58, 62, 70];
    const COLOR_SOME = [230, 190, 40];
    const COLOR_DONE = [40, 170, 80];

    function cellColor(fraction) {
        if (fraction <= 0) return `rgb(${COLOR_EMPTY.join(',')})`;
        const t = Math.min(1, fraction);
        const from = t < 1 ? COLOR_SOME : COLOR_DONE;
        const k = t < 1 ? 0.35 + 0.65 * t : 1;
        const c = from.map((v, i) => Math.round(COLOR_EMPTY[i] + (v - COLOR_EMPTY[i]) * k));
        return `rgb(${c.join(',')})`;
    }

    function format(value, digits) {
        return value === null || value === undefined || !Number.isFinite(value) ? '-' : value.toFixed(digits);
    }

    class LensFitWindow {
        constructor(options = {}) {
            this.mode = MODE_LENSFIT;
            this.options = options;
            this.info = null;
            this.collector = null;
            this.state = 'waiting';     // waiting, collecting, solving, result, unsupported
            this.result = null;
            this.closed = false;
            this.lastComplete = 0;
            this.lastSample = 0;
            this._frame = 0;
        }

        t(text) { return OF.i18n ? OF.i18n.t(text) : text; }

        open() {
            const canvas = doc.createElement('canvas');
            canvas.className = 'fullscreen-canvas';
            const overlay = doc.createElement('div');
            overlay.className = 'fullscreen-window mode-lensfit';
            overlay.tabIndex = -1;
            overlay.setAttribute('role', 'dialog');
            overlay.setAttribute('aria-label', this.t('Refine camera lens'));
            const retry = doc.createElement('button');
            retry.className = 'fullscreen-retry';
            retry.type = 'button';
            retry.hidden = true;
            retry.textContent = this.t('Click here to show this window in fullscreen');
            retry.addEventListener('click', () => this._requestFullscreen());
            const panel = doc.createElement('div');
            panel.className = 'lensfit-panel';
            panel.hidden = true;
            overlay.append(canvas, panel, retry);
            (doc.getElementById('overlay-root') || doc.body).append(overlay);
            this.overlay = overlay;
            this.canvas = canvas;
            this.panel = panel;
            this.retryButton = retry;

            const app = doc.getElementById('app');
            this._appWasInert = app ? app.inert : false;
            if (app) app.inert = true;

            this._onKey = (event) => {
                if (event.key === 'Escape') {
                    if (OF.UI && OF.UI.modalOpen && OF.UI.modalOpen()) return;
                    event.preventDefault();
                    event.stopPropagation();
                    this.shutdown(null);
                }
            };
            this._onResize = () => this.render();
            this._onFullscreen = () => {
                if (!doc.fullscreenElement && this._wasFullscreen && !this.closed) this.shutdown(null);
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
            this._infoTimer = root.setTimeout(() => {
                if (!this.info && !this.closed) { this.state = 'unsupported'; this._showPanel(); this.render(); }
            }, INFO_TIMEOUT);
            this._tick = root.setInterval(() => this.render(), 250);
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

        /** values: { k1, k2, focal } (null = camera record value) when applied, null when cancelled. */
        shutdown(values) {
            if (this.closed) return;
            this.closed = true;
            root.clearTimeout(this._infoTimer);
            root.clearTimeout(this._solveTimer);
            root.clearInterval(this._tick);
            root.cancelAnimationFrame(this._frame);
            doc.removeEventListener('keydown', this._onKey, true);
            root.removeEventListener('resize', this._onResize);
            doc.removeEventListener('fullscreenchange', this._onFullscreen);
            if (doc.fullscreenElement === this.overlay && doc.exitFullscreen) doc.exitFullscreen().catch(() => {});
            this.overlay.remove();
            const app = doc.getElementById('app');
            if (app) app.inert = this._appWasInert;
            if (this.options.onExit) this.options.onExit(this.mode, values || null);
        }

        escape() { this.shutdown(null); }

        // The IR test events are not used here (firmware without the lens flag could send them).
        drawTest() {}
        setTestBlobs() {}

        setLensInfo(payload) {
            if (this.closed || this.info) return;
            const info = OF.LensFit.parseInfo(payload);
            if (!info) return;
            this.info = info;
            this.collector = new OF.LensFit.Collector(info);
            this.state = info.record.focal > 0 ? 'collecting' : 'unsupported';
            if (this.state === 'unsupported') this._showPanel();
            this.render();
        }

        addLensSample(payload) {
            if (this.closed || this.state !== 'collecting') return;
            const sample = OF.LensFit.parseSample(payload);
            if (!sample) return;
            const now = Date.now();
            this.lastSample = now;
            if (sample.complete) this.lastComplete = now;
            this.collector.add(sample);
            if (this.collector.progress().done) this._solve();
            this._schedule();
        }

        _schedule() {
            if (this._frame || this.closed) return;
            this._frame = root.requestAnimationFrame(() => { this._frame = 0; this.render(); });
        }

        _solve() {
            this.state = 'solving';
            this.solveProgress = 0;
            this.render();
            // The fit gives the page back between its iterations: the progress is drawn and ESC
            // closes the window (the fit then stops).
            const frames = this.collector.frames.slice();
            this._solveTimer = root.setTimeout(async () => {
                if (this.closed) return;
                let result;
                try {
                    result = await OF.LensFit.solveAsync(this.info, frames, {
                        onProgress: (p) => { this.solveProgress = p; this._schedule(); },
                        isCancelled: () => this.closed || this.state !== 'solving',
                    });
                } catch (error) { result = { ok: false, reason: String(error) }; }
                if (this.closed || this.state !== 'solving') return;
                this.result = result;
                // Last run kept for diagnosis (also after the window is closed), from the browser
                // console: copy(JSON.stringify(OF.LensFit.lastRun))
                OF.LensFit.lastRun = { info: this.info, result, frames };
                this.state = 'result';
                this._showPanel();
                this.render();
            }, 60);
        }

        _restart() {
            this.collector = new OF.LensFit.Collector(this.info);
            this.result = null;
            this.state = 'collecting';
            this.panel.hidden = true;
            this.panel.replaceChildren();
            this.overlay.focus();
            this.render();
        }

        // ----- Result panel ---------------------------------------------------------------

        _showPanel() {
            const t = (s) => this.t(s);
            const el = (tag, cls, text) => {
                const node = doc.createElement(tag);
                if (cls) node.className = cls;
                if (text !== undefined) node.textContent = text;
                return node;
            };
            const button = (text, onClick, primary) => {
                const b = el('button', primary ? 'lensfit-button primary' : 'lensfit-button', text);
                b.type = 'button';
                b.addEventListener('click', onClick);
                return b;
            };
            const panel = this.panel;
            panel.replaceChildren();
            panel.append(el('h2', '', t('Refine camera lens')));
            const buttons = el('div', 'lensfit-buttons');

            if (this.state === 'unsupported') {
                panel.append(el('p', '', this.info
                    ? t('The camera record has no focal ratio: the lens cannot be refined.')
                    : t('The gun did not answer: its firmware does not support the lens refinement. Update the firmware.')));
                buttons.append(button(t('Close'), () => this.shutdown(null), true));
            } else {
                // One screen for every outcome: the values in use next to the calculated ones
                // (ERROR in their place when the calculation failed).
                const r = this.result;
                const ok = !!(r && r.ok);
                // Confirmed: a refined tier, or the default values when they beat the values in use.
                const better = ok && !r.keepActive && (r.tier !== 'record' || r.activeCustom);
                // What can be applied: the confirmed result, a refined tier that does not beat the
                // values in use, or else the closest calculated values not confirmed (near). None
                // when nothing could be determined: the default values are never shown as calculated.
                const proposal = !ok ? null : better || r.tier !== 'record' ? r : r.near || null;
                const rec = this.info ? this.info.record : null, active = this.info ? this.info.active : null;
                const value = (v, fallback) => v === null ? fallback : v;
                const next = proposal && rec
                    ? { k1: value(proposal.k1, rec.k1), k2: value(proposal.k2, rec.k2), focal: value(proposal.focal, rec.focal) } : null;
                const rows = [[t('Lens radial correction K1'), 'k1', 4], [t('Lens radial correction K2'), 'k2', 4],
                    [t('Lens focal ratio'), 'focal', 4]];
                const changed = (key, digits) => !!next && format(next[key], digits) !== format(active[key], digits);
                const anyChange = rows.some(([, key, digits]) => changed(key, digits));
                panel.append(el('p', '', !ok
                    ? t('The lens values could not be computed. Try again, moving the gun more slowly over the whole screen.')
                    : !proposal ? t('The analysis found no reliable change: keep the current values.')
                        : !anyChange ? t('The current values are already right for your camera: nothing needs to change.')
                            : better ? t('Values better suited to your camera were found: the aim should be more precise, above all towards the edges of the screen.')
                                : t('The calculated values do not clearly improve on the current ones: better to keep the current ones, so the profiles need no new calibration.')));
                if (active) {
                    const table = el('table', 'lensfit-table');
                    const head = el('tr');
                    for (const h of ['', t('Current'), t('Calculated')]) head.append(el('th', '', h));
                    table.append(head);
                    for (const [label, key, digits] of rows) {
                        const tr = el('tr');
                        tr.append(el('td', '', label), el('td', '', format(active[key], digits)),
                            !ok ? el('td', 'error', t('ERROR'))
                                : next ? el('td', changed(key, digits) ? 'changed' : '', format(next[key], digits)) : el('td', '', '\u2014'));
                        table.append(tr);
                    }
                    panel.append(table);
                }
                const apply = () => this.shutdown({ k1: proposal.k1, k2: proposal.k2, focal: proposal.focal });
                if (anyChange && better) {
                    panel.append(el('p', 'note', t('After saving, calibrate again every profile you use.')));
                    buttons.append(button(t('Apply'), apply, true),
                        button(t('Retry'), () => this._restart()), button(t('Cancel'), () => this.shutdown(null)));
                } else if (anyChange) {
                    // Not confirmed: keeping the current ones is the suggested choice (the sentence
                    // above says that applying them means calibrating again, and so does the status
                    // line after Apply).
                    buttons.append(button(t('Keep current'), () => this.shutdown(null), true), button(t('Apply'), apply),
                        button(t('Retry'), () => this._restart()));
                } else if (ok) {
                    buttons.append(button(t('Close'), () => this.shutdown(null), true), button(t('Retry'), () => this._restart()));
                } else {
                    buttons.append(button(t('Retry'), () => this._restart(), true), button(t('Cancel'), () => this.shutdown(null)));
                }
            }
            panel.append(buttons);
            panel.hidden = false;
            const first = buttons.querySelector('button');
            if (first) first.focus();
        }

        // ----- Drawing --------------------------------------------------------------------

        render() {
            if (this.closed || !this.canvas) return;
            const canvas = this.canvas;
            const dpr = root.devicePixelRatio || 1;
            const width = this.overlay.clientWidth || root.innerWidth, height = this.overlay.clientHeight || root.innerHeight;
            if (canvas.width !== Math.round(width * dpr) || canvas.height !== Math.round(height * dpr)) {
                canvas.width = Math.round(width * dpr);
                canvas.height = Math.round(height * dpr);
            }
            const ctx = canvas.getContext('2d');
            ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            ctx.fillStyle = '#15181d';
            ctx.fillRect(0, 0, width, height);
            const L = OF.LensFit;
            const c = this.collector;

            // Screen grid: painted in step 1; afterwards it stays green, dimmed, behind the guide.
            const g = c && this.state === 'collecting' ? c.guide() : null;
            const step = g ? g.step : 0;
            const gap = Math.max(2, Math.round(Math.min(width, height) * 0.004));
            const cw = width / L.SCREEN_COLS, ch = height / L.SCREEN_ROWS;
            ctx.globalAlpha = step >= 2 ? 0.35 : 1;
            for (let row = 0; row < L.SCREEN_ROWS; ++row) for (let col = 0; col < L.SCREEN_COLS; ++col) {
                const n = c ? c.screen[row * L.SCREEN_COLS + col] : 0;
                ctx.fillStyle = cellColor(n / L.SCREEN_TARGET);
                ctx.fillRect(col * cw + gap, row * ch + gap, cw - 2 * gap, ch - 2 * gap);
            }
            ctx.globalAlpha = 1;

            const now = Date.now();
            const lost = this.state === 'collecting' && now - this.lastComplete > LOST_TIME;
            const font = Math.max(14, Math.round(Math.min(width, height) * 0.024));
            const big = Math.round(font * 1.6);
            const pulse = 0.5 + 0.5 * Math.sin(now / 180);
            const aim = c && c.aim && !lost ? [c.aim[0] * width, c.aim[1] * height] : null;

            // Step 2: the yellow circle to aim at (always inside the screen), joined to the gun's aim by an arrow.
            if (step === 2 && g.target && !lost) {
                const r = Math.min(width, height) * 0.06;
                const tx = Math.min(width - r * 1.2, Math.max(r * 1.2, g.target[0] * width));
                const ty = Math.min(height - r * 1.2, Math.max(r * 1.2, g.target[1] * height));
                if (aim && Math.hypot(tx - aim[0], ty - aim[1]) > r * 1.5) {
                    const ang = Math.atan2(ty - aim[1], tx - aim[0]);
                    const ex = tx - Math.cos(ang) * r * 1.2, ey = ty - Math.sin(ang) * r * 1.2;
                    ctx.strokeStyle = 'rgba(255,215,90,0.9)'; ctx.fillStyle = 'rgba(255,215,90,0.9)';
                    ctx.lineWidth = Math.max(4, font * 0.3);
                    ctx.beginPath(); ctx.moveTo(aim[0], aim[1]); ctx.lineTo(ex, ey); ctx.stroke();
                    const h = font * 1.2;
                    ctx.beginPath(); ctx.moveTo(ex + Math.cos(ang) * h * 0.4, ey + Math.sin(ang) * h * 0.4);
                    ctx.lineTo(ex - Math.cos(ang - 0.5) * h, ey - Math.sin(ang - 0.5) * h);
                    ctx.lineTo(ex - Math.cos(ang + 0.5) * h, ey - Math.sin(ang + 0.5) * h); ctx.closePath(); ctx.fill();
                }
                ctx.fillStyle = `rgba(255,215,90,${0.25 + 0.25 * pulse})`;
                ctx.beginPath(); ctx.arc(tx, ty, r * (1 + 0.15 * pulse), 0, 2 * Math.PI); ctx.fill();
                ctx.strokeStyle = '#ffd75a'; ctx.lineWidth = Math.max(4, font * 0.3);
                ctx.beginPath(); ctx.arc(tx, ty, r, 0, 2 * Math.PI); ctx.stroke();
                ctx.beginPath(); ctx.arc(tx, ty, r * 0.15, 0, 2 * Math.PI); ctx.fillStyle = '#ffd75a'; ctx.fill();
            }

            // Step 2 without a reachable circle: closer or farther, with a big animated frame.
            if (step === 2 && !g.target && !lost) {
                const cx = width / 2, cy = height / 2, base = Math.min(width, height) * 0.16;
                const grow = g.closer ? pulse : 1 - pulse;
                ctx.strokeStyle = '#ffd75a'; ctx.lineWidth = Math.max(4, font * 0.3);
                for (let k = 0; k < 3; ++k) {
                    const s = base * (0.6 + 0.5 * ((k / 3 + grow / 3) % 1) * 2);
                    ctx.globalAlpha = 0.35 + 0.65 * (k + 1) / 3;
                    ctx.strokeRect(cx - s * 16 / 9, cy - s, 2 * s * 16 / 9, 2 * s);
                }
                ctx.globalAlpha = 1;
            }

            // Step 3: a dial; the needle follows the tilt of the gun, each side turns green.
            if (step === 3 && !lost) {
                const cx = width / 2, cy = height * 0.62, R = Math.min(width, height) * 0.26;
                const span = 40 * Math.PI / 180, lim = g.limit;
                const arc = (from, to, color, w) => {
                    ctx.strokeStyle = color; ctx.lineWidth = w;
                    ctx.beginPath(); ctx.arc(cx, cy, R, -Math.PI / 2 + from, -Math.PI / 2 + to); ctx.stroke();
                };
                const zone = (fill) => cellColor(fill);
                arc(-span, span, '#3a3e46', R * 0.22);
                arc(-span, -lim, zone(g.left), R * 0.22);
                arc(lim, span, zone(g.right), R * 0.22);
                const a = Math.max(-span, Math.min(span, g.roll));
                ctx.strokeStyle = '#ffffff'; ctx.lineWidth = Math.max(5, font * 0.35);
                ctx.beginPath(); ctx.moveTo(cx, cy);
                ctx.lineTo(cx + Math.sin(a) * R * 1.12, cy - Math.cos(a) * R * 1.12); ctx.stroke();
                ctx.fillStyle = '#ffffff'; ctx.beginPath(); ctx.arc(cx, cy, font * 0.6, 0, 2 * Math.PI); ctx.fill();
                ctx.font = `bold ${big}px "Segoe UI", Arial, sans-serif`;
                ctx.textBaseline = 'middle'; ctx.textAlign = 'center';
                ctx.fillStyle = '#ffffff';
                ctx.fillText(g.left >= 1 ? '✔' : '↶', cx - Math.sin(span) * R * 1.35, cy - Math.cos(span) * R * 1.35);
                ctx.fillText(g.right >= 1 ? '✔' : '↷', cx + Math.sin(span) * R * 1.35, cy - Math.cos(span) * R * 1.35);
                ctx.textAlign = 'start';
            }

            // Aim of the gun (uncalibrated, only to show where it points; in step 3 too, the aim has to
            // keep moving while the gun is tilted).
            if (aim && this.state === 'collecting') {
                const [x, y] = aim;
                ctx.strokeStyle = '#ffffff';
                ctx.lineWidth = 3;
                ctx.beginPath(); ctx.arc(x, y, 18, 0, 2 * Math.PI); ctx.stroke();
                ctx.beginPath(); ctx.moveTo(x - 30, y); ctx.lineTo(x + 30, y); ctx.moveTo(x, y - 30); ctx.lineTo(x, y + 30); ctx.stroke();
            }

            // Text panel: step dots, one short instruction, the overall progress.
            ctx.font = `${font}px "Segoe UI", Arial, sans-serif`;
            ctx.textBaseline = 'top';
            const boxW = Math.min(width * 0.7, font * 36), pad = font * 0.7;
            const lines = [];
            const p = c ? c.progress() : { total: 0 };
            if (this.state === 'waiting') lines.push([this.t('Waiting for the gun...'), '#ffffff', font]);
            if (step === 1) {
                lines.push([this.t('Step 1 of 3: paint the screen'), '#ffffff', big]);
                lines.push([this.t('Move the gun slowly over the whole screen until every square is green.'), '#dddddd', font]);
            } else if (step === 2) {
                lines.push([this.t('Step 2 of 3: edges of the camera view'), '#ffffff', big]);
                lines.push([g.target ? this.t('Bring the white crosshair into the yellow circle.')
                    : g.closer ? this.t('Come closer to the screen until the yellow circle appears.')
                        : this.t('Step back from the screen until the yellow circle appears.'), '#ffd75a', font]);
            } else if (step === 3) {
                lines.push([this.t('Step 3 of 3: rotate the gun'), '#ffffff', big]);
                lines.push([this.t('Rotate the gun on its own axis, like turning a key, to one side; holding it there, slowly move the white crosshair over the screen until that side turns green, then do the same on the other side.'), '#ffd75a', font]);
            } else if (step === 4 || this.state === 'solving') {
                lines.push([this.t('Done!'), '#ffffff', big]);
                lines.push([this.t('Computing the lens values... %1% (Esc cancels)').replace('%1', String(Math.floor((this.solveProgress || 0) * 100))), '#dddddd', font]);
            }
            if (lines.length && this.state !== 'result' && this.state !== 'unsupported') {
                const wrapped = [];
                for (const [text, color, size] of lines) {
                    ctx.font = `${size === big ? 'bold ' : ''}${size}px "Segoe UI", Arial, sans-serif`;
                    const maxW = boxW - 2 * pad - (size === big && step ? big * 3.4 : 0);   // the title leaves room for the step dots
                    let line = '';
                    for (const word of text.split(' ')) {
                        const test = line ? line + ' ' + word : word;
                        if (ctx.measureText(test).width > maxW && line) { wrapped.push([line, color, size]); line = word; } else line = test;
                    }
                    if (line) wrapped.push([line, color, size]);
                }
                const bar = this.state === 'collecting' || this.state === 'solving';
                const textH = wrapped.reduce((s, w) => s + w[2] * 1.3, 0);
                const boxH = textH + 2 * pad + (bar ? font * 1.2 : 0);
                const bx = (width - boxW) / 2, by = font;
                ctx.fillStyle = 'rgba(0,0,0,0.78)';
                ctx.fillRect(bx, by, boxW, boxH);
                let ty = by + pad;
                for (const [text, color, size] of wrapped) {
                    ctx.font = `${size === big ? 'bold ' : ''}${size}px "Segoe UI", Arial, sans-serif`;
                    ctx.fillStyle = color; ctx.fillText(text, bx + pad, ty); ty += size * 1.3;
                }
                if (bar) {
                    const barY = ty + font * 0.2, barW = boxW - 2 * pad;
                    const share = this.state === 'solving' ? (this.solveProgress || 0) : p.total;
                    ctx.fillStyle = '#3a3e46'; ctx.fillRect(bx + pad, barY, barW, font * 0.5);
                    ctx.fillStyle = '#28aa50'; ctx.fillRect(bx + pad, barY, barW * share, font * 0.5);
                }
                // Step dots on the right of the title: done = green, current = yellow.
                if (step) for (let k = 1; k <= 3; ++k) {
                    const dx = bx + boxW - pad - (4 - k) * big * 1.1 + big * 0.4, dy = by + pad + big * 0.55;
                    ctx.fillStyle = k < step ? '#28aa50' : k === step ? '#ffd75a' : '#3a3e46';
                    ctx.beginPath(); ctx.arc(dx, dy, big * 0.4, 0, 2 * Math.PI); ctx.fill();
                    ctx.fillStyle = k <= step ? '#15181d' : '#aaaaaa';
                    ctx.font = `bold ${Math.round(big * 0.5)}px "Segoe UI", Arial, sans-serif`;
                    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
                    ctx.fillText(k < step ? '✔' : String(k), dx, dy + 1);
                    ctx.textAlign = 'start'; ctx.textBaseline = 'top';
                }
                ctx.font = `${Math.round(font * 0.8)}px "Segoe UI", Arial, sans-serif`;
                ctx.fillStyle = '#aaaaaa';
                if (this.state === 'collecting') ctx.fillText(this.t('Esc cancels.'), bx + pad, by + boxH + font * 0.3);
            }

            // Red warning: no valid view for a while. The gun cannot say why (an LED out of view, too
            // much tilt, a view not measured exactly): the text names the usual causes.
            if (lost) {
                const text = step <= 1 ? this.t('Invalid views: step back a little and keep all four LEDs in view.')
                    : this.t('Invalid views: keep all four LEDs in view and rotate the gun less.');
                ctx.font = `bold ${font}px "Segoe UI", Arial, sans-serif`;
                const w = Math.min(width * 0.95, ctx.measureText(text).width + font * 2);
                // Wrapped on narrow screens (a phone), not squeezed.
                const rows = [];
                let row = '';
                for (const word of text.split(' ')) {
                    const test = row ? row + ' ' + word : word;
                    if (row && ctx.measureText(test).width > w - 2 * font) { rows.push(row); row = word; } else row = test;
                }
                if (row) rows.push(row);
                const h = font * (0.9 + 1.3 * rows.length);
                ctx.fillStyle = 'rgba(200,30,30,0.9)';
                ctx.fillRect((width - w) / 2, height / 2 - font, w, h);
                ctx.fillStyle = '#ffffff';
                ctx.textBaseline = 'top';
                rows.forEach((r, i) => ctx.fillText(r, (width - w) / 2 + font, height / 2 - font * 0.45 + i * font * 1.3, w - 2 * font));
            }

            // Inset: camera view, its rings painted where the LEDs landed (bigger in step 2, in the
            // bottom corner away from the yellow circle).
            if (c && this.info && this.state === 'collecting') {
                const size = Math.min(width, height) * (step === 2 ? 0.3 : 0.2);
                const left = step === 2 && g.target && g.target[0] > 0.5;
                const ix = left ? font : width - size - font, iy = height - size - font * 1.2;
                ctx.fillStyle = 'rgba(0,0,0,0.75)';
                ctx.fillRect(ix - 6, iy - font * 1.6, size + 12, size + font * 1.6 + 6);
                ctx.font = `${Math.round(font * 0.8)}px "Segoe UI", Arial, sans-serif`;
                ctx.fillStyle = '#ffffff'; ctx.textBaseline = 'top';
                ctx.fillText(this.t('Camera view'), ix, iy - font * 1.35);
                const rings = L.FIELD_RINGS, ox = ix + size / 2, oy = iy + size / 2, rr = size / 2 / 1.06;
                ctx.save();
                ctx.beginPath(); ctx.rect(ix, iy, size, size); ctx.clip();
                for (let k = rings.length - 1; k >= 0; --k) {
                    const outer = (k + 1 < rings.length ? rings[k + 1] : 1.5) * rr;
                    ctx.fillStyle = cellColor(c.field[k] / L.FIELD_TARGET);
                    ctx.beginPath(); ctx.arc(ox, oy, outer, 0, 2 * Math.PI); ctx.fill();
                    ctx.strokeStyle = '#15181d'; ctx.lineWidth = 2; ctx.stroke();
                }
                ctx.restore();
                if (c.lastPoints) {
                    ctx.fillStyle = c.lastComplete ? '#ffffff' : '#ff4040';
                    for (const pt of c.lastPoints) {
                        const px = ox + (pt[0] / this.info.resX - 0.5) * 2 * rr, py = oy + (pt[1] / this.info.resY - 0.5) * 2 * rr;
                        if (px < ix || py < iy || px > ix + size || py > iy + size) continue;
                        ctx.beginPath(); ctx.arc(px, py, 5, 0, 2 * Math.PI); ctx.fill();
                    }
                }
            }
        }
    }

    LensFitWindow.MODE_LENSFIT = MODE_LENSFIT;
    OF.LensFitWindow = LensFitWindow;
})(typeof globalThis !== 'undefined' ? globalThis : this);
