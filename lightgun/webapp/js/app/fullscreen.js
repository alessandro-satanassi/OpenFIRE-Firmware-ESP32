/*  OpenFIRE Web App - fullscreen windows: calibration, IR emitter alignment, IR camera test
    (Qt App: appcali.cpp). Drawn on a canvas with the App's 8x8 bitmap typeface.

        const win = new OF.FullscreenWindow('calibrate', { onExitRequest, onExit });
        win.open();                      // call it inside the click handler (fullscreen needs a user gesture)
        win.setStage(stage); win.setInfo(payload); win.setTestBlobs(blobs); win.drawTest(coords);
        win.showIrWarning(bits);         // calibration: the board refused a target shot (sCaliIrWarning)
        win.shutdown();
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};
    const doc = root.document;

    const MODE_CALIBRATE = 'calibrate';
    const MODE_ALIGNMENT = 'alignment';
    const MODE_IRTEST = 'irtest';

    const STAGE = { init: 0, top: 1, bottom: 2, left: 3, right: 4, center: 5, verify: 6, end: 7 };

    const CALI_PREFIXES = ['Top Offset: ', 'Bottom Offset: ', 'Left Offset: ', 'Right Offset: ', 'Top Left LED: ', 'Top Right LED: '];
    const CALI_KEYS = ['topOffset', 'bottomOffset', 'leftOffset', 'rightOffset', 'TLled', 'TRled'];
    const CALI_SHOT_TIME = 120; // ms: confirms a target acquired, not a save to flash
    const CALI_SHOT_HOLD = 40;  // ms fully white, followed by a quick fade

    const CROSSHAIR_SIZE = 61.44;
    const CROSSHAIR_SVG = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 61.44 61.44" width="61.44" height="61.44">' +
        '<g fill="none" stroke="#ff8758" stroke-linecap="square">' +
        '<circle cx="30.72" cy="30.72" r="24.42" stroke-width="2.4"/>' +
        '<path stroke-width="1.44" d="M0.72 30.72h11.16M30.72 0.73v11.16M60.71 30.72h-11.16M30.72 60.71v-11.16M26.94 30.72h7.58M30.72 26.94v7.58"/>' +
        '</g></svg>';

    // ----- IR test blobs (sTestBlobs) ---------------------------------------------------
    // Firmware 'sTestBlobs': area and brightness of the blob seen at each emitter. The camera
    // drivers already normalise the values, so every camera uses the same scale here.
    // Radius (1920x1080 test space) grows with sqrt(area): IRTEST_BLOB_RADIUS_MAX at
    // IRTEST_BLOB_AREA_MAX, the largest blob the PAJ7025 DSP reports. The firmware DFRobot
    // driver relies on these two values for its equivalent area when it reads the camera's
    // Extended format (ReadDFRobotExtended); with the Full format it reports the box area.
    const IRTEST_BLOB_AREA_MAX = 300;
    const IRTEST_BLOB_RADIUS_MAX = 60;
    const IRTEST_BLOB_RADIUS_MIN = 6;
    // Size of the circles of the emitters seen: multiplies the radius above (min and max
    // included). Change only this to draw them bigger or smaller.
    const IRTEST_BLOB_RADIUS_SCALE = 3;
    const IRTEST_EMITTER_RADIUS = 25; // circle not seen / older firmware (Qt App size), not scaled
    // Detected blobs are always above the camera brightness threshold (130-150): this
    // range becomes 0..1 of the fill (radial gradient: centre = max, edge = average).
    const IRTEST_BRIGHTNESS_MIN = 128;
    const IRTEST_BRIGHTNESS_MAX = 255;
    // Exponential average of radius and brightness (not of the position) against flicker:
    // weight of the new 20 Hz sample.
    const IRTEST_SMOOTHING = 0.5;
    // Blob data older than this (ms) is not used: classic view as with older firmware.
    const IRTEST_BLOBS_MAX_AGE = 250;
    const IRTEST_BLOBS_VERSION = 1;
    const IRTEST_BLOBS_LENGTH = 21;
    const IRTEST_NOT_SEEN_COLOR = '#ff3030';
    const IRTEST_BLOB_WEAK = 0x02;          // sTestBlobs flag: emitter seen but weak (firmware IR_WEAK_MAX_BRIGHTNESS)

    // ----- IR view in the calibration window ------------------------------------------
    // Firmware that knows caliFlagIrView sends sTestBlobs + sTestCoords during calibration too:
    // the crosshair changes colour with the weakest emitter (red when one is missing), a small panel
    // in the bottom right corner shows each emitter, and the board refuses every target shot
    // (sCaliIrWarning) while one is missing. A weak emitter (orange crosshair) is accepted.
    const CALI_IR_PANEL_WIDTH = 0.14;      // of the window width (the panel is square)
    const CALI_IR_PANEL_MAX_WIDTH = 300;   // px
    const CALI_IR_PANEL_MIN_WIDTH = 140;   // px, also on small screens
    const CALI_IR_MAX_AGE = 500;           // ms without coordinates: older firmware, no panel
    const CALI_IR_WARNING_TIME = 4000;     // ms the refusal message stays on screen
    const CALI_IR_MISSING = 0x01;          // sCaliIrWarning bit
    const CALI_IR_PANEL_COLOR = [160, 160, 164]; // panel border and completed progress dots, fixed
    // Colour of each emitter seen in the panel, and of the crosshair from the emitter with the
    // lowest brightness; missing = intense red (the crosshair only; the panel keeps the dashed
    // circle). Continuous scale of the brightness, without jumps: red at the camera threshold
    // (130), light orange at CALI_IR_BRIGHTNESS_WEAK - CALI_IR_BLEND, then through yellow-green
    // to light green at CALI_IR_BRIGHTNESS_WEAK + CALI_IR_BLEND, full green at CALI_IR_FULL_GREEN.
    // CALI_IR_BRIGHTNESS_WEAK is the firmware limit (IR_WEAK_MAX_BRIGHTNESS). The colours
    // follow brightness continuously, independently of the firmware weak flag/hysteresis.
    const CALI_IR_BRIGHTNESS_MIN = 130;
    const CALI_IR_BRIGHTNESS_WEAK = 170;
    const CALI_IR_BLEND = 20;              // half width of the orange to green blend around the limit
    const CALI_IR_FULL_GREEN = 230;
    const CALI_IR_SCALE = {
        missing: [255, 20, 20],
        weakLow: [255, 40, 40], weakHigh: [255, 190, 90],     // 130, limit - blend
        goodLow: [170, 255, 140], goodHigh: [0, 220, 0],      // limit + blend, full green
    };

    // ----- Rotating target ring in calibration ----------------------------------------
    // Only the outer dashed ring rotates. The original crosshair stays still and unchanged;
    // both use the same colour. No coloured ring while verifying; reduced motion is static.
    // Acquisition holds the crosshair at its target during a solid white flash, including the final centre.
    const CROSSHAIR_RING = 24.42 / 61.44;   // ring radius / crosshair size (CROSSHAIR_SVG)
    const CROSSHAIR_COLOR = [255, 135, 88]; // #ff8758, crosshair without IR data
    const TARGET_RING_RADIUS = 48 / 39;    // outer / original ring radius, as in the preview
    const TARGET_RING_PERIOD = 8000;       // ms of one full clockwise turn
    const TARGET_RING_DASH = 21;           // degrees; ten identical segments, 15-degree gaps
    const TARGET_RING_STEP = 36;
    const TARGET_RING_WIDTH = 2.1;         // slightly thinner than the original crosshair's 2.4 stroke
    const CALI_ANIMATION_FRAME_TIME = 33;  // ms, about 30 animation frames per second

    function reducedMotion() {
        try { return !!(root.matchMedia && root.matchMedia('(prefers-reduced-motion: reduce)').matches); } catch (e) { return false; }
    }

    /** [r, g, b] + alpha -> 'rgba(...)'. */
    function rgba(rgb, alpha) {
        return `rgba(${rgb.map(Math.round).join(',')},${alpha.toFixed(3)})`;
    }

    /** ?irdebug in the address: shows area and brightness under each emitter circle. */
    function irTestDebug() {
        try { return new URLSearchParams(root.location.search).has('irdebug'); } catch (e) { return false; }
    }

    /** Brightness 0..255 -> 0..1 fill over the range of detected blobs. */
    function brightnessLevel(value) {
        const t = (value - IRTEST_BRIGHTNESS_MIN) / (IRTEST_BRIGHTNESS_MAX - IRTEST_BRIGHTNESS_MIN);
        return Math.min(1, Math.max(0, t));
    }

    /** '#rrggbb' + alpha -> rgba() */
    function withAlpha(color, alpha) {
        const n = parseInt(color.slice(1), 16);
        return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${alpha.toFixed(3)})`;
    }

    // ----- Bitmap typeface ------------------------------------------------------------

    const GLYPH = 8;
    const glyphCache = new Map();
    let crosshairImage = null;

    function fontPixels() {
        const font = OF.TEST_FONT;
        const bytes = Uint8Array.from(atob(font.data), (c) => c.charCodeAt(0));
        const pixels = new Uint8Array(font.count * GLYPH * GLYPH);
        for (let i = 0; i < pixels.length; ++i) pixels[i] = (bytes[i >> 2] >> (6 - 2 * (i & 3))) & 3;
        return pixels;
    }

    let decodedFont = null;

    // Accents drawn over lowercase letters without ascenders (rows 0-1), '#' main, 's' shadow.
    const ACCENTS = {
        '\u0300': ['.##s', '..##s'],     // grave
        '\u0301': ['...##s', '..##s'],   // acute
        '\u0302': ['..##s', '.#s.#s'],   // circumflex
        '\u0308': ['', '##s##s'],        // diaeresis
        '\u0303': ['.##s#s', '#s##s'],   // tilde
    };
    const PUNCTUATION = { '\u2018': "'", '\u2019': "'", '\u201C': '"', '\u201D': '"', '\u00AB': '"', '\u00BB': '"',
        '\u2013': '-', '\u2014': '-', '\u00A0': ' ', '\u00B0': 'o' };

    /** 8x8 pixel map (0 none, 1 main, 2 shadow) of a character, or null when the typeface has no glyph. */
    function glyphPixels(ch) {
        const mapped = PUNCTUATION[ch] || ch;
        let code = mapped.codePointAt(0);
        let accent = null;
        if (code > 126) {
            const parts = mapped.normalize('NFD');
            code = parts.codePointAt(0);
            const mark = parts.slice(1, 2);
            accent = ACCENTS[mark] || null;
            // Cedilla/ogonek have no room below the letter: plain letter.
            if (code < 33 || code > 126 || (mark && !accent && mark !== '\u0327' && mark !== '\u0328')) return null;
        }
        if (code === 32) return new Uint8Array(GLYPH * GLYPH);
        if (code < 33 || code > 126) return null;
        if (!decodedFont) decodedFont = fontPixels();
        const offset = (code - 33) * GLYPH * GLYPH;
        const pixels = decodedFont.slice(offset, offset + GLYPH * GLYPH);
        if (accent) {
            const top = pixels.subarray(0, GLYPH * 2);
            const dotted = code === 105 || code === 106; // i, j: the dot makes room for the accent
            if (top.some((v) => v) && !dotted) return pixels; // ascender or capital: no room, plain letter
            top.fill(0);
            const shift = dotted ? 1 : 0;
            accent.forEach((row, y) => {
                for (let x = 0; x < row.length && x + shift < GLYPH; ++x)
                    if (row[x] !== '.') pixels[y * GLYPH + x + shift] = row[x] === '#' ? 1 : 2;
            });
        }
        return pixels;
    }

    /** Qt GenerateText tint: main pixels take the colour, shadows a 140 darker shade. */
    function palette(tint) {
        if (!tint) return ['#ffffff', 'rgb(115,115,115)'];
        const dark = tint.map((v) => (v > 140 ? v - 140 : 0));
        return [`rgb(${tint.join(',')})`, `rgb(${dark.join(',')})`];
    }

    function glyphCanvas(ch, tint) {
        const key = ch + '|' + (tint ? tint.join(',') : '');
        if (glyphCache.has(key)) return glyphCache.get(key);
        const pixels = glyphPixels(ch);
        let canvas = null;
        if (pixels) {
            canvas = doc.createElement('canvas');
            canvas.width = GLYPH;
            canvas.height = GLYPH;
            const ctx = canvas.getContext('2d');
            const colors = palette(tint);
            for (let i = 0; i < pixels.length; ++i) {
                if (!pixels[i]) continue;
                ctx.fillStyle = colors[pixels[i] - 1];
                ctx.fillRect(i % GLYPH, Math.floor(i / GLYPH), 1, 1);
            }
        }
        glyphCache.set(key, canvas);
        return canvas;
    }

    /** Size in scene pixels of a text block (lines are centred inside it). */
    function textSize(lines, scale) {
        const width = Math.max(0, ...lines.map((line) => [...line].length)) * GLYPH * scale;
        return { width, height: lines.length * GLYPH * scale };
    }

    /** Draws lines with top-left corner (x, y); each line is centred in the block. */
    function drawText(ctx, lines, x, y, scale, tint) {
        const block = textSize(lines, scale);
        const cell = GLYPH * scale;
        const colors = palette(tint);
        lines.forEach((line, row) => {
            const chars = [...line];
            let cx = x + (block.width - chars.length * cell) / 2;
            const cy = y + row * cell;
            for (const ch of chars) {
                const glyph = glyphCanvas(ch, tint);
                if (glyph) {
                    ctx.drawImage(glyph, cx, cy, cell, cell);
                } else {
                    // Not in the typeface (other scripts): system font in the same cell.
                    ctx.save();
                    ctx.font = `bold ${Math.round(cell * 0.95)}px monospace`;
                    ctx.textAlign = 'center';
                    ctx.textBaseline = 'middle';
                    ctx.fillStyle = colors[1];
                    ctx.fillText(ch, cx + cell / 2 + scale, cy + cell / 2 + scale);
                    ctx.fillStyle = colors[0];
                    ctx.fillText(ch, cx + cell / 2, cy + cell / 2);
                    ctx.restore();
                }
                cx += cell;
            }
        });
        return block;
    }

    /** Splits each line at the spaces so that no line is longer than maxChars (longer words stay whole). */
    function wrapLines(textLines, maxChars) {
        const out = [];
        for (const line of textLines) {
            let current = '';
            for (const word of line.split(' ')) {
                if (current && [...current].length + 1 + [...word].length > maxChars) {
                    out.push(current);
                    current = word;
                } else {
                    current = current ? current + ' ' + word : word;
                }
            }
            out.push(current);
        }
        return out;
    }

    /** Translated lines: the Qt texts are padded for a left-aligned block, here lines are centred. */
    function lines(...keys) {
        return keys.map((key) => (key ? OF.i18n.t(key) : '').replace(/\u2026/g, '...').trim());
    }

    function loadCrosshair() {
        if (!crosshairImage) {
            crosshairImage = new Image();
            crosshairImage.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(CROSSHAIR_SVG);
            crosshairImage.crosshairColor = CROSSHAIR_COLOR;
        }
        return crosshairImage;
    }

    // Crosshair tinted by the calibration IR view: one image per colour (colours are rounded,
    // so only a few dozen exist).
    const tintedCrosshairs = new Map();
    function loadTintedCrosshair(rgb) {
        const key = rgb.map((v) => Math.round(v / 8) * 8).map((v) => Math.min(255, v));
        const color = '#' + key.map((v) => v.toString(16).padStart(2, '0')).join('');
        let image = tintedCrosshairs.get(color);
        if (!image) {
            image = new Image();
            image.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(CROSSHAIR_SVG.replace('#ff8758', color));
            image.crosshairColor = key;
            tintedCrosshairs.set(color, image);
        }
        return image;
    }

    function mix(a, b, t) {
        const k = Math.min(1, Math.max(0, t));
        return a.map((v, i) => v + (b[i] - v) * k);
    }

    /** Calibration IR view: colour [r, g, b] of an emitter seen, from its brightness (continuous scale). */
    function caliIrEmitterColor(blob) {
        const orange = CALI_IR_BRIGHTNESS_WEAK - CALI_IR_BLEND;
        const green = CALI_IR_BRIGHTNESS_WEAK + CALI_IR_BLEND;
        if (blob.max < orange)
            return mix(CALI_IR_SCALE.weakLow, CALI_IR_SCALE.weakHigh,
                (blob.max - CALI_IR_BRIGHTNESS_MIN) / (orange - CALI_IR_BRIGHTNESS_MIN));
        if (blob.max < green)
            return mix(CALI_IR_SCALE.weakHigh, CALI_IR_SCALE.goodLow, (blob.max - orange) / (green - orange));
        return mix(CALI_IR_SCALE.goodLow, CALI_IR_SCALE.goodHigh,
            (blob.max - green) / (CALI_IR_FULL_GREEN - green));
    }

    /** [r, g, b] -> '#rrggbb'. */
    function rgbHex(rgb) {
        return '#' + rgb.map((v) => Math.round(v).toString(16).padStart(2, '0')).join('');
    }

    // ----- Window -------------------------------------------------------------------------

    class FullscreenWindow {
        constructor(mode, options = {}) {
            this.mode = mode;
            this.options = options;
            this.stage = STAGE.init;
            this.values = {};
            this.coords = null;
            this.coordsTime = 0;     // when the last sTestCoords arrived (calibration IR panel)
            this.irWarning = null;   // calibration: { bits, time } of the last refused shot
            this.blobs = null;       // sTestBlobs: { time, entries[4] }
            this.blobLevels = [null, null, null, null]; // smoothed { radius, center, edge } per emitter
            this.irDebug = irTestDebug();
            this.mouse = null;
            this.closed = false;
            this.targetTime = Date.now(); // calibration: start of the target ring's rotation
            this.acceptedShot = null; // { stage, positionStage, time }: acquisition and held visual target
            this.resetValues();
        }

        resetValues() {
            for (const key of CALI_KEYS) this.values[key] = -1;
            const prefixes = CALI_PREFIXES.map((key) => OF.i18n.t(key));
            const width = Math.max(...prefixes.map((text) => [...text].length));
            this.infoPrefixes = prefixes.map((text) => ' '.repeat(width - [...text].length) + text);
            this.infoText = this.infoPrefixes.slice();
            this.infoVisible = false;
        }

        open() {
            const canvas = doc.createElement('canvas');
            canvas.className = 'fullscreen-canvas';
            const overlay = doc.createElement('div');
            overlay.className = `fullscreen-window mode-${this.mode}`;
            overlay.tabIndex = -1;
            overlay.setAttribute('role', 'dialog');
            overlay.setAttribute('aria-label', this.title());
            // Shown only when the browser refuses fullscreen mode (it needs a recent click).
            const retry = doc.createElement('button');
            retry.className = 'fullscreen-retry';
            retry.type = 'button';
            retry.hidden = true;
            retry.textContent = OF.i18n.t('Click here to show this window in fullscreen');
            retry.addEventListener('click', () => this._requestFullscreen());
            overlay.append(canvas, retry);
            (doc.getElementById('overlay-root') || doc.body).append(overlay);
            this.overlay = overlay;
            this.canvas = canvas;
            this.retryButton = retry;

            // The page behind stays out of reach (keyboard focus, screen readers), like a Qt fullscreen window.
            const app = doc.getElementById('app');
            this._appWasInert = app ? app.inert : false;
            if (app) app.inert = true;

            this._onKey = (event) => {
                if (event.key === 'Escape') {
                    if (OF.UI && OF.UI.modalOpen && OF.UI.modalOpen()) return; // ESC closes the message box first
                    event.preventDefault();
                    event.stopPropagation();
                    this.escape();
                }
            };
            this._onResize = () => this.render();
            this._onFullscreen = () => {
                if (!doc.fullscreenElement && this._wasFullscreen && !this.closed) this.escape();
                this._wasFullscreen = !!doc.fullscreenElement;
                this.render();
            };
            this._onMove = (event) => {
                if (this.mode === MODE_CALIBRATE && this.stage === STAGE.verify) {
                    this.mouse = { x: event.clientX, y: event.clientY };
                    this.render();
                }
            };
            doc.addEventListener('keydown', this._onKey, true);
            root.addEventListener('resize', this._onResize);
            doc.addEventListener('fullscreenchange', this._onFullscreen);
            overlay.addEventListener('pointermove', this._onMove);
            overlay.addEventListener('contextmenu', (event) => event.preventDefault());

            this._wasFullscreen = false;
            this._requestFullscreen();
            overlay.focus();
            loadCrosshair().decode?.().then(() => this.render()).catch(() => {});
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

        title() {
            const t = OF.i18n.t.bind(OF.i18n);
            if (this.mode === MODE_ALIGNMENT) return t('Alignment Assistant');
            if (this.mode === MODE_IRTEST) return t('IR Emitters Test');
            return t('Calibration Window');
        }

        escape() {
            if (this.closed) return;
            if (this.mode === MODE_CALIBRATE) {
                if (this.options.onExitRequest) this.options.onExitRequest();
            } else {
                this.shutdown();
            }
        }

        /** Closes the window; calibration returns its values (all -1 when not completed). */
        shutdown(values) {
            if (this.closed) return;
            this.closed = true;
            clearTimeout(this._irWarningTimer);
            clearTimeout(this._irExpiryTimer);
            clearTimeout(this._targetRingTimer);
            this._targetRingTimer = 0;
            clearTimeout(this._shotTimer);
            this._shotTimer = 0;
            root.cancelAnimationFrame(this._frame);
            this.acceptedShot = null;
            doc.removeEventListener('keydown', this._onKey, true);
            root.removeEventListener('resize', this._onResize);
            doc.removeEventListener('fullscreenchange', this._onFullscreen);
            if (doc.fullscreenElement === this.overlay && doc.exitFullscreen) doc.exitFullscreen().catch(() => {});
            this.overlay.remove();
            const app = doc.getElementById('app');
            if (app) app.inert = this._appWasInert;
            if (this.options.onExit) this.options.onExit(this.mode, values || null);
        }

        // ----- Calibration events ----------------------------------------------------------

        setStage(stage) {
            if (this.closed || this.mode !== MODE_CALIBRATE) return;
            if (!(stage >= STAGE.init && stage <= STAGE.end)) return; // unknown stage: Qt CaliModeSet does nothing
            if (stage === STAGE.end) {
                this.shutdown(Object.assign({}, this.values));
                return;
            }
            this._setCalibrationStage(stage);
            this.irWarning = null; // a new stage: the refused shot is no longer current
            if (stage === STAGE.init) {
                this.resetValues();
                this.mouse = null;
            } else if (stage >= STAGE.top) {
                this.infoVisible = true;
            }
            this.render();
        }

        /** Only a consecutive forward update acknowledges a shot; duplicates and resets do not. */
        _setCalibrationStage(stage) {
            if (stage === this.stage) return;
            const now = Date.now();
            const previous = this.acceptedShot;
            const positionStage = previous && now - previous.time < CALI_SHOT_TIME ? previous.positionStage : this.stage;
            this.acceptedShot = stage === this.stage + 1 && stage <= STAGE.verify ?
                { stage: this.stage, positionStage, time: now } : null;
            clearTimeout(this._shotTimer);
            this._shotTimer = 0;
            clearTimeout(this._targetRingTimer);
            this._targetRingTimer = 0;
            this.targetTime = now;
            this.stage = stage;
            if (stage > STAGE.center) {
                clearTimeout(this._targetRingTimer);
                this._targetRingTimer = 0;
            }
        }

        /** sCaliInfoUpd payload: type (1-4 int32 offsets, 5-6 float LED positions) + 4 bytes. */
        setInfo(payload) {
            if (this.closed || this.mode !== MODE_CALIBRATE || payload.length !== 5) return;
            const type = payload[0];
            if (type < 1 || type > 6) return;
            const view = new DataView(payload.buffer, payload.byteOffset + 1, 4);
            const value = type <= 4 ? view.getInt32(0, true) : view.getFloat32(0, true);
            this.values[CALI_KEYS[type - 1]] = value;
            // Qt QString::number(int) / QString::number(float, 'g', 6)
            const shown = type <= 4 ? String(value) : (OF.UI && OF.UI.formatG ? OF.UI.formatG(value) : String(Number(value.toPrecision(6))));
            this.infoText[type - 1] = this.infoPrefixes[type - 1] + shown;
            if (type === 6) this._setCalibrationStage(STAGE.verify); // may arrive after the verify stage
            this.render();
        }

        /** sTestCoords: 12 int32 (TL, TR, BL, BR with the outside-FOV flag in bit 0 of X; mouse; D). */
        drawTest(payload) {
            if (this.closed || (this.mode !== MODE_IRTEST && this.mode !== MODE_CALIBRATE)) return;
            if (payload.length !== 48) return;
            const view = new DataView(payload.buffer, payload.byteOffset, 48);
            this.coords = Array.from({ length: 12 }, (_, i) => view.getInt32(i * 4, true));
            this.coordsTime = Date.now();
            this.render();
            this._armIrExpiry();
        }

        /**
         * Calibration IR view: redraws when the IR data gets old (blobs first, then coordinates),
         * so that a link that stops sending without closing does not leave a green crosshair.
         */
        _armIrExpiry() {
            clearTimeout(this._irExpiryTimer);
            if (this.closed || this.mode !== MODE_CALIBRATE) return;
            const now = Date.now();
            const ends = [this.coordsTime + CALI_IR_MAX_AGE, this.blobs ? this.blobs.time + IRTEST_BLOBS_MAX_AGE : 0]
                .filter((end) => end >= now);
            if (!ends.length) return;
            this._irExpiryTimer = setTimeout(() => { this.render(); this._armIrExpiry(); }, Math.min(...ends) - now + 20);
        }

        /** Calibration: the board refused a target shot (sCaliIrWarning bit 1: an emitter is missing). */
        showIrWarning(bits) {
            if (this.closed || this.mode !== MODE_CALIBRATE || !bits) return;
            this.irWarning = { bits, time: Date.now() };
            clearTimeout(this._irWarningTimer);
            this._irWarningTimer = setTimeout(() => this.render(), CALI_IR_WARNING_TIME + 50);
            this.render();
        }

        /**
         * sTestBlobs (sent just before sTestCoords): version, then for each emitter in the
         * sTestCoords order: flags (bit 0 seen), average brightness, max brightness, area (uint16 LE).
         * Unknown versions and lengths are ignored (classic view).
         */
        setTestBlobs(payload) {
            if (this.closed || (this.mode !== MODE_IRTEST && this.mode !== MODE_CALIBRATE)) return;
            if (payload.length < IRTEST_BLOBS_LENGTH || payload[0] !== IRTEST_BLOBS_VERSION) return;
            const entries = [];
            for (let i = 0; i < 4; ++i) {
                const o = 1 + i * 5;
                entries.push({
                    seen: (payload[o] & 1) !== 0,
                    weak: (payload[o] & IRTEST_BLOB_WEAK) !== 0,
                    avg: payload[o + 1],
                    max: payload[o + 2],
                    area: payload[o + 3] | (payload[o + 4] << 8),
                });
            }
            if (!this._freshBlobs()) this.blobLevels = [null, null, null, null]; // after a gap: no stale average
            this.blobs = { time: Date.now(), entries };
            this._updateBlobLevels(entries);
            this._armIrExpiry();
        }

        /** Current blob data, or null when missing or stale (older firmware: classic view). */
        _freshBlobs() {
            const b = this.blobs;
            return b && Date.now() - b.time <= IRTEST_BLOBS_MAX_AGE ? b.entries : null;
        }

        /** Circle of each emitter from the blob data, smoothed over the 20 Hz updates. */
        _updateBlobLevels(entries) {
            for (let i = 0; i < 4; ++i) {
                const blob = entries[i];
                if (!blob.seen) {
                    this.blobLevels[i] = null; // a new sighting starts from its own values
                    continue;
                }
                const k = Math.sqrt(Math.min(blob.area, IRTEST_BLOB_AREA_MAX) / IRTEST_BLOB_AREA_MAX);
                const target = {
                    radius: Math.max(IRTEST_BLOB_RADIUS_MIN, IRTEST_BLOB_RADIUS_MAX * k) * IRTEST_BLOB_RADIUS_SCALE,
                    center: brightnessLevel(blob.max),
                    edge: brightnessLevel(blob.avg),
                };
                const prev = this.blobLevels[i];
                if (!prev) {
                    this.blobLevels[i] = target;
                } else {
                    const a = IRTEST_SMOOTHING;
                    this.blobLevels[i] = {
                        radius: prev.radius + a * (target.radius - prev.radius),
                        center: prev.center + a * (target.center - prev.center),
                        edge: prev.edge + a * (target.edge - prev.edge),
                    };
                }
            }
        }

        // ----- Rendering ----------------------------------------------------------------------

        textScale(type) {
            const h = this.height;
            const table = h >= 1440 ? [5, 4, 3, 4] : h >= 1080 ? [4, 3, 2, 3] : h >= 720 ? [3, 2, 2, 3] : [2, 2, 1, 2];
            return table[{ heading: 0, sub: 1, small: 2, crosshair: 3 }[type]];
        }

        render() {
            if (this.closed || !this.canvas) return;
            if (this._frame) return;
            this._frame = root.requestAnimationFrame(() => {
                this._frame = 0;
                if (!this.closed) this._draw();
            });
        }

        _draw() {
            const rect = this.overlay.getBoundingClientRect();
            const dpr = root.devicePixelRatio || 1;
            this.width = rect.width;
            this.height = rect.height;
            const canvas = this.canvas;
            const pw = Math.round(rect.width * dpr);
            const ph = Math.round(rect.height * dpr);
            if (canvas.width !== pw || canvas.height !== ph) {
                canvas.width = pw;
                canvas.height = ph;
            }
            const ctx = canvas.getContext('2d');
            ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
            ctx.imageSmoothingEnabled = false;

            if (this.mode === MODE_CALIBRATE) this._drawCalibration(ctx);
            else if (this.mode === MODE_ALIGNMENT) this._drawAlignment(ctx);
            else this._drawIRTest(ctx);
        }

        _background(ctx, color) {
            ctx.fillStyle = color;
            ctx.fillRect(0, 0, this.width, this.height);
        }

        /** One line centred on x, with its top at y, drawn in pieces of their own colour: [[text, tint], ...]. */
        _centeredSegments(ctx, segments, y, scale) {
            const cell = GLYPH * scale;
            let x = this.width / 2 - segments.reduce((n, [text]) => n + [...text].length, 0) * cell / 2;
            for (const [text, tint] of segments) {
                drawText(ctx, [text], x, y, scale, tint);
                x += [...text].length * cell;
            }
        }

        /** Text block centred on x, with its top at y. */
        _centered(ctx, text, y, scale, tint) {
            const size = textSize(text, scale);
            drawText(ctx, text, this.width / 2 - size.width / 2, y, scale, tint);
            return size;
        }

        _polyline(ctx, points, color, width) {
            ctx.beginPath();
            points.forEach(([x, y], i) => (i ? ctx.lineTo(x, y) : ctx.moveTo(x, y)));
            ctx.closePath();
            ctx.strokeStyle = color;
            ctx.lineWidth = width;
            ctx.stroke();
        }

        _crossLines() {
            const w = this.width;
            const h = this.height;
            return [[w / 2, -10], [w / 2, h + 10], [-10, h + 10], [-10, h / 2], [w + 10, h / 2], [w + 10, -10]];
        }

        _drawCalibration(ctx) {
            const w = this.width;
            const h = this.height;
            const heading = this.textScale('heading');
            const sub = this.textScale('sub');
            const red = [225, 25, 25];
            this._background(ctx, 'dimgray');

            const irColor = this._caliIrColor();
            let image = irColor ? loadTintedCrosshair(irColor) : loadCrosshair();
            if (irColor && !(image.complete && image.naturalWidth)) {
                // A new colour is still loading: keep the previous one for this frame.
                image.onload = () => this.render();
                image = this._lastCrosshair || loadCrosshair();
            }
            this._drawTargetRing(ctx, image.crosshairColor); // below the original crosshair and guides
            if (this.stage !== STAGE.init) this._polyline(ctx, this._crossLines(), 'white', 2);

            let stageText;
            let header;
            let tutorial;
            let stageY = h * 0.25;
            let headerStageY = null;  // the header keeps the previous stage position (Qt rare verify branch)
            let tint = null;

            switch (this.stage) {
            case STAGE.init:
                stageText = lines('Initialize Calibration:');
                header = lines('Shoot at the target to start calibration.');
                tutorial = lines('Calibration can be exited without changes', '  by pressing either Button A, Button B, ', '       or Button C (if available).       ');
                break;
            case STAGE.top:
            case STAGE.bottom:
            case STAGE.left:
            case STAGE.right:
            case STAGE.center:
                stageText = lines(`Step ${this.stage}:`);
                header = lines(['', 'Shoot at the top edge of the screen.', 'Shoot at the bottom edge of the screen.',
                    'Shoot at the left edge of the screen.', 'Shoot at the right edge of the screen.',
                    'Shoot at the final target in the center.'][this.stage]);
                tutorial = lines('The calibration process can be reset by pressing', 'either Button A or Button B, and can be canceled', '      by pressing Button C (if available).      ');
                break;
            default: {
                stageY = h * 0.15;
                const v = this.values;
                const inRange = CALI_KEYS.every((key) => v[key] >= -32768 && v[key] <= 32768);
                if (inRange) {
                    stageText = lines('Verify aiming:');
                    header = lines('    Confirm that the bullseye     ', '   lines up with the gun sight.   ', 'If this calibration is acceptable,', '  confirm by pulling the trigger. ');
                    tutorial = lines("       If this target accuracy isn't desirable,      ", '          press either Button A or Button B          ',
                        '         to restart the calibration process.         ', '', '[You can also exit calibration without saving changes', '        by pressing Button C (if available).]        ');
                } else if (v.TLled === -1 || v.TRled === -1) {
                    stageText = lines('Verify aiming:');
                    header = lines('Shoot at the final target in the center.');
                    headerStageY = h * 0.25;
                    tutorial = [];
                } else {
                    tint = red;
                    stageText = lines('WARNING: Possibly Malformed Calibration!!');
                    header = lines('  The current pending values for this profile  ', 'will likely cause incorrect or broken tracking.');
                    tutorial = lines('  Press Button A or Button B to restart calibration, ', '   or pull trigger to continue with these settings.  ', '',
                        '[You can also exit calibration without saving changes', '        by pressing Button C (if available).]        ');
                }
                break;
            }
            }

            const stageSize = textSize(stageText, heading);
            const stageTop = stageY - stageSize.height / 2;
            this._drawCalibrationProgress(ctx, stageTop);
            this._centered(ctx, stageText, stageTop, heading, tint);
            const headerTop = headerStageY === null ? stageTop + stageSize.height :
                headerStageY - textSize(lines('Step 5:'), heading).height / 2 + textSize(lines('Step 5:'), heading).height;
            const headerSize = this._centered(ctx, header, headerTop, heading, tint);
            // IR view: a target shot is refused only with the crosshair red (an emitter missing).
            // The colour word (%1) is drawn in the full green of the crosshair.
            if (this._caliIrColor() && this.stage <= STAGE.center) {
                const [before, after = ''] = lines('Preferably shoot when the crosshair is %1.')[0].split('%1');
                this._centeredSegments(ctx, [[before], [OF.i18n.t('green'), CALI_IR_SCALE.goodHigh], [after]],
                    headerTop + headerSize.height + GLYPH * sub, sub);
            }
            // A refused shot (IR view): the reason takes the place of the tutorial for a few seconds.
            // At most 60% of the width, so that the IR panel at the bottom right stays clear.
            const warningChars = Math.max(16, Math.floor((w * 0.6) / (GLYPH * sub)));
            const warningLines = (bits) => wrapLines(this._caliIrWarningText(bits), warningChars);
            // The IR panel is sized for this text too, so it keeps its size when a warning appears.
            let textRight = w / 2 + textSize(warningLines(CALI_IR_MISSING), sub).width / 2 + 8 * sub;
            const irWarning = this.irWarning;
            if (irWarning && Date.now() - irWarning.time <= CALI_IR_WARNING_TIME) {
                const scale = sub;
                const warning = warningLines(irWarning.bits);
                const size = textSize(warning, scale);
                const top = h * 0.8 - size.height / 2;
                ctx.save();
                ctx.fillStyle = 'rgba(0, 0, 0, 0.6)';
                ctx.fillRect(w / 2 - size.width / 2 - 8 * scale, top - 6 * scale, size.width + 16 * scale, size.height + 12 * scale);
                ctx.restore();
                this._centered(ctx, warning, top, scale, [255, 90, 90]);
            } else if (tutorial.length) {
                const size = textSize(tutorial, sub);
                this._centered(ctx, tutorial, h * 0.8 - size.height / 2, sub, tint);
            }
            if (tutorial.length) textRight = Math.max(textRight, w / 2 + textSize(tutorial, sub).width / 2);

            if (this.infoVisible) {
                let y = h / 2 - 24 * sub;
                for (const text of this.infoText) {
                    drawText(ctx, [text], 80 * sub, y, sub);
                    y += GLYPH * sub;
                }
            }

            const [x, y] = this._crosshairPosition();
            const size = CROSSHAIR_SIZE * this.textScale('crosshair');
            if (image.complete && image.naturalWidth) {
                ctx.imageSmoothingEnabled = true;
                ctx.drawImage(image, x - size / 2, y - size / 2, size, size);
                ctx.imageSmoothingEnabled = false;
                this._lastCrosshair = image;
            }

            this._drawCaliIrPanel(ctx, textRight);
        }

        /** Six targets: centre, top, bottom, left, right, centre; green is progress, not IR quality. */
        _drawCalibrationProgress(ctx, stageTop) {
            const scale = this.textScale('sub') / 2;
            const y = stageTop - 22 * scale;
            const x = this.width / 2 - 2.5 * 22 * scale;
            ctx.save();
            for (let i = 0; i < 6; ++i) {
                const current = i === this.stage;
                const done = i < this.stage;
                ctx.beginPath();
                ctx.arc(x + i * 22 * scale, y, (current ? 8 : 6) * scale, 0, Math.PI * 2);
                ctx.fillStyle = current ? '#00dc00' : 'rgb(160,160,164)';
                ctx.strokeStyle = 'rgb(160,160,164)';
                ctx.lineWidth = 2 * scale;
                if (current || done) ctx.fill();
                else ctx.stroke();
            }
            ctx.restore();
        }

        /** Hold only the drawing during the flash; verification cursor data continues to update. */
        _crosshairPosition(stage = this.acceptedShot ? this.acceptedShot.positionStage : this.stage) {
            const w = this.width;
            const h = this.height;
            const positions = [[w / 2, h / 2], [w / 2, 0], [w / 2, h], [0, h / 2], [w, h / 2], [w / 2, h / 2]];
            if (stage === STAGE.verify && this.mouse) return [this.mouse.x, this.mouse.y];
            return positions[stage] || positions[5];
        }

        /** Outer dashed target ring; colour comes from the SVG actually drawn in this frame. */
        _drawTargetRing(ctx, color = CROSSHAIR_COLOR) {
            const still = reducedMotion();
            if (this.closed || this.stage > STAGE.center || still || this.acceptedShot) {
                clearTimeout(this._targetRingTimer);
                this._targetRingTimer = 0;
            }
            if (this.closed) return;
            if (this._drawShotConfirmation(ctx, still)) return;
            if (this.stage <= STAGE.center) {
                this._strokeTargetRing(ctx, this._crosshairPosition(), still ? 0 : Date.now() - this.targetTime, rgbHex(color));
                if (!still && !this._targetRingTimer)
                    this._targetRingTimer = setTimeout(() => { this._targetRingTimer = 0; this.render(); }, CALI_ANIMATION_FRAME_TIME);
            }
        }

        /** Solid white flash; only the crosshair's presentation waits for its expiry. */
        _drawShotConfirmation(ctx, still) {
            const shot = this.acceptedShot;
            if (!shot) return false;
            const elapsed = Math.max(0, Date.now() - shot.time);
            if (elapsed >= CALI_SHOT_TIME) {
                this.acceptedShot = null;
                clearTimeout(this._shotTimer);
                this._shotTimer = 0;
                this.targetTime = Date.now();
                return false;
            }
            const alpha = still ? 1 : Math.min(1, (CALI_SHOT_TIME - elapsed) / (CALI_SHOT_TIME - CALI_SHOT_HOLD));
            this._strokeTargetRing(ctx, this._crosshairPosition(shot.positionStage), 0,
                rgba([255, 255, 255], alpha), true);
            // The final target also confirms during verification; reduced motion uses one expiry timer.
            if (!this._shotTimer)
                this._shotTimer = setTimeout(() => {
                    this._shotTimer = 0;
                    this.render();
                }, still ? CALI_SHOT_TIME - elapsed : Math.min(CALI_ANIMATION_FRAME_TIME, CALI_SHOT_TIME - elapsed));
            return true;
        }

        /** Same radius/stroke for the dashed rotating target and its solid stationary flash. */
        _strokeTargetRing(ctx, position, elapsed, color, continuous = false) {
            const [x, y] = position;
            const scale = this.textScale('crosshair');
            const radius = CROSSHAIR_SIZE * scale * CROSSHAIR_RING * TARGET_RING_RADIUS;
            const arcUnit = 2 * Math.PI * radius / 360;
            ctx.save();
            ctx.translate(x, y);
            if (!continuous) ctx.rotate((Math.max(0, elapsed) % TARGET_RING_PERIOD) * 2 * Math.PI / TARGET_RING_PERIOD);
            ctx.strokeStyle = color;
            ctx.lineWidth = TARGET_RING_WIDTH * scale;
            ctx.lineCap = 'round';
            ctx.setLineDash(continuous ? [] : [TARGET_RING_DASH * arcUnit, (TARGET_RING_STEP - TARGET_RING_DASH) * arcUnit]);
            ctx.beginPath();
            ctx.arc(0, 0, radius, 0, Math.PI * 2);
            ctx.stroke();
            ctx.restore();
        }

        /**
         * Calibration IR view: crosshair colour [r, g, b] from the worst emitter, or null without
         * IR data (older firmware: the crosshair keeps its colour). Red when the board refuses the
         * shot (an emitter missing), otherwise the continuous scale of the weakest emitter.
         */
        _caliIrColor() {
            if (!this.coords || Date.now() - this.coordsTime > CALI_IR_MAX_AGE) return null;
            const blobs = this._freshBlobs();
            if (!blobs) return null;
            if (blobs.some((b) => !b.seen)) return CALI_IR_SCALE.missing;
            // The emitter with the lowest brightness gives the colour, so the crosshair always has
            // the colour of one of the four circles.
            const worst = blobs.reduce((a, b) => (b.max < a.max ? b : a));
            return caliIrEmitterColor(worst);
        }

        /**
         * Calibration: a square, rounded IR diagram at the bottom right, without text.
         * Drawn only while the firmware sends the coordinates (caliFlagIrView).
         */
        _drawCaliIrPanel(ctx, textRight) {
            const c = this.coords;
            if (!c || Date.now() - this.coordsTime > CALI_IR_MAX_AGE) return;
            const w = this.width;
            const h = this.height;
            const small = this.textScale('small');
            // Right of the tutorial text (bottom centre), which it must not cover on small screens.
            const free = w - textRight - 2 * 12 * small;
            // The emitters are drawn in their layout on a square of this side, not where the camera sees them.
            const side = Math.max(CALI_IR_PANEL_MIN_WIDTH, Math.min(w * CALI_IR_PANEL_WIDTH, CALI_IR_PANEL_MAX_WIDTH, free));
            const margin = 12 * small;

            const blobs = this._freshBlobs();
            const x0 = w - margin - side;
            const y0 = h - margin - side;
            const radius = side * 0.075;

            ctx.save();
            // Explicit path also works on browsers without CanvasRenderingContext2D.roundRect.
            ctx.beginPath();
            ctx.moveTo(x0 + radius, y0);
            ctx.lineTo(x0 + side - radius, y0);
            ctx.quadraticCurveTo(x0 + side, y0, x0 + side, y0 + radius);
            ctx.lineTo(x0 + side, y0 + side - radius);
            ctx.quadraticCurveTo(x0 + side, y0 + side, x0 + side - radius, y0 + side);
            ctx.lineTo(x0 + radius, y0 + side);
            ctx.quadraticCurveTo(x0, y0 + side, x0, y0 + side - radius);
            ctx.lineTo(x0, y0 + radius);
            ctx.quadraticCurveTo(x0, y0, x0 + radius, y0);
            ctx.closePath();
            ctx.fillStyle = 'rgba(0, 0, 0, 0.3)';
            ctx.fill();
            ctx.strokeStyle = `rgb(${CALI_IR_PANEL_COLOR.join(',')})`;
            ctx.lineWidth = 2;
            ctx.stroke();
            ctx.clip();

            // Each emitter at its place in the layout, on a virtual square: Square at the corners
            // (sTestCoords order TL, TR, BL, BR), Diamond at the middle of the sides (Diamond slots:
            // 0 top, 1 left, 2 bottom, 3 right on the screen). No lines between them.
            const d = side * 0.25;
            const places = (this.options.diamond ?
                [[side / 2, d], [d, side / 2], [side / 2, side - d], [side - d, side / 2]] :
                [[d, d], [side - d, d], [d, side - d], [side - d, side - d]]);
            // Same radii as before: the largest seen blob gets 18% of the side; missing gets 9%.
            const seenScale = side * 0.18 / (IRTEST_BLOB_RADIUS_MAX * IRTEST_BLOB_RADIUS_SCALE);
            const notSeenScale = side * 0.09 / IRTEST_EMITTER_RADIUS;
            const colors = ['#00ff00', '#00ff00', '#00ffff', '#00ffff'];
            places.forEach(([px, py], i) => {
                const blob = blobs ? blobs[i] : { seen: c[i * 2] % 2 === 0 };
                // Seen: the crosshair scale with its own brightness (older firmware: layout colour);
                // not seen: red dashed circle with the red X (only here, the IR test view keeps its colours).
                const emitterColor = !blob.seen ? IRTEST_NOT_SEEN_COLOR :
                    blobs ? rgbHex(caliIrEmitterColor(blob)) : colors[i];
                // Without blob data (older firmware) the classic circle size is used.
                const level = blobs ? this.blobLevels[i] : { radius: IRTEST_EMITTER_RADIUS, center: 1, edge: 0.4 };
                const scale = blob && blob.seen && level ? seenScale : notSeenScale;
                this._drawEmitter(ctx, x0 + px, y0 + py, scale, emitterColor, blob, level);
            });
            ctx.restore();

        }

        /** Calibration: lines explaining why the board refused a target shot (sCaliIrWarning bits). */
        _caliIrWarningText(bits) {
            return lines(
                bits & CALI_IR_MISSING ? 'Shot refused: the camera does not see all four IR emitters.' : '',
                'Check the emitters, your distance and the IR sensitivity, then shoot again.');
        }

        _drawAlignment(ctx) {
            const w = this.width;
            const h = this.height;
            const heading = this.textScale('heading');
            const sub = this.textScale('sub');
            const small = this.textScale('small');
            this._background(ctx, 'darkslategray');

            // Square layout: emitters at the top and bottom, base width independent from the aspect ratio.
            const offset = (h * 0.711) / 2;
            const leftX = w / 2 - offset;
            const rightX = w / 2 + offset;
            // QGraphicsRectItem: brush colour and the default pen (black, 1 px).
            const box = (x1, y1, x2, y2, color) => {
                const x = Math.min(x1, x2);
                const y = Math.min(y1, y2);
                ctx.fillStyle = color;
                ctx.fillRect(x, y, Math.abs(x2 - x1), Math.abs(y2 - y1));
                ctx.strokeStyle = '#000';
                ctx.lineWidth = 1;
                ctx.strokeRect(x, y, Math.abs(x2 - x1), Math.abs(y2 - y1));
            };
            // Same order as the Qt scene: square and diamond box of each index.
            const square = [[leftX, -10, 10 * sub], [rightX, -10, 10 * sub], [leftX, h + 10, h - 10 * sub], [rightX, h + 10, h - 10 * sub]];
            const diamond = [
                [w / 2 - 15 * sub, -10, w / 2 + 15 * sub, 10 * sub],
                [w / 2 - 15 * sub, h + 10, w / 2 + 15 * sub, h - 10 * sub],
                [-10, h / 2 - 15 * sub, 10 * sub, h / 2 + 15 * sub],
                [w + 10, h / 2 - 15 * sub, w - 10 * sub, h / 2 + 15 * sub],
            ];
            for (let i = 0; i < 4; ++i) {
                const [x, y1, y2] = square[i];
                box(x - 15 * sub, y1, x + 15 * sub, y2, 'firebrick');
                box(...diamond[i], 'olivedrab');
            }

            this._polyline(ctx, [[leftX, -10], [leftX, h + 10], [rightX, h + 10], [rightX, -10]], 'firebrick', 2);
            this._polyline(ctx, this._crossLines(), 'olivedrab', 2);

            const header = lines('       Depending on your desired layout,       ', '      your IR emitters should be aligned       ',
                '         to either one of the two sets         ', '               of colored boxes:               ');
            const headerSize = textSize(header, heading);
            this._centered(ctx, header, h * 0.15 - headerSize.height / 2, heading);

            const leftText = lines('      For Square Layout,     ', '    the emitters should be   ', '    placed at the top and    ',
                '    bottom of the display;   ', 'each one being aligned to the');
            const leftColored = lines('      Red-colored boxes.     ');
            const rightText = lines('     For Diamond Layout,     ', 'the emitters should be placed', '  at the center of the four  ',
                '    edges of the display;    ', 'each one being aligned to the');
            const rightColored = lines('     Green-colored boxes.    ');

            const drawSide = (text, colored, tint, alignRight) => {
                const size = textSize(text.concat(colored), small);
                const main = textSize(text, small);
                const x = alignRight ? w * 0.98 - size.width : w * 0.02;
                const y = h * 0.3 + main.height / 2;
                drawText(ctx, text, x + (size.width - main.width) / 2, y, small);
                const coloredSize = textSize(colored, small);
                drawText(ctx, colored, x + (size.width - coloredSize.width) / 2, y + main.height, small, tint);
            };
            drawSide(leftText, leftColored, [255, 100, 100], false);
            drawSide(rightText, rightColored, [100, 255, 100], true);

            const tutorial = lines('Press ESC to exit alignment tool.');
            const size = textSize(tutorial, sub);
            drawText(ctx, tutorial, w * 0.05, h * 0.9 - size.height / 2, sub);
        }

        _drawIRTest(ctx) {
            const w = this.width;
            const h = this.height;
            const heading = this.textScale('heading');
            const sub = this.textScale('sub');
            this._background(ctx, 'midnightblue');

            const header = lines('     The array of shapes displayed onscreen     ', 'represents the emitters that the camera can see.',
                '   The colored points should move opposite to   ', '     your aim, and the gray circle should be    ', '          lining up with your gun sight.        ');
            const headerSize = textSize(header, heading);
            this._centered(ctx, header, h * 0.1 - headerSize.height / 2, heading);
            this._centered(ctx, lines('Press ESC to exit test mode.'), h * 0.85, sub);

            const c = this.coords;
            if (!c) return;

            // Emitters and the D point use the 1920x1080 space kept at its aspect ratio;
            // the mouse position covers the whole screen.
            const scaleX = w / 1920;
            const scaleY = h / 1080;
            const scale = Math.min(scaleX, scaleY);
            const offsetX = (w - 1920 * scale) / 2;
            const offsetY = (h - 1080 * scale) / 2;
            const px = (x) => offsetX + x * scale;
            const py = (y) => offsetY + y * scale;

            const points = [];
            for (let i = 0; i < 4; ++i) {
                const encodedX = c[i * 2];
                const outside = encodedX % 2 !== 0;
                points.push({ x: (encodedX - (outside ? 1 : 0)) / 2, y: c[i * 2 + 1], outside });
            }

            // Circles (Qt::green/cyan/gray/red, 3 px) in a scaled space, like the Qt scene items:
            // the mouse circle uses the screen scale on each axis, its outline too.
            const circle = (x0, y0, sx, sy, cx, cy, color, fill) => {
                ctx.save();
                ctx.translate(x0, y0);
                ctx.scale(sx, sy);
                ctx.beginPath();
                ctx.ellipse(cx, cy, 25, 25, 0, 0, Math.PI * 2);
                if (fill) { ctx.fillStyle = color; ctx.fill(); }
                ctx.strokeStyle = color;
                ctx.lineWidth = 3;
                ctx.stroke();
                ctx.restore();
            };
            const colors = ['#00ff00', '#00ff00', '#00ffff', '#00ffff'];
            const blobs = this._freshBlobs();
            if (blobs) {
                points.forEach((p, i) => this._drawEmitter(ctx, px(p.x), py(p.y), scale, colors[i], blobs[i], this.blobLevels[i]));
            } else {
                points.forEach((p, i) => circle(offsetX, offsetY, scale, scale, p.x, p.y, colors[i], p.outside));
            }
            circle(0, 0, scaleX, scaleY, c[8], c[9], '#a0a0a4', false);
            circle(offsetX, offsetY, scale, scale, c[10], c[11], '#ff0000', false);

            // The emitters box is added last to the Qt scene: drawn over the circles.
            const [tl, tr, bl, br] = points;
            this._polyline(ctx, [tl, tr, br, bl].map((p) => [px(p.x), py(p.y)]), '#a0a0a4', 2 * scale);

            if (blobs && this.irDebug) {
                const small = this.textScale('small');
                points.forEach((p, i) => this._drawBlobInfo(ctx, px(p.x), py(p.y), scale, small, blobs[i], this.blobLevels[i]));
            }
        }

        /**
         * One emitter from the blob data. Seen: filled circle, radius from the blob area, radial
         * gradient from the max brightness (centre) to the average brightness (edge), outline in
         * the emitter colour. Not seen: dashed circle in the emitter colour with a red X.
         */
        _drawEmitter(ctx, x, y, scale, color, blob, level) {
            ctx.save();
            if (blob && blob.seen && level) {
                const r = level.radius * scale;
                const gradient = ctx.createRadialGradient(x, y, 0, x, y, r);
                gradient.addColorStop(0, withAlpha(color, level.center));
                gradient.addColorStop(1, withAlpha(color, level.edge));
                ctx.beginPath();
                ctx.arc(x, y, r, 0, Math.PI * 2);
                ctx.fillStyle = gradient;
                ctx.fill();
                ctx.strokeStyle = color;
                ctx.lineWidth = Math.max(1, 2 * scale);
                ctx.stroke();
            } else {
                const r = IRTEST_EMITTER_RADIUS * scale;
                ctx.beginPath();
                ctx.arc(x, y, r, 0, Math.PI * 2);
                ctx.setLineDash([6 * scale, 5 * scale]);
                ctx.strokeStyle = color;
                ctx.lineWidth = 3 * scale;
                ctx.stroke();
                ctx.setLineDash([]);
                const d = r * 0.62;
                ctx.beginPath();
                ctx.moveTo(x - d, y - d);
                ctx.lineTo(x + d, y + d);
                ctx.moveTo(x + d, y - d);
                ctx.lineTo(x - d, y + d);
                ctx.lineCap = 'round';
                ctx.strokeStyle = IRTEST_NOT_SEEN_COLOR;
                ctx.lineWidth = 3 * scale;
                ctx.stroke();
            }
            ctx.restore();
        }

        /** ?irdebug: area and brightness (max/average) under the emitter circle. */
        _drawBlobInfo(ctx, x, y, scale, textScale, blob, level) {
            const text = blob && blob.seen ? [`A ${blob.area}`, `${blob.max}/${blob.avg}`] : ['--'];
            const radius = (level ? level.radius : IRTEST_EMITTER_RADIUS) * scale;
            const size = textSize(text, textScale);
            drawText(ctx, text, x - size.width / 2, y + radius + 4 * scale, textScale);
        }
    }

    FullscreenWindow.MODE_CALIBRATE = MODE_CALIBRATE;
    FullscreenWindow.MODE_ALIGNMENT = MODE_ALIGNMENT;
    FullscreenWindow.MODE_IRTEST = MODE_IRTEST;
    FullscreenWindow.STAGE = STAGE;
    FullscreenWindow.glyphPixels = glyphPixels;

    OF.FullscreenWindow = FullscreenWindow;
})(typeof globalThis !== 'undefined' ? globalThis : this);
