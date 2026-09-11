/*  OpenFIRE Web App - fullscreen windows: calibration, IR emitter alignment, IR camera test
    (Qt App: appcali.cpp). Drawn on a canvas with the App's 8x8 bitmap typeface.

        const win = new OF.FullscreenWindow('calibrate', { onExitRequest, onExit });
        win.open();                      // call it inside the click handler (fullscreen needs a user gesture)
        win.setStage(stage); win.setInfo(payload); win.drawTest(coords);
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

    const CALI_PREFIXES = ['   Top Offset: ', 'Bottom Offset: ', '  Left Offset: ', ' Right Offset: ', ' Top Left LED: ', 'Top Right LED: '];
    const CALI_KEYS = ['topOffset', 'bottomOffset', 'leftOffset', 'rightOffset', 'TLled', 'TRled'];

    const CROSSHAIR_SIZE = 61.44;
    const CROSSHAIR_SVG = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 61.44 61.44" width="61.44" height="61.44">' +
        '<g fill="none" stroke="#ff8758" stroke-linecap="square">' +
        '<circle cx="30.72" cy="30.72" r="24.42" stroke-width="2.4"/>' +
        '<path stroke-width="1.44" d="M0.72 30.72h11.16M30.72 0.73v11.16M60.71 30.72h-11.16M30.72 60.71v-11.16M26.94 30.72h7.58M30.72 26.94v7.58"/>' +
        '</g></svg>';

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

    /** Translated lines: the Qt texts are padded for a left-aligned block, here lines are centred. */
    function lines(...keys) {
        return keys.map((key) => (key ? OF.i18n.t(key) : '').replace(/\u2026/g, '...').trim());
    }

    function loadCrosshair() {
        if (!crosshairImage) {
            crosshairImage = new Image();
            crosshairImage.src = 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(CROSSHAIR_SVG);
        }
        return crosshairImage;
    }

    // ----- Window -------------------------------------------------------------------------

    class FullscreenWindow {
        constructor(mode, options = {}) {
            this.mode = mode;
            this.options = options;
            this.stage = STAGE.init;
            this.values = {};
            this.coords = null;
            this.mouse = null;
            this.closed = false;
            this.resetValues();
        }

        resetValues() {
            for (const key of CALI_KEYS) this.values[key] = -1;
            this.infoText = CALI_PREFIXES.slice();
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
            this.stage = stage;
            if (stage === STAGE.init) {
                this.resetValues();
                this.mouse = null;
            } else if (stage >= STAGE.top) {
                this.infoVisible = true;
            }
            this.render();
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
            this.infoText[type - 1] = CALI_PREFIXES[type - 1] + shown;
            if (type === 6) this.stage = STAGE.verify; // may arrive after the verify stage
            this.render();
        }

        /** sTestCoords: 12 int32 (TL, TR, BL, BR with the outside-FOV flag in bit 0 of X; mouse; D). */
        drawTest(payload) {
            if (this.closed || this.mode !== MODE_IRTEST) return;
            if (payload.length !== 48) return;
            const view = new DataView(payload.buffer, payload.byteOffset, 48);
            this.coords = Array.from({ length: 12 }, (_, i) => view.getInt32(i * 4, true));
            this.render();
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
                this._draw();
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
                stageText = lines(`Cali Step ${this.stage}:`);
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
                    stageText = lines('Verify New Calibration:');
                    header = lines('    Confirm that the bullseye     ', '   lines up with the gun sight.   ', 'If this calibration is acceptable,', '  confirm by pulling the trigger. ');
                    tutorial = lines("       If this target accuracy isn't desirable,      ", '          press either Button A or Button B          ',
                        '         to restart the calibration process.         ', '', '[You can also exit calibration without saving changes', '        by pressing Button C (if available).]        ');
                } else if (v.TLled === -1 || v.TRled === -1) {
                    stageText = lines('Verify New Calibration:');
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
            this._centered(ctx, stageText, stageTop, heading, tint);
            const headerTop = headerStageY === null ? stageTop + stageSize.height :
                headerStageY - textSize(lines('Cali Step 5:'), heading).height / 2 + textSize(lines('Cali Step 5:'), heading).height;
            this._centered(ctx, header, headerTop, heading, tint);
            if (tutorial.length) {
                const size = textSize(tutorial, sub);
                this._centered(ctx, tutorial, h * 0.8 - size.height / 2, sub, tint);
            }

            if (this.infoVisible) {
                let y = h / 2 - 24 * sub;
                for (const text of this.infoText) {
                    drawText(ctx, [text], 80 * sub, y, sub);
                    y += GLYPH * sub;
                }
            }

            // Crosshair position of each stage; while verifying it follows the gun's cursor.
            const positions = [[w / 2, h / 2], [w / 2, 0], [w / 2, h], [0, h / 2], [w, h / 2], [w / 2, h / 2]];
            let [x, y] = positions[this.stage] || positions[5];
            if (this.stage === STAGE.verify && this.mouse) ({ x, y } = this.mouse);
            const image = loadCrosshair();
            const size = CROSSHAIR_SIZE * this.textScale('crosshair');
            if (image.complete && image.naturalWidth) {
                ctx.imageSmoothingEnabled = true;
                ctx.drawImage(image, x - size / 2, y - size / 2, size, size);
                ctx.imageSmoothingEnabled = false;
            }
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
            points.forEach((p, i) => circle(offsetX, offsetY, scale, scale, p.x, p.y, colors[i], p.outside));
            circle(0, 0, scaleX, scaleY, c[8], c[9], '#a0a0a4', false);
            circle(offsetX, offsetY, scale, scale, c[10], c[11], '#ff0000', false);

            // The emitters box is added last to the Qt scene: drawn over the circles.
            const [tl, tr, bl, br] = points;
            this._polyline(ctx, [tl, tr, br, bl].map((p) => [px(p.x), py(p.y)]), '#a0a0a4', 2 * scale);
        }
    }

    FullscreenWindow.MODE_CALIBRATE = MODE_CALIBRATE;
    FullscreenWindow.MODE_ALIGNMENT = MODE_ALIGNMENT;
    FullscreenWindow.MODE_IRTEST = MODE_IRTEST;
    FullscreenWindow.STAGE = STAGE;
    FullscreenWindow.glyphPixels = glyphPixels;

    OF.FullscreenWindow = FullscreenWindow;
})(typeof globalThis !== 'undefined' ? globalThis : this);
