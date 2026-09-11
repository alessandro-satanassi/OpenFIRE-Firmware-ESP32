/*  OpenFIRE Web App - configuration state and editing rules.

    Pure logic (no DOM), ported from the Qt App main window (appmainwindow.cpp):
      orig      configuration as loaded from (or last saved to) the board
      cur       configuration being edited, sent by Save
      testPins  pins map at the last load/save: the Gun Tests tab describes the board as it is

    config.pins[function] = GPIO (-1 unmapped), the Qt inputsMap. When custom pins are off
    cur.pins holds the board's default layout, like the pin boxes of the Qt App.

        const state = new OF.AppState(OpenFIREshared);
        state.load(board, config);
        const messages = state.changePin(gpio, functionValue);   // I2C clash messages
        state.isDirty();
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const FIRECON_NAME = (n) => `FIRECon P${n}`;
    const NAME_MAX = 15;

    /** Ranges of the spin boxes of the Qt App (appmainwindow.ui), by settingsTypes_e name. */
    const SETTING_RANGES = {
        rumbleStrength: [0, 255],
        rumbleInterval: [0, 9999],
        solenoidOnLength: [0, 99],
        solenoidOffLength: [0, 99],
        solenoidHoldLength: [0, 9999],
        holdToPauseLength: [0, 9999],
        customLEDcount: [1, 100],
        customLEDstatic: [0, 3],
        tempWarning: [20, 99],
        tempShutdown: [25, 60],
    };

    class AppState {
        constructor(shared) {
            this.S = shared || OF.Boards.shared;
            this.E = this.S.boardInputs_e;
            this.B = this.S.boolTypes_e;
            this.T = this.S.settingsTypes_e;
            this.cap = this.S.pinCapabilities_e;
            this.codec = OF.ProtocolCodec ? new OF.ProtocolCodec(this.S) : null;
            this.reset();
        }

        reset() {
            this.loaded = false;
            this.forceDirty = false;      // restored edits: Save stays available
            this.board = null;
            this.orig = null;
            this.cur = null;
            this.testPins = [];
            this.altPresetIndex = -1;
        }

        clone(config) {
            return this.codec ? this.codec.cloneConfig(config) : JSON.parse(JSON.stringify(config));
        }

        // ----- Board data -----------------------------------------------------

        /** Board key of the data maps: its own entry, else the generic layout. */
        get layoutBoard() {
            const type = this.board && this.board.type;
            return type && this.S.boardsBoxPositions[type] ? type : 'generic';
        }

        get boxPositions() { return this.S.boardsBoxPositions[this.layoutBoard] || []; }
        /** Default layout: only boards listed in boardsPresetsMap have one (Qt BoxesUpdate). */
        get presetPins() {
            const type = this.board && this.board.type;
            return (type && this.S.boardsPresetsMap[type]) || [];
        }
        get pinCount() { return this.boxPositions.length; }

        get altPresets() {
            const type = this.board && this.board.type;
            return (type && (this.S.boardsAltPresets || {})[type]) || [];
        }

        /** Capability bits of every GPIO: board override, else architecture map. */
        get capabilities() {
            const maps = this.S.mcuCapableMaps;
            const type = this.board && this.board.type;
            return (type && maps[type]) || maps[this.board ? this.board.arch : ''] || maps[this.S.boardArchs[0]] || [];
        }

        get isRP() {
            return !!this.board && this.board.arch === this.S.boardArchs[this.S.boardArchs_e.boardRP];
        }

        get profileCount() { return this.cur ? this.cur.profiles.length : 0; }

        /** Title of the loaded gun (Qt PrettifyName). */
        prettyName(name, unnamed = 'Unnamed Device') {
            const label = name || unnamed;
            const names = this.S.boardNames || {};
            const type = this.board && this.board.type;
            return `${label} | ${names[type] || names.generic || type || ''}`;
        }

        /** 'i2c0' | 'i2c1' | 'any' | 'none': colour class of the GPIO label. */
        pinI2CClass(gpio) {
            const bits = this.capabilities[gpio] || 0;
            if (bits & this.cap.pinAnyI2C) return 'any';
            if (bits & this.cap.pinCanI2C) return (bits & this.cap.pinIsI2C1) ? 'i2c1' : 'i2c0';
            return 'none';
        }

        /** True when a function cannot be assigned to this GPIO (Qt SetComboBoxItemEnabled). */
        functionDisabled(gpio, value) {
            const E = this.E;
            const bits = this.capabilities[gpio] || 0;
            if (!(bits & this.cap.pinHasADC) && (value === E.analogX || value === E.analogY || value === E.tempPin))
                return true;
            if (!(bits & this.cap.pinAnyI2C)) {
                if (bits & this.cap.pinCanI2C) {
                    if (bits & this.cap.pinIsI2CSCL) return value === E.camSDA || value === E.periphSDA;
                    return value === E.camSCL || value === E.periphSCL;
                }
                return value === E.camSDA || value === E.camSCL || value === E.periphSDA || value === E.periphSCL;
            }
            return false;
        }

        // ----- Load / save -------------------------------------------------------

        load(board, config) {
            this.board = Object.assign({}, board);
            this.orig = this.clone(config);
            this.cur = this.clone(config);
            this.altPresetIndex = -1;
            this.forceDirty = false;
            this.loaded = true;
            this._normalise(this.cur);

            const pins = this.cur.pins;
            pins.fill(this.E.btnUnmapped);
            if (this.cur.toggles[this.B.customPins]) {
                for (let fn = 0; fn < this.orig.pins.length; ++fn) {
                    const gpio = this.orig.pins[fn];
                    if (gpio > this.E.btnUnmapped && gpio < this.pinCount && fn < this.E.boardInputsCount)
                        this._setPin(gpio, fn + 1, []);
                }
            } else {
                this._applyLayout(this.presetPins, []);
            }
            this.testPins = pins.slice();
        }

        /** Values the Qt widgets cannot show are corrected in the edited copy (Save becomes available). */
        _normalise(config) {
            const T = this.T;
            const B = this.B;
            const settings = config.settings;
            for (const [name, [min, max]] of Object.entries(SETTING_RANGES)) {
                const i = T[name];
                if (i === undefined || i >= settings.length) continue;
                settings[i] = Math.min(max, Math.max(min, settings[i] >>> 0));
            }
            if (settings[T.customLEDstatic] > settings[T.customLEDcount])
                settings[T.customLEDstatic] = settings[T.customLEDcount];
            const toggles = config.toggles;
            if (toggles[B.solenoid] && toggles[B.rumbleFF]) toggles[B.solenoid] = false;
            if (!(toggles[B.solenoid] || (toggles[B.rumble] && toggles[B.rumbleFF]))) toggles[B.autofire] = false;
        }

        /** Edits kept from a session that was lost before saving (see main.js). */
        snapshot() {
            return this.loaded ? { type: this.board.type, cur: this.clone(this.cur), orig: this.clone(this.orig), forceDirty: this.forceDirty } : null;
        }

        /** True when the snapshot has something the board just loaded does not. */
        snapshotDiffers(snapshot) {
            if (!this.loaded || !snapshot || snapshot.type !== this.board.type) return false;
            // Against the board now (edits not sent) or against the board then (calibration kept only in RAM).
            return !!snapshot.forceDirty || this._isDirtyBetween(snapshot.cur, this.orig) ||
                this._isDirtyBetween(snapshot.cur, snapshot.orig);
        }

        _isDirtyBetween(cur, orig) {
            const saved = [this.cur, this.orig, this.forceDirty];
            this.cur = cur;
            this.orig = orig;
            this.forceDirty = false;
            try {
                return this.isDirty();
            } finally {
                [this.cur, this.orig, this.forceDirty] = saved;
            }
        }

        /** Puts back the edits of a lost session; Save stays available until they are saved. */
        restoreSnapshot(snapshot) {
            if (!this.loaded || !snapshot || snapshot.type !== this.board.type) return false;
            const cur = this.clone(snapshot.cur);
            if (cur.pins.length !== this.cur.pins.length || cur.profiles.length !== this.cur.profiles.length) return false;
            this.cur = cur;
            this.altPresetIndex = -1;
            this.forceDirty = true;
            return true;
        }

        /** Qt on_confirmButton_clicked after a successful commit. */
        commitDone() {
            this.forceDirty = false;
            const orig = this.clone(this.cur);
            if (!orig.toggles[this.B.customPins]) orig.pins.fill(this.E.btnUnmapped);
            this.orig = orig;
            this.testPins = this.cur.pins.slice();
        }

        get cameraChanged() {
            const i = this.T.cameraModel;
            return !!this.cur && this.cur.settings[i] !== this.orig.settings[i];
        }

        /** Qt DiffUpdate: anything to save. */
        isDirty() {
            if (!this.loaded) return false;
            if (this.forceDirty) return true;
            const a = this.cur;
            const b = this.orig;
            const same = (x, y) => x.length === y.length && x.every((v, i) => v === y[i]);
            if (!same(a.toggles, b.toggles)) return true;
            if (a.toggles[this.B.customPins] && !same(a.pins, b.pins)) return true;
            if (!same(a.settings, b.settings)) return true;
            if (!a.buttons.every((button, i) => same(button, b.buttons[i]))) return true;
            if (a.tinyUSB.id !== b.tinyUSB.id || a.tinyUSB.name !== b.tinyUSB.name) return true;
            if (a.selectedProfile !== b.selectedProfile) return true;
            return a.profiles.some((profile, i) =>
                Object.keys(profile).some((key) => profile[key] !== (b.profiles[i] || {})[key]));
        }

        /** Qt PixelsDiff: NeoPixel settings changed (a power cycle may be needed). */
        pixelsChanged() {
            if (!this.loaded) return false;
            const T = this.T;
            return [T.customLEDcount, T.customLEDstatic, T.customLEDcolor1, T.customLEDcolor2, T.customLEDcolor3]
                .some((i) => this.cur.settings[i] !== this.orig.settings[i]);
        }

        // ----- Pins (Qt pinBoxes_currentIndexChanged / BoxesUpdate) ----------------

        pinOf(fn) { return this.cur ? this.cur.pins[fn] : -1; }

        functionAt(gpio) {
            return this.cur ? this.cur.pins.indexOf(gpio) : -1;
        }

        pinMapped(fn) { return this.pinOf(fn) > this.E.btnUnmapped; }

        get customPins() { return !!this.cur && !!this.cur.toggles[this.B.customPins]; }

        /** A pin box changed: index = function value + 1 (0 = unmapped). Returns status messages. */
        changePin(gpio, index) {
            const messages = [];
            this.altPresetIndex = -1;
            this._setPin(gpio, index, messages);
            this._pinSideEffects();
            return messages;
        }

        _setPin(gpio, index, messages) {
            const E = this.E;
            const pins = this.cur.pins;
            const prev = pins.indexOf(gpio) + 1;

            if (index <= 0) {
                if (prev > 0) pins[prev - 1] = E.btnUnmapped;
                return;
            }
            if (prev === index) return;

            const fn = index - 1;
            if (prev > 0) pins[prev - 1] = E.btnUnmapped;

            // Remove the pin that this function was mapped to.
            if (pins[fn] > E.btnUnmapped && pins[fn] !== gpio)
                this._setPin(pins[fn], 0, messages);

            const caps = this.capabilities;
            if (!((caps[gpio] || 0) & this.cap.pinAnyI2C)) {
                const channel = (pin) => (caps[pin] || 0) & this.cap.pinIsI2C1;
                const unmap = (pin) => this._setPin(pin, 0, messages);

                if (fn === E.camSDA || fn === E.periphSDA) {
                    const pair = pins[fn + 1];
                    if (pair > E.btnUnmapped && channel(gpio) !== channel(pair)) {
                        if (fn === E.camSDA) {
                            unmap(pins[E.camSCL]);
                            messages.push("Camera pins are not on the same I2C channel! Please check camera pins' mappings.");
                        } else {
                            unmap(pins[E.periphSCL]);
                            messages.push("Peripheral pins are not on the same I2C channel! Please check peripheral pins' mappings.");
                        }
                    }
                    if (fn === E.camSDA && pins[E.periphSDA] > E.btnUnmapped && channel(gpio) === channel(pins[E.periphSDA])) {
                        unmap(pins[E.periphSDA]);
                        messages.push('Camera and Peripheral Data pins clashed! Please remap Peripheral SDA.');
                    } else if (fn === E.periphSDA && pins[E.camSDA] > E.btnUnmapped && channel(gpio) === channel(pins[E.camSDA])) {
                        unmap(pins[E.camSDA]);
                        messages.push('Camera and Peripheral Data pins clashed! Please remap Camera SDA.');
                    }
                } else if (fn === E.camSCL || fn === E.periphSCL) {
                    const pair = pins[fn - 1];
                    if (pair > E.btnUnmapped && channel(gpio) !== channel(pair)) {
                        if (fn === E.camSCL) {
                            unmap(pins[E.camSDA]);
                            messages.push("Camera pins are not on the same I2C channel! Please check camera pins' mappings.");
                        } else {
                            unmap(pins[E.periphSDA]);
                            messages.push("Peripheral pins are not on the same I2C channel! Please check peripheral pins' mappings.");
                        }
                    }
                    if (fn === E.camSCL && pins[E.periphSCL] > E.btnUnmapped && channel(gpio) === channel(pins[E.periphSCL])) {
                        unmap(pins[E.periphSCL]);
                        messages.push('Camera and Peripheral Data pins clashed! Please remap Peripheral SCL.');
                    } else if (fn === E.periphSCL && pins[E.camSCL] > E.btnUnmapped && channel(gpio) === channel(pins[E.camSCL])) {
                        unmap(pins[E.camSCL]);
                        messages.push('Camera and Peripheral Data pins clashed! Please remap Camera SCL.');
                    }
                }
            }

            pins[fn] = gpio;
        }

        /** Clears every pin, then applies a layout given as GPIO -> function value. */
        _applyLayout(layout, messages) {
            const count = this.pinCount;
            for (let gpio = 0; gpio < count; ++gpio) this._setPin(gpio, 0, messages);
            for (let gpio = 0; gpio < count && gpio < layout.length; ++gpio)
                this._setPin(gpio, layout[gpio] + 1, messages);
        }

        /** Feedback outputs without a pin cannot stay enabled. */
        _pinSideEffects() {
            const E = this.E;
            const B = this.B;
            if (!this.pinMapped(E.rumblePin)) {
                this.setToggle(B.rumble, false);
                this.setToggle(B.rumbleFF, false);
            }
            if (!this.pinMapped(E.solenoidPin)) this.setToggle(B.solenoid, false);
        }

        /** Qt on_customPinsEnabled_stateChanged + BoxesUpdate. */
        setCustomPins(on) {
            const B = this.B;
            if (!!this.cur.toggles[B.customPins] === !!on) return [];
            const messages = [];
            this.cur.toggles[B.customPins] = !!on;
            if (on) {
                if (this.orig.toggles[B.customPins]) {
                    for (let gpio = 0; gpio < this.pinCount; ++gpio) this._setPin(gpio, 0, messages);
                    for (let fn = 0; fn < this.orig.pins.length; ++fn) {
                        const gpio = this.orig.pins[fn];
                        if (gpio > this.E.btnUnmapped && gpio < this.pinCount && fn < this.E.boardInputsCount)
                            this._setPin(gpio, fn + 1, messages);
                    }
                }
                // Otherwise the default layout already shown becomes the custom one.
            } else {
                this.altPresetIndex = -1;
                this._applyLayout(this.presetPins, messages);
            }
            this._pinSideEffects();
            return messages;
        }

        /** Qt on_presetsBox_currentIndexChanged. */
        applyAltPreset(index) {
            const preset = this.altPresets[index];
            if (!preset) return [];
            const messages = [];
            this.cur.toggles[this.B.customPins] = true;
            this._applyLayout(preset.pin, messages);
            this._pinSideEffects();
            this.altPresetIndex = index;
            return messages;
        }

        /** Custom layout file (.ofl): board type + '\n', then (function name, NUL, GPIO byte) records. */
        exportLayout() {
            const bytes = [];
            const push = (text) => { for (const c of text) bytes.push(c.charCodeAt(0) & 0xFF); };
            push(this.board.type + '\n');
            for (const [name, value] of Object.entries(this.S.boardInputs_Strings)) {
                if (value > this.E.btnUnmapped && this.pinMapped(value)) {
                    push(name);
                    bytes.push(0, this.pinOf(value) & 0xFF);
                }
            }
            return Uint8Array.from(bytes);
        }

        /** Returns 'ok' | 'mismatch'. */
        importLayout(data) {
            const bytes = data instanceof Uint8Array ? data : new Uint8Array(data);
            const newline = bytes.indexOf(0x0A);
            const firstLine = String.fromCharCode(...bytes.subarray(0, newline < 0 ? bytes.length : newline)).trim();
            if (newline < 0 || firstLine !== this.board.type) return 'mismatch';

            const layout = new Array(this.pinCount).fill(this.E.btnUnmapped);
            let offset = newline + 1;
            while (offset < bytes.length) {
                const end = bytes.indexOf(0, offset);
                if (end < 0 || end + 1 >= bytes.length) break;
                const name = String.fromCharCode(...bytes.subarray(offset, end));
                const gpio = bytes[end + 1];
                offset = end + 2;
                const value = this.S.boardInputs_Strings[name];
                if (value !== undefined && gpio < layout.length) layout[gpio] = value;
            }

            const messages = [];
            this.cur.toggles[this.B.customPins] = true;
            this.altPresetIndex = -1;
            this._applyLayout(layout, messages);
            this._pinSideEffects();
            // Qt: "recover previous toggles that may have been turned off in the process" (saved values).
            if (this.pinMapped(this.E.solenoidPin)) this.setToggle(this.B.solenoid, this.orig.toggles[this.B.solenoid]);
            if (this.pinMapped(this.E.rumblePin)) this.setToggle(this.B.rumble, this.orig.toggles[this.B.rumble]);
            return 'ok';
        }

        // ----- Toggles and settings (Qt on_*Toggle_stateChanged) ------------------

        toggle(index) { return !!this.cur.toggles[index]; }

        setToggle(index, value) {
            const B = this.B;
            const toggles = this.cur.toggles;
            value = !!value;
            if (index === B.customPins) return this.setCustomPins(value);
            if (!!toggles[index] === value) return;
            toggles[index] = value;

            switch (index) {
            case B.rumble:
                if (!value) this.setToggle(B.rumbleFF, false);
                this._autofireCheck();
                break;
            case B.solenoid:
                if (value) this.setToggle(B.rumbleFF, false);
                this._autofireCheck();
                break;
            case B.rumbleFF:
                if (value) this.setToggle(B.solenoid, false);
                this._autofireCheck();
                break;
            default:
                break;
            }
        }

        /** Autofire needs the solenoid, or rumble used as force feedback. */
        get autofireAvailable() {
            const B = this.B;
            return this.toggle(B.solenoid) || (this.toggle(B.rumble) && this.toggle(B.rumbleFF));
        }

        _autofireCheck() {
            if (!this.autofireAvailable) this.cur.toggles[this.B.autofire] = false;
        }

        setting(index) { return this.cur.settings[index]; }

        setSetting(index, value) {
            const T = this.T;
            const settings = this.cur.settings;
            value = Math.max(0, Math.trunc(Number(value) || 0)) >>> 0;
            if (index === T.customLEDstatic && value > settings[T.customLEDcount]) value = settings[T.customLEDcount];
            settings[index] = value;
            if (index === T.customLEDcount && value < settings[T.customLEDstatic]) settings[T.customLEDstatic] = value;
        }

        // ----- Button mapping (Qt btnFuncTypeBox / btnFuncBox) ----------------------

        buttonType(button, slot) { return this.cur.buttons[button][slot * 2]; }
        buttonValue(button, slot) { return this.cur.buttons[button][slot * 2 + 1]; }

        /** Output type of a slot changed: the saved output is kept when going back to the saved type. */
        setButtonType(button, slot, type) {
            const row = this.cur.buttons[button];
            if (row[slot * 2] === type) return;
            row[slot * 2] = type;
            const saved = this.orig.buttons[button];
            const list = OF.Maps.outputs[type] || [];
            if (saved[slot * 2] === type && OF.Maps.outputIndex(type, saved[slot * 2 + 1]) >= 0)
                row[slot * 2 + 1] = saved[slot * 2 + 1];
            else if (list.length)
                row[slot * 2 + 1] = list[0].value;
        }

        setButtonValue(button, slot, value) {
            this.cur.buttons[button][slot * 2 + 1] = value;
        }

        // ----- TinyUSB identifier ------------------------------------------------

        /** 1..4 for the simple player presets, else 0. */
        get usbPreset() {
            const id = this.cur.tinyUSB.id;
            return id >= 1 && id <= 4 ? id : 0;
        }

        setUsbPreset(n) {
            this.cur.tinyUSB.id = n;
            this.cur.tinyUSB.name = FIRECON_NAME(n);
        }

        setUsbId(id) {
            const value = Math.min(0xFFFF, Math.max(0, Math.trunc(Number(id) || 0)));
            this.cur.tinyUSB.id = value;
            if (!this.cur.tinyUSB.name && value >= 1 && value <= 4) this.setUsbPreset(value);
        }

        /** Returns false (nothing changed) when the name has characters outside Latin-1. */
        setUsbName(text) {
            if (!AppState.isLatin1(text)) return false;
            this.cur.tinyUSB.name = String(text).slice(0, NAME_MAX);
            return true;
        }

        static isLatin1(text) {
            return ![...String(text)].some((c) => c.codePointAt(0) > 255);
        }

        // ----- Calibration profiles -------------------------------------------------

        setProfileField(profile, key, value) {
            this.cur.profiles[profile][key] = value;
        }

        /** Returns false (nothing changed) for an empty name or characters outside Latin-1. */
        renameProfile(profile, name) {
            if (!name || !AppState.isLatin1(name)) return false;
            this.cur.profiles[profile].name = String(name).slice(0, NAME_MAX);
            return true;
        }

        /** sCurrentProf from the board. */
        setSelectedProfile(profile) {
            if (profile < this.profileCount) this.cur.selectedProfile = profile;
        }

        /** Values returned by a calibration (Qt CaliWindowExiting). Returns 'ok' | 'malformed' | 'cancelled'. */
        applyCalibration(profile, values) {
            const keys = ['topOffset', 'bottomOffset', 'leftOffset', 'rightOffset', 'TLled', 'TRled'];
            if (keys.every((key) => values[key] === -1)) return 'cancelled';
            const target = this.cur.profiles[profile];
            for (const key of keys) target[key] = values[key];
            return keys.every((key) => values[key] >= -32768 && values[key] <= 32768) ? 'ok' : 'malformed';
        }

        /** The profile's IR layout differs from the saved one (the IR test uses the saved settings). */
        layoutChangedUnsaved() {
            const i = this.cur.selectedProfile;
            const a = this.cur.profiles[i];
            const b = this.orig.profiles[i];
            return !!a && !!b && a.layoutType !== b.layoutType;
        }
    }

    AppState.NAME_MAX = NAME_MAX;
    AppState.SETTING_RANGES = SETTING_RANGES;
    OF.AppState = AppState;
})(typeof globalThis !== 'undefined' ? globalThis : this);
