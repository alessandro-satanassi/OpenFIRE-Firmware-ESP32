/*  OpenFIRE Web App tests: in-memory model of the firmware side of the App protocol.

    Mirrors lightgun/src/OpenFIREserial.cpp (OF_Serial::AppSerial*), the docked
    loop of main.cpp (ExecGunModeDocked / ExecRunModeProcessing) and the desktop
    calibration flow of OpenFIREcommon.cpp (FW_Common::ExecCalMode). Data is kept
    in raw byte arrays so that the firmware's memcpy semantics are reproduced.

    Node only: never packed into the web app.
*/
'use strict';

const START_1 = 0xA5;
const START_2 = 0x5A;
const MAX_PAYLOAD = 200;
const OVERHEAD = 7;
const MAX_FRAME = MAX_PAYLOAD + OVERHEAD;
const FRAME_TIMEOUT = 250;
const ACK_TIMEOUT = 500;
const MAX_RETRIES = 3;
const ERR_COMMIT_RETRY = 0x82;

const T_REQUEST = 0x00, T_RESPONSE = 0x01, T_EVENT = 0x02, T_ACK = 0x03, T_MASK = 0x03, F_FINAL = 0x80;

const PROFILE_COUNT = 4;
const PROFILE_SIZE = 68;
const BUTTON_COUNT = 14;           // ButtonCount in OpenFIREprefs.h
const USB_SIZE = 18;

const tick = (ms = 1) => new Promise((resolve) => setTimeout(resolve, ms));

// Independent CRC implementation (same algorithm as OF_Serial::AppSerialCRC8).
function crc8(bytes) {
    let crc = 0;
    for (const b of bytes) {
        crc ^= b;
        for (let i = 0; i < 8; ++i)
            crc = (crc & 0x80) ? ((crc << 1) ^ 0x9B) & 0xFF : (crc << 1) & 0xFF;
    }
    return crc;
}

// Deterministic PRNG for reproducible fault injection.
function mulberry32(seed) {
    return function () {
        seed |= 0; seed = seed + 0x6D2B79F5 | 0;
        let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
        t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
        return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
}

// ===== In-memory link between the web transport and the mock ===============

class MemoryLink {
    constructor(options = {}) {
        this.latency = options.latency ?? 2;
        this.random = mulberry32(options.seed ?? 1234);
        this.dropRate = options.dropRate ?? 0;           // per write chunk, both directions
        this.corruptRate = options.corruptRate ?? 0;     // per write chunk, flips one byte
        this.maxChunk = options.maxChunk ?? 0;           // split writes into random chunk sizes
        this.toFirmware = null;
        this.toApp = null;
        this.filterToFirmware = null;                    // (bytes) => bytes | null
        this.filterToApp = null;
        this.connected = true;
        this.stats = { dropped: 0, corrupted: 0 };
    }

    _impair(bytes, deliver) {
        if (!this.connected)
            return;
        if (this.dropRate && this.random() < this.dropRate) {
            this.stats.dropped++;
            return;
        }
        let data = Uint8Array.from(bytes);
        if (this.corruptRate && this.random() < this.corruptRate && data.length) {
            data[Math.floor(this.random() * data.length)] ^= 0x10;
            this.stats.corrupted++;
        }
        const chunks = [];
        if (this.maxChunk > 0) {
            for (let i = 0; i < data.length;) {
                const size = 1 + Math.floor(this.random() * this.maxChunk);
                chunks.push(data.subarray(i, i + size));
                i += size;
            }
        } else {
            chunks.push(data);
        }
        for (const chunk of chunks)
            setTimeout(() => { if (this.connected) deliver(chunk); }, this.latency);
    }

    appWrite(bytes) {
        const data = this.filterToFirmware ? this.filterToFirmware(Uint8Array.from(bytes)) : bytes;
        if (data)
            this._impair(data, (chunk) => this.toFirmware && this.toFirmware(chunk));
    }

    firmwareWrite(bytes) {
        const data = this.filterToApp ? this.filterToApp(Uint8Array.from(bytes)) : bytes;
        if (data)
            this._impair(data, (chunk) => this.toApp && this.toApp(chunk));
    }
}

/// Web-side transport implementing the interface of js/core/transport.js.
function createMemoryTransport(OF, link, kind = 'webserial') {
    class MemoryTransport extends OF.Transport {
        constructor() {
            super(kind);
            this.writes = [];
        }
        async open() {
            if (!link.connected)
                throw OF.TransportError('open_failed', 'link down');
            link.toApp = (bytes) => this._emitData(bytes);
            this._closeNotified = false;
            this.isOpen = true;
        }
        async _write(bytes) {
            this.writes.push(Uint8Array.from(bytes));
            link.appWrite(bytes);
            return true;
        }
        async close() {
            await this._writeChain;
            this.isOpen = false;
            link.toApp = null;
            this._emitClose('closed');
        }
        simulateLoss() {
            link.toApp = null;
            this._emitClose('device_lost');
        }
    }
    return new MemoryTransport();
}

// ===== Mock firmware ========================================================

class MockFirmware {
    /// `link` is the serial link (USB / wireless dongle). `options.webLink` adds the
    /// WebSocket of the web configuration mode (OF_WebConfigModeActive).
    constructor(shared, link, options = {}) {
        this.S = shared;
        this.C = shared.serialCmdTypes_e;
        this.links = { serial: link, web: options.webLink || null };
        this.webMode = !!options.webLink;
        this.sessionLink = 'serial';        // WebAppSerial::link
        this.boardType = options.boardType || 'esp32-s3-devkitc-1';
        this.version = options.version || '6.2-abcdef0';
        this.log = [];

        this.rxQueues = { serial: [], web: [] };
        this.rawDock = { serial: { state: 0, timestamp: 0 }, web: { state: 0, timestamp: 0 } };
        link.toFirmware = (bytes) => { for (const b of bytes) this.rxQueues.serial.push(b); };
        if (this.links.web)
            this.links.web.toFirmware = (bytes) => { for (const b of bytes) this.rxQueues.web.push(b); };

        this.running = false;
        this.failNextSave = false;
        this.camNotAvailable = !!options.camNotAvailable;
        this.pendingEvents = [];
        this.trigger = false;
        this.buttonA = false;
        this.rebootedToBootloader = false;
        this.flashCleared = false;
        this.mouse = { x: 1200, y: 90 };
        this.clientLost = false;          // WebApp_TakeClientLost() source (WebSocket close)
        this.sessionCounter = 0;

        this.resetPrefs();
        this.gunMode = 'run';
        this.runMode = 'normal';
        this._resetSession();
        this.appSerialPinsBeforeCommitValid = false;
        this.pinsBeforeCommit = new Int8Array(this.S.boardInputs_e.boardInputsCount);
        this.dockedSaving = false;
    }

    // ----- Preferences (OF_Prefs defaults) --------------------------------

    resetPrefs() {
        const S = this.S;
        this.toggles = Uint8Array.from([0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0]);
        this.settings = new Uint8Array(S.settingsTypes_e.settingsTypesCount * 4);
        const defaults = [255, 150, 20, 80, 500, 2500, 1, 0, 0xFF0000, 0x00FF00, 0x0000FF, 38, 45, 0, 1];
        const sv = new DataView(this.settings.buffer);
        defaults.forEach((value, i) => sv.setUint32(i * 4, value, true));

        this.pins = new Int8Array(S.boardInputs_e.boardInputsCount).fill(-1);

        // ButtonDesc defaults (report type/code on, off, gamepad).
        this.buttonDesc = new Uint8Array(BUTTON_COUNT * 6);
        const desc = [
            [0, 1, 0, 1, 2, 9], [0, 2, 0, 2, 2, 8], [0, 4, 0, 4, 2, 4], [1, 0xFF, 1, 0xFF, 2, 11],
            [1, 0xFE, 1, 0xFE, 2, 10], [0, 8, 0, 8, 2, 0], [2, 15, 2, 15, 2, 15], [2, 16, 2, 16, 2, 16],
            [2, 17, 2, 17, 2, 17], [2, 18, 2, 18, 2, 18], [0, 8, 0, 8, 2, 3], [0, 16, 0, 16, 2, 1],
            [0, 2, 0, 2, 2, 8], [3, 0, 3, 0, 3, 0]
        ];
        desc.forEach((d, i) => this.buttonDesc.set(d, i * 6));
        this.backupButtonDesc = Uint8Array.from(this.buttonDesc);

        this.profiles = [];
        const names = ['Profile A', 'Profile B', 'Profile Start', 'Profile Select'];
        const colors = [0xFF0000, 0x00FF00, 0x0000FF, 0xFF00FF];
        for (let i = 0; i < PROFILE_COUNT; ++i) {
            const raw = new Uint8Array(PROFILE_SIZE);
            const v = new DataView(raw.buffer);
            v.setInt32(0, 10 + i, true); v.setInt32(4, 20 + i, true); v.setInt32(8, -30 - i, true); v.setInt32(12, 40 + i, true);
            v.setFloat32(16, 1023.5, true); v.setFloat32(20, 6656.25, true); v.setFloat32(24, 512, true); v.setFloat32(28, 384, true);
            v.setUint32(32, 0, true); v.setUint32(36, 0, true); v.setUint32(40, 0, true); v.setUint32(44, 0, true);
            v.setUint32(48, colors[i], true);
            for (let c = 0; c < names[i].length; ++c) raw[52 + c] = names[i].charCodeAt(c);
            this.profiles.push(raw);
        }
        this.currentProfile = 0;

        this.usb = new Uint8Array(USB_SIZE);
        new DataView(this.usb.buffer).setUint16(0, 1, true);
        'FIRECon P1'.split('').forEach((ch, i) => { this.usb[2 + i] = ch.charCodeAt(0); });
    }

    loadPresets() {
        this.pins.fill(-1);
        const preset = this.S.boardsPresetsMap[this.boardType];
        if (preset) {
            preset.forEach((func, gpio) => { if (func > -1) this.pins[func] = gpio; });
        }
        this.backupButtonDesc = Uint8Array.from(this.buttonDesc);
    }

    // ----- Byte stream -----------------------------------------------------

    // Stream of the current App session (WebAppSerial::ops).
    available() { return this.rxQueues[this.sessionLink].length; }
    read() { const q = this.rxQueues[this.sessionLink]; return q.length ? q.shift() : -1; }

    // Direct access to one link.
    linkAvailable(name) { return name === 'web' && !this.webMode ? 0 : this.rxQueues[name].length; }
    linkRead(name) { const q = this.rxQueues[name]; return q.length ? q.shift() : -1; }
    otherLink() { return this.sessionLink === 'serial' ? 'web' : 'serial'; }

    /// OF_Serial::AppSerialInputPending()
    inputPending() {
        return this.available() > 0 || this.linkAvailable('serial') > 0 || this.linkAvailable('web') > 0 ||
               (this.webMode && this.clientLost);
    }

    /// Client lost is meaningful only for a WebSocket session (AppSerialTakeLinkLost).
    takeLinkLost() {
        return this.takeClientLost() && this.sessionLink === 'web';
    }

    enterDocked(name) {
        this.sessionLink = name;
        this.sessionBegin();
        this.gunMode = 'docked';
    }

    /// OF_Serial::AppSerialPollDock(link): raw dock on one link while no session is
    /// active; bytes of a link that does not own the active session are discarded.
    pollDock(name) {
        if (name === 'web' && !this.webMode)
            return false;
        const dock = this.rawDock[name];
        if (this.sessionActive) {
            if (name !== this.sessionLink) {
                while (this.linkAvailable(name)) this.linkRead(name);
                dock.state = 0;
            }
            return false;
        }
        if (dock.state !== 0 && Date.now() - dock.timestamp > FRAME_TIMEOUT)
            dock.state = 0;
        while (this.linkAvailable(name)) {
            const incoming = this.linkRead(name);
            if (dock.state !== 0 && incoming === this.C.sDock2) {
                dock.state = 0;
                this.enterDocked(name);
                return true;
            }
            dock.state = incoming === this.C.sDock1 ? 1 : 0;
            if (dock.state !== 0)
                dock.timestamp = Date.now();
        }
        return false;
    }

    writeFrame(typeFlags, command, sequence, payload) {
        const length = payload ? payload.length : 0;
        if (length > MAX_PAYLOAD)
            return false;
        const frame = new Uint8Array(OVERHEAD + length);
        frame[0] = START_1; frame[1] = START_2;
        frame[2] = typeFlags; frame[3] = command; frame[4] = sequence; frame[5] = length;
        if (length) frame.set(payload, 6);
        frame[6 + length] = crc8(frame.subarray(2, 6 + length));
        this.log.push({ dir: 'fw->app', type: typeFlags & T_MASK, final: !!(typeFlags & F_FINAL), command, sequence, length });
        this.links[this.sessionLink].firmwareWrite(frame);
        return true;
    }

    _resetSession() {
        this.sessionActive = false;
        this.rxBuf = [];
        this.rxTimestamp = 0;
        this.txSequence = 0;
        if (this.rawDock) { this.rawDock.serial.state = 0; this.rawDock.web.state = 0; }
        this.lastRxValid = false;
        this.lastRx = null;
        this.commitActive = false;
        this.commitFailed = false;
        this.calibrationCancel = false;
        this.waitingForAck = false;
        this.dispatching = false;
        this.processingDeferred = false;
        this.deferredValid = false;
        this.deferredFrame = null;
        this.rxClean = true;
        this.redockRequested = false;
        this.sessionConfirmed = false;
    }

    sessionBegin() {
        this._resetSession();
        this.sessionActive = true;
        this.dockedSaving = this.appSerialPinsBeforeCommitValid;
        ++this.sessionCounter;
    }

    takeClientLost() {
        const lost = this.clientLost;
        this.clientLost = false;
        return lost;
    }

    calibrating() {
        return this.gunMode === 'calibration' || this.gunMode === 'verification';
    }

    // OF_Serial::AppSerialAbandonSession()
    abandonSession() {
        this.redockRequested = false;
        if (!this.sessionActive)
            return;
        if (this.calibrating()) {
            this.calibrationCancel = true;
            this.sessionActive = false;
            return;
        }
        this.commitActive = false;
        this.commitFailed = false;
        this.sessionEnd();
        if (!this.appSerialPinsBeforeCommitValid)
            this.restoreRunState();
    }

    sessionEnd() {
        this.sessionActive = false;
        this.rxBuf = [];
        if (this.rawDock) { this.rawDock.serial.state = 0; this.rawDock.web.state = 0; }
        this.lastRxValid = false;
        this.commitActive = false;
        this.commitFailed = false;
        this.calibrationCancel = false;
        this.dockedSaving = this.appSerialPinsBeforeCommitValid;
        this.waitingForAck = false;
        this.deferredValid = false;
        this.rxClean = true;
        this.redockRequested = false;
    }

    nextSequence() {
        this.txSequence = (this.txSequence + 1) & 0xFF;
        if (this.txSequence === 0) this.txSequence = 1;
        return this.txSequence;
    }

    readFrame() {
        if (this.rxBuf.length > 0 && Date.now() - this.rxTimestamp > FRAME_TIMEOUT) {
            this.rxBuf = [];
            this.rxClean = true;
        }

        for (;;) {
            while (this.rxBuf.length >= 2) {
                let start = 0;
                while (start + 1 < this.rxBuf.length && (this.rxBuf[start] !== START_1 || this.rxBuf[start + 1] !== START_2))
                    ++start;
                if (start + 1 >= this.rxBuf.length) {
                    if (this.sessionActive && this.sessionConfirmed && this.rxClean && this.rxBuf.length === 2 &&
                        this.rxBuf[0] === this.C.sDock1 && this.rxBuf[1] === this.C.sDock2) {
                        this.rxBuf = [];
                        this.redockRequested = true;
                        return null;
                    }
                    this.rxBuf = this.rxBuf[this.rxBuf.length - 1] === START_1 ? [START_1] : [];
                    this.rxClean = false;
                    break;
                }
                if (start > 0) {
                    this.rxBuf.splice(0, start);
                    this.rxClean = false;
                }
                if (this.rxBuf.length < 6)
                    break;

                const typeFlags = this.rxBuf[2];
                const type = typeFlags & T_MASK;
                const sequence = this.rxBuf[4];
                const length = this.rxBuf[5];
                const invalidFlags = (typeFlags & ~(T_MASK | F_FINAL) & 0xFF) !== 0;
                const invalidFinal = (typeFlags & F_FINAL) && (type !== T_RESPONSE || length !== 0);
                const invalidSeq = type === T_EVENT ? sequence !== 0 : sequence === 0;
                const invalidAck = type === T_ACK && length !== 0;
                if (invalidFlags || invalidFinal || invalidSeq || invalidAck || length > MAX_PAYLOAD) {
                    this.rxBuf.shift();
                    this.rxClean = false;
                    continue;
                }
                const frameLength = OVERHEAD + length;
                if (this.rxBuf.length < frameLength)
                    break;
                const crc = this.rxBuf[6 + length];
                if (crc !== crc8(this.rxBuf.slice(2, 6 + length))) {
                    this.rxBuf.shift();
                    this.rxClean = false;
                    continue;
                }
                const frame = {
                    typeFlags, command: this.rxBuf[3], sequence, length, crc,
                    payload: Uint8Array.from(this.rxBuf.slice(6, 6 + length))
                };
                this.rxBuf.splice(0, frameLength);
                this.rxClean = true;
                if (this.sessionActive)
                    this.sessionConfirmed = true;
                this.log.push({ dir: 'app->fw', type, command: frame.command, sequence, length });
                return frame;
            }

            if (!this.available())
                return null;
            const incoming = this.read();
            if (this.rxBuf.length >= MAX_FRAME) {
                this.rxBuf.shift();
                this.rxClean = false;
            }
            this.rxBuf.push(incoming);
            this.rxTimestamp = Date.now();
        }
    }

    requestMatchesLast(frame) {
        return this.lastRxValid && this.lastRx.command === frame.command && this.lastRx.sequence === frame.sequence &&
               this.lastRx.length === frame.length && this.lastRx.crc === frame.crc;
    }

    sendAck(command, sequence) {
        this.writeFrame(T_ACK, command, sequence, null);
    }

    async sendReliable(typeFlags, command, payload, sequence = 0) {
        if (!this.sessionActive || (payload && payload.length > MAX_PAYLOAD))
            return false;
        if (sequence === 0)
            sequence = this.nextSequence();

        let success = false;
        this.waitingForAck = true;

        for (let attempt = 0; attempt <= MAX_RETRIES && this.sessionActive && !success; ++attempt) {
            if (!this.writeFrame(typeFlags, command, sequence, payload))
                break;
            const started = Date.now();
            while (this.sessionActive && Date.now() - started < ACK_TIMEOUT && !success) {
                let frame;
                while ((frame = this.readFrame())) {
                    const type = frame.typeFlags & T_MASK;
                    if (type === T_ACK) {
                        if (frame.command === command && frame.sequence === sequence) {
                            success = true;
                            break;
                        }
                        continue;
                    }
                    if (type === T_REQUEST) {
                        if (this.requestMatchesLast(frame)) {
                            this.sendAck(frame.command, frame.sequence);
                        } else {
                            if (!this.deferredValid) {
                                this.deferredFrame = frame;
                                this.deferredValid = true;
                            }
                            success = true;
                            break;
                        }
                    }
                }
                if (!success && (this.redockRequested || this.takeLinkLost())) {
                    this.abandonSession();
                    break;
                }
                if (!success)
                    await tick();
            }
        }

        this.waitingForAck = false;
        if (!this.dispatching)
            await this.processDeferredRequest();
        return success;
    }

    sendResponse(command, payload = null, final = false) {
        if (final && payload && payload.length)
            return Promise.resolve(false);
        return this.sendReliable(T_RESPONSE | (final ? F_FINAL : 0), command, payload);
    }

    sendEvent(command, payload = null) {
        if (!this.sessionActive)
            return false;
        return this.writeFrame(T_EVENT, command, 0, payload);
    }

    async sendError(error = 0) {
        if (error) return this.sendResponse(this.C.sError, Uint8Array.of(error));
        return this.sendResponse(this.C.sError);
    }

    async sendCommitError(frame) {
        return this.sendReliable(T_RESPONSE, this.C.sError, Uint8Array.of(ERR_COMMIT_RETRY, frame.command), frame.sequence);
    }

    takeCalibrationCancel() {
        const requested = this.calibrationCancel;
        this.calibrationCancel = false;
        return requested;
    }

    async processDeferredRequest() {
        if (this.waitingForAck || this.dispatching || this.processingDeferred)
            return;
        this.processingDeferred = true;
        while (this.sessionActive && this.deferredValid && !this.waitingForAck) {
            const frame = this.deferredFrame;
            this.deferredValid = false;
            await this.handleFrame(frame);
        }
        this.processingDeferred = false;
    }

    async handleFrame(frame) {
        if ((frame.typeFlags & T_MASK) !== T_REQUEST)
            return;
        const duplicate = this.requestMatchesLast(frame);
        if (!duplicate) {
            this.lastRxValid = true;
            this.lastRx = { command: frame.command, sequence: frame.sequence, length: frame.length, crc: frame.crc };
        }
        this.sendAck(frame.command, frame.sequence);
        if (!duplicate) {
            const already = this.dispatching;
            this.dispatching = true;
            await this.dispatchRequest(frame);
            this.dispatching = already;
        }
        await this.processDeferredRequest();
    }

    async serialProcessingDocked() {
        const calibrating = this.calibrating();
        if (this.takeLinkLost())
            this.abandonSession();

        if (this.sessionActive && this.deferredValid && !this.waitingForAck && calibrating) {
            const frame = this.deferredFrame;
            this.deferredValid = false;
            await this.handleFrame(frame);
            if (!this.sessionActive)
                return;
        }

        if (!this.sessionActive && !calibrating) {
            if (!this.pollDock('serial'))
                this.pollDock('web');
        } else if (this.sessionActive) {
            this.pollDock(this.otherLink());
        }

        if (this.sessionActive) {
            let frame;
            while (this.sessionActive && (frame = this.readFrame()))
                await this.handleFrame(frame);
        }

        if (this.redockRequested)
            this.abandonSession();
    }

    // ----- Records -----------------------------------------------------------

    fieldsReversed(map) {
        // std::unordered_map order is unspecified: use a different order from the web app.
        return Object.entries(map).reverse();
    }

    async sendRecords(command, data, fields, dataSize, profile = -1, sendFinal = true) {
        const profileData = profile >= 0;
        const buttonData = data === this.backupButtonDesc;
        const profCurrent = this.S.profSyncTypes_e.profCurrent;
        const profName = this.S.profSyncTypes_e.profName;

        for (const [name, index] of this.fieldsReversed(fields)) {
            if (index < 0) continue;
            if (profileData && index === profCurrent && profile !== 0) continue;

            const nameBytes = Uint8Array.from(name, (ch) => ch.charCodeAt(0));
            let body;
            if (profileData) {
                if (index === profCurrent) {
                    body = Uint8Array.of(1, this.currentProfile);
                } else {
                    const valueSize = index === profName ? 16 : dataSize;
                    body = new Uint8Array(2 + valueSize);
                    body[0] = valueSize;
                    body[1] = profile;
                    body.set(data.subarray(dataSize * index, dataSize * index + valueSize), 2);
                }
            } else {
                if (buttonData && index >= BUTTON_COUNT - 1) continue;
                body = new Uint8Array(1 + dataSize);
                body[0] = dataSize;
                body.set(new Uint8Array(data.buffer, data.byteOffset + dataSize * index, dataSize), 1);
            }
            const payload = new Uint8Array(nameBytes.length + 1 + body.length);
            payload.set(nameBytes, 0);
            payload.set(body, nameBytes.length + 1);
            if (!await this.sendResponse(command, payload))
                return false;
        }
        if (!sendFinal)
            return true;
        return this.sendResponse(command, null, true);
    }

    receiveRecord(payload, data, fields, dataSize, profileData = false) {
        if (!payload || payload.length < 2)
            return false;
        const buttonData = data === this.backupButtonDesc;
        let nameLength = 0;
        while (nameLength < payload.length && payload[nameLength] !== 0) ++nameLength;
        if (nameLength >= payload.length || nameLength + 1 >= payload.length)
            return false;
        const name = String.fromCharCode(...payload.subarray(0, nameLength));
        let pos = nameLength + 1;
        const valueSize = payload[pos++];
        if (!Object.prototype.hasOwnProperty.call(fields, name))
            return true;
        const index = fields[name];
        if (index < 0)
            return true;

        if (profileData) {
            if (index === this.S.profSyncTypes_e.profCurrent) {
                if (valueSize !== 1 || pos + valueSize !== payload.length) return false;
                if (payload[pos] >= PROFILE_COUNT) return false;
                this.currentProfile = payload[pos];
                return true;
            }
            if (pos >= payload.length) return false;
            const profNum = payload[pos++];
            if (profNum >= PROFILE_COUNT) return false;
            const expected = index === this.S.profSyncTypes_e.profName ? 16 : dataSize;
            if (valueSize !== expected || pos + valueSize !== payload.length) return false;
            this.profiles[profNum].set(payload.subarray(pos, pos + valueSize), dataSize * index);
            return true;
        }

        if (valueSize !== dataSize || pos + valueSize !== payload.length)
            return false;
        if (buttonData && index >= BUTTON_COUNT - 1)
            return true;
        new Uint8Array(data.buffer, data.byteOffset, data.byteLength).set(payload.subarray(pos, pos + valueSize), dataSize * index);
        return true;
    }

    // ----- Dispatch ----------------------------------------------------------

    async dispatchCommit(frame, fullDisconnect) {
        const C = this.C;
        const S = this.S;
        if (frame.command === C.sCommitStart) {
            if (!this.appSerialPinsBeforeCommitValid) {
                this.pinsBeforeCommit = Int8Array.from(this.pins);
                this.appSerialPinsBeforeCommitValid = true;
            }
            this.loadPresets();
            this.dockedSaving = true;
            this.commitActive = true;
            this.commitFailed = false;
            if (!await this.sendReliable(T_RESPONSE, frame.command, null, frame.sequence))
                this.commitActive = false;
            return true;
        }

        if (!this.commitActive || frame.command === C.sClearFlash || frame.command === C.sRebootToBootloader)
            return false;

        switch (frame.command) {
        case C.sSave: {
            const complete = !this.commitFailed;
            const saved = complete && !this.failNextSave;
            this.failNextSave = false;
            if (saved) {
                this.buttonDesc = Uint8Array.from(this.backupButtonDesc);
                this.appSerialPinsBeforeCommitValid = false;
                this.dockedSaving = false;
                this.savedCount = (this.savedCount || 0) + 1;
            }
            this.commitActive = false;
            const text = saved ? ' (Successfully saved to LittleFS Storage)' : ' (Failed to save to LittleFS Storage)';
            const payload = Uint8Array.from([saved ? 1 : 0, ...Array.from(text, (ch) => ch.charCodeAt(0))]);
            await this.sendReliable(T_RESPONSE, C.sSave, payload, frame.sequence);
            break;
        }
        case C.serialTerminator:
            this.commitActive = false;
            this.commitFailed = false;
            await this.sendResponse(C.serialTerminator, null, true);
            if (fullDisconnect)
                this.sessionEnd();
            break;
        case C.sCommitID:
            if (frame.length === USB_SIZE) this.usb = Uint8Array.from(frame.payload);
            else this.commitFailed = true;
            break;
        case C.sCommitToggles:
            if (!this.receiveRecord(frame.payload, this.toggles, S.boolTypes_Strings, 1)) this.commitFailed = true;
            break;
        case C.sCommitPins:
            if (!this.receiveRecord(frame.payload, this.pins, S.boardInputs_Strings, 1)) this.commitFailed = true;
            break;
        case C.sCommitSettings:
            if (!this.receiveRecord(frame.payload, this.settings, S.settingsTypes_Strings, 4)) this.commitFailed = true;
            break;
        case C.sCommitBtns:
            if (!this.receiveRecord(frame.payload, this.backupButtonDesc, S.boardInputs_Strings, 6)) this.commitFailed = true;
            break;
        case C.sCommitProfile:
            if (!this.receiveRecord(frame.payload, null, S.profSettingTypes_Strings, 4, true)) this.commitFailed = true;
            break;
        default:
            this.commitFailed = true;
            this.commitActive = false;
            await this.sendCommitError(frame);
            break;
        }
        return true;
    }

    restoreRunState() {
        this.gunMode = 'run';
        this.runMode = 'normal';
    }

    async dispatchRequest(frame) {
        const C = this.C;
        const S = this.S;
        const fullDisconnect = frame.command === C.serialTerminator && frame.length === 2 &&
                               frame.payload[0] === C.serialTerminator && frame.payload[1] === C.serialTerminator;

        if ((this.gunMode === 'calibration' || this.gunMode === 'verification') && frame.command !== C.serialTerminator) {
            this.calibrationCancel = true;
            await this.sendError();
            return;
        }

        if (await this.dispatchCommit(frame, fullDisconnect))
            return;

        if (this.appSerialPinsBeforeCommitValid) {
            switch (frame.command) {
            case C.sGetToggles: case C.sGetPins: case C.sGetSettings: case C.sGetBtns: case C.sGetProfile:
            case C.serialTerminator: case C.sClearFlash: case C.sRebootToBootloader:
                break;
            default:
                await this.sendCommitError(frame);
                return;
            }
        }

        switch (frame.command) {
        case C.serialTerminator:
            if (this.gunMode === 'calibration' || this.gunMode === 'verification') {
                this.calibrationCancel = true;
                if (fullDisconnect) {
                    await this.sendResponse(C.serialTerminator, null, true);
                    this.sessionActive = false;
                }
                break;
            }
            if (!fullDisconnect)
                break;
            await this.sendResponse(C.serialTerminator, null, true);
            this.sessionEnd();
            if (this.appSerialPinsBeforeCommitValid)
                break;
            this.restoreRunState();
            break;
        case C.sGetToggles:
            await this.sendRecords(frame.command, this.toggles, S.boolTypes_Strings, 1);
            break;
        case C.sGetPins:
            await this.sendRecords(frame.command, this.pins, S.boardInputs_Strings, 1);
            break;
        case C.sGetSettings:
            await this.sendRecords(frame.command, this.settings, S.settingsTypes_Strings, 4);
            break;
        case C.sGetBtns:
            await this.sendRecords(frame.command, this.backupButtonDesc, S.boardInputs_Strings, 6);
            break;
        case C.sGetProfile: {
            let sent = true;
            for (let prof = 0; prof < PROFILE_COUNT && sent; ++prof)
                sent = await this.sendRecords(frame.command, this.profiles[prof], S.profSettingTypes_Strings, 4, prof, false);
            if (sent && await this.sendResponse(frame.command, null, true) && this.appSerialPinsBeforeCommitValid)
                await this.sendCommitError(frame);
            break;
        }
        case C.sIRTest:
            if (frame.payload[0] && this.camNotAvailable) await this.sendError(C.sErrCam);
            else if (frame.payload[0]) this.runMode = 'processing';
            else if (this.runMode === 'processing') this.runMode = 'normal';
            break;
        case C.sCaliProfile: {
            const operation = frame.payload[0];
            const profile = frame.payload[1];
            if (profile < PROFILE_COUNT) this.currentProfile = profile;
            if (!await this.sendResponse(C.sCurrentProf, Uint8Array.of(this.currentProfile))) {
                await this.sendError();
                break;
            }
            if (operation === C.sCaliStart) {
                if (this.camNotAvailable) {
                    await this.sendError(C.sErrCam);
                } else {
                    const settings = frame.payload[2];
                    new DataView(this.profiles[this.currentProfile].buffer).setUint32(32, settings & 0x0F, true);
                    new DataView(this.profiles[this.currentProfile].buffer).setUint32(40, settings >> 4, true);
                    this.gunMode = 'calibration';
                    await this.execCalMode(true);
                    if (!this.sessionActive) {
                        this.sessionEnd();
                        this.restoreRunState();
                    }
                }
            }
            break;
        }
        case C.sTestSolenoid:
        case C.sTestRumble:
        case C.sTestLEDR:
        case C.sTestLEDG:
        case C.sTestLEDB:
            this.lastTest = frame.command;
            await tick(20);
            break;
        case C.sClearFlash:
            this.resetPrefs();
            this.flashCleared = true;
            this.sessionEnd();
            this.gunMode = 'run';
            break;
        case C.sRebootToBootloader:
            this.rebootedToBootloader = true;
            this.sessionEnd();
            break;
        default:
            break;
        }
    }

    // ----- Calibration (FW_Common::ExecCalMode, fromDesktop) -----------------

    profileView() { return new DataView(this.profiles[this.currentProfile].buffer); }

    async execCalMode(fromDesktop) {
        const C = this.C;
        const view = this.profileView();
        const backup = Uint8Array.from(this.profiles[this.currentProfile]);
        let calStage = 0;
        let communicationFailed = false;
        let top = 0, bottom = 0, left = 0, right = 0;
        const info = (type, bytes) => { const p = new Uint8Array(5); p[0] = type; p.set(bytes, 1); return p; };
        const i32 = (v) => { const b = new Uint8Array(4); new DataView(b.buffer).setInt32(0, v, true); return b; };

        view.setInt32(0, 0, true); view.setInt32(4, 0, true); view.setInt32(8, 0, true); view.setInt32(12, 0, true);
        this.gunMode = 'calibration';

        const fail = async () => { communicationFailed = true; return this._calCancelled(backup, communicationFailed); };

        if (!await this.sendResponse(C.sCaliStageUpd, Uint8Array.of(0)))
            return fail();

        while (this.gunMode === 'calibration') {
            await this.serialProcessingDocked();
            const desktopCancel = fromDesktop && this.takeCalibrationCancel();
            if (!this.sessionActive && !desktopCancel)
                return fail();
            if (desktopCancel)
                return this._calCancelled(backup, communicationFailed);

            if (this.trigger) {
                this.trigger = false;
                ++calStage;
                if (!await this.sendResponse(C.sCaliStageUpd, Uint8Array.of(calStage)))
                    return fail();
                switch (calStage) {
                case 2: top = this.mouse.y; if (!await this.sendResponse(C.sCaliInfoUpd, info(1, i32(top)))) return fail(); break;
                case 3: bottom = 55; if (!await this.sendResponse(C.sCaliInfoUpd, info(2, i32(bottom)))) return fail(); break;
                case 4: left = this.mouse.x; if (!await this.sendResponse(C.sCaliInfoUpd, info(3, i32(left)))) return fail(); break;
                case 5:
                    right = -7;
                    if (!await this.sendResponse(C.sCaliInfoUpd, info(4, i32(right)))) return fail();
                    view.setInt32(0, top, true); view.setInt32(4, bottom, true); view.setInt32(8, left, true); view.setInt32(12, right, true);
                    break;
                case 6: {
                    if (!await this.sendResponse(C.sCaliInfoUpd, info(5, this.profiles[this.currentProfile].subarray(16, 20)))) return fail();
                    if (!await this.sendResponse(C.sCaliInfoUpd, info(6, this.profiles[this.currentProfile].subarray(20, 24)))) return fail();
                    this.gunMode = 'verification';
                    while (this.gunMode === 'verification') {
                        await this.serialProcessingDocked();
                        const cancel = fromDesktop && this.takeCalibrationCancel();
                        if (cancel)
                            return this._calCancelled(backup, communicationFailed);
                        if (this.trigger) {
                            this.trigger = false;
                            calStage++;
                            this.gunMode = 'docked-end';
                            break;
                        } else if (this.buttonA) {
                            this.buttonA = false;
                            calStage = 0;
                            if (!await this.sendResponse(C.sCaliStageUpd, Uint8Array.of(0))) return fail();
                            view.setInt32(0, 0, true); view.setInt32(4, 0, true); view.setInt32(8, 0, true); view.setInt32(12, 0, true);
                            this.gunMode = 'calibration';
                        }
                        await tick();
                    }
                    break;
                }
                default: break;
                }
            }
            await tick();
        }

        this.gunMode = 'docked';
        calStage = 7;
        if (!await this.sendResponse(C.sCaliStageUpd, Uint8Array.of(calStage)))
            await this.sendError();
    }

    async _calCancelled(backup, communicationFailed) {
        const C = this.C;
        this.profiles[this.currentProfile].set(backup);
        this.gunMode = 'docked';
        if (!communicationFailed) {
            if (!await this.sendResponse(C.sCaliStageUpd, Uint8Array.of(0))) communicationFailed = true;
            else if (!await this.sendResponse(C.sCaliStageUpd, Uint8Array.of(7))) communicationFailed = true;
        }
        if (communicationFailed)
            await this.sendError();
    }

    // ----- Main loops (main.cpp) --------------------------------------------

    boardInfo() {
        const bytes = [];
        for (const ch of this.version) bytes.push(ch.charCodeAt(0));
        bytes.push(this.C.serialTerminator);
        for (const ch of this.boardType) bytes.push(ch.charCodeAt(0));
        bytes.push(this.C.serialTerminator);
        bytes.push(...this.usb);
        if (this.camNotAvailable)
            bytes.push(this.C.serialTerminator, this.C.sError);
        return Uint8Array.from(bytes);
    }

    start() {
        if (this.running)
            return;
        this.running = true;
        this._loop = this._mainLoop();
    }

    async stop() {
        this.running = false;
        try { await this._loop; } catch (e) { /* ignored */ }
    }

    async _mainLoop() {
        while (this.running) {
            if (this.gunMode === 'run') {
                // OF_Serial::SerialProcessing(): Enter Docked Mode from the serial link.
                if (this.linkAvailable('serial')) {
                    const incoming = this.linkRead('serial');
                    if (incoming === this.C.sDock1) {
                        const started = Date.now();
                        while (!this.linkAvailable('serial') && Date.now() - started < 1000 && this.running)
                            await tick();
                        if (!this.sessionActive && this.linkAvailable('serial') && this.linkRead('serial') === this.C.sDock2)
                            this.enterDocked('serial');
                    }
                }
                // OF_Serial::SerialProcessingWebDock()
                if (this.webMode && this.gunMode === 'run') {
                    if (!this.sessionActive) this.takeClientLost();
                    this.pollDock('web');
                }
            } else if (this.gunMode === 'docked') {
                await this._execGunModeDocked();
            }
            await tick();
        }
    }

    async _execGunModeDocked() {
        let sendBoardInfo = true;
        let lastCoords = 0;
        while (this.running) {
            if (sendBoardInfo) {
                const sent = await this.sendResponse(this.C.sDock2, this.boardInfo());
                sendBoardInfo = false;
                // Nobody answered the board info: end the unanswered session, so
                // neither link stays blocked by it.
                if (!sent && this.sessionActive && !this.sessionConfirmed)
                    this.abandonSession();
                if (this.gunMode !== 'docked')
                    return;
            }

            if (!this.dockedSaving) {
                while (this.pendingEvents.length) {
                    const [command, payload] = this.pendingEvents.shift();
                    this.sendEvent(command, payload);
                }
            }

            if (this.inputPending()) {
                const sessionBefore = this.sessionCounter;
                await this.serialProcessingDocked();
                if (this.sessionActive && this.sessionCounter !== sessionBefore)
                    sendBoardInfo = true;
            }

            if (this.gunMode !== 'docked')
                return;

            if (!this.dockedSaving && this.runMode === 'processing' && Date.now() - lastCoords >= 50) {
                lastCoords = Date.now();
                const coords = new Uint8Array(48);
                const v = new DataView(coords.buffer);
                const values = [600 * 2, 300, 1320 * 2 + 1, 300, 600 * 2, 780, 1320 * 2, 780, 960, 540, 955, 545];
                values.forEach((value, i) => v.setInt32(i * 4, value, true));
                this.sendEvent(this.C.sTestCoords, coords);
            }
            await tick();
        }
    }

    queueEvent(command, payload) {
        this.pendingEvents.push([command, payload]);
    }

    // ----- Accessors for tests ------------------------------------------------

    getSetting(index) { return new DataView(this.settings.buffer).getUint32(index * 4, true); }
    getProfileName(i) {
        const raw = this.profiles[i].subarray(52, 68);
        let s = '';
        for (const b of raw) { if (!b) break; s += String.fromCharCode(b); }
        return s;
    }
    getProfileU32(i, index) { return new DataView(this.profiles[i].buffer).getUint32(index * 4, true); }
    getProfileI32(i, index) { return new DataView(this.profiles[i].buffer).getInt32(index * 4, true); }
    getUSB() {
        const v = new DataView(this.usb.buffer);
        let s = '';
        for (const b of this.usb.subarray(2)) { if (!b) break; s += String.fromCharCode(b); }
        return { id: v.getUint16(0, true), name: s };
    }
}

module.exports = { MockFirmware, MemoryLink, createMemoryTransport, crc8, mulberry32 };
