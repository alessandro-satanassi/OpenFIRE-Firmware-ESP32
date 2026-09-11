/*  OpenFIRE Web App: App <-> firmware serial protocol.

    JavaScript port of the Qt App's AppSerial (src/appserial.cpp), the peer of
    the firmware's OF_Serial (lightgun/src/OpenFIREserial.cpp). Behaviour,
    timeouts, retries and recovery rules intentionally mirror appserial.cpp;
    only blocking waits are replaced by promises.

    Wire format: A5 5A | TYPE_FLAGS | COMMAND | SEQUENCE | LENGTH | PAYLOAD | CRC8
    CRC-8 (poly 0x9B, init 0) covers TYPE_FLAGS through PAYLOAD.

    Every command code and settings field name comes from OpenFIREshared.js,
    which is generated from boards/OpenFIREshared.h: nothing is duplicated here.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/

(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    // ===== Constants (appserial.h / appcommon.h) =============================

    const APP_SERIAL_START_1 = 0xA5;
    const APP_SERIAL_START_2 = 0x5A;
    const APP_SERIAL_MAX_PAYLOAD = 200;
    const APP_SERIAL_OVERHEAD = 7;

    // Transport limits in milliseconds; retries exclude the first send.
    const APP_SERIAL_FRAME_TIMEOUT = 250;
    const APP_SERIAL_ACK_TIMEOUT = 500;
    const APP_SERIAL_MAX_RETRIES = 3;
    const APP_SERIAL_MAX_PROFILE_COUNT = 20;

    const TYPE_REQUEST = 0x00;
    const TYPE_RESPONSE = 0x01;
    const TYPE_EVENT = 0x02;
    const TYPE_ACK = 0x03;
    const TYPE_MASK = 0x03;
    const FLAG_FINAL = 0x80;

    // Commit-retry sError payload: [retry required] [originating request command].
    const APP_SERIAL_ERR_COMMIT_RETRY = 0x82;

    // Operation deadlines (appserial.h, private).
    const APP_SERIAL_SETTINGS_TIMEOUT = 2500;
    const APP_SERIAL_COMMIT_START_TIMEOUT = 2500;
    const APP_SERIAL_COMMIT_READY_TIMEOUT = 1000;
    const APP_SERIAL_SAVE_TIMEOUT = 5000;
    const APP_SERIAL_DISCONNECT_TIMEOUT = 2500;

    // Local Dock policy (AppSerial::BeginDock).
    const DOCK_TIMEOUT = 3000;
    const DOCK_RETRY_INTERVAL = 750;

    // appcommon.h: BUTTON_COUNT; the App maps BUTTON_COUNT - 1 buttons (Home is internal).
    const BUTTON_COUNT = 14;
    const BUTTON_DESC_SIZE = 6;

    const TINYUSB_TABLE_SIZE = 18;   // uint16 PID + char[16] name
    const PROFILE_NAME_SIZE = 16;
    const PROFILE_RAW_SIZE = 68;     // 13 * 4 + 16, see App_Common::profilesTable_s

    // Profile layout (App_Common::profilesTable_s / OF_Prefs::ProfileData_s).
    // [profSyncTypes_e name, JS property, encoding]
    const PROFILE_FIELDS = [
        ['profTopOffset',    'topOffset',     'i32'],
        ['profBottomOffset', 'bottomOffset',  'i32'],
        ['profLeftOffset',   'leftOffset',    'i32'],
        ['profRightOffset',  'rightOffset',   'i32'],
        ['profTLled',        'TLled',         'f32'],
        ['profTRled',        'TRled',         'f32'],
        ['profAdjX',         'adjX',          'f32'],
        ['profAdjY',         'adjY',          'f32'],
        ['profIrSens',       'irSensitivity', 'u32'],
        ['profRunMode',      'runMode',       'u32'],
        ['profIrLayout',     'layoutType',    'u32'],
        ['profAR',           'aspectRatio',   'u32'],
        ['profColor',        'color',         'u32'],
        ['profName',         'name',          'str']
    ];

    const TIMEOUT = Symbol('timeout');

    const now = (typeof performance !== 'undefined' && performance.now) ?
        () => performance.now() : () => Date.now();

    const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

    function hex(value) {
        return '0x' + (value & 0xFF).toString(16).toUpperCase().padStart(2, '0');
    }

    // ===== Byte helpers ======================================================

    function crc8(data, offset, length) {
        let crc = 0;
        for (let i = 0; i < length; ++i) {
            crc ^= data[offset + i];
            for (let bit = 0; bit < 8; ++bit)
                crc = (crc & 0x80) ? ((crc << 1) ^ 0x9B) & 0xFF : (crc << 1) & 0xFF;
        }
        return crc;
    }

    function buildFrame(typeFlags, command, sequence, payload) {
        const length = payload ? payload.length : 0;
        if (length > APP_SERIAL_MAX_PAYLOAD) {
            console.warn('[AppSerial] Cannot build frame: invalid payload for command', hex(command), 'length', length);
            return null;
        }

        const frame = new Uint8Array(APP_SERIAL_OVERHEAD + length);
        frame[0] = APP_SERIAL_START_1;
        frame[1] = APP_SERIAL_START_2;
        frame[2] = typeFlags & 0xFF;
        frame[3] = command & 0xFF;
        frame[4] = sequence & 0xFF;
        frame[5] = length;
        if (length > 0)
            frame.set(payload, 6);
        frame[6 + length] = crc8(frame, 2, 4 + length);
        return frame;
    }

    // char[] fields are C strings in the firmware: bytes up to the first NUL (Latin-1).
    function decodeCString(bytes) {
        let text = '';
        for (let i = 0; i < bytes.length && bytes[i] !== 0; ++i)
            text += String.fromCharCode(bytes[i]);
        return text;
    }

    // strncpy(dest, src, size - 1) semantics: truncated, NUL terminated, zero padded.
    function encodeCString(text, size) {
        const bytes = new Uint8Array(size);
        const value = String(text || '');
        const count = Math.min(value.length, size - 1);
        for (let i = 0; i < count; ++i) {
            const code = value.charCodeAt(i);
            bytes[i] = code > 255 ? 0x3F : code;
        }
        return bytes;
    }

    function concatBytes(parts) {
        let length = 0;
        for (const part of parts)
            length += part.length;
        const out = new Uint8Array(length);
        let offset = 0;
        for (const part of parts) {
            out.set(part, offset);
            offset += part.length;
        }
        return out;
    }

    function ascii(text) {
        const bytes = new Uint8Array(text.length);
        for (let i = 0; i < text.length; ++i)
            bytes[i] = text.charCodeAt(i) & 0xFF;
        return bytes;
    }

    // ===== Data model codecs =================================================

    class Codec {
        constructor(shared) {
            this.shared = shared;
            this.profileFieldByIndex = new Map();
            for (const [enumName, key, encoding] of PROFILE_FIELDS) {
                const index = shared.profSyncTypes_e[enumName];
                this.profileFieldByIndex.set(index, { key, encoding, index });
            }
        }

        get boolTypesCount() { return this.shared.boolTypes_e.boolTypesCount; }
        get settingsTypesCount() { return this.shared.settingsTypes_e.settingsTypesCount; }
        get boardInputsCount() { return this.shared.boardInputs_e.boardInputsCount; }

        createProfile() {
            return {
                topOffset: 0, bottomOffset: 0, leftOffset: 0, rightOffset: 0,
                TLled: 0, TRled: 0, adjX: 0, adjY: 0,
                irSensitivity: 0, runMode: 0, layoutType: 0, aspectRatio: 0, color: 0,
                name: ''
            };
        }

        /// Empty configuration, equivalent to App_Common's zeroed tables.
        createConfig() {
            const buttons = [];
            for (let i = 0; i < BUTTON_COUNT; ++i)
                buttons.push(new Array(BUTTON_DESC_SIZE).fill(0));

            return {
                toggles: new Array(this.boolTypesCount).fill(false),
                pins: new Array(this.boardInputsCount).fill(this.shared.boardInputs_e.btnUnmapped),
                settings: new Array(this.settingsTypesCount).fill(0),
                buttons,
                profiles: [],
                selectedProfile: -1,
                tinyUSB: { id: 0, name: '' }
            };
        }

        cloneConfig(config) {
            return {
                toggles: config.toggles.slice(),
                pins: config.pins.slice(),
                settings: config.settings.slice(),
                buttons: config.buttons.map((button) => button.slice()),
                profiles: config.profiles.map((profile) => Object.assign({}, profile)),
                selectedProfile: config.selectedProfile,
                tinyUSB: Object.assign({}, config.tinyUSB)
            };
        }

        encodeProfileField(profile, index) {
            const field = this.profileFieldByIndex.get(index);
            if (!field)
                return null;
            if (field.encoding === 'str')
                return encodeCString(profile[field.key], PROFILE_NAME_SIZE);

            const bytes = new Uint8Array(4);
            const view = new DataView(bytes.buffer);
            const value = Number(profile[field.key]) || 0;
            if (field.encoding === 'i32') view.setInt32(0, value, true);
            else if (field.encoding === 'f32') view.setFloat32(0, value, true);
            else view.setUint32(0, value >>> 0, true);
            return bytes;
        }

        decodeProfileField(profile, index, bytes) {
            const field = this.profileFieldByIndex.get(index);
            if (!field)
                return;
            if (field.encoding === 'str') {
                profile[field.key] = decodeCString(bytes);
                return;
            }
            const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
            if (field.encoding === 'i32') profile[field.key] = view.getInt32(0, true);
            else if (field.encoding === 'f32') profile[field.key] = view.getFloat32(0, true);
            else profile[field.key] = view.getUint32(0, true);
        }

        /// Raw 68-byte image of a profile (memcmp-equivalent comparisons).
        encodeProfile(profile) {
            const raw = new Uint8Array(PROFILE_RAW_SIZE);
            for (const index of this.profileFieldByIndex.keys()) {
                const bytes = this.encodeProfileField(profile, index);
                raw.set(bytes, 4 * index);
            }
            return raw;
        }

        encodeTinyUSB(tinyUSB) {
            const bytes = new Uint8Array(TINYUSB_TABLE_SIZE);
            new DataView(bytes.buffer).setUint16(0, (tinyUSB.id >>> 0) & 0xFFFF, true);
            bytes.set(encodeCString(tinyUSB.name, 16), 2);
            return bytes;
        }

        decodeTinyUSB(bytes) {
            const view = new DataView(bytes.buffer, bytes.byteOffset, TINYUSB_TABLE_SIZE);
            return { id: view.getUint16(0, true), name: decodeCString(bytes.subarray(2, TINYUSB_TABLE_SIZE)) };
        }
    }

    // ===== Protocol ==========================================================

    class OpenFIREProtocol {
        constructor(options = {}) {
            // OpenFIREshared.js declares a top-level const: it is a global lexical
            // binding shared by classic scripts, not a property of window.
            this.shared = options.shared || root.OpenFIREshared ||
                (typeof OpenFIREshared !== 'undefined' ? OpenFIREshared : null); // eslint-disable-line no-undef
            if (!this.shared)
                throw new Error('OpenFIREshared data is not loaded.');

            this.cmd = this.shared.serialCmdTypes_e;
            this.codec = new Codec(this.shared);
            this.debug = !!options.debug;

            this.transport = null;
            this._portOpen = false;
            this._closingTransport = false;

            this._listeners = new Map();
            this._waiters = [];
            this._lock = Promise.resolve();

            this._resetProtocolState();
        }

        // ----- Events --------------------------------------------------------
        // 'event'    (command, payload)  decoded device events and local errors
        //                                (Qt: Serial_EventReceived, queued)
        // 'progress' ({ range } | { value, text })   (Qt: Serial_SetProgressRange / Serial_ProgressUpdate)
        // 'closed'   (reason)            the transport was lost unexpectedly

        on(name, listener) {
            if (!this._listeners.has(name))
                this._listeners.set(name, new Set());
            this._listeners.get(name).add(listener);
            return () => this.off(name, listener);
        }

        off(name, listener) {
            const set = this._listeners.get(name);
            if (set)
                set.delete(listener);
        }

        _emit(name, ...args) {
            const set = this._listeners.get(name);
            if (!set)
                return;
            for (const listener of Array.from(set)) {
                try {
                    listener(...args);
                } catch (error) {
                    console.error(`[AppSerial] '${name}' listener failed:`, error);
                }
            }
        }

        // Qt delivers Serial_EventReceived through a queued connection.
        _emitQueued(name, ...args) {
            const deliver = () => this._emit(name, ...args);
            if (typeof queueMicrotask === 'function') queueMicrotask(deliver);
            else Promise.resolve().then(deliver);
        }

        _progressRange(range) { this._emit('progress', { range }); }
        _progress(value, text) { this._emit('progress', { value, text: text || null }); }

        // ----- State ---------------------------------------------------------

        get isOpen() { return this._portOpen && !!this.transport && this.transport.isOpen; }
        get sessionActive() { return this.appSerialSessionActive; }
        commitNeedsRetry() { return this.appSerialCommitNeedsRetry; }

        _resetProtocolState() {
            this.appSerialRxBuffer = new Uint8Array(0);
            this.appSerialResponses = [];
            this.appSerialRxTimestamp = -1;
            this.appSerialTxSequence = 0;
            this.appSerialLastRequest = 0;
            this.appSerialAckReceived = false;
            this.appSerialAckCommand = 0;
            this.appSerialAckSequence = 0;
            this.appSerialLastRxValid = false;
            this.appSerialLastRxTypeFlags = 0;
            this.appSerialLastRxCommand = 0;
            this.appSerialLastRxSequence = 0;
            this.appSerialLastRxLength = 0;
            this.appSerialLastRxCRC = 0;

            this.appSerialSessionActive = false;
            this.appSerialCommitInProgress = false;
            this.appSerialCommitNeedsRetry = false;
            this.appSerialCommitError = false;
            this.appSerialLastStartSequence = 0;
            this.appSerialLastSaveSequence = 0;
        }

        _nextSequence() {
            this.appSerialTxSequence = (this.appSerialTxSequence + 1) & 0xFF;
            if (this.appSerialTxSequence === 0)
                this.appSerialTxSequence = 1;
            return this.appSerialTxSequence;
        }

        // Serialise public operations: the Qt GUI thread is blocked while a
        // transaction runs, so two transactions can never overlap.
        _exclusive(operation) {
            const run = this._lock.then(operation, operation);
            this._lock = run.then(() => undefined, () => undefined);
            return run;
        }

        // ----- Waiting -------------------------------------------------------

        /// Resolves with the first non-undefined value returned by check(),
        /// re-evaluated whenever frames arrive or the state changes,
        /// or with TIMEOUT once timeoutMs expires.
        _waitFor(check, timeoutMs) {
            const immediate = check();
            if (immediate !== undefined)
                return Promise.resolve(immediate);
            if (timeoutMs <= 0)
                return Promise.resolve(TIMEOUT);

            return new Promise((resolve) => {
                const waiter = { check, done: false };
                waiter.finish = (value) => {
                    if (waiter.done)
                        return;
                    waiter.done = true;
                    clearTimeout(waiter.timer);
                    const index = this._waiters.indexOf(waiter);
                    if (index >= 0)
                        this._waiters.splice(index, 1);
                    resolve(value);
                };
                waiter.timer = setTimeout(() => waiter.finish(TIMEOUT), timeoutMs);
                this._waiters.push(waiter);
            });
        }

        _notifyWaiters() {
            for (const waiter of this._waiters.slice()) {
                if (waiter.done)
                    continue;
                const value = waiter.check();
                if (value !== undefined)
                    waiter.finish(value);
            }
        }

        _commitAborted() {
            return this.appSerialCommitInProgress && this.appSerialCommitError;
        }

        // ----- Connection ----------------------------------------------------

        /// Opens a transport and binds it to this protocol instance.
        async connect(transport) {
            return this._exclusive(async () => {
                if (this.transport && this.transport.isOpen)
                    await this._closeTransport();

                this._resetProtocolState();
                this.transport = transport;
                transport.onData = (bytes) => this._processIncoming(bytes);
                transport.onClose = (reason) => this._onTransportClose(transport, reason);

                await transport.open();
                this._portOpen = true;
                return true;
            });
        }

        _onTransportClose(transport, reason) {
            if (transport !== this.transport)
                return;
            const unexpected = this._portOpen && !this._closingTransport;
            this._portOpen = false;
            this._notifyWaiters();
            if (unexpected) {
                console.warn('[AppSerial] Transport closed:', reason);
                this._resetProtocolState();
                this._emitQueued('closed', reason);
            }
        }

        async _closeTransport() {
            const transport = this.transport;
            this._portOpen = false;
            this._notifyWaiters();
            if (!transport)
                return;
            this._closingTransport = true;
            try {
                await transport.close();
            } catch (error) {
                console.warn('[AppSerial] Transport close failed:', error);
            } finally {
                this._closingTransport = false;
            }
        }

        // ===== Frame encoding and incoming stream ============================

        _sendAck(command, sequence) {
            const frame = buildFrame(TYPE_ACK, command, sequence, null);
            if (!frame || !this.isOpen) {
                console.warn('[AppSerial] Cannot send ACK for command', hex(command), 'sequence', sequence,
                             'because the serial port is not available.');
                return false;
            }
            // Queued in order; never awaited by the parser (equivalent to the
            // Qt App flushing the tiny ACK immediately).
            this.transport.write(frame);
            return true;
        }

        _processIncoming(incoming) {
            if (!this._portOpen)
                return;

            const t = now();
            if (this.appSerialRxBuffer.length > 0 && this.appSerialRxTimestamp >= 0 &&
                t - this.appSerialRxTimestamp > APP_SERIAL_FRAME_TIMEOUT) {
                console.warn('[AppSerial] Partial frame timed out; discarded', this.appSerialRxBuffer.length, 'byte(s).');
                this.appSerialRxBuffer = new Uint8Array(0);
            }

            if (incoming && incoming.length) {
                this.appSerialRxBuffer = this.appSerialRxBuffer.length ?
                    concatBytes([this.appSerialRxBuffer, incoming]) : new Uint8Array(incoming);
                this.appSerialRxTimestamp = t;
            }

            this._parseIncoming();
            this._notifyWaiters();
        }

        _parseIncoming() {
            for (;;) {
                let buffer = this.appSerialRxBuffer;
                if (buffer.length < 2)
                    return;

                let start = -1;
                for (let i = 0; i + 1 < buffer.length; ++i) {
                    if (buffer[i] === APP_SERIAL_START_1 && buffer[i + 1] === APP_SERIAL_START_2) {
                        start = i;
                        break;
                    }
                }

                if (start < 0) {
                    if (this.debug)
                        console.debug('[AppSerial] Discarding', buffer.length, 'byte(s) without a complete frame start marker.');
                    this.appSerialRxBuffer = buffer[buffer.length - 1] === APP_SERIAL_START_1 ?
                        Uint8Array.of(APP_SERIAL_START_1) : new Uint8Array(0);
                    return;
                }

                if (start > 0) {
                    if (this.debug)
                        console.debug('[AppSerial] Resynchronizing after', start, 'unexpected byte(s).');
                    buffer = this.appSerialRxBuffer = buffer.subarray(start);
                }

                if (buffer.length < 6)
                    return;

                const typeFlags = buffer[2];
                const type = typeFlags & TYPE_MASK;
                const command = buffer[3];
                const sequence = buffer[4];
                const length = buffer[5];
                const invalidFlags = (typeFlags & ~(TYPE_MASK | FLAG_FINAL) & 0xFF) !== 0;
                const invalidFinal = (typeFlags & FLAG_FINAL) !== 0 && (type !== TYPE_RESPONSE || length !== 0);
                const invalidSequence = type === TYPE_EVENT ? sequence !== 0 : sequence === 0;
                const invalidAck = type === TYPE_ACK && length !== 0;

                if (invalidFlags || invalidFinal || invalidSequence || invalidAck || length > APP_SERIAL_MAX_PAYLOAD) {
                    console.warn('[AppSerial] Invalid frame header: type/flags', hex(typeFlags), 'command', hex(command),
                                 'sequence', sequence, 'length', length, '; resynchronizing.');
                    this.appSerialRxBuffer = buffer.subarray(1);
                    continue;
                }

                const frameLength = APP_SERIAL_OVERHEAD + length;
                if (buffer.length < frameLength)
                    return;

                const receivedCRC = buffer[6 + length];
                const calculatedCRC = crc8(buffer, 2, 4 + length);
                if (receivedCRC !== calculatedCRC) {
                    console.warn('[AppSerial] CRC mismatch for command', hex(command), 'sequence', sequence,
                                 'received', hex(receivedCRC), 'expected', hex(calculatedCRC), '; frame discarded.');
                    this.appSerialRxBuffer = buffer.subarray(1);
                    continue;
                }

                const frame = {
                    typeFlags,
                    command,
                    sequence,
                    payload: buffer.slice(6, 6 + length),
                    crc: receivedCRC
                };
                this.appSerialRxBuffer = buffer.subarray(frameLength);

                this._handleFrame(frame);
                if (!this._portOpen)
                    return;
            }
        }

        // ===== Response routing and recovery =================================

        _isTransactionResponse(command) {
            const c = this.cmd;
            switch (command) {
            case c.sDock2:
            case c.sGetPins:
            case c.sGetToggles:
            case c.sGetSettings:
            case c.sGetProfile:
            case c.sGetBtns:
            case c.sCommitStart:
            case c.sSave:
            case c.serialTerminator:
                return true;
            default:
                return false;
            }
        }

        _handleFrame(frame) {
            const c = this.cmd;
            const type = frame.typeFlags & TYPE_MASK;
            let deliverEvent = type === TYPE_EVENT;

            if (type === TYPE_ACK) {
                this.appSerialAckCommand = frame.command;
                this.appSerialAckSequence = frame.sequence;
                this.appSerialAckReceived = true;
                return;
            }

            // Reliable responses are ACKed even when duplicate or stale. Only an
            // accepted response reaches the transaction queue or the event path.
            if (type === TYPE_RESPONSE) {
                const length = frame.payload.length;
                const duplicate =
                    this.appSerialLastRxValid &&
                    this.appSerialLastRxTypeFlags === frame.typeFlags &&
                    this.appSerialLastRxCommand === frame.command &&
                    this.appSerialLastRxSequence === frame.sequence &&
                    this.appSerialLastRxLength === length &&
                    this.appSerialLastRxCRC === frame.crc;

                if (!duplicate) {
                    this.appSerialLastRxValid = true;
                    this.appSerialLastRxTypeFlags = frame.typeFlags;
                    this.appSerialLastRxCommand = frame.command;
                    this.appSerialLastRxSequence = frame.sequence;
                    this.appSerialLastRxLength = length;
                    this.appSerialLastRxCRC = frame.crc;
                }

                this._sendAck(frame.command, frame.sequence);

                if (duplicate)
                    return;

                // Commit results are tied to the request, not merely its command.
                if (frame.command === c.sError &&
                    frame.payload.length === 2 &&
                    frame.payload[0] === APP_SERIAL_ERR_COMMIT_RETRY) {
                    if (frame.payload[1] !== this.appSerialLastRequest ||
                        frame.sequence !== this.appSerialTxSequence)
                        return; // ACKed, but belongs to an older attempt.
                    this.appSerialCommitNeedsRetry = true;
                    this.appSerialCommitError = this.appSerialCommitInProgress;
                    if (!this.appSerialCommitInProgress)
                        this._emitQueued('event', frame.command, frame.payload);
                    return;
                }

                const transactionResponse = this._isTransactionResponse(frame.command);
                const expectedResponse =
                    transactionResponse &&
                    (frame.command === this.appSerialLastRequest ||
                     (frame.command === c.sDock2 && this.appSerialLastRequest === 0)) &&
                    ((frame.command !== c.sCommitStart && frame.command !== c.sSave) ||
                     frame.sequence === this.appSerialTxSequence);

                if (expectedResponse) {
                    if (frame.command === c.sDock2)
                        this.appSerialSessionActive = true;
                    this.appSerialResponses.push(frame);
                } else if (transactionResponse) {
                    if (this.debug)
                        console.debug('[AppSerial] Ignoring stale transaction response', hex(frame.command),
                                      'while waiting for', hex(this.appSerialLastRequest), '.');
                } else {
                    deliverEvent = true;
                }
            }

            // Both wire events and non-transaction responses use this delivery path.
            if (!deliverEvent)
                return;

            if (frame.command === c.sError) {
                console.warn(type === TYPE_EVENT ?
                    '[AppSerial] Device reported an asynchronous error; payload:' :
                    '[AppSerial] Device reported an error; payload:', Array.from(frame.payload));
            }

            if (frame.command !== c.sError || frame.payload.length > 0) {
                this._emitQueued('event', frame.command, frame.payload);
                return;
            }

            if (this.appSerialCommitInProgress)
                this.appSerialCommitError = true;
            else if (this.appSerialCommitNeedsRetry)
                this._emitQueued('event', c.sError, Uint8Array.of(APP_SERIAL_ERR_COMMIT_RETRY));
            else
                this.failOperation();
        }

        // Commit failure is recoverable while both peers retain their RAM.
        _commitFailed() {
            this.appSerialCommitInProgress = false;
            this.appSerialCommitNeedsRetry = true;
            console.warn('[AppSerial] Save failed or was not confirmed. RAM is retained; retry Save without restarting the board.');
            this._emitQueued('event', this.cmd.sError, Uint8Array.of(APP_SERIAL_ERR_COMMIT_RETRY));
            return false;
        }

        /// Close on a fatal error; the board must be restarted. Returns false.
        failOperation() {
            const notifyError = this._portOpen;
            this._resetProtocolState();
            if (this._portOpen)
                this._closeTransport();
            this._notifyWaiters();

            if (notifyError) {
                console.warn('[AppSerial] Operation failed or not confirmed. Restart the board and reconnect before continuing.');
                this._emitQueued('event', this.cmd.sError, new Uint8Array(0));
            }
            return false;
        }

        // ===== Reliable command delivery and response waits ==================

        async _waitForAck(command, sequence) {
            const result = await this._waitFor(() => {
                if (!this.isOpen || this._commitAborted())
                    return false;
                if (this.appSerialAckReceived &&
                    this.appSerialAckCommand === command &&
                    this.appSerialAckSequence === sequence) {
                    this.appSerialAckReceived = false;
                    return true;
                }
                // A valid response proves the request was delivered even if its ACK was lost.
                if (this.appSerialResponses.some((response) => response.command === command))
                    return true;
                return undefined;
            }, APP_SERIAL_ACK_TIMEOUT);

            return result === true;
        }

        /// Sends one framed command and waits for delivery acknowledgement
        /// (not operation completion). Public entry point, serialised.
        sendCommand(command, payload = null) {
            return this._exclusive(() => this._sendCommand(command, payload));
        }

        async _sendCommand(command, payload = null) {
            const length = payload ? payload.length : 0;
            if (length > APP_SERIAL_MAX_PAYLOAD) {
                console.warn('[AppSerial] SendCommand rejected invalid payload for command', hex(command), 'length', length, '.');
                return false;
            }

            if (!this.isOpen) {
                console.warn('[AppSerial] Cannot send command', hex(command), ': serial port is closed.');
                return false;
            }

            let sequence = this._nextSequence();
            // Even if record traffic wraps the byte counter, do not reuse the last
            // start/save ID for the very next attempt of the same command.
            if (command === this.cmd.sCommitStart) {
                if (sequence === this.appSerialLastStartSequence)
                    sequence = this._nextSequence();
                this.appSerialLastStartSequence = sequence;
            } else if (command === this.cmd.sSave) {
                if (sequence === this.appSerialLastSaveSequence)
                    sequence = this._nextSequence();
                this.appSerialLastSaveSequence = sequence;
            }

            const frame = buildFrame(TYPE_REQUEST, command, sequence, payload);
            if (!frame) {
                console.warn('[AppSerial] Frame creation failed for command', hex(command), 'sequence', sequence, '.');
                return false;
            }

            this.appSerialResponses = this.appSerialResponses.filter((response) => response.command !== command);
            this.appSerialLastRequest = command;

            for (let attempt = 0; attempt <= APP_SERIAL_MAX_RETRIES; ++attempt) {
                if (this._commitAborted())
                    return false;
                this.appSerialAckReceived = false;

                const written = await this.transport.write(frame);
                if (!written) {
                    console.warn('[AppSerial] Incomplete frame write for command', hex(command), 'sequence', sequence, '.');
                    return false;
                }

                if (await this._waitForAck(command, sequence))
                    return true;
                if (!this.isOpen || this._commitAborted())
                    return false;

                console.warn('[AppSerial] ACK timeout for command', hex(command), 'sequence', sequence,
                             'attempt', attempt + 1, 'of', APP_SERIAL_MAX_RETRIES + 1,
                             attempt < APP_SERIAL_MAX_RETRIES ? '; retrying.' : '.');
            }

            console.warn('[AppSerial] Command', hex(command), 'sequence', sequence, 'failed after',
                         APP_SERIAL_MAX_RETRIES + 1, 'transmission attempts.');
            return false;
        }

        /// take = true removes and returns the response frame; false only waits.
        /// Resolves to the frame, true (take = false), or null on failure/timeout.
        async _waitForResponse(command, take, timeoutMs) {
            const result = await this._waitFor(() => {
                if (!this.isOpen || this._commitAborted())
                    return null;
                const index = this.appSerialResponses.findIndex((response) => response.command === command);
                if (index < 0)
                    return undefined;
                if (!take)
                    return true;
                return this.appSerialResponses.splice(index, 1)[0];
            }, timeoutMs);

            if (result === TIMEOUT) {
                console.warn(take ? '[AppSerial] Timed out while retrieving response for command' :
                                    '[AppSerial] Response timeout for command', hex(command), 'after', timeoutMs, 'ms.');
                return null;
            }
            return result;
        }

        // ===== Dock handshake: raw request, framed response ==================

        async _beginDock() {
            const c = this.cmd;
            this._resetProtocolState();

            const dock = Uint8Array.of(c.sDock1, c.sDock2);
            const started = now();
            let nextAttempt = 0;

            for (;;) {
                if (!this.isOpen)
                    return null;

                const index = this.appSerialResponses.findIndex((response) => response.command === c.sDock2);
                if (index >= 0) {
                    const frame = this.appSerialResponses.splice(index, 1)[0];
                    if ((frame.typeFlags & FLAG_FINAL) || frame.payload.length === 0) {
                        console.warn('[AppSerial] Invalid board-information response: unexpected final flag or empty payload.');
                        return null;
                    }
                    return frame.payload;
                }

                const elapsed = now() - started;
                if (elapsed >= DOCK_TIMEOUT)
                    break;

                if (elapsed >= nextAttempt) {
                    if (!await this.transport.write(dock))
                        console.warn('[AppSerial] Incomplete initial Dock write; the request will be retried.');
                    nextAttempt = (now() - started) + DOCK_RETRY_INTERVAL;
                }

                const current = now() - started;
                const waitTime = Math.max(1, Math.min(DOCK_TIMEOUT - current, nextAttempt - current));
                await this._waitFor(() => {
                    if (!this.isOpen)
                        return false;
                    return this.appSerialResponses.some((response) => response.command === c.sDock2) ? true : undefined;
                }, waitTime);
            }

            console.warn('[AppSerial] No valid board-information frame received during Dock.');
            return null;
        }

        // ===== Initial synchronization =======================================

        /// Docks and loads the device settings (AppSerial::GetSettings).
        /// Resolves to { ok: true, board, config } or { ok: false, error, rebootSuggested }.
        getSettings() {
            return this._exclusive(() => this._getSettings());
        }

        async _getSettings() {
            const c = this.cmd;
            const shared = this.shared;
            const codec = this.codec;

            if (!this.isOpen)
                return { ok: false, error: 'not_open', rebootSuggested: false };

            const boardInfo = await this._beginDock();
            if (!boardInfo)
                return { ok: false, error: 'dock_timeout', rebootSuggested: this.isOpen };

            const separator = c.serialTerminator;
            const firstSeparator = boardInfo.indexOf(separator);
            const secondSeparator = firstSeparator >= 0 ? boardInfo.indexOf(separator, firstSeparator + 1) : -1;
            const usbOffset = secondSeparator + 1;

            if (firstSeparator <= 0 || secondSeparator <= firstSeparator + 1 ||
                boardInfo.length < usbOffset + TINYUSB_TABLE_SIZE) {
                console.warn('[AppSerial] Port did not respond with a valid board-information payload.');
                return { ok: false, error: 'bad_board_info', rebootSuggested: true };
            }

            this._progressRange(6);
            this._progress(1, 'Getting Board Info');

            const board = {
                version: decodeCString(boardInfo.subarray(0, firstSeparator)),
                type: decodeCString(boardInfo.subarray(firstSeparator + 1, secondSeparator)),
                arch: '',
                cameraError: false
            };
            board.arch = boardArch(shared, board.type);

            const config = codec.createConfig();
            config.tinyUSB = codec.decodeTinyUSB(boardInfo.subarray(usbOffset, usbOffset + TINYUSB_TABLE_SIZE));

            const tailOffset = usbOffset + TINYUSB_TABLE_SIZE;
            if (boardInfo.length >= tailOffset + 2 &&
                boardInfo[tailOffset] === separator &&
                boardInfo[tailOffset + 1] === c.sError)
                board.cameraError = true;

            const fail = (error) => ({ ok: false, error, rebootSuggested: false, board });
            const boolTypes = shared.boolTypes_e;

            // For the fixed tables, wait for the first record before decoding.
            if (!await this._sendCommand(c.sGetToggles) ||
                !await this._waitForResponse(c.sGetToggles, false, APP_SERIAL_SETTINGS_TIMEOUT))
                return fail('toggles');
            if (!await this._receiveSettingsRecords(c.sGetToggles, config))
                return fail('toggles');

            this._progress(2, 'Getting Settings (1)');

            if (config.toggles[boolTypes.customPins]) {
                if (!await this._sendCommand(c.sGetPins) ||
                    !await this._waitForResponse(c.sGetPins, false, APP_SERIAL_SETTINGS_TIMEOUT))
                    return fail('pins');
                if (!await this._receiveSettingsRecords(c.sGetPins, config))
                    return fail('pins');
            }

            this._progress(3, 'Getting Settings (2)');

            if (!await this._sendCommand(c.sGetSettings) ||
                !await this._waitForResponse(c.sGetSettings, false, APP_SERIAL_SETTINGS_TIMEOUT))
                return fail('settings');
            if (!await this._receiveSettingsRecords(c.sGetSettings, config))
                return fail('settings');

            this._progress(4, 'Getting Button Mappings');

            if (!await this._sendCommand(c.sGetBtns) ||
                !await this._waitForResponse(c.sGetBtns, false, APP_SERIAL_SETTINGS_TIMEOUT))
                return fail('buttons');
            if (!await this._receiveSettingsRecords(c.sGetBtns, config))
                return fail('buttons');

            this._progress(5, 'Getting Profiles Data');

            if (!await this._sendCommand(c.sGetProfile) ||
                !await this._waitForResponse(c.sGetProfile, false, APP_SERIAL_SETTINGS_TIMEOUT))
                return fail('profiles');
            if (!await this._receiveSettingsRecords(c.sGetProfile, config))
                return fail('profiles');

            this._progress(6, 'Successfully synced data!');
            return { ok: true, board, config };
        }

        // ===== Settings record decoding ======================================
        // Payload: field name + NUL, value size, optional profile number, value bytes.
        // CurrentProf has no profile-number byte; an empty FINAL ends the record stream.

        async _receiveSettingsRecords(command, config) {
            const c = this.cmd;
            const shared = this.shared;
            const codec = this.codec;

            let fieldMap;
            let dataSize;
            switch (command) {
            case c.sGetToggles:  fieldMap = shared.boolTypes_Strings;        dataSize = 1; break;
            case c.sGetPins:     fieldMap = shared.boardInputs_Strings;      dataSize = 1; break;
            case c.sGetSettings: fieldMap = shared.settingsTypes_Strings;    dataSize = 4; break;
            case c.sGetBtns:     fieldMap = shared.boardInputs_Strings;      dataSize = BUTTON_DESC_SIZE; break;
            case c.sGetProfile:  fieldMap = shared.profSettingTypes_Strings; dataSize = 4; break;
            default: return false;
            }

            const profCurrent = shared.profSyncTypes_e.profCurrent;
            const profName = shared.profSyncTypes_e.profName;

            for (;;) {
                const frame = await this._waitForResponse(command, true, APP_SERIAL_SETTINGS_TIMEOUT);
                if (!frame)
                    return false;

                if (frame.typeFlags & FLAG_FINAL) {
                    if (frame.payload.length !== 0) {
                        console.warn('[AppSerial] Malformed final record for command', hex(command), ': final payload must be empty.');
                        return false;
                    }

                    if (command === c.sGetProfile) {
                        if (config.profiles.length === 0 ||
                            config.profiles.length > APP_SERIAL_MAX_PROFILE_COUNT ||
                            config.selectedProfile < 0 ||
                            config.selectedProfile >= config.profiles.length) {
                            console.warn('[AppSerial] Incomplete or invalid profile set received.');
                            return false;
                        }
                    }
                    return true;
                }

                const payload = frame.payload;
                const nullPos = payload.indexOf(0);
                if (nullPos < 0 || nullPos + 1 >= payload.length) {
                    console.warn('[AppSerial] Malformed settings record for command', hex(command),
                                 ': missing field terminator or value length; payload size', payload.length, '.');
                    return false;
                }

                const fieldName = decodeCString(payload.subarray(0, nullPos));
                let pos = nullPos + 1;
                const valueSize = payload[pos++];

                if (!Object.prototype.hasOwnProperty.call(fieldMap, fieldName)) {
                    if (this.debug)
                        console.debug('[AppSerial] Ignoring unknown settings field', fieldName, 'for forward compatibility.');
                    continue;
                }

                const index = fieldMap[fieldName];
                if (index < 0) {
                    if (this.debug)
                        console.debug('[AppSerial] Ignoring disabled settings field', fieldName, '.');
                    continue;
                }

                if (command === c.sGetPins) {
                    if (valueSize !== 1 || pos + valueSize !== payload.length) {
                        console.warn('[AppSerial] Malformed pin record', fieldName, ': declared value size', valueSize,
                                     'payload size', payload.length, '.');
                        return false;
                    }
                    if (index < config.pins.length)
                        config.pins[index] = (payload[pos] << 24) >> 24;
                } else if (command === c.sGetProfile) {
                    if (index === profCurrent) {
                        // CurrentProf is intentionally a one-byte wire value.
                        if (valueSize !== 1 || pos + valueSize !== payload.length) {
                            console.warn('[AppSerial] Malformed current-profile record: declared value size', valueSize,
                                         'payload size', payload.length, '.');
                            return false;
                        }
                        const selectedProfile = payload[pos];
                        if (selectedProfile >= APP_SERIAL_MAX_PROFILE_COUNT) {
                            console.warn('[AppSerial] Invalid selected profile:', selectedProfile, '.');
                            return false;
                        }
                        config.selectedProfile = selectedProfile;
                    } else {
                        if (pos >= payload.length) {
                            console.warn('[AppSerial] Malformed profile record', fieldName, ': missing profile number.');
                            return false;
                        }
                        const profNum = payload[pos++];
                        if (profNum >= APP_SERIAL_MAX_PROFILE_COUNT) {
                            console.warn('[AppSerial] Invalid profile number:', profNum, '.');
                            return false;
                        }

                        while (profNum >= config.profiles.length)
                            config.profiles.push(codec.createProfile());

                        const expectedSize = index === profName ? PROFILE_NAME_SIZE : dataSize;
                        if (valueSize !== expectedSize || pos + valueSize !== payload.length) {
                            console.warn('[AppSerial] Malformed profile record', fieldName, 'profile', profNum,
                                         ': declared value size', valueSize, 'expected', expectedSize,
                                         'payload size', payload.length, '.');
                            return false;
                        }

                        codec.decodeProfileField(config.profiles[profNum], index, payload.subarray(pos, pos + valueSize));
                    }
                } else {
                    if (valueSize !== dataSize || pos + valueSize !== payload.length) {
                        console.warn('[AppSerial] Malformed settings record', fieldName, ': declared value size', valueSize,
                                     'expected', dataSize, 'payload size', payload.length, '.');
                        return false;
                    }

                    const value = payload.subarray(pos, pos + valueSize);
                    if (command === c.sGetToggles) {
                        if (index < config.toggles.length)
                            config.toggles[index] = value[0] !== 0;
                    } else if (command === c.sGetSettings) {
                        if (index < config.settings.length)
                            config.settings[index] = new DataView(value.buffer, value.byteOffset, 4).getUint32(0, true);
                    } else if (command === c.sGetBtns) {
                        // Bounded: only button descriptors exist on the firmware side.
                        if (index < BUTTON_COUNT)
                            config.buttons[index] = Array.from(value);
                    }
                }
            }
        }

        // ===== Settings commit ===============================================

        /// Sends the current settings and waits for the device's save result
        /// (AppSerial::CommitSettings). original is the configuration last
        /// loaded from (or saved to) the board.
        /// Resolves to { ok: true, message } or { ok: false }.
        commitSettings(current, original) {
            return this._exclusive(async () => {
                // Nothing to send (e.g. the page lost its board before Save started): no transaction.
                if (!current || !current.toggles || !current.settings || !current.buttons || !current.profiles || !current.tinyUSB) {
                    console.warn('[AppSerial] Cannot commit settings: no configuration.');
                    return { ok: false };
                }
                try {
                    return await this._commitSettings(current, original);
                } catch (error) {
                    console.error('[AppSerial] Commit failed:', error);
                    return { ok: this._commitFailed() };
                }
            });
        }

        async _commitSettings(current, original) {
            const c = this.cmd;
            const shared = this.shared;

            this.appSerialCommitInProgress = true;
            this.appSerialCommitNeedsRetry = true;
            this.appSerialCommitError = false;

            if (!this.isOpen) {
                console.warn('[AppSerial] Cannot commit settings: serial port is closed.');
                return { ok: this._commitFailed() };
            }

            this._progressRange(8);
            this._progress(0, 'Waiting for board...');

            // Preserve both phases: wait without consuming, then retrieve and validate.
            if (!await this._sendCommand(c.sCommitStart) ||
                !await this._waitForResponse(c.sCommitStart, false, APP_SERIAL_COMMIT_START_TIMEOUT)) {
                console.warn('[AppSerial] Commit failed while starting the transaction.');
                return { ok: this._commitFailed() };
            }

            const ready = await this._waitForResponse(c.sCommitStart, true, APP_SERIAL_COMMIT_READY_TIMEOUT);
            if (!ready) {
                console.warn('[AppSerial] Commit-start confirmation was not available.');
                return { ok: this._commitFailed() };
            }
            if (ready.payload.length !== 0 || (ready.typeFlags & FLAG_FINAL))
                return { ok: this._commitFailed() };

            this._progress(1, 'Sending Toggles...');
            if (!await this._sendSettingsRecords(c.sCommitToggles, shared.boolTypes_Strings, 1,
                    (index) => Uint8Array.of(current.toggles[index] ? 1 : 0)))
                return { ok: this._commitFailed() };

            // Include the pin group when custom pins are being disabled as well.
            const customPins = shared.boolTypes_e.customPins;
            if (current.toggles[customPins] || (original && original.toggles[customPins])) {
                this._progress(2, 'Sending Pins Map...');
                if (!await this._sendSettingsRecords(c.sCommitPins, shared.boardInputs_Strings, 1,
                        (index) => Uint8Array.of(index < current.pins.length ? current.pins[index] & 0xFF : 0xFF)))
                    return { ok: this._commitFailed() };
            }

            this._progress(3, 'Sending Settings...');
            if (!await this._sendSettingsRecords(c.sCommitSettings, shared.settingsTypes_Strings, 4, (index) => {
                const bytes = new Uint8Array(4);
                new DataView(bytes.buffer).setUint32(0, (current.settings[index] >>> 0), true);
                return bytes;
            }))
                return { ok: this._commitFailed() };

            this._progress(4, 'Sending Buttons...');
            if (!await this._sendSettingsRecords(c.sCommitBtns, shared.boardInputs_Strings, BUTTON_DESC_SIZE, (index) => {
                if (index >= BUTTON_COUNT - 1)
                    return null; // Home (and non-button inputs) are not managed by the App.
                const bytes = new Uint8Array(BUTTON_DESC_SIZE);
                bytes.set((current.buttons[index] || []).slice(0, BUTTON_DESC_SIZE).map((v) => v & 0xFF));
                return bytes;
            }))
                return { ok: this._commitFailed() };

            this._progress(5, 'Sending Profile Data...');
            for (let i = 0; i < current.profiles.length; ++i) {
                if (!await this._sendProfileRecords(current, i))
                    return { ok: this._commitFailed() };
            }

            this._progress(6, 'Sending TinyUSB ID Data...');
            if (!await this._sendCommand(c.sCommitID, this.codec.encodeTinyUSB(current.tinyUSB)))
                return { ok: this._commitFailed() };

            this._progress(7, 'Saving...');
            if (!await this._sendCommand(c.sSave))
                return { ok: this._commitFailed() };

            // A delivery ACK is not a save result: wait for the device's actual outcome.
            const saved = await this._waitForResponse(c.sSave, true, APP_SERIAL_SAVE_TIMEOUT);
            if (!saved) {
                console.warn('[AppSerial] No final save result received from the device.');
                return { ok: this._commitFailed() };
            }
            if (saved.payload.length === 0 || saved.payload[0] > 1) {
                console.warn('[AppSerial] Malformed save result.');
                return { ok: this._commitFailed() };
            }
            const message = decodeCString(saved.payload.subarray(1));
            if (saved.payload[0] === 0) {
                console.warn('[AppSerial] Device reported that settings could not be saved:', message);
                return { ok: this._commitFailed(), message };
            }

            this.appSerialCommitInProgress = false;
            this.appSerialCommitNeedsRetry = false;
            this.appSerialCommitError = false;
            this._progress(8);
            return { ok: true, message };
        }

        // ===== Settings record encoding ======================================

        async _sendSettingsRecords(command, fieldMap, dataSize, valueFor) {
            for (const [name, index] of Object.entries(fieldMap)) {
                if (index < 0)
                    continue;

                const value = valueFor(index);
                if (!value)
                    continue;

                const nameBytes = ascii(name);
                if (nameBytes.length + 2 + dataSize > APP_SERIAL_MAX_PAYLOAD) {
                    console.warn('[AppSerial] Settings record exceeds the payload limit:', name, '.');
                    return false;
                }

                const payload = concatBytes([nameBytes, Uint8Array.of(0, dataSize), value]);
                if (!await this._sendCommand(command, payload)) {
                    console.warn('[AppSerial] Failed to send settings field', name, 'using command', hex(command), '.');
                    return false;
                }
            }
            return true;
        }

        async _sendProfileRecords(current, profNum) {
            const c = this.cmd;
            const shared = this.shared;
            const profIrSens = shared.profSyncTypes_e.profIrSens;
            const profCurrent = shared.profSyncTypes_e.profCurrent;
            const profName = shared.profSyncTypes_e.profName;
            let profCurrentMarked = false;

            for (const [name, index] of Object.entries(shared.profSettingTypes_Strings)) {
                if (index < 0)
                    continue;
                // Calibration geometry stays in the firmware; send App-writable fields.
                if (index < profIrSens)
                    continue;

                const nameBytes = ascii(name);
                let payload;

                if (index === profCurrent) {
                    // Preserve the existing one-per-batch record and ignore aliases.
                    if (profCurrentMarked)
                        continue;
                    payload = concatBytes([nameBytes, Uint8Array.of(0, 1, current.selectedProfile & 0xFF)]);
                    profCurrentMarked = true;
                } else {
                    const valueSize = index === profName ? PROFILE_NAME_SIZE : 4;
                    const value = this.codec.encodeProfileField(current.profiles[profNum], index);
                    if (!value)
                        continue;
                    payload = concatBytes([nameBytes, Uint8Array.of(0, valueSize, profNum & 0xFF), value]);
                }

                if (payload.length > APP_SERIAL_MAX_PAYLOAD) {
                    console.warn('[AppSerial] Profile record exceeds the payload limit:', name, 'profile', profNum, '.');
                    return false;
                }

                if (!await this._sendCommand(c.sCommitProfile, payload)) {
                    console.warn('[AppSerial] Failed to send settings field', name, 'using command', hex(c.sCommitProfile), '.');
                    return false;
                }
            }
            return true;
        }

        // ===== Disconnect, bootloader and flash reset ========================

        /// Undocks (FE FE FE) when a session is active, then closes the transport.
        disconnect() {
            return this._exclusive(async () => {
                const c = this.cmd;
                if (this.isOpen && this.appSerialSessionActive) {
                    // Together with the command byte this preserves the FE FE FE
                    // full-disconnect marker used by the firmware.
                    if (await this._sendCommand(c.serialTerminator, Uint8Array.of(c.serialTerminator, c.serialTerminator))) {
                        const finalResponse = await this._waitForResponse(c.serialTerminator, true, APP_SERIAL_DISCONNECT_TIMEOUT);
                        if (finalResponse && (!(finalResponse.typeFlags & FLAG_FINAL) || finalResponse.payload.length !== 0))
                            console.warn('[AppSerial] Invalid final disconnect response.');
                    }
                }

                await this._closeTransport();
                this._resetProtocolState();
            });
        }

        /// Enters the bootloader: RP2040/RP235X use the 1200-baud touch when the
        /// transport supports it (Web Serial); otherwise the framed command is
        /// sent without waiting for an ACK (AppSerial::RebootToBootldr).
        rebootToBootloader(arch) {
            return this._exclusive(async () => {
                const shared = this.shared;
                const transport = this.transport;

                if (transport && this.isOpen) {
                    const isRP = arch === shared.boardArchs[shared.boardArchs_e.boardRP];
                    if (isRP && typeof transport.touch1200 === 'function') {
                        this._closingTransport = true;
                        this._portOpen = false;
                        try {
                            await transport.touch1200();
                        } catch (error) {
                            console.warn('[AppSerial] 1200-baud bootloader touch failed:', error);
                        } finally {
                            this._closingTransport = false;
                        }
                    } else {
                        const frame = buildFrame(TYPE_REQUEST, this.cmd.sRebootToBootloader, this._nextSequence(), null);
                        if (!frame || !await transport.write(frame))
                            console.warn('[AppSerial] Incomplete reboot command write.');
                        await this._closeTransport();
                    }
                }

                this._resetProtocolState();
                this._notifyWaiters();
            });
        }

        /// Clears the board's saved data (Qt: on_clearEepromBtn_clicked).
        /// The firmware formats its storage and restarts right after the ACK.
        clearSaveMemory() {
            return this._exclusive(async () => {
                const sent = await this._sendCommand(this.cmd.sClearFlash);
                await sleep(50);
                await this._closeTransport();
                this._resetProtocolState();
                return sent;
            });
        }
    }

    /// Architecture whose name is part of the board name ("esp32-s3"), else the first
    /// one (RP2040/235X): the Qt App rule, extended to architectures added to boardArchs.
    function boardArch(shared, boardType) {
        const archs = shared.boardArchs;
        let found = archs[0];
        for (const arch of archs.slice(1)) {
            if (arch && boardType.includes(arch) && (found === archs[0] || arch.length > found.length))
                found = arch;
        }
        return found;
    }

    OF.Protocol = OpenFIREProtocol;
    OF.ProtocolCodec = Codec;
    OF.ProtocolConstants = Object.freeze({
        APP_SERIAL_START_1, APP_SERIAL_START_2, APP_SERIAL_MAX_PAYLOAD, APP_SERIAL_OVERHEAD,
        APP_SERIAL_FRAME_TIMEOUT, APP_SERIAL_ACK_TIMEOUT, APP_SERIAL_MAX_RETRIES, APP_SERIAL_MAX_PROFILE_COUNT,
        TYPE_REQUEST, TYPE_RESPONSE, TYPE_EVENT, TYPE_ACK, TYPE_MASK, FLAG_FINAL,
        APP_SERIAL_ERR_COMMIT_RETRY, BUTTON_COUNT, BUTTON_DESC_SIZE,
        TINYUSB_TABLE_SIZE, PROFILE_NAME_SIZE, PROFILE_RAW_SIZE
    });
    OF.ProtocolUtils = Object.freeze({ crc8, buildFrame, decodeCString, encodeCString, concatBytes, boardArch });

})(typeof globalThis !== 'undefined' ? globalThis : window);
