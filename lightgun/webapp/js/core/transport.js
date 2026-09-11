/*  OpenFIRE Web App: byte transports for the App serial protocol.

    A transport only moves raw bytes. Framing, acknowledgements and every
    other protocol rule live in protocol.js, so the same protocol code runs
    unchanged over Web Serial (GitHub Pages / local file / Tauri) and over the
    WebSocket served by the lightgun itself (web configuration mode).

    Common interface:
        kind                    'webserial' | 'websocket'
        isOpen                  boolean
        onData(Uint8Array)      callback, set by the owner
        onClose(reason)         callback, invoked once when the link is lost or closed
        open()                  Promise<void>, rejects with an Error carrying .code
        write(Uint8Array)       Promise<boolean>, writes are strictly ordered
        close()                 Promise<void>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
*/

(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    function TransportError(code, message, cause) {
        const error = new Error(message);
        error.code = code;
        if (cause !== undefined)
            error.cause = cause;
        return error;
    }

    class Transport {
        constructor(kind) {
            this.kind = kind;
            this.isOpen = false;
            this.onData = null;
            this.onClose = null;
            this._writeChain = Promise.resolve();
            this._closeNotified = false;
        }

        // Writes are chained so that an ACK queued by the parser can never
        // overtake, or interleave with, a request written just before it.
        write(bytes) {
            const result = this._writeChain.then(() => {
                if (!this.isOpen)
                    return false;
                return this._write(bytes);
            }).catch((error) => {
                console.warn(`[Transport:${this.kind}] Write failed:`, error);
                return false;
            });

            this._writeChain = result.then(() => undefined, () => undefined);
            return result;
        }

        _emitData(bytes) {
            if (this.onData && bytes && bytes.length)
                this.onData(bytes);
        }

        _emitClose(reason) {
            this.isOpen = false;
            if (this._closeNotified)
                return;
            this._closeNotified = true;
            if (this.onClose)
                this.onClose(reason);
        }
    }

    // ===== Web Serial ========================================================

    class WebSerialTransport extends Transport {
        // For compatibility with whitelists, OpenFIRE keeps the VID fixed.
        static get VENDOR_ID() { return 0xF143; }

        static isSupported() {
            return typeof navigator !== 'undefined' && !!navigator.serial;
        }

        static filters() {
            return [{ usbVendorId: WebSerialTransport.VENDOR_ID }];
        }

        /// Opens the browser port chooser (must be called from a user gesture).
        /// Resolves to null when the user dismisses the chooser.
        static async requestPort() {
            if (!WebSerialTransport.isSupported())
                throw TransportError('unsupported', 'Web Serial API is not available in this browser.');
            try {
                return await navigator.serial.requestPort({ filters: WebSerialTransport.filters() });
            } catch (error) {
                if (error && error.name === 'NotFoundError')
                    return null;
                throw TransportError('request_failed', 'Serial port request failed.', error);
            }
        }

        /// Ports already authorised by the user on this origin (no prompt).
        static async getKnownPorts() {
            if (!WebSerialTransport.isSupported())
                return [];
            const ports = await navigator.serial.getPorts();
            return ports.filter((port) => {
                const info = port.getInfo ? port.getInfo() : {};
                return info.usbVendorId === WebSerialTransport.VENDOR_ID;
            });
        }

        static describePort(port) {
            const info = port && port.getInfo ? port.getInfo() : {};
            const hex = (value) => '0x' + (value >>> 0).toString(16).toUpperCase().padStart(4, '0');
            return {
                vendorId: info.usbVendorId,
                productId: info.usbProductId,
                label: info.usbProductId !== undefined ?
                    `OpenFIRE (PID ${hex(info.usbProductId)})` : 'OpenFIRE'
            };
        }

        constructor(port, options = {}) {
            super('webserial');
            this.port = port;
            this.baudRate = options.baudRate || 9600;
            this._reader = null;
            this._writer = null;
            this._readLoop = null;
            this._closing = false;
            this._onDisconnect = (event) => {
                if (event.target === this.port || event.port === this.port)
                    this._lost('device_lost');
            };
        }

        async open() {
            if (this.isOpen)
                return;

            this._closing = false;
            this._closeNotified = false;

            try {
                await this.port.open({ baudRate: this.baudRate, bufferSize: 4096 });
            } catch (error) {
                // Chromium reports a port held by another program (Arduino IDE,
                // the Qt App, a terminal) as a NetworkError / InvalidStateError.
                const busy = error && (error.name === 'NetworkError' || error.name === 'InvalidStateError');
                throw TransportError(busy ? 'port_busy' : 'open_failed',
                    busy ? 'The serial port is already in use by another program.' :
                           'The serial port could not be opened.', error);
            }

            // Same as QSerialPort::setDataTerminalReady(true) in the Qt App.
            try {
                await this.port.setSignals({ dataTerminalReady: true });
            } catch (error) {
                console.warn('[Transport:webserial] Could not assert DTR:', error);
            }

            this._writer = this.port.writable.getWriter();
            this.isOpen = true;

            if (navigator.serial && navigator.serial.addEventListener)
                navigator.serial.addEventListener('disconnect', this._onDisconnect);

            this._readLoop = this._runReadLoop();
        }

        async _runReadLoop() {
            while (this.isOpen && !this._closing && this.port.readable) {
                this._reader = this.port.readable.getReader();
                try {
                    for (;;) {
                        const { value, done } = await this._reader.read();
                        if (done)
                            break;
                        if (value)
                            this._emitData(value);
                    }
                } catch (error) {
                    // Framing/parity/overrun/break errors are recoverable:
                    // a new reader is created and the protocol resynchronises.
                    const recoverable = error && ['BufferOverrunError', 'BreakError', 'FramingError', 'ParityError']
                        .includes(error.name);
                    if (!recoverable && !this._closing) {
                        console.warn('[Transport:webserial] Read loop stopped:', error);
                        try { this._reader.releaseLock(); } catch (e) { /* already released */ }
                        this._reader = null;
                        this._lost('device_lost');
                        return;
                    }
                } finally {
                    if (this._reader) {
                        try { this._reader.releaseLock(); } catch (e) { /* already released */ }
                        this._reader = null;
                    }
                }

                if (this._closing)
                    break;
            }

            if (!this._closing && this.isOpen)
                this._lost('device_lost');
        }

        async _write(bytes) {
            if (!this._writer)
                return false;
            await this._writer.write(bytes);
            return true;
        }

        async _releaseAndClose() {
            this._closing = true;
            this.isOpen = false;

            if (navigator.serial && navigator.serial.removeEventListener)
                navigator.serial.removeEventListener('disconnect', this._onDisconnect);

            if (this._reader) {
                try { await this._reader.cancel(); } catch (e) { /* port already gone */ }
            }
            if (this._readLoop) {
                try { await this._readLoop; } catch (e) { /* ignored */ }
                this._readLoop = null;
            }
            if (this._writer) {
                try { this._writer.releaseLock(); } catch (e) { /* ignored */ }
                this._writer = null;
            }
            try { await this.port.close(); } catch (e) { /* already closed or unplugged */ }
        }

        _lost(reason) {
            if (this._closing)
                return;
            this._releaseAndClose().finally(() => this._emitClose(reason));
        }

        async close() {
            if (this._closing && !this.isOpen)
                return;
            // Let an already queued write (e.g. the final disconnect ACK) leave first.
            try { await this._writeChain; } catch (e) { /* ignored */ }
            await this._releaseAndClose();
            this._emitClose('closed');
        }

        /// RP2040 "magic baud" bootloader entry, equivalent to the Qt App:
        /// set 1200 baud and drop DTR. Web Serial cannot change the baud rate
        /// of an open port, so the port is reopened at 1200 baud.
        async touch1200() {
            if (this.isOpen)
                await this.close();

            try {
                await this.port.open({ baudRate: 1200 });
                try {
                    await this.port.setSignals({ dataTerminalReady: false });
                } catch (e) { /* some platforms reset on open already */ }
            } catch (error) {
                throw TransportError('open_failed', 'Could not reopen the port at 1200 baud.', error);
            } finally {
                try { await this.port.close(); } catch (e) { /* the board may already be gone */ }
            }
        }
    }

    // ===== WebSocket (lightgun web configuration mode) =======================

    class WebSocketTransport extends Transport {
        static isSupported() {
            return typeof WebSocket !== 'undefined';
        }

        /// ws://<host>/ws for pages served by the lightgun itself.
        static defaultUrl(location) {
            const loc = location || (typeof window !== 'undefined' ? window.location : null);
            if (!loc || !loc.host)
                return null;
            return (loc.protocol === 'https:' ? 'wss://' : 'ws://') + loc.host + '/ws';
        }

        constructor(url, options = {}) {
            super('websocket');
            this.url = url;
            this.connectTimeoutMs = options.connectTimeoutMs || 4000;
            this.socket = null;
        }

        open() {
            if (this.isOpen)
                return Promise.resolve();

            this._closeNotified = false;

            return new Promise((resolve, reject) => {
                let settled = false;
                let socket;

                try {
                    socket = new WebSocket(this.url);
                } catch (error) {
                    reject(TransportError('open_failed', 'Invalid WebSocket address.', error));
                    return;
                }

                socket.binaryType = 'arraybuffer';
                this.socket = socket;

                const timer = setTimeout(() => {
                    if (settled)
                        return;
                    settled = true;
                    try { socket.close(); } catch (e) { /* ignored */ }
                    reject(TransportError('timeout', 'WebSocket connection timed out.'));
                }, this.connectTimeoutMs);

                socket.onopen = () => {
                    if (settled)
                        return;
                    settled = true;
                    clearTimeout(timer);
                    this.isOpen = true;
                    resolve();
                };

                socket.onerror = (event) => {
                    if (settled)
                        return;
                    settled = true;
                    clearTimeout(timer);
                    reject(TransportError('open_failed', 'WebSocket connection failed.', event));
                };

                socket.onmessage = (event) => {
                    if (event.data instanceof ArrayBuffer) {
                        this._emitData(new Uint8Array(event.data));
                    } else if (typeof Blob !== 'undefined' && event.data instanceof Blob) {
                        event.data.arrayBuffer().then((buffer) => this._emitData(new Uint8Array(buffer)));
                    }
                };

                socket.onclose = () => {
                    clearTimeout(timer);
                    if (!settled) {
                        settled = true;
                        reject(TransportError('open_failed', 'WebSocket closed before opening.'));
                        return;
                    }
                    this.socket = null;
                    this._emitClose(this.isOpen ? 'connection_lost' : 'closed');
                };
            });
        }

        async _write(bytes) {
            if (!this.socket || this.socket.readyState !== WebSocket.OPEN)
                return false;
            this.socket.send(bytes);
            return true;
        }

        async close() {
            try { await this._writeChain; } catch (e) { /* ignored */ }
            const socket = this.socket;
            this.isOpen = false;
            if (socket) {
                socket.onclose = null;
                try { socket.close(); } catch (e) { /* ignored */ }
                this.socket = null;
            }
            this._emitClose('closed');
        }
    }

    OF.TransportError = TransportError;
    OF.Transport = Transport;
    OF.WebSerialTransport = WebSerialTransport;
    OF.WebSocketTransport = WebSocketTransport;

})(typeof globalThis !== 'undefined' ? globalThis : window);
