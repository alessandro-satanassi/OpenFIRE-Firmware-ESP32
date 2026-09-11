/*  OpenFIRE Web App - connection controller.

    Chooses the transport and keeps the App session alive:
      - on the lightgun (OF.BUILD.target === 'device') the page connects by itself to
        ws://<host>/ws and connects again when the WebSocket drops;
      - on the site / Tauri the user picks the port (Web Serial needs a click).

        const connection = new OF.Connection();
        connection.on('state', (state, detail) => ...);   // idle | connecting | docked | lost | error
        connection.on('docked', ({ board, config }) => ...);
        connection.on('event', (command, payload) => ...); // firmware events (calibration, tests)
        connection.on('progress', (progress) => ...);
        connection.start();                                // device: auto connect
        button.onclick = () => connection.connectSerial(); // site

    options.onDockFailure(result, protocol, kind): awaited before the port of a failed dock is
    closed (the Qt App offers a reboot to bootloader while the port is still open).
    notifyClosed(reason): the page closed the link itself (bootloader, clear save memory, fatal error).

    Lightgun page: the firmware serves one App page at a time and a new page takes the link over.
    A page retries only while it is visible, and stops retrying (state 'lost', detail.stopped)
    after losing a docked session LOSS_LIMIT times within LOSS_WINDOW_MS: two open pages would
    otherwise take the gun from each other forever. reconnect() starts again.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const RECONNECT_DELAY_MS = 1500;
    const DOCK_RETRY_DELAY_MS = 2000;
    const LOSS_LIMIT = 3;
    const LOSS_WINDOW_MS = 30000;

    class Connection {
        constructor(options = {}) {
            this.protocol = options.protocol || new OF.Protocol({ debug: !!options.debug });
            this.onDockFailure = options.onDockFailure || null;
            this.state = 'idle';
            this.board = null;
            this.config = null;
            this.transportKind = null;
            this._listeners = new Map();
            this._retryTimer = null;
            this._stopped = false;
            this._busy = false;
            this._losses = [];
            this._retryWhenVisible = false;
            this.autoStopped = false;     // lightgun page: retries stopped, see reconnect()

            this.protocol.on('event', (command, payload) => this._emit('event', command, payload));
            this.protocol.on('progress', (progress) => this._emit('progress', progress));
            this.protocol.on('closed', (reason) => this._onClosed(reason, true));

            const doc = root.document;
            if (doc && doc.addEventListener) {
                doc.addEventListener('visibilitychange', () => {
                    if (doc.hidden || !this._retryWhenVisible) return;
                    this._retryWhenVisible = false;
                    if (!this._stopped && !this.autoStopped && !this._busy && this.state !== 'docked')
                        this.connectWebSocket();
                });
            }
        }

        get _pageHidden() {
            const doc = root.document;
            return !!(doc && doc.hidden);
        }

        get isDevice() { return OF.Boards ? OF.Boards.isDevice : false; }
        get serialSupported() { return OF.WebSerialTransport.isSupported(); }
        get isDocked() { return this.state === 'docked'; }

        on(name, listener) {
            if (!this._listeners.has(name)) this._listeners.set(name, new Set());
            this._listeners.get(name).add(listener);
        }

        off(name, listener) {
            const set = this._listeners.get(name);
            if (set) set.delete(listener);
        }

        _emit(name, ...args) {
            const set = this._listeners.get(name);
            if (!set) return;
            for (const listener of [...set]) {
                try {
                    listener(...args);
                } catch (error) {
                    console.error(`[Connection] '${name}' listener failed:`, error);
                }
            }
        }

        _setState(state, detail = null) {
            this.state = state;
            this._emit('state', state, detail);
        }

        /** WebSocket URL to use automatically: the lightgun page, or ?ws[=url] (development). */
        get autoWebSocketUrl() {
            const location = root.location;
            if (this.isDevice) return location ? OF.WebSocketTransport.defaultUrl(location) : null;
            const query = location && location.search ? new URLSearchParams(location.search) : null;
            if (!query || !query.has('ws')) return null;
            return query.get('ws') || OF.WebSocketTransport.defaultUrl(location);
        }

        /** Connects by itself when a WebSocket is expected (lightgun page) and keeps retrying. */
        start() {
            this._stopped = false;
            if (!this.autoWebSocketUrl || !OF.WebSocketTransport.isSupported()) return;
            // A page opened in a background tab does not take the gun from the page in use.
            if (this._pageHidden) this._retryWhenVisible = true;
            else this.connectWebSocket();
        }

        /** Lightgun page: starts again after the retries stopped (or at any time). */
        reconnect() {
            this._losses = [];
            this.autoStopped = false;
            this._stopped = false;
            if (this.state === 'docked' || this._busy) return Promise.resolve(false);
            return this.connectWebSocket();
        }

        _scheduleRetry(delay) {
            clearTimeout(this._retryTimer);
            if (this._stopped || this.autoStopped || !this.autoWebSocketUrl) return;
            this._retryTimer = setTimeout(() => {
                if (this._stopped || this.autoStopped) return;
                if (this._pageHidden) this._retryWhenVisible = true;
                else this.connectWebSocket();
            }, delay);
        }

        async connectWebSocket(url) {
            clearTimeout(this._retryTimer);
            this._retryWhenVisible = false;
            if (this._busy || this.state === 'docked') return false;
            const transport = new OF.WebSocketTransport(url || this.autoWebSocketUrl || OF.WebSocketTransport.defaultUrl(root.location));
            const ok = await this._connect(transport, 'websocket');
            if (!ok) this._scheduleRetry(this.state === 'error' ? DOCK_RETRY_DELAY_MS : RECONNECT_DELAY_MS);
            return ok;
        }

        /** Site: asks the browser for an OpenFIRE serial port (call it from a click). */
        async connectSerial(port) {
            if (this._busy) return false;
            if (!this.serialSupported) {
                this._setState('error', { error: 'serial_unsupported' });
                return false;
            }
            let selected = port;
            if (!selected) {
                try {
                    selected = await OF.WebSerialTransport.requestPort();
                } catch (error) {
                    this._setState('error', { error: 'port_request_failed', message: error.message });
                    return false;
                }
                if (!selected) return false; // chooser closed
            }
            return this._connect(new OF.WebSerialTransport(selected), 'serial');
        }

        async _connect(transport, kind) {
            this._busy = true;
            this._stopped = false;
            try {
                return await this._connectAndDock(transport, kind);
            } finally {
                this._busy = false;
            }
        }

        async _connectAndDock(transport, kind) {
            if (kind === 'serial') {
                this._stopped = true; // no automatic WebSocket retries while the user picked a port
                clearTimeout(this._retryTimer);
            }
            this._setState('connecting', { transport: kind });
            this.transportKind = kind;
            try {
                await this.protocol.connect(transport);
            } catch (error) {
                this._setState(kind === 'websocket' ? 'lost' : 'error', { error: error.code || 'open_failed', message: error.message });
                return false;
            }

            const result = await this.protocol.getSettings();
            if (!result.ok) {
                if (this.onDockFailure) {
                    try {
                        await this.onDockFailure(result, this.protocol, kind);
                    } catch (error) {
                        console.error('[Connection] onDockFailure failed:', error);
                    }
                }
                await this.protocol.disconnect().catch(() => {});
                this._setState('error', { error: result.error, rebootSuggested: result.rebootSuggested });
                return false;
            }

            this.board = result.board;
            this.config = result.config;
            this._setState('docked', { board: result.board });
            this._emit('docked', { board: result.board, config: result.config });
            return true;
        }

        /** The link was closed on purpose by the page: same as a lost link (the lightgun page reconnects). */
        notifyClosed(reason = 'closed') {
            this._onClosed(reason, false);
        }

        _onClosed(reason, unexpected) {
            if (this._busy) return; // reported by the connect flow
            const wasDocked = this.state === 'docked';
            this.board = null;
            this.config = null;
            let stopped = false;
            if (unexpected && wasDocked && this.autoWebSocketUrl && !this._stopped) {
                const now = Date.now();
                this._losses = this._losses.filter((time) => now - time < LOSS_WINDOW_MS);
                this._losses.push(now);
                if (this._losses.length >= LOSS_LIMIT) {
                    this.autoStopped = true;
                    stopped = true;
                    clearTimeout(this._retryTimer);
                }
            }
            this._setState('lost', { reason, unexpected: !!unexpected, wasDocked, stopped });
            this._scheduleRetry(RECONNECT_DELAY_MS);
        }

        /** Undocks the gun (it goes back to Run mode) and closes the link. */
        async disconnect() {
            this._stopped = true;
            clearTimeout(this._retryTimer);
            await this.protocol.disconnect().catch(() => {});
            this.board = null;
            this.config = null;
            this._setState('idle');
        }
    }

    OF.Connection = Connection;
})(typeof globalThis !== 'undefined' ? globalThis : this);
