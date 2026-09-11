/*  OpenFIRE Web App - main window (Qt: appmainwindow.cpp).

    Menu bar, device selector (site), title, the five tabs, status bar with Save; connection
    and board events. The lightgun page (OF.BUILD.target 'device') connects by itself over
    WebSocket and has no device selector or board previews.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};
    const doc = root.document;

    const THEME_KEY = 'of_theme';
    const TAB_ORDER = ['pins', 'buttons', 'settings', 'profiles', 'tests'];
    const DOCS_URL = 'https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/blob/OpenFIRE-dev/OpenFIREmain/README.md';
    const WIKI_URL = 'https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/wiki';

    const storage = {
        get(key) { try { return root.localStorage.getItem(key); } catch (e) { return null; } },
        set(key, value) { try { root.localStorage.setItem(key, value); } catch (e) { /* private mode */ } },
    };

    /** Menu texts of the Qt App carry '&' accelerators. */
    const menuText = (key) => OF.UI.t(key).replace(/&(?=[^&\s])/, '');

    class App {
        constructor() {
            this.S = OF.Boards.shared;
            this.C = this.S.serialCmdTypes_e;
            this.state = new OF.AppState(this.S);
            this.connection = new OF.Connection({ onDockFailure: (result, protocol, kind) => this.onDockFailure(result, protocol, kind) });
            this.protocol = this.connection.protocol;
            this.busy = false;            // Qt serialActive: a transaction is running
            this.irTestActive = false;
            this.window = null;           // fullscreen window
            this.caliProfile = -1;
            this.expectClose = false;     // the page closes the link itself (bootloader, clear)
            this.session = 0;             // increases at every dock and undock: stale awaits check it
            this.lostEdits = null;        // unsaved edits of a session lost by accident
            this.portArch = new Map();    // Web Serial port -> architecture of its last dock
            this.irTestStarting = false;
            this.showUnsafe = false;      // like the Qt App, not remembered
            this.theme = storage.get(THEME_KEY) || 'system';
            this.currentTab = 'pins';
            this.ports = [];
            this.activePort = null;
            this.tabs = {};
            this.titleName = '';
        }

        get isDevice() { return OF.Boards.isDevice; }
        get commitNeedsRetry() { return this.protocol.commitNeedsRetry(); }
        get connected() { return this.state.loaded && this.connection.isDocked; }

        start() {
            this.applyTheme();
            this.render();
            OF.i18n.on('change', () => this.render());

            this.connection.on('state', (state, detail) => this.onConnectionState(state, detail || {}));
            this.connection.on('docked', ({ board, config }) => this.onDocked(board, config));
            this.connection.on('event', (command, payload) => this.onEvent(command, payload));
            this.connection.on('progress', (progress) => this.onProgress(progress));

            if (root.matchMedia) {
                root.matchMedia('(prefers-color-scheme: dark)').addEventListener?.('change', () => this.applyTheme());
            }

            if (!this.isDevice && OF.WebSerialTransport.isSupported()) {
                navigator.serial.addEventListener('connect', () => this.refreshPorts());
                navigator.serial.addEventListener('disconnect', () => this.refreshPorts());
                this.refreshPorts();
            }
            this.status.show(OF.UI.t('Welcome to the OpenFIRE app!'), 3000);
            this.connection.start();
            // Unsaved edits: the browser asks before leaving the page.
            root.addEventListener('beforeunload', (event) => {
                if (this.state.isDirty() || this.commitNeedsRetry && this.state.loaded) {
                    event.preventDefault();
                    event.returnValue = '';
                }
            });
            // The gun goes back to Run mode when the page really goes away.
            root.addEventListener('pagehide', () => { if (this.connection.isDocked) this.protocol.disconnect(); });
            // Back from the browser cache the link is closed: start again.
            root.addEventListener('pageshow', (event) => { if (event.persisted) root.location.reload(); });
            // Qt App shortcuts of the Help menu.
            doc.addEventListener('keydown', (event) => {
                if (!event.altKey || event.ctrlKey || event.metaKey || event.shiftKey) return;
                if (this.window) return; // the fullscreen windows have no menu
                const url = event.code === 'KeyD' ? DOCS_URL : event.code === 'KeyS' ? WIKI_URL : null;
                if (!url) return;
                event.preventDefault();
                root.open(url, '_blank', 'noopener');
            });
        }

        // ----- Theme ----------------------------------------------------------------------

        applyTheme() {
            const html = doc.documentElement;
            if (this.theme === 'light' || this.theme === 'dark') html.dataset.theme = this.theme;
            else delete html.dataset.theme;
        }

        setTheme(theme) {
            this.theme = theme;
            storage.set(THEME_KEY, theme);
            this.applyTheme();
        }

        // ----- Page structure -----------------------------------------------------------------

        render() {
            const { el, t, icon } = OF.UI;
            const app = doc.getElementById('app');
            const settingsView = this.tabs.settings && this.tabs.settings.saveViewState();
            app.textContent = '';
            OF.UI.closeMenus();

            // Menu bar
            const menubar = el('nav', { class: 'menubar', attrs: { 'aria-label': 'Menu' } },
                el('span', { class: 'menubar-logo', html: OF.LOGO_SVG, attrs: { 'aria-hidden': 'true' } }),
                this.isDevice ? null : el('button', { class: 'menu-button', text: t('Board Previews'), on: { click: () => OF.Windows.openPreviewer(this.state.board && this.state.board.type) } }),
                el('button', { class: 'menu-button', text: t('Emitter Alignment'), on: { click: () => this.openAlignment() } }),
                OF.UI.menu(t('View'), () => [
                    { label: t('Show Unsafe Settings'), checked: this.showUnsafe, action: () => this.setUnsafe(!this.showUnsafe) },
                    { label: t('Debug Window'), action: () => this.openDebugWindow(), hidden: !this.debugEnabled },
                    { separator: true },
                    { label: t('System Theme'), checked: this.theme === 'system', action: () => this.setTheme('system') },
                    { label: t('Light Theme'), checked: this.theme === 'light', action: () => this.setTheme('light') },
                    { label: t('Dark Theme'), checked: this.theme === 'dark', action: () => this.setTheme('dark') },
                ]),
                OF.UI.menu(menuText('&Help'), () => [
                    { label: t('View Compatible Boards'), action: () => OF.Windows.openPreviewer(this.state.board && this.state.board.type), hidden: this.isDevice },
                    { label: menuText('Open &IR Emitter Alignment Assistant'), action: () => this.openAlignment() },
                    { separator: true },
                    { label: menuText('&OpenFIRE Documentation on the Repo...'), href: DOCS_URL, shortcut: 'Alt+D' },
                    { label: menuText('OpenFIRE &Serial Usage Docs on the Wiki...'), href: WIKI_URL, shortcut: 'Alt+S' },
                ]),
                el('button', { class: 'menu-button', text: t('About'), on: { click: () => OF.Windows.openAbout() } }),
                el('span', { class: 'spacer' }),
                this.languageSelector());

            // Device selector (site)
            let deviceBar = null;
            if (!this.isDevice) {
                this.portSelector = el('select', { class: 'port-selector', attrs: { 'aria-label': t('COM Port:') } });
                this.portSelector.addEventListener('change', () => this.onPortSelected());
                this.addDeviceButton = el('button', { class: 'icon-text-button', title: t('Add a Device...'), on: { click: () => this.addDevice() } },
                    icon('usb'), el('span', { text: t('Add a Device...') }));
                deviceBar = el('div', { class: 'device-bar' },
                    el('label', { class: 'device-label', text: t('COM Port:') }), this.portSelector, this.addDeviceButton);
            }

            this.titleNode = el('h1', { class: 'board-title' });
            this.versionNode = el('div', { class: 'fw-version' });
            const header = el('header', { class: 'app-header' }, menubar, deviceBar,
                el('div', { class: 'title-block' }, this.titleNode, this.versionNode));

            // Welcome (not connected)
            this.welcome = this.buildWelcome();

            // Tabs
            this.panels = {};
            const main = el('main', { class: 'content' }, this.welcome);
            const tabbar = el('nav', { class: 'tabbar', attrs: { role: 'tablist' } });
            this.tabButtons = {};
            for (const id of TAB_ORDER) {
                const module = OF.Tabs[id];
                const tab = module.build(this);
                this.tabs[id] = tab;
                const panel = el('section', { class: 'tab-panel', id: `tab-${id}`, attrs: { role: 'tabpanel' } }, tab.root);
                this.panels[id] = panel;
                main.append(panel);
                const button = el('button', { class: 'tab-button', attrs: { role: 'tab', 'aria-controls': `tab-${id}` },
                    on: { click: () => this.selectTab(id) } }, icon(module.icon), el('span', { class: 'tab-label', text: t(module.label) }));
                this.tabButtons[id] = button;
                tabbar.append(button);
            }
            this.main = main;
            this.tabbar = tabbar;

            // Status bar
            const statusText = el('span', { class: 'status-text', attrs: { role: 'status' } });
            const progress = el('progress', { class: 'status-progress', hidden: true });
            this.saveButton = el('button', { class: 'primary save-button', on: { click: () => this.save() } },
                icon('save'), el('span', { class: 'save-label' }));
            const previousStatus = this.status;
            this.status = new OF.UI.StatusBar(statusText, progress);
            this.status.adopt(previousStatus);
            const footer = el('footer', { class: 'app-footer' }, tabbar,
                el('div', { class: 'statusbar' }, statusText, progress, this.saveButton));

            app.append(header, main, footer);

            if (this.state.loaded) {
                if (settingsView) this.tabs.settings.restoreViewState(settingsView);
                this.tabs.tests.resetReadings();
                this.updateHeader();
            }
            this.renderPortSelector();
            this.selectTab(this.currentTab, true);
            this.refresh();
        }

        languageSelector() {
            const { el } = OF.UI;
            const select = el('select', { class: 'lang-selector', attrs: { 'aria-label': 'Language' } });
            for (const code of OF.i18n.languages())
                select.append(el('option', { value: code, text: OF.i18n.languageName(code) }));
            select.value = OF.i18n.currentLang;
            select.addEventListener('change', () => OF.i18n.setLanguage(select.value));
            return select;
        }

        buildWelcome() {
            const { el, t, icon } = OF.UI;
            this.welcomeText = el('p', { class: 'welcome-text' });
            this.welcomeDetail = el('p', { class: 'welcome-detail', hidden: true });
            const logo = el('div', { class: 'welcome-logo', html: OF.LOGO_SVG });
            const children = [logo, el('div', { class: 'wordmark big', text: 'OpenFIRE' }), this.welcomeText, this.welcomeDetail];
            if (this.isDevice) {
                this.spinner = el('div', { class: 'spinner', attrs: { 'aria-hidden': 'true' } });
                this.reconnectButton = el('button', { class: 'primary big-button', hidden: true, on: { click: () => this.connection.reconnect() } },
                    icon('wifi'), el('span', { text: t('Reconnect') }));
                children.push(this.spinner, this.reconnectButton);
            } else if (OF.WebSerialTransport.isSupported()) {
                this.connectButton = el('button', { class: 'primary big-button', on: { click: () => this.addDevice() } },
                    icon('usb'), el('span', { text: t('Connect') }));
                children.push(this.connectButton);
            }
            return el('section', { class: 'welcome' }, ...children);
        }

        updateWelcome() {
            const t = OF.UI.t;
            let text;
            if (this.isDevice) {
                const stopped = this.connection.autoStopped;
                text = stopped ? t('The lightgun was taken over by another page. Close the other page, then press Reconnect.') :
                    this.connection.state === 'error' ? t('Waiting for the lightgun... another App may be using it.') : t('Connecting to the lightgun...');
                this.spinner.hidden = stopped;
                this.reconnectButton.hidden = !stopped;
            } else if (!OF.WebSerialTransport.isSupported()) {
                text = t('This browser does not support Web Serial: use Chrome, Edge or Opera on a computer.');
            } else if (this.connection.state === 'connecting') {
                text = t('Connecting...');
            } else {
                text = t('Connect the lightgun with its USB cable, then select it.');
            }
            this.welcomeText.textContent = text;
            // Why the last attempt failed: useful where there is no console (phone).
            const detail = this.connection.state === 'error' ? this.lastConnectionError : null;
            this.welcomeDetail.textContent = detail ? this.connectionErrorText(detail) : '';
            this.welcomeDetail.hidden = !detail;
            if (this.connectButton) this.connectButton.disabled = this.connection.state === 'connecting';
        }

        selectTab(id, silent) {
            this.currentTab = id;
            for (const tabId of TAB_ORDER) {
                const active = tabId === id;
                this.panels[tabId].hidden = !active;
                this.tabButtons[tabId].classList.toggle('active', active);
                this.tabButtons[tabId].setAttribute('aria-selected', String(active));
            }
            const tab = this.tabs[id];
            if (tab && tab.onShow && !silent) tab.onShow();
        }

        setUnsafe(on) {
            this.showUnsafe = on;
            this.refresh();
        }

        /** Qt DiffUpdate + enable/disable rules, applied to every tab. */
        refresh() {
            const t = OF.UI.t;
            const loaded = this.state.loaded;
            this.welcome.hidden = loaded;
            this.updateWelcome();
            for (const id of TAB_ORDER) {
                const panel = this.panels[id];
                const lockedByTest = this.irTestActive && id !== 'buttons' && id !== 'tests';
                panel.classList.toggle('not-loaded', !loaded);
                panel.inert = !loaded || this.busy || lockedByTest;
                panel.classList.toggle('locked', !!lockedByTest);
                if (!loaded) panel.hidden = true;
                else panel.hidden = id !== this.currentTab;
                this.tabButtons[id].disabled = !loaded || this.busy;
                if (loaded) this.tabs[id].update();
            }

            const label = this.saveButton.querySelector('.save-label');
            let text;
            let enabled = false;
            if (!loaded) text = t('[Currently Not Connected]');
            else if (this.irTestActive) text = t('[Disabled while in Test Mode]');
            else if (this.state.isDirty() || this.commitNeedsRetry) { text = t('Save and Send Settings'); enabled = !this.busy; }
            else text = t('[Nothing To Save]');
            label.textContent = text;
            this.saveButton.disabled = !enabled;
            this.saveButton.classList.toggle('has-changes', enabled);

            if (this.portSelector) {
                this.portSelector.disabled = this.busy || this.connection.state === 'connecting' || !OF.WebSerialTransport.isSupported();
                this.addDeviceButton.disabled = this.portSelector.disabled;
            }
            if (!loaded) {
                this.titleNode.textContent = 'OpenFIRE';
                this.versionNode.textContent = '';
            }
        }

        updateHeader() {
            const board = this.state.board;
            this.titleNode.textContent = this.state.prettyName(this.titleName, OF.UI.t('Unnamed Device'));
            const version = board.version || '';
            const dash = version.indexOf('-');
            this.versionNode.textContent = '';
            if (dash > -1) {
                const org = this.state.isRP ? 'TeamOpenFIRE' : 'alessandro-satanassi';
                const hash = version.slice(dash + 1);
                this.versionNode.append('FW ', OF.UI.el('tt', null, `v${version.slice(0, dash + 1)}`,
                    OF.UI.el('a', { href: `https://github.com/${org}/OpenFIRE-Firmware/commit/${hash}`, target: '_blank', rel: 'noopener', text: hash })));
            } else if (version) {
                this.versionNode.append('FW ', OF.UI.el('tt', { text: `v${version}` }));
            }
        }

        // ----- Device selector (Web Serial) ------------------------------------------------------

        async refreshPorts() {
            if (this.isDevice || !OF.WebSerialTransport.isSupported()) return;
            try {
                this.ports = await OF.WebSerialTransport.getKnownPorts();
            } catch (error) {
                this.ports = [];
            }
            if (this.activePort && !this.ports.includes(this.activePort) && this.connection.isDocked)
                this.status.show(OF.UI.t('Current board has been disconnected.'));
            this.renderPortSelector();
        }

        renderPortSelector() {
            const select = this.portSelector;
            if (!select) return;
            const { el, t } = OF.UI;
            select.textContent = '';
            const docked = this.state.loaded;
            const first = docked ? t('[Disconnect Current Device]') :
                this.ports.length ? t('[Select a Device to Configure]') : t('[No devices currently available]');
            select.append(el('option', { value: 'none', text: first }));
            const seen = new Map();
            this.ports.forEach((port, index) => {
                let { label } = OF.WebSerialTransport.describePort(port);
                const count = (seen.get(label) || 0) + 1;
                seen.set(label, count);
                if (count > 1) label += ` #${count}`;
                select.append(el('option', { value: String(index), text: label }));
            });
            const activeIndex = this.ports.indexOf(this.activePort);
            select.value = docked && activeIndex >= 0 ? String(activeIndex) : 'none';
        }

        async onPortSelected() {
            const value = this.portSelector.value;
            if (value === 'none') {
                if (this.connection.isDocked) await this.disconnect();
                return;
            }
            const port = this.ports[Number(value)];
            if (!port || port === this.activePort && this.connection.isDocked) return;
            await this.connectPort(port);
        }

        async addDevice() {
            if (!OF.WebSerialTransport.isSupported()) {
                await OF.UI.alert(OF.UI.t('Connection failed'), OF.UI.escape(OF.UI.t('This browser does not support Web Serial: use Chrome, Edge or Opera on a computer.')), 'error');
                return;
            }
            let port;
            try {
                port = await OF.WebSerialTransport.requestPort();
            } catch (error) {
                return;
            }
            if (!port) return;
            await this.refreshPorts();
            await this.connectPort(port);
        }

        async connectPort(port) {
            if (this.connection.isDocked) await this.disconnect();
            this.activePort = port;
            this.renderPortSelector();
            const ok = await this.connection.connectSerial(port);
            if (!ok) this.activePort = null;
            this.renderPortSelector();
        }

        async disconnect() {
            this.closeWindows(false);
            await this.connection.disconnect();
            this.activePort = null;
            this.renderPortSelector();
        }

        // ----- Connection ----------------------------------------------------------------------------

        onConnectionState(state, detail) {
            const t = OF.UI.t;
            this.logNote(`state: ${state} ${JSON.stringify(detail || {})}`);
            switch (state) {
            case 'connecting':
                this.status.show(t('Connecting...'));
                this.lastConnectionError = null;
                break;
            case 'lost':
            case 'idle': {
                const wasLoaded = this.state.loaded;
                const accidental = state === 'lost' && wasLoaded && !this.expectClose;
                // Kept to be offered back when the same gun docks again.
                if (accidental && (this.state.isDirty() || this.state.forceDirty))
                    this.lostEdits = Object.assign(this.state.snapshot(), { port: this.activePort });
                else if (wasLoaded)
                    this.lostEdits = null;
                if (this.expectClose) {
                    this.expectClose = false;
                } else if (detail.stopped) {
                    this.status.show(t('The lightgun was taken over by another page. Close the other page, then press Reconnect.'));
                } else if (accidental) {
                    this.status.show(this.isDevice ? t('Connection lost: reconnecting...') : t('Current board has been disconnected.'));
                }
                // After the message: a calibration cancelled by the loss reports last, like the Qt App.
                this.onDisconnected();
                break;
            }
            case 'error':
                this.lastConnectionError = detail;
                this.expectClose = false;
                this.onDisconnected();
                this.onConnectError(detail);
                break;
            default:
                break;
            }
            this.refresh();
        }

        onProgress(progress) {
            if ('range' in progress) {
                this.status.progressRange(progress.range);
            } else {
                this.status.progress(progress.value);
                if (progress.text) this.status.show(OF.UI.t(progress.text), 5000);
            }
        }

        onDocked(board, config) {
            this.closeWindows(false);
            OF.UI.closeDialogs('session');
            this.session++;
            this.expectClose = false;
            this.irTestActive = false;
            this.irTestStarting = false;
            this.state.load(board, config);
            if (this.activePort) this.portArch.set(this.activePort, board.arch);
            this.titleName = config.tinyUSB.name;
            for (const id of TAB_ORDER) {
                const tab = this.tabs[id];
                if (tab.reload) tab.reload();
                if (tab.onLoad) tab.onLoad();
            }
            this.updateHeader();
            this.renderPortSelector();
            this.selectTab('pins');
            this.status.progressRange(0);
            this.refresh();
            if (board.cameraError) this.showCameraError('warning', false);
            if (this.retryNoticePending) {
                this.retryNoticePending = false;
                if (this.commitNeedsRetry) this.showSaveNotConfirmed();
            }
            this.offerLostEdits();
        }

        /** The gun docked again after a lost session with unsaved edits: offer them back. */
        async offerLostEdits() {
            const snapshot = this.lostEdits;
            this.lostEdits = null;
            if (!snapshot || (this.portSelector && snapshot.port !== this.activePort)) return;
            if (!this.state.snapshotDiffers(snapshot)) return;
            const t = OF.UI.t;
            const session = this.session;
            const yes = await OF.UI.confirm(t('Unsaved Changes'), t('The connection was lost before your changes were saved.'),
                t('Restore the changes made before the connection was lost?\n\nCalibrations stay in the board only until it restarts: save them to keep them.'),
                { icon: 'warning', scope: 'session' });
            if (!yes || session !== this.session || !this.state.loaded) return;
            if (!this.state.restoreSnapshot(snapshot)) return;
            for (const id of TAB_ORDER) {
                const tab = this.tabs[id];
                if (tab.reload) tab.reload();
            }
            this.status.show(t('Changes restored: press Save and Send Settings to keep them.'), 5000);
            this.refresh();
        }

        onDisconnected() {
            this.closeWindows(false); // the alignment assistant does not need the board
            if (this.state.loaded) {
                OF.UI.closeDialogs('session');
                this.session++;
            }
            this.irTestActive = false;
            this.irTestStarting = false;
            this.retryNoticePending = false;
            this.busy = false;
            this.state.reset();
            this.status.progressRange(0);
            if (!this.connection.isDocked && this.portSelector) {
                if (this.connection.state !== 'connecting') this.activePort = null;
                this.renderPortSelector();
            }
        }

        /** What went wrong, in the words of the App protocol (js/core/protocol.js results). */
        connectionErrorText(detail) {
            const t = OF.UI.t;
            switch (detail.error) {
            case 'dock_timeout':
                return t('The lightgun did not answer the connection request: it may be in use by an App on the USB cable or by another page.');
            case 'bad_board_info':
                return t('The lightgun answered with data this App does not understand: check that its firmware and this App belong to the same version.');
            case 'not_open':
            case 'open_failed':
            case 'closed':
                return t('The link to the lightgun closed before the settings were read.');
            default:
                return t('The settings could not be read (%1): the connection dropped during the transfer.', String(detail.error || '?'));
            }
        }

        onConnectError(detail) {
            const t = OF.UI.t;
            const esc = OF.UI.escape;
            if (this.isDevice) {
                // The lightgun page keeps trying: another App (USB) may hold the session.
                this.status.show(t('Waiting for the lightgun... another App may be using it.'));
                this.logNote('connect error: ' + JSON.stringify(detail));
                return;
            }
            this.activePort = null;
            this.renderPortSelector();
            switch (detail.error) {
            case 'port_busy':
            case 'open_failed':
                OF.UI.alert(t('Serial port is already in use!'), esc(t("This usually indicates that the port is being used by something else, e.g. Arduino IDE's serial monitor, or another command line app (stty, screen).\n\nPlease close the offending application and try selecting this port again.")), 'warning');
                break;
            case 'serial_unsupported':
            case 'port_request_failed':
                OF.UI.alert(t('Connection failed'), esc(t('This browser does not support Web Serial: use Chrome, Edge or Opera on a computer.')), 'error');
                break;
            case 'dock_timeout':
            case 'bad_board_info':
                break; // handled by onDockFailure while the port was open
            default:
                this.status.show(t('Connection failed'), 5000);
                break;
            }
        }

        /** Qt AppSerial::GetSettings failure paths: stale dock warning and RequestToReboot. */
        async onDockFailure(result, protocol, kind) {
            if (kind !== 'serial' || !protocol.isOpen) return; // unplugged: nothing to ask
            const t = OF.UI.t;
            if (result.error === 'dock_timeout') {
                await OF.UI.alert(t("Data hasn't arrived! (Stale state?)"), OF.UI.escape(t("Device was detected, but initial settings request wasn't received in time!\nThis can happen if the app was unexpectedly closed and the gun is in a stale docked state.\n\nTry selecting the device again.")), 'warning');
            }
            if ((result.error === 'dock_timeout' || result.error === 'bad_board_info') && protocol.isOpen) {
                const arch = result.board ? result.board.arch : (this.portArch.get(this.activePort) || '');
                const isRP = arch === this.S.boardArchs[this.S.boardArchs_e.boardRP];
                const title = isRP ? t('Reset Board to Bootloader?') : t('Reboot Microcontroller?');
                const text = isRP ?
                    t("<p>The board you selected did not respond to the app properly.</p><p>This can usually be resolved by rebooting the microcontroller to its bootloader, and then updating the board to the latest firmware, which can be found at:</p><p><a href='https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/releases/latest'><span style=' text-decoration: underline; color:#8ab4f8;'>https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/releases/latest</span></a></p><p>Would you like to reboot this board to apply an update?</p>") :
                    t("<p>The board you selected did not respond to the app properly.</p><p>This can usually be resolved by rebooting the microcontroller, and then updating the board to the latest firmware, which can be found at:</p><p><a href='https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/releases/latest'><span style=' text-decoration: underline; color:#8ab4f8;'>https://github.com/TeamOpenFIRE/OpenFIRE-Firmware/releases/latest</span></a></p><p>Would you like to reboot this board to apply an update?</p>");
                if (await OF.UI.confirm(title, null, null, { html: text, icon: 'error' }) && protocol.isOpen)
                    await protocol.rebootToBootloader(arch);
            }
        }

        // ----- Board events (Qt serialPort_readyRead) ----------------------------------------------------

        onEvent(command, payload) {
            const C = this.C;
            const t = OF.UI.t;
            this.logDebug(command, payload);
            if (this.busy && command !== C.sError) return;
            if (!this.state.loaded && command !== C.sError) return;

            if (this.tabs.tests && this.tabs.tests.onEvent(command, payload)) return;

            switch (command) {
            case C.sCurrentProf:
                if (payload.length === 1) {
                    const selection = payload[0];
                    if (selection !== this.state.cur.selectedProfile && selection < this.state.profileCount)
                        this.state.setSelectedProfile(selection);
                    this.refresh();
                }
                break;

            case C.sError:
                this.onError(payload);
                break;

            case C.sCaliStageUpd:
                if (payload.length === 1 && this.window && this.window.mode === OF.FullscreenWindow.MODE_CALIBRATE)
                    this.window.setStage(payload[0]);
                break;

            case C.sCaliInfoUpd:
                if (payload.length === 5 && this.window && this.window.mode === OF.FullscreenWindow.MODE_CALIBRATE)
                    this.window.setInfo(payload);
                break;

            case C.sTestCoords:
                if (this.window) this.window.drawTest(payload);
                break;

            case C.sClearFlash:
                this.showClearFlashDone();
                break;

            default:
                break;
            }
            void t;
        }

        onError(payload) {
            const C = this.C;
            const t = OF.UI.t;
            const retry = OF.ProtocolConstants.APP_SERIAL_ERR_COMMIT_RETRY;
            if (payload.length && payload[0] === retry) {
                // A queued error must not undo a later successful retry; a link that is gone offers
                // the edits back when the board docks again.
                if (!this.commitNeedsRetry || !this.protocol.isOpen) return;
                if (!this.state.loaded) {
                    this.retryNoticePending = true; // the board reported a pending save while docking
                    return;
                }
                this.showSaveNotConfirmed();
                return;
            }
            if (!payload.length) {
                this.closeWindows(true);
                this.expectClose = true;
                this.connection.notifyClosed('operation_failed');
                this.status.show(t('Operation not confirmed: restart the microcontroller and reconnect.'));
                this.showError(t('Operation not confirmed'), OF.UI.escape(t('The operation failed or its result could not be confirmed. Restart the microcontroller and reconnect before continuing. If saving had already started, verify the settings and save them again after resolving the error.')), 'error');
                return;
            }
            switch (payload[0]) {
            case C.sErrCam:
                this.closeWindows(true);
                this.showCameraError('error', true);
                break;
            case C.sErrPeriphGeneric:
                this.showError(t('Peripheral Device Error!'), t("<p>Data received from the board indicates that an I2C peripheral device failed to initialize.<br>" +
                    "This can happen if, for example, the peripheral's wires are crossed<br>" +
                    '(data wire to clock pin, clock wire to data pin),<br>' +
                    'or the pins for the peripheral are set to a different component,<br>' +
                    'such as a button or Force Feedback output.</p>' +
                    '<p>Confirm that the wires for the peripheral are connected to the correct <i>Peripheral I2C</i> pins<br>' +
                    'in the <i>Boards Layout</i> tab.</p>'), 'error');
                break;
            default:
                break;
            }
        }

        showSaveNotConfirmed() {
            const t = OF.UI.t;
            this.refresh();
            this.status.show(t('Save not confirmed. Keep the board powered and retry Save.'));
            this.showError(t('Save not confirmed'), OF.UI.escape(t('Saving failed or its result could not be confirmed. Your edits remain in this App; completed calibration remains in the powered board.\n\n' +
                'Check the connection or the reported storage error, then press Save and Send Settings again. No restart is required for a recoverable save error. Hardware tests stay disabled until a save is confirmed.\n\n' +
                'Disconnecting/reloading can replace unsent App edits with the board\'s partly updated settings. Restarting loses unsaved board RAM.')), 'error', 'session');
        }

        /** Qt AppSerial::ShowError: one board error box at a time, repeated reports are ignored. */
        showError(title, html, iconName, scope) {
            if (this.errorShown) return;
            this.errorShown = true;
            OF.UI.alert(title, html, iconName, { scope }).finally(() => { this.errorShown = false; });
        }

        showCameraError(iconName, testsUnavailable) {
            const t = OF.UI.t;
            const text = testsUnavailable ?
                t('<p>Data received from the board indicates that the camera is in a bad state.<br>' +
                    'This can happen if, for example, the camera wires are crossed<br>' +
                    '(data wire to clock pin, clock wire to data pin),<br>' +
                    'or the camera pins are wired to a different component,<br>' +
                    'such as a button or Force Feedback output.</p>' +
                    '<p>You are able to change the camera pins in the <i>Boards Layout</i> tab<br>' +
                    'if they should be mapped different GPIO;<br>' +
                    'Otherwise, the camera wires must be resoldered to resolve this error.</p>' +
                    '<p>IR Testing and Calibration will not be available while in this state.</p>') :
                t('<p>Data received from the board indicates that the camera is in a bad state.<br>' +
                    'This can happen if, for example, the camera wires are crossed<br>' +
                    '(data wire to clock pin, clock wire to data pin),<br>' +
                    'or the camera pins are wired to a different component,<br>' +
                    'such as a button or Force Feedback output.</p>' +
                    '<p>You are able to change the camera pins in the <i>Boards Layout</i> tab<br>' +
                    'if they should be mapped different GPIO;<br>' +
                    'Otherwise, the camera wires must be resoldered to resolve this error.</p>');
            this.showError(t('Device Error: Camera not available!'), text, iconName);
        }

        // ----- Debug window (Qt AppDebugWindow, non-Release builds) -------------------------------------

        /** Available from the unbundled webapp folder or with ?debug in the address. */
        get debugEnabled() {
            if (!OF.BUILD) return true;
            try { return new URLSearchParams(root.location.search).has('debug'); } catch (e) { return false; }
        }

        openDebugWindow() {
            const { el, t, icon } = OF.UI;
            if (this.debugWindow) {
                this.debugWindow.root.hidden = false;
                return;
            }
            const text = el('textarea', { class: 'debug-text mono', readOnly: true, attrs: { 'aria-label': 'Text' } });
            const hex = el('textarea', { class: 'debug-text mono', readOnly: true, attrs: { 'aria-label': 'Hex' } });
            const panel = el('section', { class: 'debug-window', attrs: { role: 'dialog', 'aria-label': t('Serial Debug') } },
                el('div', { class: 'debug-head' },
                    el('strong', { text: t('Serial Debug') }),
                    el('span', { class: 'spacer' }),
                    el('button', { text: t('Clear'), on: { click: () => { text.value = ''; hex.value = ''; } } }),
                    el('button', { class: 'icon-button', title: t('Close'), on: { click: () => { panel.hidden = true; } } }, icon('close'))),
                text, hex);
            doc.body.append(panel);
            this.debugWindow = { root: panel, text, hex };
        }

        /** Line written in the debug window (unbundled folder or ?debug), without a board event. */
        logNote(text) {
            const win = this.debugWindow;
            if (!win) return;
            win.text.value = (win.text.value + `[${new Date().toLocaleTimeString()}] ${text}\n`).slice(-40000);
            win.text.scrollTop = win.text.scrollHeight;
        }

        logDebug(command, payload) {
            const win = this.debugWindow;
            if (!win) return;
            const bytes = Array.from(payload || []);
            const append = (area, value) => {
                area.value = (area.value + value).slice(-40000);
                area.scrollTop = area.scrollHeight;
            };
            append(win.text, String.fromCharCode(...bytes) + '\n');
            append(win.hex, ` [${command.toString(16).toUpperCase().padStart(2, '0')}]` +
                bytes.map((b) => ' ' + b.toString(16).toUpperCase().padStart(2, '0')).join(''));
        }

        // ----- Commands ----------------------------------------------------------------------------------

        async sendCommand(command, payload) {
            if (!this.connection.isDocked) return false;
            return this.protocol.sendCommand(command, payload ? Uint8Array.from(payload) : null);
        }

        selectProfile(profile) {
            this.sendCommand(this.C.sCaliProfile, [this.C.sCaliProfile, profile]);
        }

        calibrate(profile) {
            const C = this.C;
            if (!this.connected || this.window && this.window.mode !== OF.FullscreenWindow.MODE_ALIGNMENT) return;
            const settings = this.state.cur.profiles[profile];
            if (this.window) this.window.shutdown();
            this.caliProfile = profile;
            const win = new OF.FullscreenWindow(OF.FullscreenWindow.MODE_CALIBRATE, {
                onExitRequest: () => this.sendCommand(C.serialTerminator).then((ok) => { if (!ok) this.failOperation(); }),
                onExit: (mode, values) => this.onWindowExit(win, mode, values),
            });
            this.window = win;
            win.open();
            this.sendCommand(C.sCaliProfile, [C.sCaliStart, profile, (settings.irSensitivity + (settings.layoutType << 4)) & 0xFF])
                .then((ok) => { if (!ok) this.failOperation(); });
        }

        failOperation() {
            if (this.connection.isDocked) this.protocol.failOperation();
        }

        openAlignment() {
            if (this.window && this.window.mode !== OF.FullscreenWindow.MODE_ALIGNMENT) return;
            if (this.window) this.window.shutdown();
            const win = new OF.FullscreenWindow(OF.FullscreenWindow.MODE_ALIGNMENT, { onExit: (mode) => this.onWindowExit(win, mode) });
            this.window = win;
            win.open();
        }

        async openIRTest() {
            const t = OF.UI.t;
            if (!this.connected || this.irTestStarting || this.irTestActive) return; // double click
            if (this.window && this.window.mode !== OF.FullscreenWindow.MODE_ALIGNMENT) return;
            const session = this.session;
            this.irTestStarting = true;
            try {
                if (this.state.layoutChangedUnsaved()) {
                    const proceed = await OF.UI.confirm(t('Warning: Unsaved Changes'),
                        t("Your current calibration profile has an unsaved IR Layout type change.\nTest Mode relies on the profile's settings since the last save, and may not work as expected.\nIt is recommended that you save your changes before continuing.\n\nContinue to Test Mode?"),
                        null, { icon: 'info', scope: 'session' });
                    if (!proceed || session !== this.session) return;
                }
                if (!await this.sendCommand(this.C.sIRTest, [1]) || session !== this.session) return;
            } finally {
                if (session === this.session) this.irTestStarting = false;
            }
            if (this.window) this.window.shutdown();
            this.irTestActive = true;
            const win = new OF.FullscreenWindow(OF.FullscreenWindow.MODE_IRTEST, { onExit: (mode) => this.onWindowExit(win, mode) });
            this.window = win;
            win.open();
            this.refresh();
        }

        /** Qt CaliWindowExiting. */
        async onWindowExit(win, mode, values) {
            const t = OF.UI.t;
            if (this.window === win) this.window = null;
            if (mode === OF.FullscreenWindow.MODE_CALIBRATE) {
                const profile = this.caliProfile;
                if (!this.state.loaded || profile < 0) return;
                const empty = { topOffset: -1, bottomOffset: -1, leftOffset: -1, rightOffset: -1, TLled: -1, TRled: -1 };
                const result = this.state.applyCalibration(profile, values || empty);
                const n = profile + 1;
                if (result === 'ok') this.status.show(t('Calibration for Profile ') + n + t(' successful'), 5000);
                else if (result === 'malformed') this.status.show(t('Calibration for Profile ') + n + t(' returned with malformed values.'), 10000);
                else this.status.show(t('Cancelled Calibration for Profile ') + n, 10000);
                this.refresh();
            } else if (mode === OF.FullscreenWindow.MODE_IRTEST) {
                if (!this.connection.isDocked) {
                    this.irTestActive = false;
                    this.refresh();
                } else if (await this.sendCommand(this.C.sIRTest, [0])) {
                    this.irTestActive = false;
                    this.refresh();
                } else {
                    // Test mode could not be left: the board must be restarted (fatal error dialog).
                    this.failOperation();
                }
            }
        }

        /** Closes the fullscreen windows (the alignment assistant stays unless everything closes). */
        closeWindows(all) {
            // Like Qt Shutdown(): the window reports its exit (cancelled calibration, IR test mode off).
            if (this.window && (all || this.window.mode !== OF.FullscreenWindow.MODE_ALIGNMENT))
                this.window.shutdown();
            if (all) this.irTestActive = false;
        }

        async rebootToBootloader() {
            const t = OF.UI.t;
            if (!this.connected) return;
            const rp = this.state.isRP;
            this.expectClose = true;
            await this.protocol.rebootToBootloader(this.state.board.arch);
            this.connection.notifyClosed('bootloader');
            this.status.show(rp ? t('Board reset to bootloader.') : t('Board restarted.'), 5000);
        }

        async clearSaveMemory() {
            const t = OF.UI.t;
            if (!this.connected) return;
            const session = this.session;
            const yes = await OF.UI.confirm(t('Delete Confirmation'), t('Really delete saved data?'),
                t('This operation will delete all saved data, including:\n\n - Calibration Profiles\n - Toggles\n - Settings\n - Custom Identifiers\n\nAre you sure about this?'),
                { icon: 'warning', defaultNo: true, danger: true, scope: 'session' });
            if (session !== this.session || !this.connected) return;
            if (!yes) {
                this.status.show(t('Clear operation canceled.'), 3000);
                return;
            }
            this.status.show(t('Board reset to initial settings.'));
            this.expectClose = true;
            this.clearFlashShown = false;
            const sent = await this.protocol.clearSaveMemory();
            this.connection.notifyClosed('clear_flash');
            if (sent) this.showClearFlashDone();
        }

        showClearFlashDone() {
            if (this.clearFlashShown) return;
            this.clearFlashShown = true;
            const t = OF.UI.t;
            OF.UI.alert(t('Successfully reset board settings'), OF.UI.escape(t('Please unplug the board and reinsert it into the PC.')), 'info');
        }

        // ----- Save (Qt on_confirmButton_clicked) ------------------------------------------------------

        async save() {
            const t = OF.UI.t;
            if (!this.connected || this.busy || this.irTestActive) return;
            const session = this.session;
            const yes = await OF.UI.confirm(t('Commit Confirmation'), t('Are these settings okay?'),
                t('These settings will be committed to your lightgun. Is that okay?'), { icon: 'info', scope: 'session' });
            if (session !== this.session || !this.connected || this.busy) return;
            if (!yes) {
                this.status.show(t('Save operation canceled.'), 3000);
                return;
            }
            const cameraChanged = this.state.cameraChanged;
            // The edits being sent: later changes (none while busy) are not marked as saved.
            const sent = this.state.clone(this.state.cur);
            this.busy = true;
            this.refresh();
            let result;
            try {
                result = await this.protocol.commitSettings(sent, this.state.orig);
            } catch (error) {
                console.error(error);
                result = { ok: false };
            }
            if (session !== this.session) return; // the session ended during the save (onDisconnected reset the UI)
            this.busy = false;
            this.status.progressRange(0);
            if (result.ok && this.state.loaded) {
                this.status.show(t('Sent settings successfully!'), 5000);
                this.state.commitDone();
                this.titleName = this.state.cur.tinyUSB.name;
                this.updateHeader();
                this.tabs.tests.resetReadings();
                if (cameraChanged) {
                    OF.UI.dialog({
                        title: t('Camera Changed'), text: t('Camera Model Changed Successfully'), icon: 'info',
                        info: t('The lightgun has automatically swapped the active sensors.\n\nIf you experience any unexpected issues or the new camera fails to respond, we recommend physically restarting the microcontroller.'),
                        buttons: [{ label: t('OK'), value: true, primary: true }],
                    });
                }
            }
            this.refresh();
        }
    }

    OF.App = App;

    if (doc && !root.OF_NO_AUTOSTART) {
        const boot = () => {
            OF.app = new App();
            OF.app.start();
        };
        if (doc.readyState === 'loading') doc.addEventListener('DOMContentLoaded', boot);
        else boot();
    }
})(typeof globalThis !== 'undefined' ? globalThis : this);
