/*  OpenFIRE Web App - the version check of a published App.

    Every version of the App stays published, each in its own folder (v/<version>/ of the
    site), and the home page of the site is the one that chooses: it reads the version of
    the firmware and opens the App published for it (webapp/launcher/). An App, once open,
    never goes looking for another one - it is the App of its version and it stays that.

    All it does is this: if the lightgun runs a firmware of another version, say so and let
    the user decide whether to carry on or to go back to the home page, which opens the
    right one. That is the whole of this file.

      afterDock(app, board)  once docked: the question, when the versions differ.
      resume(app)            arrived from the home page: take the port back, say why.

    The App embedded in the lightgun does not include this file at all
    (scripts/webapp_build.py leaves it out of the device build): the gun serves the App of
    its own firmware, so there is nothing to check, and every call to OF.Version elsewhere
    is guarded.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const JUMP_KEY = 'of_version_jump';    // left by the home page for the App it opens

    const store = {
        get(key) { try { return root.sessionStorage.getItem(key); } catch (e) { return null; } },
        drop(key) { try { root.sessionStorage.removeItem(key); } catch (e) { /* private mode */ } },
    };

    /** The folder of this page, without the file name. */
    const pageFolder = (path) => String(path || '').replace(/[^/]*$/, '');

    /** True when this page is one of the published versions, v/<version>/ of a site: then
        there is a home page to go back to. A copy opened from anywhere else (a folder on
        the disk, the desktop app) has nowhere to send anybody. */
    const isPublished = (path) => /\/v\/[^/]+\/$/.test(pageFolder(path));

    /** The home page of the site this version belongs to. */
    const siteHome = (path) => {
        const folder = pageFolder(path);
        const match = folder.match(/^(.*\/)v\/[^/]+\/$/);
        return match ? match[1] : folder;
    };

    const here = () => (root.location && root.location.pathname) || '';
    const isSite = () => !!(OF.BUILD && OF.BUILD.target === 'site');
    const myVersion = () => String((OF.BUILD && OF.BUILD.version) || '');
    const myLabel = () => String((OF.BUILD && (OF.BUILD.versionLabel || OF.BUILD.version)) || '');
    const langSuffix = () => {
        const lang = OF.i18n ? OF.i18n.currentLang : '';
        return lang ? '?lang=' + encodeURIComponent(lang) : '';
    };

    /** The short version of a firmware: "6.2" out of "6.2-<git hash>". It is all a
        firmware before 7.0 says about itself. */
    const shortOf = (board) => String((board && board.version) || '').split('-')[0].trim();

    /** How this firmware names itself: the complete version, 7.0.0 or 7.0.0-beta1. */
    const idOf = (board) => String((board && board.versionFull) || '').trim() || shortOf(board);

    /** The names under which this same firmware may have been published: its complete
        version, and the short one - all a firmware before 7.0 could say about itself, and
        how those apps were archived (v/6.2/). Both mean "this very firmware", so neither
        is a mismatch. What is NOT in here is the version without the suffix: 7.0.0 and
        7.0.0-beta1 are two different firmwares, and the app of one says so to the other. */
    const namesOf = (board) => [idOf(board), shortOf(board)].filter(Boolean);

    OF.Version = {
        isPublished, siteHome,

        /** Once docked: when the firmware is of another version, ask what to do. Returns
            true when the question was asked. */
        async afterDock(app, board) {
            if (!isSite()) return false;
            const mine = myVersion();
            const theirs = idOf(board);
            if (!mine || !theirs || namesOf(board).indexOf(mine) >= 0) return false;
            if (this.asked === theirs) return false;  // once per version, not at every dock
            this.asked = theirs;

            const t = OF.UI.t;
            const text = t('This is the App of version %1, while the lightgun runs firmware %2.', myLabel(), theirs);

            // Nowhere to go back to: this copy is not part of a published site.
            if (!isPublished(here())) {
                await OF.UI.alert(t('Versions do not match'), OF.UI.escape(text +
                    '\n\n' + t('They can still work together, but something may not match.')), 'warning');
                return true;
            }

            const carryOn = await OF.UI.dialog({
                title: t('Versions do not match'),
                icon: 'warning',
                text,
                info: t('They can still work together, but something may not match.\n\nCarry on anyway, or go back to the home page, which opens the App published for this firmware?'),
                buttons: [
                    { label: t('Carry on'), value: true },
                    { label: t('Go back'), value: false, primary: true },
                ],
                cancelValue: false,          // closing the window is going back
            });
            if (!carryOn) await this.goHome(app);
            return true;
        },

        /** Undocks the lightgun for good and goes back to the home page of the site. */
        async goHome(app) {
            try {
                if (app && app.connection && app.connection.isDocked) await app.disconnect();
            } catch (error) { /* going away anyway */ }
            root.location.href = siteHome(here()) + langSuffix();
        },

        /** Opened by the home page: take the port back if there is one matching candidate
            among the ports already allowed on this site. Otherwise the user selects it
            with Connect. Say why we are here once the messages of the docking are over. */
        async resume(app) {
            let jump = null;
            try {
                const raw = store.get(JUMP_KEY);
                store.drop(JUMP_KEY);
                jump = raw ? JSON.parse(raw) : null;
            } catch (error) {
                return;
            }
            if (!jump) return;

            if (jump.latest && String(jump.latest) !== myVersion()) {
                app.versionNotice = OF.UI.t('App %1, published for the firmware of this lightgun. A newer firmware is available: %2.',
                                            myLabel(), String(jump.latestLabel || jump.latest));
                app.status.show(app.versionNotice, 12000);
            }

            if (!OF.WebSerialTransport.isSupported()) return;
            let ports = [];
            try {
                ports = await OF.WebSerialTransport.getKnownPorts();
            } catch (error) {
                return;
            }
            if (!ports.length) return;
            const wanted = jump.port || {};
            const matches = ports.filter((item) => {
                const info = OF.WebSerialTransport.describePort(item);
                return info.productId === wanted.productId && info.vendorId === wanted.vendorId;
            });
            // VID/PID identify a USB type, not an individual lightgun: do not guess.
            if (matches.length === 1) await app.connectPort(matches[0]);
        },
    };
})(typeof globalThis !== 'undefined' ? globalThis : this);
