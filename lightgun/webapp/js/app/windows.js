/*  OpenFIRE Web App - secondary windows: Boards Previewer and About (Qt: apppreviewer, appabout).
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    OF.LOGO_SVG = '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 90.691 90.691"><rect width="90.691" height="90.691" ry="9.6" fill="#2a2a2a"/>' +
        '<g transform="translate(2,2)"><path fill="#fff" d="M86.6 40H81.3C79.6 22.2 65.4 7.9 47.6 6.3V0h-7V6.1C22.8 7.9 8.5 22.1 6.9 40h-7v7h7c1.7 17.8 15.9 32.1 33.7 33.7v5.6h7V80.7C65.4 79 79.7 64.8 81.3 47h5.3zm-72.9 3.5c0-16.8 13.6-30.4 30.4-30.4 16.8 0 28.4 11.8 30.2 26.9H40.6V73.7C25.5 72 13.7 59.1 13.7 43.5ZM47.6 73.6V60.9H59.8V55.2H47.6V46.9H74.3C72.7 60.9 61.6 72 47.6 73.6Z"/>' +
        '<g fill="#ed1a3b"><circle cx="24.9" cy="49" r="3.1"/><circle cx="34" cy="49" r="3.1"/><circle cx="24.9" cy="57.8" r="3.1"/><circle cx="34" cy="57.8" r="3.1"/></g>' +
        '<path fill="#00b3f0" d="M44.3 20.2C33 20.2 23.6 28.7 22.2 39.7h4.7c1.3-8.4 8.6-14.9 17.4-14.9 8.8 0 9.1 1.8 12.5 5.2l3.3-3.3C55.9 22.5 50.3 20.2 44.4 20.2Z"/></g></svg>';

    const GPIO_TOOLTIP = 'GPIO Pin No. %1.\n\n' +
        'Blue pin numbers are members of I2C0.\n' +
        'Orange are members of I2C1.\n' +
        'Purple pin numbers can automatically select any I2C channel in software.\n' +
        'Gray cannot use I2C devices.';
    const CAPS_TOOLTIP = 'ADC indicates whether this pin can read Analog Inputs.\n' +
        'I2C indicates if this pin can interact with I2C devices, and what channel and type it uses.\n' +
        'SPI indicates if this pin can interact with SPI devices, and what channel and type it uses.\n' +
        '(*) means pin can use any type via automated software selectable channels/type.';
    const ESP_FORK = "<p>Compatible with the " +
        "<a href='https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32'><span style=' text-decoration: underline; color:#8ab4f8;'>ESP-IDF fork of the OpenFIRE Firmware</span></a> by <i>Alessandro Satanassi.</i><br>" +
        "Any issues should be reported <b><a href='https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32/issues'><span style=' text-decoration: underline; color:#8ab4f8;'>here!</span></a></b></p>";
    const UPSTREAM = "<p>Compatible with " +
        "<a href='https://github.com/TeamOpenFIRE/OpenFIRE-Firmware'><span style=' text-decoration: underline; color:#8ab4f8;'>upstream OpenFIRE Firmware</span></a> by <i>Team OpenFIRE.</i></p>";

    /** Capability marks of a GPIO (Qt pinCapabilityMarks). */
    function capabilityMarks(S, bits) {
        const { el } = OF.UI;
        const C = S.pinCapabilities_e;
        const mark = (text, cls, strong) => el(strong ? 'b' : 'span', { class: 'cap ' + cls, text });
        const marks = [];
        marks.push(bits & C.pinHasADC ? mark('ADC', 'adc', true) : mark('ADC', 'off'));
        if (bits & C.pinAnyI2C) marks.push(mark('I2C(*)', 'any', true));
        else if (bits & C.pinCanI2C) {
            const channel = (bits & C.pinIsI2C1) >> 3;
            marks.push(mark(`I2C${channel}${OF.Maps.i2cTypeLabels[(bits & C.pinIsI2CSCL) >> 2]}`, channel ? 'i2c1' : 'i2c0', true));
        } else marks.push(mark('I2C', 'off'));
        if (bits & C.pinAnySPI) marks.push(mark('SPI(*)', 'spi-any'));
        else if (bits & C.pinCanSPI) {
            marks.push(mark(`SPI${(bits & C.pinIsSPI1) >> 4}${OF.Maps.spiTypeLabels[((bits & C.pinCanSPI) >> 5) - 1]}`, 'spi', true));
        } else marks.push(mark('SPI', 'off'));
        const node = el('span', { class: 'cap-marks mono' });
        marks.forEach((m, i) => { if (i) node.append(' '); node.append(m); });
        return node;
    }

    let previewer = null;

    /** Boards Previewer (Qt AppBoardsPreviewer, a window that stays open next to the main one):
        default layout and pin capabilities of every compatible board. */
    function openPreviewer(currentBoard) {
        if (previewer) {
            previewer.select(currentBoard);
            previewer.node.querySelector('.preview-select').focus();
            return previewer.done;
        }
        const { el, t } = OF.UI;
        const S = OF.Boards.shared;
        const P = S.boardBoxPositions_e;
        const names = OF.Maps.functionNames(S);
        // std::map order of the Qt boardNames (sorted by board name).
        const boards = OF.Boards.list().slice().sort((a, b) => (a < b ? -1 : a > b ? 1 : 0));

        const selector = el('select', { class: 'preview-select', attrs: { 'aria-label': t('Boards Previewer') } },
            el('option', { value: '', text: t('Click Here to Select a Board'), disabled: true, hidden: true }),
            ...boards.map((board) => el('option', { value: board, text: OF.Boards.displayName(board) })));
        const subtext = el('div', { class: 'preview-subtext rich' });
        const line = el('hr', { class: 'preview-line', hidden: true });
        const left = el('div', { class: 'pins-col pins-left' });
        const right = el('div', { class: 'pins-col pins-right' });
        const middle = el('div', { class: 'pins-middle' });
        const picture = el('div', { class: 'board-picture' });
        const layout = el('div', { class: 'pins-layout preview' }, left, el('div', { class: 'pins-center' }, picture, middle), right);
        const scroll = el('div', { class: 'preview-scroll' }, layout);
        const content = el('div', { class: 'previewer' }, el('div', { class: 'preview-selector-row' }, selector),
            scroll, line, subtext);
        // Like the board layout tab: the board is shown whole instead of being scrolled.
        const wide = root.matchMedia ? root.matchMedia('(min-width: 761px)') : null;
        const fit = () => OF.UI.fitToHeight(scroll, layout, 0.6, !wide || wide.matches);

        function show(board) {
            left.textContent = '';
            right.textContent = '';
            middle.textContent = '';
            picture.hidden = !board;
            subtext.innerHTML = '';
            line.hidden = !board;
            if (!board) return;
            subtext.innerHTML = t(OF.Boards.isEsp32(board) ? ESP_FORK : UPSTREAM);

            const presets = S.boardsPresetsMap[board] || [];
            const positions = S.boardsBoxPositions[board] || [];
            const caps = OF.Boards.capabilities(board);
            const placed = [];
            let maxLeft = 0;
            let maxRight = 0;

            for (let gpio = 0; gpio < positions.length; ++gpio) {
                const side = positions[gpio] & P.posCheck;
                if (side === P.posNothing) continue;
                const slot = positions[gpio] ^ side;
                const func = el('span', { class: 'preview-func', text: t(names[(presets[gpio] ?? -1) + 1] || '') });
                const label = el('span', { class: 'gpio-label ' + i2cClass(S, caps[gpio] || 0), text: `«GPIO${gpio}»`, title: t(GPIO_TOOLTIP, gpio) });
                const marks = capabilityMarks(S, caps[gpio] || 0);
                marks.title = t(CAPS_TOOLTIP);
                let item;
                if (side === P.posLeft) {
                    item = el('div', { class: 'pin-item', style: { '--row': slot } }, func, label, marks);
                    maxLeft = Math.max(maxLeft, slot);
                    placed.push([left, slot, item]);
                } else if (side === P.posRight) {
                    item = el('div', { class: 'pin-item', style: { '--row': slot } }, marks, label, func);
                    maxRight = Math.max(maxRight, slot);
                    placed.push([right, slot, item]);
                } else {
                    item = el('div', { class: 'pin-item vertical', style: { order: slot } }, marks, label, func);
                    placed.push([middle, slot, item]);
                }
                item.addEventListener('mouseenter', () => { OF.Boards.highlightPin(picture, gpio, true); item.classList.add('hover'); });
                item.addEventListener('mouseleave', () => { OF.Boards.highlightPin(picture, gpio, false); item.classList.remove('hover'); });
            }
            placed.sort((a, b) => a[1] - b[1]).forEach(([column, , item]) => column.append(item));
            left.style.setProperty('--rows', Math.max(1, maxLeft));
            right.style.setProperty('--rows', Math.max(1, maxRight));
            middle.hidden = !middle.children.length;
            fit();
            OF.Boards.showPicture(picture, board);
        }

        selector.addEventListener('change', () => show(selector.value));
        const select = (board) => {
            const value = boards.includes(board) ? board : '';
            if (selector.value === value && value) return;
            selector.value = value;
            show(value);
        };
        select(currentBoard);

        OF.UI.fitOnResize(scroll, fit);
        previewer = { node: content, select };
        previewer.done = OF.UI.dialog({
            title: t('Boards Previewer'), content, wide: true, className: 'preview-dialog', cancelValue: null,
            modeless: true, buttons: [], onClose: () => { previewer = null; },
        });
        return previewer.done;
    }

    function i2cClass(S, bits) {
        const C = S.pinCapabilities_e;
        if (bits & C.pinAnyI2C) return 'any';
        if (bits & C.pinCanI2C) return (bits & C.pinIsI2C1) ? 'i2c1' : 'i2c0';
        return 'none';
    }

    function openAbout() {
        const { el, t } = OF.UI;
        const logo = el('div', { class: 'about-logo', html: OF.LOGO_SVG });
        const content = el('div', { class: 'about' },
            el('div', { class: 'about-head' }, logo,
                el('div', null, el('div', { class: 'wordmark', text: OF.APP_NAME || 'OpenFIRE Esp32' }), el('div', { class: 'about-sub', text: t('Web App') }))),
            el('div', { class: 'rich center', html: t('<html><head/><body><p>Primary OpenFIRE Firmware maintainership and Desktop App developed by <a href="https://github.com/SeongGino"><span style=" text-decoration: underline; color:#8ab4f8;">That One Seong</span></a>.<br/>If this software has given you any value,</p></body></html>') }),
            el('p', { class: 'center' }, el('a', { class: 'kofi', href: 'https://ko-fi.com/thatoneseong', target: '_blank', rel: 'noopener', text: t('Support on Ko-fi') })),
            el('p', { class: 'rich center', html: t('ESP-IDF fork and this WebApp by <b>Alessandro Satanassi</b>.') }),
            el('div', { class: 'rich center', html: t('<html><head/><body><p align="center"><br/>Special Thanks:<br/>Samuel Ballantyne, for their work on the original <a href="https://github.com/samuelballantyne/IR-Light-Gun"><span style=" text-decoration: underline; color:#8ab4f8;">SAMCO lightgun system</span></a>.<br/>Mike Lynch, AKA Prow7, for their enhanced fork and libraries.<br/><a href="https://forum.arcadecontrols.com/index.php/board,20.0.html"><span style=" text-decoration: underline; color:#8ab4f8;">The ArcadeControls forums</span></a> for their support.<br/>All the early GitHub testers for their feedback.<br/>And You, for enjoying this software!</p><p align="center">OpenFIRE, as the name implies, is FREE SOFTWARE;<br/>if you have paid for the firmware or this utility, you have been scammed and should demand your money back!</p></body></html>') }),
            el('div', { class: 'rich center', html: t('<html><head/><body><p><a href="https://github.com/TeamOpenFIRE/OpenFIRE-Firmware"><span style=" text-decoration: underline; color:#8ab4f8;">Firmware</span></a> | <a href="https://github.com/TeamOpenFIRE/OpenFIRE-App"><span style=" text-decoration: underline; color:#8ab4f8;">Desktop Source</span></a> | <a href="discord.com/invite/dFw5z6PBQv"><span style=" text-decoration: underline; color:#8ab4f8;">Discord</span></a></p></body></html>')
                .replace('href="discord.com/', 'href="https://discord.com/') }),
            el('p', { class: 'rich center', html: '<a href="https://github.com/alessandro-satanassi/OpenFIRE-Firmware-ESP32">OpenFIRE Firmware ESP32</a>' }));
        for (const link of content.querySelectorAll('a')) { link.target = '_blank'; link.rel = 'noopener'; }
        return OF.UI.dialog({ title: t('About %1', OF.APP_NAME || 'OpenFIRE Esp32'), content, className: 'about-dialog', cancelValue: null, buttons: [{ label: t('Close'), value: null, primary: true }] });
    }

    OF.Windows = { openPreviewer, openAbout, capabilityMarks };
})(typeof globalThis !== 'undefined' ? globalThis : this);
