/*  OpenFIRE Web App - button output maps and fixed lists (Qt App: appcommon.h).

    Button mapping (config.buttons[button]) is [type, value] for three slots:
    0 on-screen, 1 off-screen, 2 gamepad mode. Types: 0 Mouse, 1 Keyboard, 2 Gamepad.
    Each list below is in display order: { name, value } where value is the code used by the firmware.
*/
(function (root) {
    'use strict';

    const OF = root.OF = root.OF || {};

    const ch = (c) => c.charCodeAt(0);

    const MOUSE = [
        ['Left Click', 0b00000001],
        ['Right Click', 0b00000010],
        ['Middle Click', 0b00000100],
        ['Side Button Back', 0b00001000],
        ['Side Button Forward', 0b00010000],
    ];

    const KEYBOARD = [
        ['Player-relative Start Key', 0xFF],
        ['Player-relative Coin Key', 0xFE],
        ['Up Arrow', 0xDA],
        ['Down Arrow', 0xD9],
        ['Left Arrow', 0xD8],
        ['Right Arrow', 0xD7],
        ['Enter/Return', 0xB0],
        ['Backspace', 0xB2],
        ['Escape', 0xB1],
        ['Left Ctrl', 0x80],
        ['Right Ctrl', 0x84],
        ['Left Alt', 0x82],
        ['Right Alt', 0x86],
        ['Left Shift', 0x81],
        ['Right Shift', 0x85],
        ['Tab', 0xB3],
        ...'ABCDEFGHIJKLMNOPQRSTUVWXYZ'.split('').map((letter) => [letter, ch(letter.toLowerCase())]),
        ...Array.from({ length: 12 }, (_, i) => [`F${i + 1}`, 0xC2 + i]),
    ];

    const GAMEPAD = [
        ['A Button', 0],
        ['B Button', 1],
        ['X Button', 3],
        ['Y Button', 4],
        ['Left Shoulder', 6],
        ['Right Shoulder', 7],
        ['Left Trigger', 8],
        ['Right Trigger', 9],
        ['Select Button', 10],
        ['Start Button', 11],
        ['Left Stick Click', 13],
        ['Right Stick Click', 14],
        ['D-Pad Up', 15],
        ['D-Pad Down', 16],
        ['D-Pad Left', 17],
        ['D-Pad Right', 18],
    ];

    const toList = (entries) => Object.freeze(entries.map(([name, value]) => Object.freeze({ name, value })));

    const Maps = {
        BUTTON_COUNT: 14,
        INPUT_MOUSE: 0,
        INPUT_KEYBOARD: 1,
        INPUT_GAMEPAD: 2,
        inputTypes: Object.freeze(['Mouse', 'Keyboard', 'Gamepad']),
        outputs: Object.freeze([toList(MOUSE), toList(KEYBOARD), toList(GAMEPAD)]),

        irSensitivity: Object.freeze(['Default', 'Higher', 'Highest']),
        runModes: Object.freeze(['Normal', '1-Frame Avg', '2-Frame Avg']),
        layouts: Object.freeze(['Square', 'Diamond']),
        aspectRatios: Object.freeze(['16:9', '16:10', '3:2', '5:4', '4:3']),
        analogModes: Object.freeze(['Gamepad Analog Stick (Left/Right)', 'Gamepad D-Pad', 'Keyboard Arrows']),
        cameraModels: Object.freeze(['DFRobot SEN0158 / WiiCam', 'PixArt PAJ7025 R2', 'PixArt PAJ7025 R3']),
        i2cTypeLabels: Object.freeze(['SDA', 'SCL']),
        spiTypeLabels: Object.freeze(['RX', 'TX', 'SCK', 'CSn']),

        /** Index in the display list of an output value, or -1. */
        outputIndex(type, value) {
            const list = this.outputs[type];
            return list ? list.findIndex((entry) => entry.value === value) : -1;
        },

        /**
         * Names of the board functions ordered by value + 1 (Qt: boardInputs_sortedStr):
         * index 0 "Unmapped", 1 "Trigger", ...
         */
        functionNames(shared) {
            const S = shared || OF.Boards.shared;
            const names = new Array(S.boardInputs_e.boardInputsCount + 1).fill('');
            for (const [name, value] of Object.entries(S.boardInputs_Strings))
                if (value + 1 >= 0 && value + 1 < names.length) names[value + 1] = name;
            return names;
        },
    };

    OF.Maps = Maps;
})(typeof globalThis !== 'undefined' ? globalThis : this);
