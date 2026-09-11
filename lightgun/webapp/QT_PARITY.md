# Qt App -> Web App correspondence

Every slot and function of the Qt App (`OpenFIRE-App/src`) and where the web app implements it.
Tests: `tests/*.test.js` (Node) and `tests/browser-e2e.js` (browser, simulated lightgun).

## Main window (`appmainwindow.cpp`)

| Qt | Web app |
| --- | --- |
| `aliveTimer_timeout`, `serialPort_SearchFinished` (ports with VID 0xF143) | `main.js` `refreshPorts` / `renderPortSelector` (Web Serial `getPorts`, connect/disconnect events); the lightgun page connects by itself |
| `on_comPortSelector_currentTextChanged` (load board, build tabs, disconnect) | `main.js` `onPortSelected`, `connectPort`, `onDocked`, `onDisconnected` |
| `eventFilter` (pin circle on hover, description boxes, no wheel on combo boxes) | `tab-pins.js` hover/focus, `ui.js` `DescriptionBox`, `select()` wheel guard |
| `BoxesUpdate`, `pinBoxes_currentIndexChanged` (dedupe, I2C rules, enables) | `state.js` `setCustomPins`, `changePin`, `_setPin`; enables in `tab-settings.js` / `tab-buttons.js` `update` |
| `on_customPinsEnabled_stateChanged`, `on_presetsBox_currentIndexChanged` | `state.js` `setCustomPins`, `applyAltPreset` |
| `on_actionImport_Custom_Layout_triggered`, `on_actionExport_Custom_Layout_triggered` (.ofl) | `tab-pins.js` User Layouts menu, `state.js` `importLayout` / `exportLayout` |
| `btnFuncTypeBox_currentIndexChanged`, `btnFuncBox_currentTextChanged`, `on_aStickModeBox_currentIndexChanged` | `tab-buttons.js`, `state.js` `setButtonType` / `setButtonValue` |
| `on_*Toggle_stateChanged` (rumble, solenoid, rumble FF, autofire, pause, anode, low buttons, OLED...) | `state.js` `setToggle` (same cascade), `tab-settings.js` |
| `on_*Box_valueChanged` (lengths, intensity, temperatures, NeoPixels) | `tab-settings.js` spin boxes (same ranges), `state.js` `setSetting` |
| `on_customLEDstaticBtn1..3_clicked`, `on_invertStaticPixelsBox_stateChanged`, `PixelsDiff` | `tab-settings.js` colour buttons, First/Last prefix, power cycle notice |
| `on_cameraModelUnlockToggle_stateChanged`, `on_camModel_*_clicked` | `tab-settings.js` camera group |
| `on_tinyUSBLayoutToggle_stateChanged`, `on_tUSB_p1..4_toggled`, `on_productIdInput_valueChanged`, `on_productNameInput_textEdited` | `tab-settings.js` TinyUSB group, `state.js` `setUsbPreset` / `setUsbId` / `setUsbName` |
| `renameBoxes_clicked`, `colorBoxes_clicked`, `profileBoxes_activated`, `selectedProfile_isChecked` | `tab-profiles.js`, `main.js` `selectProfile` |
| `caliBtns_clicked`, `NewCaliWindow`, `CaliWindowExiting`, `CaliWindowRequestedExit` | `main.js` `calibrate`, `openIRTest`, `openAlignment`, `onWindowExit`; `fullscreen.js` |
| `on_rumbleTestBtn_clicked`, `on_solenoidTestBtn_clicked`, `on_red/green/blueLedTestBtn_clicked` | `tab-tests.js` feedback tests |
| `on_testBtn_clicked` (unsaved IR layout warning) | `main.js` `openIRTest` |
| `on_baudResetBtn_clicked`, `on_clearEepromBtn_clicked` | `main.js` `rebootToBootloader`, `clearSaveMemory` |
| `LabelsUpdate` (button labels, N/C, temperature, LED tests, title) | `tab-tests.js` `resetReadings`, `main.js` `updateHeader` |
| `serialPort_readyRead` (buttons, temperature, analog, current profile, errors, calibration, test coords, clear flash) | `main.js` `onEvent` / `onError`, `tab-tests.js` `onEvent` |
| `serialPort_progressSet`, `serialPort_progressUpdate` | `main.js` `onProgress`, status bar progress |
| `on_confirmButton_clicked`, `DiffUpdate`, `CommitRecoveryUiUpdate` | `main.js` `save`, `refresh`; `state.js` `isDirty` / `commitDone` |
| `on_tabWidget_currentChanged` (reset description boxes) | `main.js` `selectTab` -> tab `onShow` |
| `on_actionShow_Unsafe_Settings_toggled` | View > Show Unsafe Settings (not remembered, like Qt) |
| `on_actionCompatible_Boards_triggered`, `on_actionAbout_UI_triggered` | `windows.js` `openPreviewer`, `openAbout` |
| `on_actionOpenFIRE_Documentation_triggered` (Alt+D), `on_actionOpenFIRE_Serial_Usage_triggered` (Alt+S) | Help menu and the same shortcuts |
| `on_actionOpen_IR_Emitter_Alignment_Assistant_triggered` | Emitter Alignment / Help menu |
| `on_actionDebug_Window_triggered` (non-Release builds) | View > Debug Window, unbundled folder or `?debug` |

## Other windows

| Qt | Web app |
| --- | --- |
| `appcali.cpp` (calibration stages, info texts, verify, malformed warning, alignment boxes, IR test points, bitmap text, text scales) | `fullscreen.js` |
| `apppreviewer.cpp` (board list, default functions, GPIO colours, ADC/I2C/SPI marks, fork/upstream link) | `windows.js` `openPreviewer` |
| `appabout.ui` | `windows.js` `openAbout` |
| `appdebug.cpp` (text and hex of received payloads) | `main.js` `openDebugWindow` / `logDebug` |
| `appserial.cpp` (framing, ACK, dock, settings sync, commit, recovery, reboot, errors, RequestToReboot, port busy / stale dock dialogs) | `js/core/protocol.js`, `js/app/connection.js`, `main.js` `onDockFailure` / `onConnectError` |

## Intentional differences

- Linux checks of the Qt App (user not in `dialout`, running as root): not applicable to a browser.
- File dialogs: the browser file chooser and download replace "Open New Layout" / "Save New Layout"; a download cannot be cancelled, so "Canceled custom layout save operation." does not exist.
- Device list: Web Serial only lists ports the user allowed, hence "Add a Device...".
- Menu bar: View, Help and About with the Qt entries; "Board Previews" and "Emitter Alignment" (also in Help) are added as buttons, with the language and theme selectors. "OpenFIRE on the Web" (a disabled placeholder in the Qt menu) is not shown; About opens directly (the Qt About menu has only "About OpenFIRE...").
- Turning custom pins off/on or choosing an alternative preset keeps the rumble/solenoid toggles when their pins are mapped again: the Qt App turns them off as a side effect of clearing the pin boxes one by one. Importing a layout gives rumble and solenoid back their saved values, exactly like the Qt App.
- Changing the output type of a button selects the first output of the new list (the Qt App showed it but kept the previous value).
- Renaming a profile starts from the current name (selected, typing replaces it); the Qt box starts empty. Characters outside Latin-1 are refused like in the product name.
- After Clear Save Memory the confirmation is shown as soon as the command is delivered: the firmware restarts without answering.
- Fullscreen screens: the browser owns fullscreen mode, so leaving it (ESC, F11, browser gesture) closes the screen like ESC; the page behind cannot be used meanwhile. If the browser refuses fullscreen mode, a button asks for it again.
- Leaving IR test mode without the board confirming it is reported as an unconfirmed operation (reconnect), instead of leaving the tabs locked.
- The Boards Previewer and the colour chooser are page windows; the previewer, like the Qt window, stays open while the main window is used.
- About also names the ESP-IDF fork and the web app; the Qt version line is not shown.
- Colours of the GPIO labels: exact Qt colours in the dark theme, slightly darker ones in the light theme for contrast.
- The rumble and solenoid tests follow the unsaved toggles like the Qt App; the firmware ignores the rumble test while the saved map has no rumble pin.

## Added for the browser

- Lightgun page: automatic reconnection; a page in a background tab does not take the gun from the page in use; after three takeovers in 30 seconds a page stops and offers `Reconnect` (two open pages would otherwise take the gun from each other forever).
- Unsaved edits of a lost link (calibrations kept only in the board RAM included) are offered back when the same gun docks again; the browser asks before leaving a page with unsaved edits.
- Questions left open when the link is lost (Save, Clear Save Memory, rename, colours, layout import, IR test warning) close by themselves and do nothing.
- Links of the Qt texts open in a new tab, so the page keeps its board.
