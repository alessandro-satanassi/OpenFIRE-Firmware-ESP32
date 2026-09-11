# OpenFIRE Web App

The web version of the OpenFIRE App. One source folder, three outputs:

| Output | Boards | Connection | Built by |
| --- | --- | --- | --- |
| Lightgun (web configuration mode) | only its own board (picture inside `app.js`) | WebSocket `ws://<gun>/ws` | `scripts/pack_webapp.py` on every PlatformIO build of an ESP32 environment -> `include/web_assets.h` |
| Site (GitHub Pages) | all boards | Web Serial | `python scripts/webapp_build.py site` -> `dist/site` |
| Tauri app | all boards | Web Serial | same files as the site |

## Folder

```
index.html              page skeleton; the scripts between the <!-- OF:SCRIPTS --> markers are joined into app.js
style.css               styles, light and dark theme (CSS variables)
js/core/transport.js    Web Serial and WebSocket byte streams
js/core/protocol.js     App protocol (peer of src/OpenFIREserial.cpp in the firmware)
js/core/i18n.js         translations (OF.i18n.t)
js/core/boards.js       board data helpers (OF.Boards)
js/app/connection.js    connection controller (auto WebSocket on the gun, Web Serial on the site)
js/app/maps.js          mouse / keyboard / gamepad outputs, profile option lists
js/app/state.js         configuration being edited and the editing rules of the Qt App (no DOM)
js/app/ui.js            elements, dialogs, spin boxes, menus, status bar
js/app/tab-*.js         Board Layout, Button Mapping, Gun Settings, Calibration Profiles, Gun Tests
js/app/fullscreen.js    calibration, IR emitter alignment and IR camera test screens
js/app/testfont.js      8x8 bitmap typeface of those screens (from the Qt App)
js/app/windows.js       Boards Previewer and About
js/app/main.js          main window: menus, device selector, tabs, save, board events
lang/<code>.json        translations, key = English text, value = translation
boards/pics/*.js        board pictures, one script per picture of src/boards/boardPics with the same
                        name (rpipico.svg -> rpipico.js), lighter SVG inside (generated, keep them in git)
tests/                  Node tests, simulated lightgun, browser checks

generated (not in git):
boards/OpenFIREshared.js   from src/boards/OpenFIREshared.h
lang/translations.js       from lang/*.json
```

## Translations

`lang/en.json` lists every text of the app (key and value are the English text).
To add a language copy it to `lang/<code>.json` (`de.json`, `fr.json`, `pt-BR.json`...) and
translate the values; texts left identical or empty are shown in English. The language appears
in the selector at the next build. Keep `%1`, `%2`... and the HTML tags of the texts that have them.
Texts of the calibration screens are drawn with a bitmap typeface: accented lowercase letters are
supported, other scripts use the system font.

## User interface

The same page is used by the gun, the site and Tauri; differences:

| | lightgun page | site / Tauri |
| --- | --- | --- |
| connection | automatic WebSocket, reconnects by itself | device selector (`COM Port`) and `Add a Device...` (Web Serial) |
| Board Previews / View Compatible Boards | hidden (only its own board is included) | all boards |

View menu: `Show Unsafe Settings` (solenoid temperature thresholds) and the theme (system, light, dark).
The theme is remembered by the browser; `Show Unsafe Settings` starts off every time, like in the Qt App.
The fullscreen screens use the browser fullscreen mode: ESC works like in the Qt App (leaving the
browser fullscreen mode in any other way also closes the screen). When the browser refuses fullscreen
mode a button at the top asks for it again.

Lightgun page: the gun serves one page at a time and a newly opened page takes it over. A page in a
background tab waits until it is shown; a page that loses the gun three times within 30 seconds stops
retrying and shows `Reconnect`. When the link drops with unsaved edits (calibrations included), they
are offered back as soon as the same gun docks again; the browser also asks before leaving a page
with unsaved edits.
`QT_PARITY.md` maps every function of the Qt App to the web app and lists the intentional differences.

## Commands (from the `lightgun` folder)

```
python scripts/webapp_build.py sizes                 # refresh generated files, print gzipped sizes
python scripts/webapp_build.py site                  # dist/site
python scripts/webapp_build.py device --board waveshare-esp32-s3-zero --out dist/device
python scripts/build_board_pics.py                   # after changing a board picture (needs Pillow)
```

## Boards: what updates by itself

Every PlatformIO build (and every `webapp_build.py` command) re-reads:

- `src/boards/OpenFIREshared.h`: boards, names, default pins, box positions, pin capabilities,
  alternative layouts (`boardsAltPresets`), architectures, commands and setting names;
- `src/boards/boardPics/boards.qrc` and the board SVGs it lists;
- `webapp/lang/*.json`.

Rules:

- all board data and pictures come only from `src/boards`;
- the picture of a board is the file whose alias in `boards.qrc` is the board name;
- the circle drawn over a pin has id `OF_pin<gpio>` and style `opacity:0`;
- a board name containing an architecture name of `boardArchs` (e.g. `esp32-s3`) uses that
  architecture, otherwise RP2040/235X;
- the firmware board comes from the `OPENFIRE_BOARD` `#ifdef`/`#elifdef`/`#elif defined(...)` chain.

The build prints `[WebApp] WARNING` for header values it cannot read, boards missing from a map,
missing pictures and environments whose board is not recognised.
A changed SVG with embedded photos needs Pillow to be compressed: the PlatformIO build installs
it by itself in its Python the first time it is needed (internet required once); without it the
picture is used uncompressed and the next build tries again.

## Opening the page as a local file

`dist/site/index.html` (or `webapp/index.html` after a build) also works opened directly from the
disk (`file://`) in Chrome/Edge: Web Serial is available there, and the board pictures are scripts
(`boards/pics/<picture>.js`), so the pin circles light up as on the site. Copy the whole folder.

## Development without hardware (from `lightgun/webapp`)

```
node --test --test-concurrency=1 tests/*.test.js     # protocol, sessions, transports, i18n, boards, editing rules
node tests/sim-server.js --root .                    # open http://localhost:8080/?ws
node tests/sim-server.js --root ../dist/device       # the page exactly as served by the gun
node tests/browser-e2e.js                            # browser checks of the three builds (needs Playwright)
```

`sim-server.js` serves the files and connects the page WebSocket to the mock firmware
(`tests/mock-firmware.js`), which follows `src/OpenFIREserial.cpp`.
Opening the unbundled folder with `?ws` makes the site version use the WebSocket.
