# OpenFIRE Web App

The web version of the OpenFIRE App. One source folder, three outputs:

| Output | Boards | Connection | Built by |
| --- | --- | --- | --- |
| Lightgun (web configuration mode) | only its own board (picture inside `app.js`) | WebSocket `ws://<gun>/ws` | `scripts/pack_webapp.py` on every PlatformIO build of an ESP32 environment -> `include/web_assets.h` |
| Published App (GitHub Pages) | all boards | Web Serial | `python scripts/webapp_build.py site` -> `dist/site`, the App of this firmware: it gets published as `v/<version>/` of the site |
| Tauri app | all boards | Web Serial | same files as the published App |
| Home page of the site | - | Web Serial, only to read the version | `python scripts/webapp_build.py launcher` -> `dist/launcher` |

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
js/app/version.js       says so when the lightgun runs a firmware of another version
                        (site only: left out of the build embedded in the lightgun)
js/app/main.js          main window: menus, device selector, tabs, save, board events
lang/<code>.json        translations, key = English text, value = translation
boards/pics/*.js        board pictures, one script per picture of src/boards/boardPics with the same
                        name (rpipico.svg -> rpipico.js), lighter SVG inside (generated, keep them in git)
launcher/               home page of the published site: the Connect button, which reads the version
                        of the lightgun and opens the App published for it (its own small bundle)
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
python scripts/webapp_build.py site                  # dist/site: the App of this firmware
python scripts/webapp_build.py site --out <folder>   # into any folder
python scripts/webapp_build.py launcher              # dist/launcher: the home page of the site
python scripts/webapp_build.py device --board waveshare-esp32-s3-zero --out dist/device
python scripts/build_board_pics.py                   # after changing a board picture (needs Pillow)
```

## Publishing a version

Every PlatformIO build writes two folders under `dist` (generated, not in git):

```
dist/site       the App of the firmware just built - one version of it
dist/launcher   the home page of the site, with the Connect button
```

Publishing a version means copying `dist/site` into `v/<version>/` of the site repository,
putting `dist/launcher` at its root, and adding the version to `versions.json`. Nothing of that
happens by itself: the build only prepares the two folders, and what reaches the published site
stays a deliberate step.

`dist/site` can be written somewhere else - a folder you keep for testing - by setting its path in
`platformio.ini` (under `[env]`, or in a single environment):

```
custom_webapp_site_dir = ../../OpenFIRE-WebApp
```

It then receives `index.html`, `style.css`, `app.js` and `boards/pics/*.js`, leaving everything
else alone (`.git`, `README`, `LICENSE`, `CNAME`), and removes the pictures of boards that no
longer exist. Files are rewritten only when their content changes, so `git status` in that folder
shows exactly what the build changed. A relative path starts from the `lightgun` folder;
`OPENFIRE_WEBAPP_SITE_DIR` in the environment has priority over `platformio.ini`, and the value
`off` writes no App folder at all.

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

## The published site: one folder per version

Every firmware has its own App, built at the same time, and every version of the App stays
published. What the build writes is not a site: `dist/site` is **the App of the firmware being
built**, one version of it. Publishing a version means copying that folder into the site:

```
/                    the home page: the Connect button, and nothing else (dist/launcher)
/versions.json       the versions published so far, written when one is published
/v/6.2/              the App of firmware 6.2, a copy of dist/site as it was then
/v/6.1/              the App of firmware 6.1, untouched since the day it was published
/.nojekyll           GitHub Pages serves the files as they are
```

The name of the folder is the number the firmware sends when the App docks (`"%.1f"` of
`OPENFIRE_VERSION` in `src/OpenFIREversion.h`, so `6.2`): it is the only version an already
installed lightgun can tell the App about. Each folder is a complete copy, board pictures
included, so a version keeps working on its own for as long as it is there.

```json
{ "latest": "6.2",
  "versions": [ { "id": "6.2", "label": "6.2.0", "type": "stable" },
                { "id": "6.1", "label": "6.1.0", "type": "stable" } ] }
```

### The home page chooses, the App does not

`webapp/launcher/` is the home page, built into `dist/launcher` by
`python scripts/webapp_build.py launcher`. It is not the App: one button. It opens the port, asks
the lightgun who it is with `Protocol.getBoardInfo()` - which docks for the first answer alone and
reads no setting - undocks, looks the version up in `versions.json` and opens `v/<version>/`. The
port is handed over in `sessionStorage`, so the App that opens takes it back with
`navigator.serial.getPorts()` and docks by itself: the lightgun is docked twice, but the user only
sees the page change and clicks once.

It uses the App's own `transport.js` and `protocol.js`, joined into `launcher.js` the same way the
App's scripts are joined into `app.js`: there is one implementation of the protocol, not two, and
the home page cannot drift away from it.

When the firmware's version is not published, **nothing is opened by itself**: the versions that
are there are offered, so a firmware nobody made an App for can still be tried with a neighbouring
one. Since a version of the App is published together with every firmware, that should not happen;
it is there so that it fails politely if it ever does.

### What a published App checks

An App, once open, never goes looking for another one: it is the App of its version and it stays
that. All it does is `js/app/version.js` (`OF.Version`):

```
OF.Version.afterDock(app, board) once docked: when the firmware is of another version, asks
                                 whether to carry on or to go back to the home page
OF.Version.resume(app)           opened by the home page: takes the port back, and says in the
                                 status bar which App this is and that a newer firmware exists
```

The question is asked once per version, not at every reconnection. *Carry on* keeps the App and
the connection; *Go back* - and closing the window, which is the same thing - undocks the lightgun
for good and opens the home page again, which then opens the right App. A copy that is not part of
a published site (a folder on the disk, the desktop app) has nowhere to go back to, so it only
says that the versions differ.

**The App embedded in the lightgun does not include this file.** `scripts/webapp_build.py` leaves
it out of the device build (`SITE_ONLY_SCRIPTS`): the gun serves the App of its own firmware, so
there is nothing to check, and every call to `OF.Version` in the rest of the App is guarded.

`?lang=it` in the address opens the App, and the home page, in that language (the pages of the
project pass the language on to each other); it is not remembered, so a link does not change what
this browser usually shows.

### Reading the version, for the years to come

The version is the first field of the first answer the lightgun gives, and the App reads it
**before it checks anything else**: it is returned even when the rest of the payload cannot be
read (`bad_board_info`), which is exactly the case a firmware of another generation would produce,
and the one where opening the right App matters most.

After the USB table the answer carries items, each a separator and a marker. The rule is fixed:
the two markers that are only a flag (`sError`, `sPedalWireless`) are two bytes; **every other
marker carries a length byte and then its data**. An App that does not know a marker therefore
skips it whole and still reads the ones after it, however many are added in the years to come.
`sVersionFull` is the first of those: the complete `6.2.0-stable`, which the home page uses to
look for an exact folder before falling back to the `6.2` one. Both Apps parse the trailer this
way (web: `protocol.js`; Qt: `appserial.cpp`), and the firmware writes it in `main.cpp`.

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
python tests/site-build.test.py                      # what the build writes: App and home page
```

`sim-server.js` serves the files and connects the page WebSocket to the mock firmware
(`tests/mock-firmware.js`), which follows `src/OpenFIREserial.cpp`.
Opening the unbundled folder with `?ws` makes the site version use the WebSocket.
