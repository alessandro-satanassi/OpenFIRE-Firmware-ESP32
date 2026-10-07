<a id="english-version"></a>

[English](#english-version) · [Italiano](#versione-italiana)

# Building the OpenFIRE ESP32 firmware

> **Just want to use the gun?** You do not need to build anything: install the firmware with the [Web Flasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=en) and follow the [lightgun guide](../README.md#english-version).

> [!NOTE]
> If you report a problem with a firmware you built yourself, **say what you changed in the code**. If it is a general firmware problem, check first whether it also happens with the official release files.

## What you need
- [Visual Studio Code](https://code.visualstudio.com/) with the **PlatformIO IDE** extension, or [PlatformIO Core](https://platformio.org/install/cli) for the command line.
- Git, to clone the repository.
- An Internet connection for the first build: PlatformIO downloads the ESP32 platform and the libraries listed in `platformio.ini`.

## Project layout
The repository contains three separate PlatformIO projects, one per device, each with its own `platformio.ini`:

| Folder | Device |
| --- | --- |
| `lightgun/` | the lightgun (it also builds and embeds the WebApp) |
| `dongle/` | the USB receiver |
| `pedal/` | the wireless pedal |

`shared_boards/` holds the board definitions and partition tables, `shared_lib/` the libraries shared by the three projects. Opening `OpenFIRE.code-workspace` in VS Code shows all of them in one window.

## Choosing the board
Each project has one environment per supported board. Pick it in the PlatformIO toolbar, or set `default_envs` at the top of `platformio.ini`.

| Lightgun environment | Dongle environment | Pedal environment |
| --- | --- | --- |
| `ESP32_S3_WROOM1_DevKitC_1_N16R8` | `LILYGO_T_DONGLE_S3` | `ESP32_S3_WROOM1_DevKitC_1_N16R8` |
| `ESP32_S3_WROOM1_DevKitC_1_N8R2` | `GNPE_POCKET_DONGLE_S3_N16R8` | `ESP32_S3_WROOM1_DevKitC_1_N8R2` |
| `WAVESHARE_ESP32_S3_PICO` | `ESP32_S3_WROOM1_DevKitC_1_N16R8` | `WAVESHARE_ESP32_S3_PICO` |
| `WAVESHARE_ESP32_S3_ZERO_N8R8` | `ESP32_S3_WROOM1_DevKitC_1_N8R2` | `WAVESHARE_ESP32_S3_ZERO_N8R8` |
| `WAVESHARE_ESP32_S3_ZERO_N4R2` | `WAVESHARE_ESP32_S3_PICO` | `WAVESHARE_ESP32_S3_ZERO_N4R2` |
| | `WAVESHARE_ESP32_S3_ZERO_N8R8` | |
| | `WAVESHARE_ESP32_S3_ZERO_N4R2` | |

Choose the exact flash/PSRAM variant of your board, as explained in [Which board variant do I have?](../README.md#board-variant). The lightgun project also has environments for the Raspberry Pi Pico family (`rpipico`, `rpipicow`, `rpipico2`, `rpipico2w`): they build a wired-only firmware, without the wireless features.

## Building and uploading
In VS Code use the PlatformIO **Build** and **Upload** buttons. From the command line, inside the device folder:

```bash
cd lightgun
pio run -e WAVESHARE_ESP32_S3_PICO              # build
pio run -e WAVESHARE_ESP32_S3_PICO -t upload    # build and upload through the board's USB port
pio run -e WAVESHARE_ESP32_S3_PICO -t erase     # erase the whole flash (clean installation)
```

Uploading keeps the settings saved in the gun; erase first for a clean installation. If the board does not enter update mode by itself, use its BOOT/RESET buttons (on a lightgun already running 7.0.0 you can also start it holding **Trigger + A**).

## The WebApp inside the lightgun
Every lightgun build runs `scripts/pack_webapp.py`, which rebuilds the WebApp from `webapp/` and embeds it in the firmware (`include/web_assets.h`), the version served by the offline WebApp mode. The same build also writes the published WebApp of this firmware to `dist/site` and the home page of the site to `dist/launcher`. The first time a board picture needs compressing, the build installs Pillow in PlatformIO's Python. Details are in [webapp/README.md](../webapp/README.md).

## Build options
- **Features:** the optional parts (OLED display, solenoid, rumble, temperature sensor, analog stick, NeoPixels, RGB LED, hardware switches, MAMEHOOKER support...) are enabled by the `-D USES_...` and similar entries of `build_flags` in the `[common]` section of `lightgun/platformio.ini`. The official files have them all enabled except hardware switches (`USES_SWITCHES`); everything can still be turned off at runtime in the WebApp.
- **Default camera:** `CAMERA_DEFAULT` in the `[camera]` section sets the camera selected after a clean installation (DFRobot/Wii in the official files). Both camera drivers are always built in.
- **Fixed player number:** uncomment `-D PLAYER_NUMBER=1` (1 to 4) in `build_flags` to tie the Start/Select keys to that player, regardless of the player chosen in the WebApp. Without it, the player can be changed in the WebApp or at any time with the serial command `XR#`.
- **USB identity:** `MANUFACTURER_NAME`, `DEVICE_NAME` (at most 15 characters) and `DEVICE_VID` are defined in `src/OpenFIREDefines.h`. For multiplayer there is no need to change them: the player number in the WebApp sets a different Product ID for each gun. Keep `DEVICE_VID` at `0xF143` and `MANUFACTURER_NAME` at `OpenFIRE`: the Apps and several programs recognise OpenFIRE guns by them.
- **Boards and default pins:** board names, default pin layouts and the pictures used by the Apps are defined in `src/boards/OpenFIREshared.h` (see [src/boards/README.md](../src/boards/README.md)). The current default pins of the ESP32 boards are also listed in [BOARDS.md](BOARDS.md).

---

<a id="versione-italiana"></a>

[English](#english-version) · [Italiano](#versione-italiana)

# Compilare il firmware OpenFIRE ESP32

> **Vuoi solo usare la pistola?** Non serve compilare nulla: installa il firmware con il [Web Flasher](https://alessandro-satanassi.github.io/OpenFIRE-ESP32-WebFlasher/?lang=it) e segui la [guida lightgun](../README.md#versione-italiana).

> [!NOTE]
> Se segnali un problema con un firmware compilato da te, **indica cosa hai modificato nel codice**. Se è un problema generale del firmware, verifica prima se si presenta anche con i file ufficiali della release.

## Cosa serve
- [Visual Studio Code](https://code.visualstudio.com/) con l'estensione **PlatformIO IDE**, oppure [PlatformIO Core](https://platformio.org/install/cli) per la riga di comando.
- Git, per clonare il repository.
- Una connessione Internet per la prima compilazione: PlatformIO scarica la piattaforma ESP32 e le librerie elencate in `platformio.ini`.

## Struttura del progetto
Il repository contiene tre progetti PlatformIO separati, uno per dispositivo, ognuno con il proprio `platformio.ini`:

| Cartella | Dispositivo |
| --- | --- |
| `lightgun/` | la lightgun (compila e incorpora anche la WebApp) |
| `dongle/` | il ricevitore USB |
| `pedal/` | il pedale wireless |

`shared_boards/` contiene le definizioni delle schede e le tabelle delle partizioni, `shared_lib/` le librerie comuni ai tre progetti. Aprendo `OpenFIRE.code-workspace` in VS Code li vedi tutti in un'unica finestra.

## Scegliere la scheda
Ogni progetto ha un ambiente (environment) per ogni scheda supportata. Sceglilo nella barra di PlatformIO, oppure imposta `default_envs` all'inizio di `platformio.ini`.

| Ambiente lightgun | Ambiente dongle | Ambiente pedale |
| --- | --- | --- |
| `ESP32_S3_WROOM1_DevKitC_1_N16R8` | `LILYGO_T_DONGLE_S3` | `ESP32_S3_WROOM1_DevKitC_1_N16R8` |
| `ESP32_S3_WROOM1_DevKitC_1_N8R2` | `GNPE_POCKET_DONGLE_S3_N16R8` | `ESP32_S3_WROOM1_DevKitC_1_N8R2` |
| `WAVESHARE_ESP32_S3_PICO` | `ESP32_S3_WROOM1_DevKitC_1_N16R8` | `WAVESHARE_ESP32_S3_PICO` |
| `WAVESHARE_ESP32_S3_ZERO_N8R8` | `ESP32_S3_WROOM1_DevKitC_1_N8R2` | `WAVESHARE_ESP32_S3_ZERO_N8R8` |
| `WAVESHARE_ESP32_S3_ZERO_N4R2` | `WAVESHARE_ESP32_S3_PICO` | `WAVESHARE_ESP32_S3_ZERO_N4R2` |
| | `WAVESHARE_ESP32_S3_ZERO_N8R8` | |
| | `WAVESHARE_ESP32_S3_ZERO_N4R2` | |

Scegli la variante flash/PSRAM esatta della tua scheda, come spiegato in [Quale variante di scheda ho?](../README.md#variante-scheda). Il progetto lightgun ha anche ambienti per la famiglia Raspberry Pi Pico (`rpipico`, `rpipicow`, `rpipico2`, `rpipico2w`): producono un firmware solo via cavo, senza le funzioni wireless.

## Compilare e caricare
In VS Code usa i pulsanti **Build** e **Upload** di PlatformIO. Da riga di comando, dentro la cartella del dispositivo:

```bash
cd lightgun
pio run -e WAVESHARE_ESP32_S3_PICO              # compila
pio run -e WAVESHARE_ESP32_S3_PICO -t upload    # compila e carica dalla porta USB della scheda
pio run -e WAVESHARE_ESP32_S3_PICO -t erase     # cancella tutta la flash (installazione pulita)
```

Il caricamento conserva le impostazioni salvate nella pistola; per un'installazione pulita cancella prima la flash. Se la scheda non entra da sola in modalità aggiornamento, usa i suoi pulsanti BOOT/RESET (su una lightgun che esegue già la 7.0.0 puoi anche avviarla tenendo premuti **Grilletto + A**).

## La WebApp dentro la lightgun
Ogni compilazione della lightgun esegue `scripts/pack_webapp.py`, che ricostruisce la WebApp da `webapp/` e la incorpora nel firmware (`include/web_assets.h`): è la versione servita dalla modalità WebApp offline. La stessa compilazione scrive anche la WebApp pubblicata di questo firmware in `dist/site` e la pagina iniziale del sito in `dist/launcher`. La prima volta che un'immagine di scheda va compressa, la compilazione installa Pillow nel Python di PlatformIO. I dettagli sono in [webapp/README.md](../webapp/README.md).

## Opzioni di compilazione
- **Funzioni:** le parti opzionali (display OLED, solenoide, rumble, sensore di temperatura, stick analogico, NeoPixel, LED RGB, interruttori fisici, supporto MAMEHOOKER...) si abilitano con le voci `-D USES_...` e simili di `build_flags` nella sezione `[common]` di `lightgun/platformio.ini`. I file ufficiali le hanno tutte attive tranne gli interruttori fisici (`USES_SWITCHES`); tutto può comunque essere disattivato durante l'uso dalla WebApp.
- **Telecamera predefinita:** `CAMERA_DEFAULT` nella sezione `[camera]` imposta la telecamera selezionata dopo un'installazione pulita (DFRobot/Wii nei file ufficiali). Entrambi i driver delle telecamere sono sempre inclusi.
- **Numero di giocatore fisso:** togli il commento a `-D PLAYER_NUMBER=1` (da 1 a 4) in `build_flags` per legare i tasti Start/Select a quel giocatore, indipendentemente dal giocatore scelto nella WebApp. Senza, il giocatore si cambia nella WebApp o in qualsiasi momento con il comando seriale `XR#`.
- **Identità USB:** `MANUFACTURER_NAME`, `DEVICE_NAME` (al massimo 15 caratteri) e `DEVICE_VID` sono definiti in `src/OpenFIREDefines.h`. Per il multigiocatore non serve cambiarli: il numero del giocatore nella WebApp imposta un Product ID diverso per ogni pistola. Lascia `DEVICE_VID` a `0xF143` e `MANUFACTURER_NAME` a `OpenFIRE`: le App e diversi programmi riconoscono le pistole OpenFIRE da questi valori.
- **Schede e pin predefiniti:** nomi delle schede, pin predefiniti e immagini usate dalle App sono definiti in `src/boards/OpenFIREshared.h` (vedi [src/boards/README.md](../src/boards/README.md)). I pin predefiniti attuali delle schede ESP32 sono elencati anche in [BOARDS.md](BOARDS.md).
