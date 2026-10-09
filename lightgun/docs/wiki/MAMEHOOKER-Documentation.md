### Table of Contents:
 - [MAMEHOOKER Setup Guide](#mamehooker-setup-guide)
 - [How does MAMEHOOKER work?](#how-does-mamehooker-work)
 - [How do configuration files work?](#how-do-configuration-files-work)
 - [Opening, Writing to, and Closing COM Ports](#opening-writing-to-and-closing-com-ports)
 - [Light Gun Serial Commands](#light-gun-serial-commands)
 - [MAMEHOOKER's Pipe & State Modifiers](#mamehookers-pipe--state-modifiers)

> **OpenFIRE ESP32 7.0.0:** these commands work exactly as in the original OpenFIRE. They reach the gun through its serial port: the gun's own USB port, or its dongle's port when it plays wirelessly. Only one program can use the port at a time, so close the WebApp first. Start the gun normally: in the offline WebApp mode (B held at startup) its USB serial port is not available.

> **OpenFIRE ESP32 7.0.0:** questi comandi funzionano esattamente come nell'OpenFIRE originale. Arrivano alla pistola attraverso la sua porta seriale: la porta USB della pistola, oppure quella del suo dongle quando gioca senza fili. La porta può essere usata da un solo programma alla volta, quindi chiudi prima la WebApp. Avvia la pistola normalmente: nella modalità WebApp offline (B premuto all'avvio) la sua porta seriale USB non è disponibile.


## MAMEHOOKER Setup Guide

[MAMEHOOKER](http://dragonking.arcadecontrols.com/static.php?page=aboutmamehooker) is... *complicated,* as the force feedback capabilities are all controlled from ini files (herein referred to as "configuration files") in `MAMEHOOKER/ini/MAME`. And rather than keeping plaintext codes a *secret*, we maintain extensive documentation here about all the odds and ends of getting game-controlled force feedback for the currently compatible gun systems--primarily focusing on OpenFIRE, though most of this still applies to certain other systems, unless otherwise noted.

###### ([Direct DL link for MameHooker 5.1 in case of issues](http://dragonking.arcadecontrols.com/files/mamehooker5.1.rar))

### Setup:
* On Windows, MAMEHOOKER requires running with Administrative privileges - otherwise, you'll have an abstract error regarding a DLL file.
* On Linux, MAMEHOOKER requires a wineprefix with the `vbrun6` (Visual Basic 6 Runtime) winetricks verb. This will only work with games running *in the same prefix* and, obviously, only works with Windows games run through Wine (afaik, there is no such outputs support in any native Linux program or emulator... for the two programs that would use it, MAME & RetroArch, anyways).

To have MAMEHOOKER "communicate" with the gun, you'll have to make sure the microcontroller's COM port is mapped relative to its intended player number - so if your gun is intended to be Player 1, it has to be COM1, Player 2 is COM2, etc. Mamehook will only see up to COM4, and it's expected that your gun is dedicated to the appropriate slot.

* **To set a COM port on Windows:**
  * While the device is plugged in, go to Device Manager and under the *"Ports (COM & LPT)"* device category, right click on the one that's likely to be your gun (for me, it mapped by default to COM5/6), and click on **Properties**
![Win1](https://puu.sh/JUWfh/1543d0b2fe.png)
  * Navigate to the **Port Settings** tab on the window that pops up, and click on **Advanced...** 
![Win2](https://puu.sh/JUWfq/73de58f681.png)
  * Finally, click on the **COM Port Number** dropdown list, and set it to the player port you'd like your gun to be.
![Win3](https://puu.sh/JUWfu/22c9ddf96c.png)
  * When finished, click **OK** on the Advanced settings window, and **OK** on the driver properties window.

* **To set a COM port on Linux:**
  * The Wineprefix that Mamehook's running in just needs to have these assigned registry keys in the Wine Registry, in `HKEY_LOCAL_MACHINE\Software\Wine\Ports` - where `ttyACM#` is the serial device of the Arduino, and `ttyS0` is built into Linux by default (and normally takes the place of COM1, hence the need to redefine them here).
![WineReg](https://puu.sh/JUWiE/133f88b9ed.png)

If you're using this with outputs provided by DemulShooter (make sure *Outputs* are enabled in `DemulShooter_GUI.exe`), Mamehook *needs to start and have its icon in the system tray first* before DS starts, followed by the game.

On a sidenote: if you happen to only need feedback for lightguns, consider using @SeongGino's own [QMamehook](https://github.com/SeongGino/QMamehook) utility; a minified, modernized background server made to be easier to use and cross-platform. None of the above setup applies to QMH, but the rest of the document does!
___
## How does MAMEHOOKER work?

Games to be hooked into by MAMEHOOKER have an .ini file named after their MAME rom name - or failing that, the shorthand name *of the `-rom`* used by DemulShooter for that game. So if you're configuring for *Terminator 2 - Judgement Day* in MAME, the file used will be `term2.ini`  in the `MAMEHOOKER/ini/MAME` directory - likewise, if you're configuring for *Sega Golden Gun* being run in JConfig w/ DemulShooter, the file used will be `sgg.ini`, to match the DS commandline `demulshooter.exe -target=ringwide -rom=sgg`.

 * MAME provides outputs natively, and are what MAMEHOOKER routes to PC devices, but it has to be manually enabled by changing the following setting in `mame.ini` to:
```
#
# OSD OUTPUT OPTIONS
#
output                    windows
```
 * DemulShooter is needed to provide outputs for PC/Windows-based arcade games and other compatible things. Games in DemulShooter will usually provide a `PX_CtmRecoil` & `PX_Damaged` (PX = player number, i.e. P1 or P2); the former intended for solenoid, the latter intended for rumble.
 * Supermodel & native Flycast *technically* provides outputs natively, but it only provides lamp control and no other significant peripheral functionality.

<a id="how-does-configuration-files-work"></a>

## How do configuration files work?

This is the basic anatomy of a MAMEHOOKER-compatible ini file - games can have either all of or less outputs available, and depends on the game whether they provide that output channel or not, but these are what you'll want to focus on mainly:

```
[General]
MameStart=
MameStop=
StateChange=
OnRotate=
OnPause=
[Output]
P1_LmpStart=
P2_LmpStart=
P1_LmpPanel=
P2_LmpPanel=
P1_GunRecoil=
P2_GunRecoil=
P1_GunMotor=
P2_GunMotor=
P1_Ammo=
P2_Ammo=
P1_Clip=
P2_Clip=
P1_CtmRecoil=
P2_CtmRecoil=
P1_Life=
P2_Life=
P1_Damaged=
P2_Damaged=
```

## `[General]`
We only care about two pieces here:
### `MameStart`
The commands sent when the program is initially hooked and the game is started, and where we establish com ports used, open them, and send the codes to signal that the gun is in serial handoff mode - in that order.

### `MameStop`
The commands sent when the program is closed, where the signals are sent to end serial communication with the gun and closes the COM port - in that order.

The others in this section don't matter, for Lightguns at least.

## `[Output]`
The meat and potatoes of the configuration, this is where we associate commands with the output 'channels' the game uses here.
### `PX_LmpStart/LmpPanel`, `LmpBillboard`
These correlate to the lamps/LEDs of the respective Start Buttons or the Panel of the original arcade cabinet corresponding to the player side, and the marquee lighting. For guns, you'll probably like to set the *LmpStart* for the LEDs.

### `PX_GunRecoil`
Correlates to the respective player's gun's **recoil feedback,** (usually) in the form of a solenoid, *provided natively by the game.* Note that for some games with DemulShooter-provided feedback, this output channel might not have the effect you'll anticipate (*Operation G.H.O.S.T.* being an example where this doesn't work well).

### `PX_GunMotor`
Correlates to the respective player's gun's **rumble motor,** *provided natively by the game.* This can either be the primary form of force feedback in some games (like *Kidou Senshi Gundam: Spirits of Zeon*), or used more conventionally as a subtle feedback during certain game events (as is the case in the *Haunted Museum* games).

### `PX_CtmRecoil`
Correlates to the respective player's gun's ***custom solenoid feedback***, *provided by DemulShooter.* This is the feedback channel for solenoids you'll want to be using if you're playing a native Windows game/arcade dump.

### `PX_Damaged`
Correlates to the respective player's gun's **vibration feedback when they're damaged,** *provided by DemulShooter.* If your game uses `GunMotor` already, you'll probably want to pick either one or the other, or else you won't be able to tell the difference!

### `PX_Ammo/Life`
Correlates to **the number of the respective player's status,** *provided by DemulShooter.* Ammo is for current ammo on the player's weapon, and Life for the player's lives count (0-9+) or life percentage (0-100), depending on the game's life system - i.e. those with *life icons* vs. life bars/percentage out of 100%. This is used mainly for representing these values on OLED displays, if available.

## Opening, Writing to, and Closing COM Ports

Now that you know what the outputs are connected to, how do we use them?

If you'd like to send one-off commands *without* using any utility, this can be done with the `echo x > y` command - where `x` is the string to send, and `y` is the device path to send it to; for device paths, Windows uses `COM#` (where `#` is the gun's COM port, as denoted in Device Manager), and Linux *normally* uses a path like `/dev/ttyACM#` (where `#` is an ascending number assigned by the kernel, from first to last plugged in device).

For MAMEHOOKER, any output we want will be sent to a COM port, and doing that is pretty simple: every command needs to start with either a COM port **OPEN** (`cmo #`), **WRITE** (`cmw #`), or **CLOSE** (`cmc #`) - where # is the port that it's directed to. Commands need to be done one at a time, with `,` (comma) telling MAMEHOOKER to do these in succession.

COM port **OPEN** commands start with `cmo`, followed by the port number, then followed by a string that initializes the com port environment. For simplicity sake,
```
cmo 1 baud=9600_parity=N_data=8_stop=1, cmo 2 baud=9600_parity=N_data=8_stop=1,
```
...is what you'll be using in `MameStart` to open COM1 and COM2. But this won't do anything on its own - we need to actually be able to send commands.

COM port **WRITE** commands start with `cmw`, followed by the port number, then followed by the "message" you'll be sending for the gun to process. The gun needs to be signaled that the PC wants to let it hand control over, so after the `cmo` COM open commands, we write:
```
cmw 1 S, cmw 2 S,
```
...this sends a basic start command to the gun, which signals it to enter serial handoff mode. We do the same, only with a different command, for the `MameStop` line for when we're *ending* communication...

COM port **CLOSE** commands, therefore, does exactly that; starting with `cmc`, followed by the port number, this closes the communication channel for the port. So for `MameStop`, we'll instead use:
```
cmw 1 E, cmw 2 E, cmc 1, cmc 2
```
...this sends the serial END command to the gun, telling it that the PC is done handling force feedback and relents control over FF back to the gun's own routines, and then properly closes the channel thereafter.

## Light Gun Serial Commands

Other gun systems have already established a common 'language', or a standardized set of shorthand strings that a gun looks for and responds to. *OpenFIRE* follows the preexisting syntax, and adds a few extra commands specific to this system; for public documentation purposes, the entire 'language' will be inscribed here, with OpenFIRE-specific notes where needed.

### `S` - Start Commands
Does what it says on the tin, *starts* the serial mode on the gun - OpenFIRE will detect this, flip the `serialMode` switch on, and light the board's builtin &/or external LED (if any) with a mid-intensity gray to indicate this switchover. The START command technically has "modes" that it defines:
 * `S0` - Start with solenoid enabled
 * `S1` - Start with rumble enabled
 * `S2` - Start with the RED LED enabled
 * `S3` - Start with the GREEN LED enabled
 * `S4` - Start with the BLUE LED enabled
 * `S6` - Start with everything enabled

...but it seems more like a *suggestion* than anything else as, again, it's the configuration's job to determine what devices it uses. The only thing OpenFIRE looks for is the `S`; everything else is irrelevant.

### `E` - End Command
Again, does what it says on the tin; signals the gun to exit serial mode and go back to handling its force feedback internally. Any force feedback ongoing will get shut off, and the LED will likewise turn off. There's only one kind of END command:
 * `E` - Ends serial communication

### `M` - Mode Commands
These set different parameters on the gun to tell it how it should behave in the game. The different "modes" are:
 * `M0` - Device Output Mode
   * `x0` - Mouse & Keyboard (on OpenFIRE ESP32 an analog stick used as a gamepad stick also goes back to the left stick, as at startup)
   * `x1` - Gamepad, w/ Camera mapped to Right Stick
     * `L` - Maps Camera to Left Stick instead (**OpenFIRE exclusive**)
   * `x2` - 'Hybrid'
 * `M1` - Offscreen Firing Mode
   * `x0` - Disabled (*not used in OpenFIRE*)
   * `x1` - Fire in bottom-left corner (*not used in OpenFIRE*)
   * `x2` - Offscreen Button Mode enabled (i.e. offscreen trigger pulls generates a Right Click instead of a Left Click)
   * `x3` - Normal shot (*always on when Offscreen Button Mode isn't set in OpenFIRE*)
 * `M2` - Pedal Mapping
   * `x0` - Separate Button (as mapped)
   * `x1` - As Right Mouse
   * `x2` - As Middle Mouse (**OpenFIRE exclusive**)
 * `M3` - Aspect Ratio correction
   * `x0` - Fullscreen
   * `x1` - 4:3 Correction
 * `M4` - Temp Sensor Control (*not used in OpenFIRE*)
   * `x0` - Disabled
   * `x1` - Enabled
 * `M5` - Auto Reload (*not used in OpenFIRE*)
   * `x0` - Disabled
   * `x1` - Enabled
 * `M6` - Rumble Only Mode
   * `x0` - Disabled (Solenoid allowed)
   * `x1` - Enabled (Solenoid disabled, Rumble enabled)
 * `M8` - Autofire Mode
   * `x0` - Disabled (*sustained fire is kept enabled in OpenFIRE*) 
   * `x1` - Auto fire on (*enables Burst Fire in OpenFIRE*)
   * `x2` - Auto fire always on rapid fire
 * `MD` - Display Mode (**OpenFIRE exclusive**)
   * `x1` - Life Only
   * `x2` - Ammo Only
   * `x3` - Life & Ammo Splitscreen
     * `B` - Life Bar (Life Glyphs otherwise)

OpenFIRE detects all of these commands except for Auto Reload and any Offscreen mode other than `M1x2` (Offscreen Button Mode). The `B` in `MDx#` is for determining whether the Life element is a Life Bar (a filled display element from 0-100%) or Life Glyphs (icons that display, from 0-10)

### `R` - Pulse Override Commands
When using force feedback devices in Pulse mode (more in the next section), these set custom overrides on how the pulses behave for the current serial session:
 * `R0` - Solenoid
 * `R1` - Rumble Motor
 * `R2` - RGB Red LED
 * `R3` - RGB Green LED
 * `R4` - RGB Blue LED
   * (For ANY of the above devices:) 
   * `x0x#` - Pulse "ON" length, where # is the amount of time (in ms) the FFB device is engaged ON (draws power) in a pulse event.
   * `x1x#` - Pulse "OFF" length, where # is the amount of time (in ms) the FFB device is engaged OFF (doesn't draw power) between pulses.
   * `x2x#` - "Analog" output, which seems to go unused as far as anyone can tell.

### `F` - Force Feedback Commands
These are the commands that will directly activate/deactivate the gun's peripherals. They are as follows:
 * `F0` - Solenoid
   * `x0` - Off
   * `x1` - On
   * `x2x#` - Pulse (gun handles on/off), where # is the number of pulses to queue (from 0-255). NOTE: a count of 0 is the same as requesting a count of 1, so long as no pulses are currently queued.
 * `F1` - Rumble
   * `x0` - Off
   * `x1` - On
   * `x2x#` - Pulse (gun handles on/off), where # is the number of pulses to queue (from 0-255). NOTE: a count of 0 is the same as requesting a count of 1, so long as no pulses are currently queued.
 * `F2` - RGB Red Channel
 * `F3` - RGB Green Channel
 * `F4` - RGB Blue Channel
   * (For ALL RGB commands:)
   * `x0` - Off
   * `x1x#` - On, where # determines the strength of the color (from 0=off to 255=max)
   * `x2x#` - Pulse (gun handles on/off), where # is the number of pulses to queue (from 0-255). NOTE: a count of 0 is the same as requesting a count of 1, so long as no pulses are currently queued.
 * `FD` - Display Event (**OpenFIRE exclusive**)
   * `Ax#` - New Ammo Count
   * `Lx#` - New Life Count

The peripheral on/off signals are straightforward (keeping in mind that Solenoids will forcefully shut off after a set time of being enabled), but the pulses bit is the amount that the gun itself is supposed to "pulse" the specific peripheral for on its own - or an automated ramping/descending action of each bit, where appropriate. OpenFIRE supports all of these commands, with LED events being reflected in 4-pin LEDs, onboard NeoPixels, and external NeoPixels *that aren't marked as static*.

These commands can be combined - so if you want to have a command both pull the solenoid and start the motor, you'll make that command:
```
cmw # F0x1xF1x1
```

However, these are just for set states - many of these games have different states that each output channel can be in. This is where you'll want to use...

## MAMEHOOKER's Pipe & State Modifiers

For example, you want the `P1_CtmRecoil` to flip the solenoid on and off, depending on if a shot is happening. It's a very straightforward output, but it has two states - *on* and *off*. MAMEHOOKER takes care of both *engaging* and *disabling* peripherals, but there's only one line for each channel, right?

This is where the `%s%` modifier comes in handy; when used, MAMEHOOKER will automatically substitute `%s%` for the on/off bit of that output channel. So our basic solenoid channel would look something like:
```
cmw 1 F0x%s%
```
This will essentially act as *both* `F0x0` (solenoid disable) *and* `F0x1` (solenoid enable), determined by the state it should be in as determined by the game.

That's great and all for the FF, but what about the LEDs? What if you want them to *pulse*? Well, a game doesn't normally provide a `2` state bit (as far as I'm aware), but we *can* substitute what the "on" and "off" functions do without necessarily having to match the state bit sent by the game!

This is why `|` exists; when used, it will indicate to MAMEHOOKER that these are different commands used for different states of this output channel, and we can individually define what these states do. So if you want the green LED to only pulse when the start button lamp is set on, you would do:
```
cmw 1 F3x0|cmw 1 F3x2x1
```
The left side in this case corresponds to "off", the right side to "on". So when the Start button lamp would turn on, it sends a pulse signal to the green RGB emitter channel - otherwise, sends a green RGB off command. Or, if you want the LED to be used to correspond to the player's health, like for DemulShooter's `P1_Health` output channel, you would do:
```
cmw 1 F2x1x255|cmw 1 F2x1x150|cmw 1 F2x1x50|cmw 1 F3x1x100|cmw 1 F3x1x255
```
This means, from left to right, set the Red LED all the way on, semi-on, dimly lit, then the Green LED goes semi-on, then all the way on to correspond to the player having 1 to 5 hitpoints left in that order.

How do you find out what feedback channels are available? DemulShooter will expose all of these when it initially hooks into a compatible game, and MAMEHOOKER/QMamehook will autofill these lines in the configuration file for that game. Failing this, for MAME, feedback channels will be exposed as the game is played, and any new channels found will be recorded into the configuration file.
