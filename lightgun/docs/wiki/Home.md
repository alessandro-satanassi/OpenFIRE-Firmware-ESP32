# OpenFIRE-Firmware Wiki Main Page

> **OpenFIRE ESP32 7.0.0:** these commands work exactly as in the original OpenFIRE. They reach the gun through its serial port: the gun's own USB port, or its dongle's port when it plays wirelessly. Only one program can use the port at a time, so close the WebApp first. Start the gun normally: in the offline WebApp mode (B held at startup) its USB serial port is not available.

> **OpenFIRE ESP32 7.0.0:** questi comandi funzionano esattamente come nell'OpenFIRE originale. Arrivano alla pistola attraverso la sua porta seriale: la porta USB della pistola, oppure quella del suo dongle quando gioca senza fili. La porta può essere usata da un solo programma alla volta, quindi chiudi prima la WebApp. Avvia la pistola normalmente: nella modalità WebApp offline (B premuto all'avvio) la sua porta seriale USB non è disponibile.


If all you're looking for is information about MAMEHOOKER integration, [find it here!](MAMEHOOKER-Documentation.md)

Use the commands below when a game launcher or script needs to change the player number or autofire interval. For normal setup, use the WebApp.

## OpenFIRE Commands - making the most of your guns!

Listed below is a chart of all the common serial commands OpenFIRE recognizes that can be used outside of MAMEHOOKER-type integration. These can be fired one-off through either the command line on Windows or Linux, or in a batch script of some kind, to change some aspects of any lightgun.

To send these commands, use:
 - **Windows:** `echo commandsHere > COM#`, where # is the assigned COM port of the lightgun, according to Windows' device manager - these are *fixed values*.
 - **Linux:** `echo commandsHere > /dev/ttyACM#`, where # is the number assigned to the lightgun's serial port object, ordered from first to last microcontroller plugged in - these are *dynamic values* and can change depending on the order the guns have been plugged in/recognized by the kernel.

OpenFIRE commands are as follows:
 * `XR#` - Remap to (#) Player, where # is the desired player number. Guns default to the keyboard binds of the player number chosen in the WebApp (**Gun Settings → TinyUSB Identifier**), falling back to the Player 1 binds `1` & `5` for Start & Select with a custom Product ID.
 * `XI#` - Sets Autofire interval (either set via hardware switch, or from an `M8x2`/`M8x1` command); 0 for normal OFF wait lengths, 1 for doubled OFF wait lengths.
 
