# OpenFIRE-Firmware Wiki Main Page

**OpenFIRE ESP32 7.0.0:** these game/serial commands remain separate from the new configuration App protocol. Disconnect the configuration App before using the same serial port. Start the lightgun normally for direct USB serial: the special offline WebApp USB NCM mode replaces that port with a network interface. The paired dongle retains its serial port.

**OpenFIRE ESP32 7.0.0:** questi comandi di gioco/seriali restano separati dal nuovo protocollo dell'App di configurazione. Disconnetti l'App prima di usare la stessa porta seriale. Per la seriale USB diretta avvia normalmente la lightgun: la modalità speciale WebApp offline USB NCM sostituisce quella porta con una rete. Il dongle associato mantiene la sua porta seriale.


If all you're looking for is information about MAMEHOOKER integration, [find it here!](MAMEHOOKER-Documentation.md)

Here is where you'll find related miscellanea regarding OpenFIRE (as it's needed).

## OpenFIRE Commands - making the most of your guns!

Listed below is a chart of all the common serial commands OpenFIRE recognizes that can be used outside of MAMEHOOKER-type integration. These can be fired one-off through either the command line on Windows or Linux, or in a batch script of some kind, to change some aspects of any lightgun.

To send these commands, use:
 - **Windows:** `echo commandsHere > COM#`, where # is the assigned COM port of the lightgun, according to Windows' device manager - these are *fixed values*.
 - **Linux:** `echo commandsHere > /dev/ttyACM#`, where # is the number assigned to the lightgun's serial port object, ordered from first to last microcontroller plugged in - these are *dynamic values* and can change depending on the order the guns have been plugged in/recognized by the kernel.

OpenFIRE commands are as follows:
 * `XR#` - Remap to (#) Player, where # is the desired player number. Guns default to the respective keyboard binds for the Basic TinyUSB Identifier selected in the OpenFIRE App, or else falling back to Player 1 binds `1` & `5` for Start & Select for any custom Product ID.
 * `XI#` - Sets Autofire interval (either set via hardware switch, or from an `M8x2`/`M8x1` command); 0 for normal OFF wait lengths, 1 for doubled OFF wait lengths.
 
