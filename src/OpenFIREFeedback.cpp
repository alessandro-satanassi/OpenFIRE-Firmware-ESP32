 /*!
 * @file OpenFIREFeedback.cpp
 * @brief Force feedback subsystems.
 *
 * @copyright That One Seong, 2025
 * @copyright GNU Lesser General Public License
 */ 

#include <Arduino.h>

#include "OpenFIREFeedback.h"

// ============ 696969 ========== redifinizione di Serial per gestire le connessione wireless seriali ========
#ifdef OPENFIRE_WIRELESS_ENABLE
    extern Stream* Serial_OpenFIRE_Stream;
    #ifdef Serial
        #define AUX_SERIAL Serial
        #undef Serial
    #endif
    #define Serial (*Serial_OpenFIRE_Stream)
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ 696969 ===== fine redifinizione di Serial per gestire le connessione wireless seriali ========


void OF_FFB::FFBOnScreen()
{
    if(OF_Prefs::toggles[OF_Const::solenoid]) {                             // (Only activate when the solenoid switch is on!)
        // If this is the first time we're firing,
        if(!triggerHeld) {
            if(burstFireActive && !burstFiring) {  // Are we in burst firing mode?
                solenoidFirstShot = true;               // Set this so we use the instant solenoid fire path,
                SolenoidActivation(0);                  // Engage it,
                solenoidFirstShot = false;              // And disable the flag to mitigate confusion.
                burstFiring = true;                     // Set that we're in a burst fire event.
                burstFireCount = 1;                     // Set this as the first shot in a burst fire sequence,
                burstFireCountLast = 1;                 // And reset the stored counter,
            } else if(!burstFireActive) {
                solenoidFirstShot = true;
                SolenoidActivation(0);
                if(OF_Prefs::toggles[OF_Const::autofire])
                    solenoidFirstShot = false;
            }
        // Else, these below are all if we've been holding the trigger.
        } else if(burstFiring) { BurstFire();
        } else if(OF_Prefs::toggles[OF_Const::autofire] && !burstFireActive) {
            SolenoidActivation(digitalRead(OF_Prefs::pins[OF_Const::solenoidPin]) ? OF_Prefs::settings[OF_Const::solenoidOnLength] : OF_Prefs::settings[OF_Const::solenoidOffLength] << autofireDoubleLengthWait);
        } else if(solenoidFirstShot) {                  // If we aren't in autofire mode, are we waiting for the initial shot timer still?
            if(digitalRead(OF_Prefs::pins[OF_Const::solenoidPin])) {              // If so, are we still engaged? We need to let it go normally, but maintain the single shot flag.
                currentMillis = millis();
                if(currentMillis - previousMillisSol >= OF_Prefs::settings[OF_Const::solenoidOnLength])
                    digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], LOW);
            } else if(millis() - previousMillisSol >= OF_Prefs::settings[OF_Const::solenoidHoldLength]) { // We're waiting on the extended wait before repeating in single shot mode.
                SolenoidActivation(0); // Process it now.
                    solenoidFirstShot = false;          // We're gonna turn this off so we don't have to pass through this check anymore.
            }
        } else if(!burstFireActive) {                   // if we don't have the single shot wait flag on (holding the trigger w/out autofire)
            SolenoidActivation(digitalRead(OF_Prefs::pins[OF_Const::solenoidPin]) ? OF_Prefs::settings[OF_Const::solenoidOnLength] : OF_Prefs::settings[OF_Const::solenoidOffLength]);
        }
    // only activate rumbleFF as a fallback if Solenoid is explicitly disabled
    } else if(OF_Prefs::toggles[OF_Const::rumble] && OF_Prefs::toggles[OF_Const::rumbleFF] && !rumbleHappened && !triggerHeld)
        RumbleActivation();

    if(OF_Prefs::toggles[OF_Const::rumble] && rumbleHappening && triggerHeld)
        RumbleActivation();                    // Continue processing the rumble command, to prevent infinite rumble while going from on-screen to off mid-command.
}

void OF_FFB::FFBOffScreen()
{
    if(OF_Prefs::toggles[OF_Const::rumble]) {
        // Is this the first time we're rumbling AND only started pulling the trigger (to prevent starting a rumble w/ trigger hold)?
        if(!OF_Prefs::toggles[OF_Const::rumbleFF] && !rumbleHappened && !triggerHeld) {
            RumbleActivation();                        // Start a rumble command.
        } else if(rumbleHappening) {  // We are currently processing a rumble command.
            RumbleActivation();                        // Keep processing that command then.
        }  // Else, we rumbled already, so don't do anything to prevent infinite rumbling.
    }

    if(burstFiring) BurstFire();
    else if(digitalRead(OF_Prefs::pins[OF_Const::solenoidPin]) && !burstFireActive) { // If the solenoid is engaged, since we're not shooting the screen, shut off the solenoid a'la an idle cycle
        SolenoidActivation(OF_Prefs::settings[OF_Const::solenoidOnLength]);
    }
}

void OF_FFB::FFBRelease()
{
    if(OF_Prefs::toggles[OF_Const::solenoid]) {
        if(burstFiring) BurstFire();
        else if(!burstFireActive && digitalRead(OF_Prefs::pins[OF_Const::solenoidPin])) {
            solenoidFirstShot = false;                      // Make sure this is unset to prevent "sticking" in single shot mode!
            currentMillis = millis();
            if(currentMillis - previousMillisSol >= OF_Prefs::settings[OF_Const::solenoidOnLength]) { // I guess if we're not firing, may as well use the fastest shutoff.
                previousMillisSol = currentMillis;
                digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], LOW);             // Make sure to turn it off.
            }
        }
    }
    
    // If Rumble FF is enabled and Autofire is enabled, the motor needs to be disabled when the trigger is released. Otherwise, allow RumbleActivation to deal with the activation timer
    if(OF_Prefs::toggles[OF_Const::rumbleFF] && OF_Prefs::toggles[OF_Const::autofire]) {
        if(rumbleHappening || rumbleHappened) {
            #ifdef ARDUINO_ARCH_ESP32
                analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0);      // Make sure the rumble is OFF. // 696969 per esp32
            #else // rp2040
                digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);      // Make sure the rumble is OFF. // 696969 per esp32
            #endif
            rumbleHappening = false;                                // This rumble command is done now.
            rumbleHappened = false;                                 // Make it clear we've stopped holding.
        }
    } else {
        if(rumbleHappening) {                                   // Are we currently in a rumble command? (Implicitly needs OF_Prefs::toggles[OF_Const::rumble])
            RumbleActivation();                                 // Continue processing our rumble command.
            // (This is to prevent making the lack of trigger pull actually activate a rumble command instead of skipping it like we should.)
        } else if(rumbleHappened) {                             // If rumble has happened,
            rumbleHappened = false;                             // well we're clear now that we've stopped holding.
        }
    }
}

void OF_FFB::SolenoidActivation(const int &solenoidFinalInterval)
{
    if(solenoidFirstShot) {                                       // If this is the first time we're shooting, it's probably safe to shoot regardless of temps.
        previousMillisSol = millis();                             // Calibrate the timer for future calcs.
        digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], HIGH);     // Since we're shooting the first time, just turn it on aaaaand fire.
    } else {
        if(OF_Prefs::pins[OF_Const::tempPin] >= 0) { // If a temp sensor is installed and enabled,
            TemperatureUpdate();

            if(tempStatus < Temp_Fatal) {
                if(tempStatus == Temp_Warning) {
                    if(digitalRead(OF_Prefs::pins[OF_Const::solenoidPin])) {    // Is the valve being pulled now?
                        if(currentMillis - previousMillisSol >= solenoidFinalInterval) {
                            previousMillisSol = currentMillis;
                            digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], !digitalRead(OF_Prefs::pins[OF_Const::solenoidPin])); // Flip, flop.
                        }
                    } else { // The solenoid's probably off right now, so that means we should wait a bit longer to fire again.
                        if(currentMillis - previousMillisSol >= (OF_Prefs::settings[OF_Const::solenoidOffLength] << 2)) { // We're keeping it low for a bit longer, to keep temps stable. Try to give it a bit of time to cool down before we go again.
                            previousMillisSol = currentMillis;
                            digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], !digitalRead(OF_Prefs::pins[OF_Const::solenoidPin]));
                        }
                    }
                } else {
                    if(currentMillis - previousMillisSol >= solenoidFinalInterval) {
                        previousMillisSol = currentMillis;
                        digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], !digitalRead(OF_Prefs::pins[OF_Const::solenoidPin])); // run the solenoid into the state we've just inverted it to.
                    }
                }
            } else {
                #ifdef PRINT_VERBOSE
                    Serial.println("Solenoid over safety threshold; not activating!");
                #endif
                digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], LOW);                       // Make sure it's off if we're this dangerously close to the sun.
            }
        } else { // No temp sensor, so just go ahead.
            currentMillis = millis();
            if(currentMillis - previousMillisSol >= solenoidFinalInterval) { // If we've waited long enough for this interval,
                previousMillisSol = currentMillis;                    // Since we've waited long enough, calibrate the timer
                digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], !digitalRead(OF_Prefs::pins[OF_Const::solenoidPin])); // run the solenoid into the state we've just inverted it to.
            }
        }
    }
}

void OF_FFB::TemperatureUpdate()
{
    currentMillis = millis();
    if(currentMillis - previousMillisTemp > TEMP_UPDATE_INTERVAL) {
        previousMillisTemp = currentMillis;
        temperatureGraph[tempGraphIndex] = (((analogRead(OF_Prefs::pins[OF_Const::tempPin]) * 3.3) / 4096) - 0.5) * 100; // Convert reading from mV->3.3->12-bit->Celsius
        if(tempGraphIndex < 3) {
            tempGraphIndex++;
        } else {
            // average out temperature from four samples taken 3ms apart from each other
            tempGraphIndex = 0;
            temperatureCurrent = (temperatureGraph[0] +
                                  temperatureGraph[1] +
                                  temperatureGraph[2] +
                                  temperatureGraph[3]) >> 2;
                                  
            if(tempStatus == Temp_Fatal) {
                if(temperatureCurrent < OF_Prefs::settings[OF_Const::tempShutdown]-5)
                    tempStatus = Temp_Warning;
            } else {
                if(temperatureCurrent >= OF_Prefs::settings[OF_Const::tempShutdown])
                    tempStatus = Temp_Fatal;
                else if(tempStatus == Temp_Warning) {
                    if(temperatureCurrent < OF_Prefs::settings[OF_Const::tempWarning]-5)
                        tempStatus = Temp_Safe;
                } else if(temperatureCurrent >= OF_Prefs::settings[OF_Const::tempWarning])
                        tempStatus = Temp_Warning;
            }
        }
    }
}

void OF_FFB::RumbleActivation()
{
    if(rumbleHappening) {                                         // Are we in a rumble command rn?
        currentMillis = millis();                                 // Calibrate a timer to set how long we've been rumbling.
        if(OF_Prefs::toggles[OF_Const::rumbleFF]) {
            if(!OF_Prefs::toggles[OF_Const::autofire]) {       // We only want to use the rumble timer if Autofire is not active. Otherwise, keep it going
                if(currentMillis - previousMillisRumble >= (OF_Prefs::settings[OF_Const::rumbleInterval] >> 1)) { // If we've been waiting long enough for this whole rumble command,
                    #ifdef ARDUINO_ARCH_ESP32
                        analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0);                         // Make sure the rumble is OFF. // 696969 per ESP32
                    #else // rp2040
                        digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);                         // Make sure the rumble is OFF. // 696969 per ESP32
                    #endif
                    rumbleHappening = false;                              // This rumble command is done now.
                    rumbleHappened = true;                                // And just to make sure, to prevent holding == repeat rumble commands.
                }
            }
        } else {
            if(currentMillis - previousMillisRumble >= OF_Prefs::settings[OF_Const::rumbleInterval]) { // If we've been waiting long enough for this whole rumble command,
                #ifdef ARDUINO_ARCH_ESP32
                    analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0);                         // Make sure the rumble is OFF. // 696969 per ESP32
                #else // rp2040
                    digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW);                         // Make sure the rumble is OFF. // 696969 per ESP32
                #endif
                rumbleHappening = false;                              // This rumble command is done now.
                rumbleHappened = true;                                // And just to make sure, to prevent holding == repeat rumble commands.
            }
        }
    } else {                                                      // OR, we're rumbling for the first time.
        previousMillisRumble = millis();                          // Mark this as the start of this rumble command.
        analogWrite(OF_Prefs::pins[OF_Const::rumblePin], OF_Prefs::settings[OF_Const::rumbleStrength]);
        rumbleHappening = true;                                   // Mark that we're in a rumble command rn.
    }
}

void OF_FFB::BurstFire()
{
    if(burstFireCount < 4) {  // Are we within the three shots alotted to a burst fire command?
        #ifdef USES_SOLENOID
            if(!digitalRead(OF_Prefs::pins[OF_Const::solenoidPin]) &&  // Is the solenoid NOT on right now, and the counter hasn't matched?
            (burstFireCount == burstFireCountLast)) {
                burstFireCount++;                                 // Increment the counter.
            }
            if(!digitalRead(OF_Prefs::pins[OF_Const::solenoidPin])) {  // Now, is the solenoid NOT on right now?
                SolenoidActivation(OF_Prefs::settings[OF_Const::solenoidOffLength] << autofireDoubleLengthWait ? 1 : 0);     // Hold it off a bit longer,
            } else {                         // or if it IS on,
                burstFireCountLast = burstFireCount;              // sync the counters since we completed one bullet cycle,
                SolenoidActivation(OF_Prefs::settings[OF_Const::solenoidOnLength]);         // And start trying to activate the dingus.
            }
        #endif // USES_SOLENOID
        return;
    } else {  // If we're at three bullets fired,
        burstFiring = false;                                      // Disable the currently firing tag,
        burstFireCount = 0;                                       // And set the count off.
        return;                                                   // Let's go back.
    }
}

void OF_FFB::FFBShutdown()
{
    digitalWrite(OF_Prefs::pins[OF_Const::solenoidPin], LOW);
    #ifdef ARDUINO_ARCH_ESP32
        analogWrite(OF_Prefs::pins[OF_Const::rumblePin], 0); // 696969
    #else // rp2040
        digitalWrite(OF_Prefs::pins[OF_Const::rumblePin], LOW); //696969
    #endif
    solenoidFirstShot = false;
    rumbleHappening = false;
    rumbleHappened = false;
    triggerHeld = false;
    burstFiring = false;
    burstFireCount = 0;
}

// ============ 696969 ========== ripristino di Serial dopo definizione per connessione seriali ==============
#ifdef OPENFIRE_WIRELESS_ENABLE
    #undef Serial
    #ifdef AUX_SERIAL
        #define Serial AUX_SERIAL
        #undef AuxSerial
    #endif
#endif // OPENFIRE_WIRELESS_ENABLE
// ============ 696969 ===== fine ripristino di Serial dopo definizione per connessione seriali ==============