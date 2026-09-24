/*!
 * @file OpenFIREversion.h
 * @brief OpenFIREversion.h
 * @n CPP OpenFIREversion.h
 *
 * @copyright alessandro-satanassi, https://github.com/alessandro-satanassi, 2026
 * @copyright GNU Lesser General Public License
 *
 * @author [Alessandro Satanassi](alessandro@cittini.it)
 * @version V1.0
 * @date 2026
 */

#ifndef OPENFIRE_VERSION_H
#define OPENFIRE_VERSION_H

// Current version / Versione attuale
//
// Questi quattro sono l'unica fonte della versione: tutto il resto e' derivato o
// verificato - la stringa che la lightgun manda quando si aggancia, il nome della
// cartella pubblicata della Web App, il tag di Git e il titolo della release.
// Il suffisso e' facoltativo: vuoto = versione definitiva; altrimenti "beta1",
// "rc1", "prerelease1"... una lettera seguita da lettere e cifre, perche' finisce
// in un nome di cartella, in un indirizzo del sito e nel tag di Git.
// Esempi: 7.0.0  |  7.0.0-beta1  |  7.0.1-rc2
#define OPENFIRE_VERSION_MAJOR 7
#define OPENFIRE_VERSION_MINOR 0
#define OPENFIRE_VERSION_PATCH 0
#define OPENFIRE_VERSION_SUFFIX ""

// Deve valere MAJOR.MINOR. E' il primo campo della prima risposta della lightgun,
// l'unico che la App desktop sa leggere, e per questo resta: non identifica piu'
// la versione, la Web App usa la stringa completa qui sotto.
// scripts/webapp_build.py verifica ad ogni compilazione che i due numeri coincidano.
#define OPENFIRE_VERSION 7.0
/////// #define OPENFIRE_CODENAME "Dawn Sigma rc2"
// #define GIT_HASH

// Con "%.1f" il campo qui sopra non distingue 7.10 da 7.1: MINOR si ferma a 9.
#if OPENFIRE_VERSION_MINOR > 9
#error "OPENFIRE_VERSION_MINOR oltre 9: il campo %.1f mandato alla App desktop non lo distinguerebbe (7.10 diventerebbe 7.1)."
#endif

// Coda della versione: "" oppure "-beta1". Il trattino esiste solo col suffisso,
// cosi' una versione definitiva e' "7.0.0" e non "7.0.0-".
#define OPENFIRE_VERSION_TAIL (OPENFIRE_VERSION_SUFFIX[0] ? "-" OPENFIRE_VERSION_SUFFIX : "")

// Current version string / Stringa di versione attuale: 7.0.0 oppure 7.0.0-beta1
#define OPENFIRE_VERSION_STRING (String(OPENFIRE_VERSION_MAJOR) + "." + String(OPENFIRE_VERSION_MINOR) + "." + String(OPENFIRE_VERSION_PATCH) + OPENFIRE_VERSION_TAIL)
// #define OPENFIRE_VERSION_STRING TU_STRING(OPENFIRE_VERSION_MAJOR) "." TU_STRING(OPENFIRE_VERSION_MINOR) "." TU_STRING(OPENFIRE_VERSION_PATCH)

// Combined version number of the current version / Numero combinato della versione attuale
#define OPENFIRE_VERSION_NUMBER (OPENFIRE_VERSION_MAJOR * 10000 + OPENFIRE_VERSION_MINOR * 100 + OPENFIRE_VERSION_PATCH)

// =========== latest compatible version / ultima versione compatibile ===========

// Last compatible previous version / Ultima versione precedente compatibile
#define OPENFIRE_COMPATIBLE_MAJOR 4
#define OPENFIRE_COMPATIBLE_MINOR 3
#define OPENFIRE_COMPATIBLE_PATCH 2

// Previous compatible version string / Stringa di versione compatibile precedente
#define OPENFIRE_COMPATIBLE_STRING (String(OPENFIRE_COMPATIBLE_MAJOR) + "." + String(OPENFIRE_COMPATIBLE_MINOR) + "." + String(OPENFIRE_COMPATIBLE_PATCH))
// #define OPENFIRE_COMPATIBLE_STRING TU_STRING(OPENFIRE_COMPATIBLE_MAJOR) "." TU_STRING(OPENFIRE_COMPATIBLE_MINOR) "." TU_STRING(OPENFIRE_COMPATIBLE_PATCH)

// Combined version number of the previous compatible version / Numero combinato della versione precedente compatibile
#define OPENFIRE_COMPATIBLE_NUMBER (OPENFIRE_COMPATIBLE_MAJOR * 10000 + OPENFIRE_COMPATIBLE_MINOR * 100 + OPENFIRE_COMPATIBLE_PATCH)

#endif  // OPENFIRE_VERSION_H
