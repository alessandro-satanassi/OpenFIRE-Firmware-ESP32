#ifndef OPENFIRE_VERSION_H
#define OPENFIRE_VERSION_H

// Versione attuale
#define OPENFIRE_VERSION 6.0
///////#define OPENFIRE_CODENAME "Dawn Sigma rc2"
//#define GIT_HASH
#define OPENFIRE_VERSION_MAJOR 6
#define OPENFIRE_VERSION_MINOR 0
#define OPENFIRE_VERSION_PATCH 0
#define OPENFIRE_VERSION_TYPE "rc2"  // Per indicare pre-release (alpha, beta, RC(Release Candidate), stable)

// Stringa di versione attuale
#define OPENFIRE_VERSION_STRING (String(OPENFIRE_VERSION_MAJOR) + "." + String(OPENFIRE_VERSION_MINOR) + "." + String(OPENFIRE_VERSION_PATCH)+ "-" + OPENFIRE_VERSION_TYPE)
//#define OPENFIRE_VERSION_STRING TU_STRING(OPENFIRE_VERSION_MAJOR) "." TU_STRING(OPENFIRE_VERSION_MINOR) "." TU_STRING(OPENFIRE_VERSION_PATCH)

// Numero combinato della versione attuale
#define OPENFIRE_VERSION_NUMBER (OPENFIRE_VERSION_MAJOR * 10000 + OPENFIRE_VERSION_MINOR * 100 + OPENFIRE_VERSION_PATCH)

// =============== ultima versione compatibile =================================

// Ultima versione precedente compatibile
#define OPENFIRE_COMPATIBLE_MAJOR 4
#define OPENFIRE_COMPATIBLE_MINOR 3
#define OPENFIRE_COMPATIBLE_PATCH 2

// Stringa di versione compatibile precedente
#define OPENFIRE_COMPATIBLE_STRING (String(OPENFIRE_COMPATIBLE_MAJOR) + "." + String(OPENFIRE_COMPATIBLE_MINOR) + "." + String(OPENFIRE_COMPATIBLE_PATCH))
//#define OPENFIRE_COMPATIBLE_STRING TU_STRING(OPENFIRE_COMPATIBLE_MAJOR) "." TU_STRING(OPENFIRE_COMPATIBLE_MINOR) "." TU_STRING(OPENFIRE_COMPATIBLE_PATCH)

// Numero combinato della versione precedente compatibile
#define OPENFIRE_COMPATIBLE_NUMBER (OPENFIRE_COMPATIBLE_MAJOR * 10000 + OPENFIRE_COMPATIBLE_MINOR * 100 + OPENFIRE_COMPATIBLE_PATCH)

#endif  // OPENFIRE_VERSION_H
