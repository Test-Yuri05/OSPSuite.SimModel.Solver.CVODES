#pragma once
#pragma once
#define MAJOR 0
#define MINOR 0
#define RELEASE 0
#define BUILD 0

#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)

#define VER_FILE_VERSION_STR        STRINGIZE(MAJOR)        \
                                    "." STRINGIZE(MINOR)    \
                                    "." STRINGIZE(RELEASE) \
                                    "." STRINGIZE(BUILD) \

#define COMMENTS "© Open Systems Pharmacology Community 2002-2018"
#define COMPANY_NAME "Open Systems Pharmacology Community"
#define FILE_DESCRIPTION "SimModel Solver based on CVODES"
#define FILE_VERSION VER_FILE_VERSION_STR
#define INTERNAL_NAME "SimModel CVODES Solver"
#define LEGAL_COPYRIGHT "Copyright 2002-2018"
#define LEGAL_TRADEMARKS "© Open-Systems-Pharmacology"
#define ORIGINAL_FILENAME "OSPSuite.SimModelSolver_CVODES.dll"
#define PRODUCT_NAME  "Open Systems Pharmacology SimModel CVODES Solver"
#define PRODUCT_VERSION VER_FILE_VERSION_STR
