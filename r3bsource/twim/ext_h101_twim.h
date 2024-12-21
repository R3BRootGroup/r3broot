/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2025 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

/********************************************************
 *
 * Structure for ext_data_fetch_event() filling.
 *
 * Do not edit - automatically generated.
 */

#ifndef __GUARD_H101_SOFTWIM_EXT_H101_SOFTWIM_H__
#define __GUARD_H101_SOFTWIM_EXT_H101_SOFTWIM_H__

#ifndef __CINT__
#include <stdint.h>
#else
/* For CINT (old version trouble with stdint.h): */
#ifndef uint32_t
typedef unsigned int uint32_t;
typedef int int32_t;
#endif
#endif
#ifndef EXT_STRUCT_CTRL
#define EXT_STRUCT_CTRL(x)
#endif

/********************************************************
 *
 * Plain structure (layout as ntuple/root file):
 */

typedef struct EXT_STR_h101_SOFTWIM_t
{
    /* RAW */
    uint32_t SOFTWIM_S1EM /* [1,16] */;
    uint32_t SOFTWIM_S1EMI[16 EXT_STRUCT_CTRL(SOFTWIM_S1EM)] /* [1,16] */;
    uint32_t SOFTWIM_S1EME[16 EXT_STRUCT_CTRL(SOFTWIM_S1EM)] /* [1,1600] */;
    uint32_t SOFTWIM_S1E /* [0,1600] */;
    uint32_t SOFTWIM_S1Ev[1600 EXT_STRUCT_CTRL(SOFTWIM_S1E)] /* [0,1048575] */;
    uint32_t SOFTWIM_S1TM /* [1,16] */;
    uint32_t SOFTWIM_S1TMI[16 EXT_STRUCT_CTRL(SOFTWIM_S1TM)] /* [1,16] */;
    uint32_t SOFTWIM_S1TME[16 EXT_STRUCT_CTRL(SOFTWIM_S1TM)] /* [1,1600] */;
    uint32_t SOFTWIM_S1T /* [0,1600] */;
    uint32_t SOFTWIM_S1Tv[1600 EXT_STRUCT_CTRL(SOFTWIM_S1T)] /* [0,1048575] */;
    uint32_t SOFTWIM_S1TREFM /* [1,1] */;
    uint32_t SOFTWIM_S1TREFMI[1 EXT_STRUCT_CTRL(SOFTWIM_S1TREFM)] /* [1,1] */;
    uint32_t SOFTWIM_S1TREFME[1 EXT_STRUCT_CTRL(SOFTWIM_S1TREFM)] /* [1,100] */;
    uint32_t SOFTWIM_S1TREF /* [0,100] */;
    uint32_t SOFTWIM_S1TREFv[100 EXT_STRUCT_CTRL(SOFTWIM_S1TREF)] /* [0,1048575] */;
    uint32_t SOFTWIM_S1TTRIGM /* [1,1] */;
    uint32_t SOFTWIM_S1TTRIGMI[1 EXT_STRUCT_CTRL(SOFTWIM_S1TTRIGM)] /* [1,1] */;
    uint32_t SOFTWIM_S1TTRIGME[1 EXT_STRUCT_CTRL(SOFTWIM_S1TTRIGM)] /* [1,100] */;
    uint32_t SOFTWIM_S1TTRIG /* [0,100] */;
    uint32_t SOFTWIM_S1TTRIGv[100 EXT_STRUCT_CTRL(SOFTWIM_S1TTRIG)] /* [0,1048575] */;
    uint32_t SOFTWIM_S2EM /* [1,16] */;
    uint32_t SOFTWIM_S2EMI[16 EXT_STRUCT_CTRL(SOFTWIM_S2EM)] /* [1,16] */;
    uint32_t SOFTWIM_S2EME[16 EXT_STRUCT_CTRL(SOFTWIM_S2EM)] /* [1,1600] */;
    uint32_t SOFTWIM_S2E /* [0,1600] */;
    uint32_t SOFTWIM_S2Ev[1600 EXT_STRUCT_CTRL(SOFTWIM_S2E)] /* [0,1048575] */;
    uint32_t SOFTWIM_S2TM /* [1,16] */;
    uint32_t SOFTWIM_S2TMI[16 EXT_STRUCT_CTRL(SOFTWIM_S2TM)] /* [1,16] */;
    uint32_t SOFTWIM_S2TME[16 EXT_STRUCT_CTRL(SOFTWIM_S2TM)] /* [1,1600] */;
    uint32_t SOFTWIM_S2T /* [0,1600] */;
    uint32_t SOFTWIM_S2Tv[1600 EXT_STRUCT_CTRL(SOFTWIM_S2T)] /* [0,1048575] */;
    uint32_t SOFTWIM_S2TREFM /* [1,1] */;
    uint32_t SOFTWIM_S2TREFMI[1 EXT_STRUCT_CTRL(SOFTWIM_S2TREFM)] /* [1,1] */;
    uint32_t SOFTWIM_S2TREFME[1 EXT_STRUCT_CTRL(SOFTWIM_S2TREFM)] /* [1,100] */;
    uint32_t SOFTWIM_S2TREF /* [0,100] */;
    uint32_t SOFTWIM_S2TREFv[100 EXT_STRUCT_CTRL(SOFTWIM_S2TREF)] /* [0,1048575] */;
    uint32_t SOFTWIM_S2TTRIGM /* [1,1] */;
    uint32_t SOFTWIM_S2TTRIGMI[1 EXT_STRUCT_CTRL(SOFTWIM_S2TTRIGM)] /* [1,1] */;
    uint32_t SOFTWIM_S2TTRIGME[1 EXT_STRUCT_CTRL(SOFTWIM_S2TTRIGM)] /* [1,100] */;
    uint32_t SOFTWIM_S2TTRIG /* [0,100] */;
    uint32_t SOFTWIM_S2TTRIGv[100 EXT_STRUCT_CTRL(SOFTWIM_S2TTRIG)] /* [0,1048575] */;
    uint32_t SOFTWIM_S3EM /* [1,16] */;
    uint32_t SOFTWIM_S3EMI[16 EXT_STRUCT_CTRL(SOFTWIM_S3EM)] /* [1,16] */;
    uint32_t SOFTWIM_S3EME[16 EXT_STRUCT_CTRL(SOFTWIM_S3EM)] /* [1,1600] */;
    uint32_t SOFTWIM_S3E /* [0,1600] */;
    uint32_t SOFTWIM_S3Ev[1600 EXT_STRUCT_CTRL(SOFTWIM_S3E)] /* [0,1048575] */;
    uint32_t SOFTWIM_S3TM /* [1,16] */;
    uint32_t SOFTWIM_S3TMI[16 EXT_STRUCT_CTRL(SOFTWIM_S3TM)] /* [1,16] */;
    uint32_t SOFTWIM_S3TME[16 EXT_STRUCT_CTRL(SOFTWIM_S3TM)] /* [1,1600] */;
    uint32_t SOFTWIM_S3T /* [0,1600] */;
    uint32_t SOFTWIM_S3Tv[1600 EXT_STRUCT_CTRL(SOFTWIM_S3T)] /* [0,1048575] */;
    uint32_t SOFTWIM_S3TREFM /* [1,1] */;
    uint32_t SOFTWIM_S3TREFMI[1 EXT_STRUCT_CTRL(SOFTWIM_S3TREFM)] /* [1,1] */;
    uint32_t SOFTWIM_S3TREFME[1 EXT_STRUCT_CTRL(SOFTWIM_S3TREFM)] /* [1,100] */;
    uint32_t SOFTWIM_S3TREF /* [0,100] */;
    uint32_t SOFTWIM_S3TREFv[100 EXT_STRUCT_CTRL(SOFTWIM_S3TREF)] /* [0,1048575] */;
    uint32_t SOFTWIM_S3TTRIGM /* [1,1] */;
    uint32_t SOFTWIM_S3TTRIGMI[1 EXT_STRUCT_CTRL(SOFTWIM_S3TTRIGM)] /* [1,1] */;
    uint32_t SOFTWIM_S3TTRIGME[1 EXT_STRUCT_CTRL(SOFTWIM_S3TTRIGM)] /* [1,100] */;
    uint32_t SOFTWIM_S3TTRIG /* [0,100] */;
    uint32_t SOFTWIM_S3TTRIGv[100 EXT_STRUCT_CTRL(SOFTWIM_S3TTRIG)] /* [0,1048575] */;
    uint32_t SOFTWIM_S4EM /* [1,16] */;
    uint32_t SOFTWIM_S4EMI[16 EXT_STRUCT_CTRL(SOFTWIM_S4EM)] /* [1,16] */;
    uint32_t SOFTWIM_S4EME[16 EXT_STRUCT_CTRL(SOFTWIM_S4EM)] /* [1,1600] */;
    uint32_t SOFTWIM_S4E /* [0,1600] */;
    uint32_t SOFTWIM_S4Ev[1600 EXT_STRUCT_CTRL(SOFTWIM_S4E)] /* [0,1048575] */;
    uint32_t SOFTWIM_S4TM /* [1,16] */;
    uint32_t SOFTWIM_S4TMI[16 EXT_STRUCT_CTRL(SOFTWIM_S4TM)] /* [1,16] */;
    uint32_t SOFTWIM_S4TME[16 EXT_STRUCT_CTRL(SOFTWIM_S4TM)] /* [1,1600] */;
    uint32_t SOFTWIM_S4T /* [0,1600] */;
    uint32_t SOFTWIM_S4Tv[1600 EXT_STRUCT_CTRL(SOFTWIM_S4T)] /* [0,1048575] */;
    uint32_t SOFTWIM_S4TREFM /* [1,1] */;
    uint32_t SOFTWIM_S4TREFMI[1 EXT_STRUCT_CTRL(SOFTWIM_S4TREFM)] /* [1,1] */;
    uint32_t SOFTWIM_S4TREFME[1 EXT_STRUCT_CTRL(SOFTWIM_S4TREFM)] /* [1,100] */;
    uint32_t SOFTWIM_S4TREF /* [0,100] */;
    uint32_t SOFTWIM_S4TREFv[100 EXT_STRUCT_CTRL(SOFTWIM_S4TREF)] /* [0,1048575] */;
    uint32_t SOFTWIM_S4TTRIGM /* [1,1] */;
    uint32_t SOFTWIM_S4TTRIGMI[1 EXT_STRUCT_CTRL(SOFTWIM_S4TTRIGM)] /* [1,1] */;
    uint32_t SOFTWIM_S4TTRIGME[1 EXT_STRUCT_CTRL(SOFTWIM_S4TTRIGM)] /* [1,100] */;
    uint32_t SOFTWIM_S4TTRIG /* [0,100] */;
    uint32_t SOFTWIM_S4TTRIGv[100 EXT_STRUCT_CTRL(SOFTWIM_S4TTRIG)] /* [0,1048575] */;

} EXT_STR_h101_SOFTWIM;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_SOFTWIM_onion_t
{
    /* RAW */
    struct
    {
        uint32_t EM;
        uint32_t EMI[16 /* EM */];
        uint32_t EME[16 /* EM */];
        uint32_t E;
        uint32_t Ev[1600 /* E */];
        uint32_t TM;
        uint32_t TMI[16 /* TM */];
        uint32_t TME[16 /* TM */];
        uint32_t T;
        uint32_t Tv[1600 /* T */];
        uint32_t TREFM;
        uint32_t TREFMI[1 /* TREFM */];
        uint32_t TREFME[1 /* TREFM */];
        uint32_t TREF;
        uint32_t TREFv[100 /* TREF */];
        uint32_t TTRIGM;
        uint32_t TTRIGMI[1 /* TTRIGM */];
        uint32_t TTRIGME[1 /* TTRIGM */];
        uint32_t TTRIG;
        uint32_t TTRIGv[100 /* TTRIG */];
    } SOFTWIM_S[4];

} EXT_STR_h101_SOFTWIM_onion;

/*******************************************************/

#define EXT_STR_h101_SOFTWIM_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                          \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S1EM, UINT32, "SOFTWIM_S1EM", 16);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1EMI, UINT32, "SOFTWIM_S1EMI", "SOFTWIM_S1EM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1EME, UINT32, "SOFTWIM_S1EME", "SOFTWIM_S1EM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S1E, UINT32, "SOFTWIM_S1E", 1600);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1Ev, UINT32, "SOFTWIM_S1Ev", "SOFTWIM_S1E");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S1TM, UINT32, "SOFTWIM_S1TM", 16);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1TMI, UINT32, "SOFTWIM_S1TMI", "SOFTWIM_S1TM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1TME, UINT32, "SOFTWIM_S1TME", "SOFTWIM_S1TM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S1T, UINT32, "SOFTWIM_S1T", 1600);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1Tv, UINT32, "SOFTWIM_S1Tv", "SOFTWIM_S1T");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S1TREFM, UINT32, "SOFTWIM_S1TREFM", 1);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1TREFMI, UINT32, "SOFTWIM_S1TREFMI", "SOFTWIM_S1TREFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1TREFME, UINT32, "SOFTWIM_S1TREFME", "SOFTWIM_S1TREFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S1TREF, UINT32, "SOFTWIM_S1TREF", 100);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1TREFv, UINT32, "SOFTWIM_S1TREFv", "SOFTWIM_S1TREF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S1TTRIGM, UINT32, "SOFTWIM_S1TTRIGM", 1);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1TTRIGMI, UINT32, "SOFTWIM_S1TTRIGMI", "SOFTWIM_S1TTRIGM"); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1TTRIGME, UINT32, "SOFTWIM_S1TTRIGME", "SOFTWIM_S1TTRIGM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S1TTRIG, UINT32, "SOFTWIM_S1TTRIG", 100);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S1TTRIGv, UINT32, "SOFTWIM_S1TTRIGv", "SOFTWIM_S1TTRIG");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S2EM, UINT32, "SOFTWIM_S2EM", 16);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2EMI, UINT32, "SOFTWIM_S2EMI", "SOFTWIM_S2EM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2EME, UINT32, "SOFTWIM_S2EME", "SOFTWIM_S2EM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S2E, UINT32, "SOFTWIM_S2E", 1600);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2Ev, UINT32, "SOFTWIM_S2Ev", "SOFTWIM_S2E");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S2TM, UINT32, "SOFTWIM_S2TM", 16);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2TMI, UINT32, "SOFTWIM_S2TMI", "SOFTWIM_S2TM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2TME, UINT32, "SOFTWIM_S2TME", "SOFTWIM_S2TM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S2T, UINT32, "SOFTWIM_S2T", 1600);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2Tv, UINT32, "SOFTWIM_S2Tv", "SOFTWIM_S2T");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S2TREFM, UINT32, "SOFTWIM_S2TREFM", 1);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2TREFMI, UINT32, "SOFTWIM_S2TREFMI", "SOFTWIM_S2TREFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2TREFME, UINT32, "SOFTWIM_S2TREFME", "SOFTWIM_S2TREFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S2TREF, UINT32, "SOFTWIM_S2TREF", 100);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2TREFv, UINT32, "SOFTWIM_S2TREFv", "SOFTWIM_S2TREF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S2TTRIGM, UINT32, "SOFTWIM_S2TTRIGM", 1);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2TTRIGMI, UINT32, "SOFTWIM_S2TTRIGMI", "SOFTWIM_S2TTRIGM"); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2TTRIGME, UINT32, "SOFTWIM_S2TTRIGME", "SOFTWIM_S2TTRIGM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S2TTRIG, UINT32, "SOFTWIM_S2TTRIG", 100);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S2TTRIGv, UINT32, "SOFTWIM_S2TTRIGv", "SOFTWIM_S2TTRIG");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S3EM, UINT32, "SOFTWIM_S3EM", 16);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3EMI, UINT32, "SOFTWIM_S3EMI", "SOFTWIM_S3EM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3EME, UINT32, "SOFTWIM_S3EME", "SOFTWIM_S3EM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S3E, UINT32, "SOFTWIM_S3E", 1600);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3Ev, UINT32, "SOFTWIM_S3Ev", "SOFTWIM_S3E");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S3TM, UINT32, "SOFTWIM_S3TM", 16);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3TMI, UINT32, "SOFTWIM_S3TMI", "SOFTWIM_S3TM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3TME, UINT32, "SOFTWIM_S3TME", "SOFTWIM_S3TM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S3T, UINT32, "SOFTWIM_S3T", 1600);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3Tv, UINT32, "SOFTWIM_S3Tv", "SOFTWIM_S3T");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S3TREFM, UINT32, "SOFTWIM_S3TREFM", 1);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3TREFMI, UINT32, "SOFTWIM_S3TREFMI", "SOFTWIM_S3TREFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3TREFME, UINT32, "SOFTWIM_S3TREFME", "SOFTWIM_S3TREFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S3TREF, UINT32, "SOFTWIM_S3TREF", 100);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3TREFv, UINT32, "SOFTWIM_S3TREFv", "SOFTWIM_S3TREF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S3TTRIGM, UINT32, "SOFTWIM_S3TTRIGM", 1);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3TTRIGMI, UINT32, "SOFTWIM_S3TTRIGMI", "SOFTWIM_S3TTRIGM"); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3TTRIGME, UINT32, "SOFTWIM_S3TTRIGME", "SOFTWIM_S3TTRIGM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S3TTRIG, UINT32, "SOFTWIM_S3TTRIG", 100);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S3TTRIGv, UINT32, "SOFTWIM_S3TTRIGv", "SOFTWIM_S3TTRIG");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S4EM, UINT32, "SOFTWIM_S4EM", 16);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4EMI, UINT32, "SOFTWIM_S4EMI", "SOFTWIM_S4EM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4EME, UINT32, "SOFTWIM_S4EME", "SOFTWIM_S4EM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S4E, UINT32, "SOFTWIM_S4E", 1600);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4Ev, UINT32, "SOFTWIM_S4Ev", "SOFTWIM_S4E");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S4TM, UINT32, "SOFTWIM_S4TM", 16);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4TMI, UINT32, "SOFTWIM_S4TMI", "SOFTWIM_S4TM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4TME, UINT32, "SOFTWIM_S4TME", "SOFTWIM_S4TM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S4T, UINT32, "SOFTWIM_S4T", 1600);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4Tv, UINT32, "SOFTWIM_S4Tv", "SOFTWIM_S4T");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S4TREFM, UINT32, "SOFTWIM_S4TREFM", 1);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4TREFMI, UINT32, "SOFTWIM_S4TREFMI", "SOFTWIM_S4TREFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4TREFME, UINT32, "SOFTWIM_S4TREFME", "SOFTWIM_S4TREFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S4TREF, UINT32, "SOFTWIM_S4TREF", 100);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4TREFv, UINT32, "SOFTWIM_S4TREFv", "SOFTWIM_S4TREF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S4TTRIGM, UINT32, "SOFTWIM_S4TTRIGM", 1);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4TTRIGMI, UINT32, "SOFTWIM_S4TTRIGMI", "SOFTWIM_S4TTRIGM"); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4TTRIGME, UINT32, "SOFTWIM_S4TTRIGME", "SOFTWIM_S4TTRIGM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFTWIM_S4TTRIG, UINT32, "SOFTWIM_S4TTRIG", 100);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SOFTWIM_S4TTRIGv, UINT32, "SOFTWIM_S4TTRIGv", "SOFTWIM_S4TTRIG");    \
                                                                                                                     \
    } while (0);
#endif /*__GUARD_H101_SOFTWIM_EXT_H101_SOFTWIM_H__*/

/*******************************************************/
