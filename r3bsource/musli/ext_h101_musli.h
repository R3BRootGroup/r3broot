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

#ifndef __GUARD_H101_EXT_H101_MUSLI_H__
#define __GUARD_H101_EXT_H101_MUSLI_H__

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

typedef struct EXT_STR_h101_MUSLI_t
{
    /* RAW */
    uint32_t MUSLI_EM /* [1,15] */;
    uint32_t MUSLI_EMI[15 EXT_STRUCT_CTRL(MUSLI_EM)] /* [1,15] */;
    uint32_t MUSLI_EME[15 EXT_STRUCT_CTRL(MUSLI_EM)] /* [1,1500] */;
    uint32_t MUSLI_E /* [0,1500] */;
    uint32_t MUSLI_Ev[1500 EXT_STRUCT_CTRL(MUSLI_E)] /* [0,1048575] */;
    uint32_t MUSLI_TM /* [1,15] */;
    uint32_t MUSLI_TMI[15 EXT_STRUCT_CTRL(MUSLI_TM)] /* [1,15] */;
    uint32_t MUSLI_TME[15 EXT_STRUCT_CTRL(MUSLI_TM)] /* [1,1500] */;
    uint32_t MUSLI_T /* [0,1500] */;
    uint32_t MUSLI_Tv[1500 EXT_STRUCT_CTRL(MUSLI_T)] /* [0,1048575] */;
    uint32_t MUSLI_TREFM /* [1,1] */;
    uint32_t MUSLI_TREFMI[1 EXT_STRUCT_CTRL(MUSLI_TREFM)] /* [1,1] */;
    uint32_t MUSLI_TREFME[1 EXT_STRUCT_CTRL(MUSLI_TREFM)] /* [1,100] */;
    uint32_t MUSLI_TREF /* [0,100] */;
    uint32_t MUSLI_TREFv[100 EXT_STRUCT_CTRL(MUSLI_TREF)] /* [0,1048575] */;
    uint32_t MUSLI_TTRIGM /* [1,1] */;
    uint32_t MUSLI_TTRIGMI[1 EXT_STRUCT_CTRL(MUSLI_TTRIGM)] /* [1,1] */;
    uint32_t MUSLI_TTRIGME[1 EXT_STRUCT_CTRL(MUSLI_TTRIGM)] /* [1,100] */;
    uint32_t MUSLI_TTRIG /* [0,100] */;
    uint32_t MUSLI_TTRIGv[100 EXT_STRUCT_CTRL(MUSLI_TTRIG)] /* [0,1048575] */;

} EXT_STR_h101_MUSLI;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_MUSLI_onion_t
{
    /* RAW */
    uint32_t MUSLI_EM;
    uint32_t MUSLI_EMI[15 /* MUSLI_EM */];
    uint32_t MUSLI_EME[15 /* MUSLI_EM */];
    uint32_t MUSLI_E;
    uint32_t MUSLI_Ev[1500 /* MUSLI_E */];
    uint32_t MUSLI_TM;
    uint32_t MUSLI_TMI[15 /* MUSLI_TM */];
    uint32_t MUSLI_TME[15 /* MUSLI_TM */];
    uint32_t MUSLI_T;
    uint32_t MUSLI_Tv[1500 /* MUSLI_T */];
    uint32_t MUSLI_TREFM;
    uint32_t MUSLI_TREFMI[1 /* MUSLI_TREFM */];
    uint32_t MUSLI_TREFME[1 /* MUSLI_TREFM */];
    uint32_t MUSLI_TREF;
    uint32_t MUSLI_TREFv[100 /* MUSLI_TREF */];
    uint32_t MUSLI_TTRIGM;
    uint32_t MUSLI_TTRIGMI[1 /* MUSLI_TTRIGM */];
    uint32_t MUSLI_TTRIGME[1 /* MUSLI_TTRIGM */];
    uint32_t MUSLI_TTRIG;
    uint32_t MUSLI_TTRIGv[100 /* MUSLI_TTRIG */];

} EXT_STR_h101_MUSLI_onion;

/*******************************************************/

#define EXT_STR_h101_MUSLI_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                            \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, MUSLI_EM, UINT32, "MUSLI_EM", 15);                 \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, MUSLI_EMI, UINT32, "MUSLI_EMI", "MUSLI_EM");       \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, MUSLI_EME, UINT32, "MUSLI_EME", "MUSLI_EM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, MUSLI_E, UINT32, "MUSLI_E", 1500);                 \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, MUSLI_Ev, UINT32, "MUSLI_Ev", "MUSLI_E");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, MUSLI_TM, UINT32, "MUSLI_TM", 15);                 \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, MUSLI_TMI, UINT32, "MUSLI_TMI", "MUSLI_TM");       \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, MUSLI_TME, UINT32, "MUSLI_TME", "MUSLI_TM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, MUSLI_T, UINT32, "MUSLI_T", 1500);                 \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, MUSLI_Tv, UINT32, "MUSLI_Tv", "MUSLI_T");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, MUSLI_TREFM, UINT32, "MUSLI_TREFM", 1);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, MUSLI_TREFMI, UINT32, "MUSLI_TREFMI", "MUSLI_TREFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, MUSLI_TREFME, UINT32, "MUSLI_TREFME", "MUSLI_TREFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, MUSLI_TREF, UINT32, "MUSLI_TREF", 100);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, MUSLI_TREFv, UINT32, "MUSLI_TREFv", "MUSLI_TREF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, MUSLI_TTRIGM, UINT32, "MUSLI_TTRIGM", 1);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, MUSLI_TTRIGMI, UINT32, "MUSLI_TTRIGMI", "MUSLI_TTRIGM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, MUSLI_TTRIGME, UINT32, "MUSLI_TTRIGME", "MUSLI_TTRIGM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, MUSLI_TTRIG, UINT32, "MUSLI_TTRIG", 100);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, MUSLI_TTRIGv, UINT32, "MUSLI_TTRIGv", "MUSLI_TTRIG");                \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_EXT_H101_MUSLI_H__*/

/*******************************************************/
