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

#ifndef __GUARD_H101_FIBZERO_EXT_H101_FIBZERO_H__
#define __GUARD_H101_FIBZERO_EXT_H101_FIBZERO_H__

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

typedef struct EXT_STR_h101_FIBZERO_t
{
    /* RAW */
    uint32_t FIBZERO_TMLCM /* [1,256] */;
    uint32_t FIBZERO_TMLCMI[256 EXT_STRUCT_CTRL(FIBZERO_TMLCM)] /* [1,256] */;
    uint32_t FIBZERO_TMLCME[256 EXT_STRUCT_CTRL(FIBZERO_TMLCM)] /* [1,5120] */;
    uint32_t FIBZERO_TMLC /* [0,5120] */;
    uint32_t FIBZERO_TMLCv[5120 EXT_STRUCT_CTRL(FIBZERO_TMLC)] /* [0,65535] */;
    uint32_t FIBZERO_TMLFM /* [1,256] */;
    uint32_t FIBZERO_TMLFMI[256 EXT_STRUCT_CTRL(FIBZERO_TMLFM)] /* [1,256] */;
    uint32_t FIBZERO_TMLFME[256 EXT_STRUCT_CTRL(FIBZERO_TMLFM)] /* [1,5120] */;
    uint32_t FIBZERO_TMLF /* [0,5120] */;
    uint32_t FIBZERO_TMLFv[5120 EXT_STRUCT_CTRL(FIBZERO_TMLF)] /* [0,65535] */;
    uint32_t FIBZERO_TMTCM /* [1,256] */;
    uint32_t FIBZERO_TMTCMI[256 EXT_STRUCT_CTRL(FIBZERO_TMTCM)] /* [1,256] */;
    uint32_t FIBZERO_TMTCME[256 EXT_STRUCT_CTRL(FIBZERO_TMTCM)] /* [1,5120] */;
    uint32_t FIBZERO_TMTC /* [0,5120] */;
    uint32_t FIBZERO_TMTCv[5120 EXT_STRUCT_CTRL(FIBZERO_TMTC)] /* [0,65535] */;
    uint32_t FIBZERO_TMTFM /* [1,256] */;
    uint32_t FIBZERO_TMTFMI[256 EXT_STRUCT_CTRL(FIBZERO_TMTFM)] /* [1,256] */;
    uint32_t FIBZERO_TMTFME[256 EXT_STRUCT_CTRL(FIBZERO_TMTFM)] /* [1,5120] */;
    uint32_t FIBZERO_TMTF /* [0,5120] */;
    uint32_t FIBZERO_TMTFv[5120 EXT_STRUCT_CTRL(FIBZERO_TMTF)] /* [0,65535] */;

} EXT_STR_h101_FIBZERO;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBZERO_onion_t
{
    /* RAW */
    uint32_t FIBZERO_TMLCM;
    uint32_t FIBZERO_TMLCMI[256 /* FIBZERO_TMLCM */];
    uint32_t FIBZERO_TMLCME[256 /* FIBZERO_TMLCM */];
    uint32_t FIBZERO_TMLC;
    uint32_t FIBZERO_TMLCv[5120 /* FIBZERO_TMLC */];
    uint32_t FIBZERO_TMLFM;
    uint32_t FIBZERO_TMLFMI[256 /* FIBZERO_TMLFM */];
    uint32_t FIBZERO_TMLFME[256 /* FIBZERO_TMLFM */];
    uint32_t FIBZERO_TMLF;
    uint32_t FIBZERO_TMLFv[5120 /* FIBZERO_TMLF */];
    uint32_t FIBZERO_TMTCM;
    uint32_t FIBZERO_TMTCMI[256 /* FIBZERO_TMTCM */];
    uint32_t FIBZERO_TMTCME[256 /* FIBZERO_TMTCM */];
    uint32_t FIBZERO_TMTC;
    uint32_t FIBZERO_TMTCv[5120 /* FIBZERO_TMTC */];
    uint32_t FIBZERO_TMTFM;
    uint32_t FIBZERO_TMTFMI[256 /* FIBZERO_TMTFM */];
    uint32_t FIBZERO_TMTFME[256 /* FIBZERO_TMTFM */];
    uint32_t FIBZERO_TMTF;
    uint32_t FIBZERO_TMTFv[5120 /* FIBZERO_TMTF */];

} EXT_STR_h101_FIBZERO_onion;

/*******************************************************/

#define EXT_STR_h101_FIBZERO_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                     \
    do                                                                                                          \
    {                                                                                                           \
        ok = 1;                                                                                                 \
        /* RAW */                                                                                               \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZERO_TMLCM, UINT32, "FIBZERO_TMLCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMLCMI, UINT32, "FIBZERO_TMLCMI", "FIBZERO_TMLCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMLCME, UINT32, "FIBZERO_TMLCME", "FIBZERO_TMLCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZERO_TMLC, UINT32, "FIBZERO_TMLC", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMLCv, UINT32, "FIBZERO_TMLCv", "FIBZERO_TMLC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZERO_TMLFM, UINT32, "FIBZERO_TMLFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMLFMI, UINT32, "FIBZERO_TMLFMI", "FIBZERO_TMLFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMLFME, UINT32, "FIBZERO_TMLFME", "FIBZERO_TMLFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZERO_TMLF, UINT32, "FIBZERO_TMLF", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMLFv, UINT32, "FIBZERO_TMLFv", "FIBZERO_TMLF");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZERO_TMTCM, UINT32, "FIBZERO_TMTCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMTCMI, UINT32, "FIBZERO_TMTCMI", "FIBZERO_TMTCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMTCME, UINT32, "FIBZERO_TMTCME", "FIBZERO_TMTCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZERO_TMTC, UINT32, "FIBZERO_TMTC", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMTCv, UINT32, "FIBZERO_TMTCv", "FIBZERO_TMTC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZERO_TMTFM, UINT32, "FIBZERO_TMTFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMTFMI, UINT32, "FIBZERO_TMTFMI", "FIBZERO_TMTFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMTFME, UINT32, "FIBZERO_TMTFME", "FIBZERO_TMTFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZERO_TMTF, UINT32, "FIBZERO_TMTF", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBZERO_TMTFv, UINT32, "FIBZERO_TMTFv", "FIBZERO_TMTF");        \
                                                                                                                \
    } while (0);

#endif /*__GUARD_H101_FIBZERO_EXT_H101_FIBZERO_H__*/

/*******************************************************/
