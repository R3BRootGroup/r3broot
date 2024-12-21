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

#ifndef __GUARD_H101_FIBFOUR_EXT_H101_FIBFOUR_H__
#define __GUARD_H101_FIBFOUR_EXT_H101_FIBFOUR_H__

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

typedef struct EXT_STR_h101_FIBFOUR_t
{
    /* RAW */
    uint32_t FIBFOUR_TMLCM /* [1,512] */;
    uint32_t FIBFOUR_TMLCMI[512 EXT_STRUCT_CTRL(FIBFOUR_TMLCM)] /* [1,512] */;
    uint32_t FIBFOUR_TMLCME[512 EXT_STRUCT_CTRL(FIBFOUR_TMLCM)] /* [1,10240] */;
    uint32_t FIBFOUR_TMLC /* [0,10240] */;
    uint32_t FIBFOUR_TMLCv[10240 EXT_STRUCT_CTRL(FIBFOUR_TMLC)] /* [0,65535] */;
    uint32_t FIBFOUR_TMLFM /* [1,512] */;
    uint32_t FIBFOUR_TMLFMI[512 EXT_STRUCT_CTRL(FIBFOUR_TMLFM)] /* [1,512] */;
    uint32_t FIBFOUR_TMLFME[512 EXT_STRUCT_CTRL(FIBFOUR_TMLFM)] /* [1,10240] */;
    uint32_t FIBFOUR_TMLF /* [0,10240] */;
    uint32_t FIBFOUR_TMLFv[10240 EXT_STRUCT_CTRL(FIBFOUR_TMLF)] /* [0,65535] */;
    uint32_t FIBFOUR_TMTCM /* [1,512] */;
    uint32_t FIBFOUR_TMTCMI[512 EXT_STRUCT_CTRL(FIBFOUR_TMTCM)] /* [1,512] */;
    uint32_t FIBFOUR_TMTCME[512 EXT_STRUCT_CTRL(FIBFOUR_TMTCM)] /* [1,10240] */;
    uint32_t FIBFOUR_TMTC /* [0,10240] */;
    uint32_t FIBFOUR_TMTCv[10240 EXT_STRUCT_CTRL(FIBFOUR_TMTC)] /* [0,65535] */;
    uint32_t FIBFOUR_TMTFM /* [1,512] */;
    uint32_t FIBFOUR_TMTFMI[512 EXT_STRUCT_CTRL(FIBFOUR_TMTFM)] /* [1,512] */;
    uint32_t FIBFOUR_TMTFME[512 EXT_STRUCT_CTRL(FIBFOUR_TMTFM)] /* [1,10240] */;
    uint32_t FIBFOUR_TMTF /* [0,10240] */;
    uint32_t FIBFOUR_TMTFv[10240 EXT_STRUCT_CTRL(FIBFOUR_TMTF)] /* [0,65535] */;
    uint32_t FIBFOUR_TSLCM /* [1,8] */;
    uint32_t FIBFOUR_TSLCMI[8 EXT_STRUCT_CTRL(FIBFOUR_TSLCM)] /* [1,8] */;
    uint32_t FIBFOUR_TSLCME[8 EXT_STRUCT_CTRL(FIBFOUR_TSLCM)] /* [1,160] */;
    uint32_t FIBFOUR_TSLC /* [0,160] */;
    uint32_t FIBFOUR_TSLCv[160 EXT_STRUCT_CTRL(FIBFOUR_TSLC)] /* [0,65535] */;
    uint32_t FIBFOUR_TSLFM /* [1,8] */;
    uint32_t FIBFOUR_TSLFMI[8 EXT_STRUCT_CTRL(FIBFOUR_TSLFM)] /* [1,8] */;
    uint32_t FIBFOUR_TSLFME[8 EXT_STRUCT_CTRL(FIBFOUR_TSLFM)] /* [1,160] */;
    uint32_t FIBFOUR_TSLF /* [0,160] */;
    uint32_t FIBFOUR_TSLFv[160 EXT_STRUCT_CTRL(FIBFOUR_TSLF)] /* [0,65535] */;
    uint32_t FIBFOUR_TSTCM /* [1,8] */;
    uint32_t FIBFOUR_TSTCMI[8 EXT_STRUCT_CTRL(FIBFOUR_TSTCM)] /* [1,8] */;
    uint32_t FIBFOUR_TSTCME[8 EXT_STRUCT_CTRL(FIBFOUR_TSTCM)] /* [1,160] */;
    uint32_t FIBFOUR_TSTC /* [0,160] */;
    uint32_t FIBFOUR_TSTCv[160 EXT_STRUCT_CTRL(FIBFOUR_TSTC)] /* [0,65535] */;
    uint32_t FIBFOUR_TSTFM /* [1,8] */;
    uint32_t FIBFOUR_TSTFMI[8 EXT_STRUCT_CTRL(FIBFOUR_TSTFM)] /* [1,8] */;
    uint32_t FIBFOUR_TSTFME[8 EXT_STRUCT_CTRL(FIBFOUR_TSTFM)] /* [1,160] */;
    uint32_t FIBFOUR_TSTF /* [0,160] */;
    uint32_t FIBFOUR_TSTFv[160 EXT_STRUCT_CTRL(FIBFOUR_TSTF)] /* [0,65535] */;

} EXT_STR_h101_FIBFOUR;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBFOUR_onion_t
{
    /* RAW */
    uint32_t FIBFOUR_TMLCM;
    uint32_t FIBFOUR_TMLCMI[512 /* FIBFOUR_TMLCM */];
    uint32_t FIBFOUR_TMLCME[512 /* FIBFOUR_TMLCM */];
    uint32_t FIBFOUR_TMLC;
    uint32_t FIBFOUR_TMLCv[10240 /* FIBFOUR_TMLC */];
    uint32_t FIBFOUR_TMLFM;
    uint32_t FIBFOUR_TMLFMI[512 /* FIBFOUR_TMLFM */];
    uint32_t FIBFOUR_TMLFME[512 /* FIBFOUR_TMLFM */];
    uint32_t FIBFOUR_TMLF;
    uint32_t FIBFOUR_TMLFv[10240 /* FIBFOUR_TMLF */];
    uint32_t FIBFOUR_TMTCM;
    uint32_t FIBFOUR_TMTCMI[512 /* FIBFOUR_TMTCM */];
    uint32_t FIBFOUR_TMTCME[512 /* FIBFOUR_TMTCM */];
    uint32_t FIBFOUR_TMTC;
    uint32_t FIBFOUR_TMTCv[10240 /* FIBFOUR_TMTC */];
    uint32_t FIBFOUR_TMTFM;
    uint32_t FIBFOUR_TMTFMI[512 /* FIBFOUR_TMTFM */];
    uint32_t FIBFOUR_TMTFME[512 /* FIBFOUR_TMTFM */];
    uint32_t FIBFOUR_TMTF;
    uint32_t FIBFOUR_TMTFv[10240 /* FIBFOUR_TMTF */];
    uint32_t FIBFOUR_TSLCM;
    uint32_t FIBFOUR_TSLCMI[8 /* FIBFOUR_TSLCM */];
    uint32_t FIBFOUR_TSLCME[8 /* FIBFOUR_TSLCM */];
    uint32_t FIBFOUR_TSLC;
    uint32_t FIBFOUR_TSLCv[160 /* FIBFOUR_TSLC */];
    uint32_t FIBFOUR_TSLFM;
    uint32_t FIBFOUR_TSLFMI[8 /* FIBFOUR_TSLFM */];
    uint32_t FIBFOUR_TSLFME[8 /* FIBFOUR_TSLFM */];
    uint32_t FIBFOUR_TSLF;
    uint32_t FIBFOUR_TSLFv[160 /* FIBFOUR_TSLF */];
    uint32_t FIBFOUR_TSTCM;
    uint32_t FIBFOUR_TSTCMI[8 /* FIBFOUR_TSTCM */];
    uint32_t FIBFOUR_TSTCME[8 /* FIBFOUR_TSTCM */];
    uint32_t FIBFOUR_TSTC;
    uint32_t FIBFOUR_TSTCv[160 /* FIBFOUR_TSTC */];
    uint32_t FIBFOUR_TSTFM;
    uint32_t FIBFOUR_TSTFMI[8 /* FIBFOUR_TSTFM */];
    uint32_t FIBFOUR_TSTFME[8 /* FIBFOUR_TSTFM */];
    uint32_t FIBFOUR_TSTF;
    uint32_t FIBFOUR_TSTFv[160 /* FIBFOUR_TSTF */];

} EXT_STR_h101_FIBFOUR_onion;

/*******************************************************/

#define EXT_STR_h101_FIBFOUR_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                     \
    do                                                                                                          \
    {                                                                                                           \
        ok = 1;                                                                                                 \
        /* RAW */                                                                                               \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TMLCM, UINT32, "FIBFOUR_TMLCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMLCMI, UINT32, "FIBFOUR_TMLCMI", "FIBFOUR_TMLCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMLCME, UINT32, "FIBFOUR_TMLCME", "FIBFOUR_TMLCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TMLC, UINT32, "FIBFOUR_TMLC", 10240); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMLCv, UINT32, "FIBFOUR_TMLCv", "FIBFOUR_TMLC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TMLFM, UINT32, "FIBFOUR_TMLFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMLFMI, UINT32, "FIBFOUR_TMLFMI", "FIBFOUR_TMLFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMLFME, UINT32, "FIBFOUR_TMLFME", "FIBFOUR_TMLFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TMLF, UINT32, "FIBFOUR_TMLF", 10240); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMLFv, UINT32, "FIBFOUR_TMLFv", "FIBFOUR_TMLF");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TMTCM, UINT32, "FIBFOUR_TMTCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMTCMI, UINT32, "FIBFOUR_TMTCMI", "FIBFOUR_TMTCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMTCME, UINT32, "FIBFOUR_TMTCME", "FIBFOUR_TMTCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TMTC, UINT32, "FIBFOUR_TMTC", 10240); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMTCv, UINT32, "FIBFOUR_TMTCv", "FIBFOUR_TMTC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TMTFM, UINT32, "FIBFOUR_TMTFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMTFMI, UINT32, "FIBFOUR_TMTFMI", "FIBFOUR_TMTFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMTFME, UINT32, "FIBFOUR_TMTFME", "FIBFOUR_TMTFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TMTF, UINT32, "FIBFOUR_TMTF", 10240); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TMTFv, UINT32, "FIBFOUR_TMTFv", "FIBFOUR_TMTF");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TSLCM, UINT32, "FIBFOUR_TSLCM", 8);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSLCMI, UINT32, "FIBFOUR_TSLCMI", "FIBFOUR_TSLCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSLCME, UINT32, "FIBFOUR_TSLCME", "FIBFOUR_TSLCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TSLC, UINT32, "FIBFOUR_TSLC", 160);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSLCv, UINT32, "FIBFOUR_TSLCv", "FIBFOUR_TSLC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TSLFM, UINT32, "FIBFOUR_TSLFM", 8);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSLFMI, UINT32, "FIBFOUR_TSLFMI", "FIBFOUR_TSLFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSLFME, UINT32, "FIBFOUR_TSLFME", "FIBFOUR_TSLFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TSLF, UINT32, "FIBFOUR_TSLF", 160);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSLFv, UINT32, "FIBFOUR_TSLFv", "FIBFOUR_TSLF");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TSTCM, UINT32, "FIBFOUR_TSTCM", 8);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSTCMI, UINT32, "FIBFOUR_TSTCMI", "FIBFOUR_TSTCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSTCME, UINT32, "FIBFOUR_TSTCME", "FIBFOUR_TSTCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TSTC, UINT32, "FIBFOUR_TSTC", 160);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSTCv, UINT32, "FIBFOUR_TSTCv", "FIBFOUR_TSTC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TSTFM, UINT32, "FIBFOUR_TSTFM", 8);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSTFMI, UINT32, "FIBFOUR_TSTFMI", "FIBFOUR_TSTFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSTFME, UINT32, "FIBFOUR_TSTFME", "FIBFOUR_TSTFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFOUR_TSTF, UINT32, "FIBFOUR_TSTF", 160);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFOUR_TSTFv, UINT32, "FIBFOUR_TSTFv", "FIBFOUR_TSTF");        \
                                                                                                                \
    } while (0);

#endif /*__GUARD_H101_FIBFOUR_EXT_H101_FIBFOUR_H__*/

/*******************************************************/
