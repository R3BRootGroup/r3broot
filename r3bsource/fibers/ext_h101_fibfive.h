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

#ifndef __GUARD_H101_FIBFIVE_EXT_H101_FIBFIVE_H__
#define __GUARD_H101_FIBFIVE_EXT_H101_FIBFIVE_H__

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

typedef struct EXT_STR_h101_FIBFIVE_t
{
    /* RAW */
    uint32_t FIBFIVE_TMLCM /* [1,256] */;
    uint32_t FIBFIVE_TMLCMI[256 EXT_STRUCT_CTRL(FIBFIVE_TMLCM)] /* [1,256] */;
    uint32_t FIBFIVE_TMLCME[256 EXT_STRUCT_CTRL(FIBFIVE_TMLCM)] /* [1,5120] */;
    uint32_t FIBFIVE_TMLC /* [0,5120] */;
    uint32_t FIBFIVE_TMLCv[5120 EXT_STRUCT_CTRL(FIBFIVE_TMLC)] /* [0,65535] */;
    uint32_t FIBFIVE_TMLFM /* [1,256] */;
    uint32_t FIBFIVE_TMLFMI[256 EXT_STRUCT_CTRL(FIBFIVE_TMLFM)] /* [1,256] */;
    uint32_t FIBFIVE_TMLFME[256 EXT_STRUCT_CTRL(FIBFIVE_TMLFM)] /* [1,5120] */;
    uint32_t FIBFIVE_TMLF /* [0,5120] */;
    uint32_t FIBFIVE_TMLFv[5120 EXT_STRUCT_CTRL(FIBFIVE_TMLF)] /* [0,65535] */;
    uint32_t FIBFIVE_TMTCM /* [1,256] */;
    uint32_t FIBFIVE_TMTCMI[256 EXT_STRUCT_CTRL(FIBFIVE_TMTCM)] /* [1,256] */;
    uint32_t FIBFIVE_TMTCME[256 EXT_STRUCT_CTRL(FIBFIVE_TMTCM)] /* [1,5120] */;
    uint32_t FIBFIVE_TMTC /* [0,5120] */;
    uint32_t FIBFIVE_TMTCv[5120 EXT_STRUCT_CTRL(FIBFIVE_TMTC)] /* [0,65535] */;
    uint32_t FIBFIVE_TMTFM /* [1,256] */;
    uint32_t FIBFIVE_TMTFMI[256 EXT_STRUCT_CTRL(FIBFIVE_TMTFM)] /* [1,256] */;
    uint32_t FIBFIVE_TMTFME[256 EXT_STRUCT_CTRL(FIBFIVE_TMTFM)] /* [1,5120] */;
    uint32_t FIBFIVE_TMTF /* [0,5120] */;
    uint32_t FIBFIVE_TMTFv[5120 EXT_STRUCT_CTRL(FIBFIVE_TMTF)] /* [0,65535] */;
    uint32_t FIBFIVE_TSLCM /* [1,4] */;
    uint32_t FIBFIVE_TSLCMI[4 EXT_STRUCT_CTRL(FIBFIVE_TSLCM)] /* [1,4] */;
    uint32_t FIBFIVE_TSLCME[4 EXT_STRUCT_CTRL(FIBFIVE_TSLCM)] /* [1,80] */;
    uint32_t FIBFIVE_TSLC /* [0,80] */;
    uint32_t FIBFIVE_TSLCv[80 EXT_STRUCT_CTRL(FIBFIVE_TSLC)] /* [0,65535] */;
    uint32_t FIBFIVE_TSLFM /* [1,4] */;
    uint32_t FIBFIVE_TSLFMI[4 EXT_STRUCT_CTRL(FIBFIVE_TSLFM)] /* [1,4] */;
    uint32_t FIBFIVE_TSLFME[4 EXT_STRUCT_CTRL(FIBFIVE_TSLFM)] /* [1,80] */;
    uint32_t FIBFIVE_TSLF /* [0,80] */;
    uint32_t FIBFIVE_TSLFv[80 EXT_STRUCT_CTRL(FIBFIVE_TSLF)] /* [0,65535] */;
    uint32_t FIBFIVE_TSTCM /* [1,4] */;
    uint32_t FIBFIVE_TSTCMI[4 EXT_STRUCT_CTRL(FIBFIVE_TSTCM)] /* [1,4] */;
    uint32_t FIBFIVE_TSTCME[4 EXT_STRUCT_CTRL(FIBFIVE_TSTCM)] /* [1,80] */;
    uint32_t FIBFIVE_TSTC /* [0,80] */;
    uint32_t FIBFIVE_TSTCv[80 EXT_STRUCT_CTRL(FIBFIVE_TSTC)] /* [0,65535] */;
    uint32_t FIBFIVE_TSTFM /* [1,4] */;
    uint32_t FIBFIVE_TSTFMI[4 EXT_STRUCT_CTRL(FIBFIVE_TSTFM)] /* [1,4] */;
    uint32_t FIBFIVE_TSTFME[4 EXT_STRUCT_CTRL(FIBFIVE_TSTFM)] /* [1,80] */;
    uint32_t FIBFIVE_TSTF /* [0,80] */;
    uint32_t FIBFIVE_TSTFv[80 EXT_STRUCT_CTRL(FIBFIVE_TSTF)] /* [0,65535] */;

} EXT_STR_h101_FIBFIVE;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBFIVE_onion_t
{
    /* RAW */
    uint32_t FIBFIVE_TMLCM;
    uint32_t FIBFIVE_TMLCMI[256 /* FIBFIVE_TMLCM */];
    uint32_t FIBFIVE_TMLCME[256 /* FIBFIVE_TMLCM */];
    uint32_t FIBFIVE_TMLC;
    uint32_t FIBFIVE_TMLCv[5120 /* FIBFIVE_TMLC */];
    uint32_t FIBFIVE_TMLFM;
    uint32_t FIBFIVE_TMLFMI[256 /* FIBFIVE_TMLFM */];
    uint32_t FIBFIVE_TMLFME[256 /* FIBFIVE_TMLFM */];
    uint32_t FIBFIVE_TMLF;
    uint32_t FIBFIVE_TMLFv[5120 /* FIBFIVE_TMLF */];
    uint32_t FIBFIVE_TMTCM;
    uint32_t FIBFIVE_TMTCMI[256 /* FIBFIVE_TMTCM */];
    uint32_t FIBFIVE_TMTCME[256 /* FIBFIVE_TMTCM */];
    uint32_t FIBFIVE_TMTC;
    uint32_t FIBFIVE_TMTCv[5120 /* FIBFIVE_TMTC */];
    uint32_t FIBFIVE_TMTFM;
    uint32_t FIBFIVE_TMTFMI[256 /* FIBFIVE_TMTFM */];
    uint32_t FIBFIVE_TMTFME[256 /* FIBFIVE_TMTFM */];
    uint32_t FIBFIVE_TMTF;
    uint32_t FIBFIVE_TMTFv[5120 /* FIBFIVE_TMTF */];
    uint32_t FIBFIVE_TSLCM;
    uint32_t FIBFIVE_TSLCMI[4 /* FIBFIVE_TSLCM */];
    uint32_t FIBFIVE_TSLCME[4 /* FIBFIVE_TSLCM */];
    uint32_t FIBFIVE_TSLC;
    uint32_t FIBFIVE_TSLCv[80 /* FIBFIVE_TSLC */];
    uint32_t FIBFIVE_TSLFM;
    uint32_t FIBFIVE_TSLFMI[4 /* FIBFIVE_TSLFM */];
    uint32_t FIBFIVE_TSLFME[4 /* FIBFIVE_TSLFM */];
    uint32_t FIBFIVE_TSLF;
    uint32_t FIBFIVE_TSLFv[80 /* FIBFIVE_TSLF */];
    uint32_t FIBFIVE_TSTCM;
    uint32_t FIBFIVE_TSTCMI[4 /* FIBFIVE_TSTCM */];
    uint32_t FIBFIVE_TSTCME[4 /* FIBFIVE_TSTCM */];
    uint32_t FIBFIVE_TSTC;
    uint32_t FIBFIVE_TSTCv[80 /* FIBFIVE_TSTC */];
    uint32_t FIBFIVE_TSTFM;
    uint32_t FIBFIVE_TSTFMI[4 /* FIBFIVE_TSTFM */];
    uint32_t FIBFIVE_TSTFME[4 /* FIBFIVE_TSTFM */];
    uint32_t FIBFIVE_TSTF;
    uint32_t FIBFIVE_TSTFv[80 /* FIBFIVE_TSTF */];

} EXT_STR_h101_FIBFIVE_onion;

/*******************************************************/

#define EXT_STR_h101_FIBFIVE_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                     \
    do                                                                                                          \
    {                                                                                                           \
        ok = 1;                                                                                                 \
        /* RAW */                                                                                               \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TMLCM, UINT32, "FIBFIVE_TMLCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMLCMI, UINT32, "FIBFIVE_TMLCMI", "FIBFIVE_TMLCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMLCME, UINT32, "FIBFIVE_TMLCME", "FIBFIVE_TMLCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TMLC, UINT32, "FIBFIVE_TMLC", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMLCv, UINT32, "FIBFIVE_TMLCv", "FIBFIVE_TMLC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TMLFM, UINT32, "FIBFIVE_TMLFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMLFMI, UINT32, "FIBFIVE_TMLFMI", "FIBFIVE_TMLFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMLFME, UINT32, "FIBFIVE_TMLFME", "FIBFIVE_TMLFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TMLF, UINT32, "FIBFIVE_TMLF", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMLFv, UINT32, "FIBFIVE_TMLFv", "FIBFIVE_TMLF");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TMTCM, UINT32, "FIBFIVE_TMTCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMTCMI, UINT32, "FIBFIVE_TMTCMI", "FIBFIVE_TMTCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMTCME, UINT32, "FIBFIVE_TMTCME", "FIBFIVE_TMTCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TMTC, UINT32, "FIBFIVE_TMTC", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMTCv, UINT32, "FIBFIVE_TMTCv", "FIBFIVE_TMTC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TMTFM, UINT32, "FIBFIVE_TMTFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMTFMI, UINT32, "FIBFIVE_TMTFMI", "FIBFIVE_TMTFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMTFME, UINT32, "FIBFIVE_TMTFME", "FIBFIVE_TMTFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TMTF, UINT32, "FIBFIVE_TMTF", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TMTFv, UINT32, "FIBFIVE_TMTFv", "FIBFIVE_TMTF");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TSLCM, UINT32, "FIBFIVE_TSLCM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSLCMI, UINT32, "FIBFIVE_TSLCMI", "FIBFIVE_TSLCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSLCME, UINT32, "FIBFIVE_TSLCME", "FIBFIVE_TSLCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TSLC, UINT32, "FIBFIVE_TSLC", 80);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSLCv, UINT32, "FIBFIVE_TSLCv", "FIBFIVE_TSLC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TSLFM, UINT32, "FIBFIVE_TSLFM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSLFMI, UINT32, "FIBFIVE_TSLFMI", "FIBFIVE_TSLFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSLFME, UINT32, "FIBFIVE_TSLFME", "FIBFIVE_TSLFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TSLF, UINT32, "FIBFIVE_TSLF", 80);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSLFv, UINT32, "FIBFIVE_TSLFv", "FIBFIVE_TSLF");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TSTCM, UINT32, "FIBFIVE_TSTCM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSTCMI, UINT32, "FIBFIVE_TSTCMI", "FIBFIVE_TSTCM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSTCME, UINT32, "FIBFIVE_TSTCME", "FIBFIVE_TSTCM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TSTC, UINT32, "FIBFIVE_TSTC", 80);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSTCv, UINT32, "FIBFIVE_TSTCv", "FIBFIVE_TSTC");        \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TSTFM, UINT32, "FIBFIVE_TSTFM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSTFMI, UINT32, "FIBFIVE_TSTFMI", "FIBFIVE_TSTFM");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSTFME, UINT32, "FIBFIVE_TSTFME", "FIBFIVE_TSTFM");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBFIVE_TSTF, UINT32, "FIBFIVE_TSTF", 80);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                  \
            ok, si, offset, struct_t, printerr, FIBFIVE_TSTFv, UINT32, "FIBFIVE_TSTFv", "FIBFIVE_TSTF");        \
                                                                                                                \
    } while (0);

#endif /*__GUARD_H101_FIBFIVE_EXT_H101_FIBFIVE_H__*/

/*******************************************************/
