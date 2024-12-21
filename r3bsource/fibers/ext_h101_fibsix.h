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

#ifndef __GUARD_H101_FIBSIX_EXT_H101_FIBSIX_H__
#define __GUARD_H101_FIBSIX_EXT_H101_FIBSIX_H__

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

typedef struct EXT_STR_h101_FIBSIX_t
{
    /* RAW */
    uint32_t FIBSIX_TMLCM /* [1,256] */;
    uint32_t FIBSIX_TMLCMI[256 EXT_STRUCT_CTRL(FIBSIX_TMLCM)] /* [1,256] */;
    uint32_t FIBSIX_TMLCME[256 EXT_STRUCT_CTRL(FIBSIX_TMLCM)] /* [1,5120] */;
    uint32_t FIBSIX_TMLC /* [0,5120] */;
    uint32_t FIBSIX_TMLCv[5120 EXT_STRUCT_CTRL(FIBSIX_TMLC)] /* [0,65535] */;
    uint32_t FIBSIX_TMLFM /* [1,256] */;
    uint32_t FIBSIX_TMLFMI[256 EXT_STRUCT_CTRL(FIBSIX_TMLFM)] /* [1,256] */;
    uint32_t FIBSIX_TMLFME[256 EXT_STRUCT_CTRL(FIBSIX_TMLFM)] /* [1,5120] */;
    uint32_t FIBSIX_TMLF /* [0,5120] */;
    uint32_t FIBSIX_TMLFv[5120 EXT_STRUCT_CTRL(FIBSIX_TMLF)] /* [0,65535] */;
    uint32_t FIBSIX_TMTCM /* [1,256] */;
    uint32_t FIBSIX_TMTCMI[256 EXT_STRUCT_CTRL(FIBSIX_TMTCM)] /* [1,256] */;
    uint32_t FIBSIX_TMTCME[256 EXT_STRUCT_CTRL(FIBSIX_TMTCM)] /* [1,5120] */;
    uint32_t FIBSIX_TMTC /* [0,5120] */;
    uint32_t FIBSIX_TMTCv[5120 EXT_STRUCT_CTRL(FIBSIX_TMTC)] /* [0,65535] */;
    uint32_t FIBSIX_TMTFM /* [1,256] */;
    uint32_t FIBSIX_TMTFMI[256 EXT_STRUCT_CTRL(FIBSIX_TMTFM)] /* [1,256] */;
    uint32_t FIBSIX_TMTFME[256 EXT_STRUCT_CTRL(FIBSIX_TMTFM)] /* [1,5120] */;
    uint32_t FIBSIX_TMTF /* [0,5120] */;
    uint32_t FIBSIX_TMTFv[5120 EXT_STRUCT_CTRL(FIBSIX_TMTF)] /* [0,65535] */;
    uint32_t FIBSIX_TSLCM /* [1,2] */;
    uint32_t FIBSIX_TSLCMI[2 EXT_STRUCT_CTRL(FIBSIX_TSLCM)] /* [1,2] */;
    uint32_t FIBSIX_TSLCME[2 EXT_STRUCT_CTRL(FIBSIX_TSLCM)] /* [1,40] */;
    uint32_t FIBSIX_TSLC /* [0,40] */;
    uint32_t FIBSIX_TSLCv[40 EXT_STRUCT_CTRL(FIBSIX_TSLC)] /* [0,65535] */;
    uint32_t FIBSIX_TSLFM /* [1,2] */;
    uint32_t FIBSIX_TSLFMI[2 EXT_STRUCT_CTRL(FIBSIX_TSLFM)] /* [1,2] */;
    uint32_t FIBSIX_TSLFME[2 EXT_STRUCT_CTRL(FIBSIX_TSLFM)] /* [1,40] */;
    uint32_t FIBSIX_TSLF /* [0,40] */;
    uint32_t FIBSIX_TSLFv[40 EXT_STRUCT_CTRL(FIBSIX_TSLF)] /* [0,65535] */;
    uint32_t FIBSIX_TSTCM /* [1,2] */;
    uint32_t FIBSIX_TSTCMI[2 EXT_STRUCT_CTRL(FIBSIX_TSTCM)] /* [1,2] */;
    uint32_t FIBSIX_TSTCME[2 EXT_STRUCT_CTRL(FIBSIX_TSTCM)] /* [1,40] */;
    uint32_t FIBSIX_TSTC /* [0,40] */;
    uint32_t FIBSIX_TSTCv[40 EXT_STRUCT_CTRL(FIBSIX_TSTC)] /* [0,65535] */;
    uint32_t FIBSIX_TSTFM /* [1,2] */;
    uint32_t FIBSIX_TSTFMI[2 EXT_STRUCT_CTRL(FIBSIX_TSTFM)] /* [1,2] */;
    uint32_t FIBSIX_TSTFME[2 EXT_STRUCT_CTRL(FIBSIX_TSTFM)] /* [1,40] */;
    uint32_t FIBSIX_TSTF /* [0,40] */;
    uint32_t FIBSIX_TSTFv[40 EXT_STRUCT_CTRL(FIBSIX_TSTF)] /* [0,65535] */;

} EXT_STR_h101_FIBSIX;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBSIX_onion_t
{
    /* RAW */
    uint32_t FIBSIX_TMLCM;
    uint32_t FIBSIX_TMLCMI[256 /* FIBSIX_TMLCM */];
    uint32_t FIBSIX_TMLCME[256 /* FIBSIX_TMLCM */];
    uint32_t FIBSIX_TMLC;
    uint32_t FIBSIX_TMLCv[5120 /* FIBSIX_TMLC */];
    uint32_t FIBSIX_TMLFM;
    uint32_t FIBSIX_TMLFMI[256 /* FIBSIX_TMLFM */];
    uint32_t FIBSIX_TMLFME[256 /* FIBSIX_TMLFM */];
    uint32_t FIBSIX_TMLF;
    uint32_t FIBSIX_TMLFv[5120 /* FIBSIX_TMLF */];
    uint32_t FIBSIX_TMTCM;
    uint32_t FIBSIX_TMTCMI[256 /* FIBSIX_TMTCM */];
    uint32_t FIBSIX_TMTCME[256 /* FIBSIX_TMTCM */];
    uint32_t FIBSIX_TMTC;
    uint32_t FIBSIX_TMTCv[5120 /* FIBSIX_TMTC */];
    uint32_t FIBSIX_TMTFM;
    uint32_t FIBSIX_TMTFMI[256 /* FIBSIX_TMTFM */];
    uint32_t FIBSIX_TMTFME[256 /* FIBSIX_TMTFM */];
    uint32_t FIBSIX_TMTF;
    uint32_t FIBSIX_TMTFv[5120 /* FIBSIX_TMTF */];
    uint32_t FIBSIX_TSLCM;
    uint32_t FIBSIX_TSLCMI[2 /* FIBSIX_TSLCM */];
    uint32_t FIBSIX_TSLCME[2 /* FIBSIX_TSLCM */];
    uint32_t FIBSIX_TSLC;
    uint32_t FIBSIX_TSLCv[40 /* FIBSIX_TSLC */];
    uint32_t FIBSIX_TSLFM;
    uint32_t FIBSIX_TSLFMI[2 /* FIBSIX_TSLFM */];
    uint32_t FIBSIX_TSLFME[2 /* FIBSIX_TSLFM */];
    uint32_t FIBSIX_TSLF;
    uint32_t FIBSIX_TSLFv[40 /* FIBSIX_TSLF */];
    uint32_t FIBSIX_TSTCM;
    uint32_t FIBSIX_TSTCMI[2 /* FIBSIX_TSTCM */];
    uint32_t FIBSIX_TSTCME[2 /* FIBSIX_TSTCM */];
    uint32_t FIBSIX_TSTC;
    uint32_t FIBSIX_TSTCv[40 /* FIBSIX_TSTC */];
    uint32_t FIBSIX_TSTFM;
    uint32_t FIBSIX_TSTFMI[2 /* FIBSIX_TSTFM */];
    uint32_t FIBSIX_TSTFME[2 /* FIBSIX_TSTFM */];
    uint32_t FIBSIX_TSTF;
    uint32_t FIBSIX_TSTFv[40 /* FIBSIX_TSTF */];

} EXT_STR_h101_FIBSIX_onion;

/*******************************************************/

#define EXT_STR_h101_FIBSIX_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                    \
    do                                                                                                        \
    {                                                                                                         \
        ok = 1;                                                                                               \
        /* RAW */                                                                                             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TMLCM, UINT32, "FIBSIX_TMLCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMLCMI, UINT32, "FIBSIX_TMLCMI", "FIBSIX_TMLCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMLCME, UINT32, "FIBSIX_TMLCME", "FIBSIX_TMLCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TMLC, UINT32, "FIBSIX_TMLC", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMLCv, UINT32, "FIBSIX_TMLCv", "FIBSIX_TMLC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TMLFM, UINT32, "FIBSIX_TMLFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMLFMI, UINT32, "FIBSIX_TMLFMI", "FIBSIX_TMLFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMLFME, UINT32, "FIBSIX_TMLFME", "FIBSIX_TMLFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TMLF, UINT32, "FIBSIX_TMLF", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMLFv, UINT32, "FIBSIX_TMLFv", "FIBSIX_TMLF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TMTCM, UINT32, "FIBSIX_TMTCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMTCMI, UINT32, "FIBSIX_TMTCMI", "FIBSIX_TMTCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMTCME, UINT32, "FIBSIX_TMTCME", "FIBSIX_TMTCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TMTC, UINT32, "FIBSIX_TMTC", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMTCv, UINT32, "FIBSIX_TMTCv", "FIBSIX_TMTC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TMTFM, UINT32, "FIBSIX_TMTFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMTFMI, UINT32, "FIBSIX_TMTFMI", "FIBSIX_TMTFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMTFME, UINT32, "FIBSIX_TMTFME", "FIBSIX_TMTFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TMTF, UINT32, "FIBSIX_TMTF", 5120);  \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TMTFv, UINT32, "FIBSIX_TMTFv", "FIBSIX_TMTF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TSLCM, UINT32, "FIBSIX_TSLCM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSLCMI, UINT32, "FIBSIX_TSLCMI", "FIBSIX_TSLCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSLCME, UINT32, "FIBSIX_TSLCME", "FIBSIX_TSLCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TSLC, UINT32, "FIBSIX_TSLC", 40);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSLCv, UINT32, "FIBSIX_TSLCv", "FIBSIX_TSLC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TSLFM, UINT32, "FIBSIX_TSLFM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSLFMI, UINT32, "FIBSIX_TSLFMI", "FIBSIX_TSLFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSLFME, UINT32, "FIBSIX_TSLFME", "FIBSIX_TSLFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TSLF, UINT32, "FIBSIX_TSLF", 40);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSLFv, UINT32, "FIBSIX_TSLFv", "FIBSIX_TSLF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TSTCM, UINT32, "FIBSIX_TSTCM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSTCMI, UINT32, "FIBSIX_TSTCMI", "FIBSIX_TSTCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSTCME, UINT32, "FIBSIX_TSTCME", "FIBSIX_TSTCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TSTC, UINT32, "FIBSIX_TSTC", 40);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSTCv, UINT32, "FIBSIX_TSTCv", "FIBSIX_TSTC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TSTFM, UINT32, "FIBSIX_TSTFM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSTFMI, UINT32, "FIBSIX_TSTFMI", "FIBSIX_TSTFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSTFME, UINT32, "FIBSIX_TSTFME", "FIBSIX_TSTFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSIX_TSTF, UINT32, "FIBSIX_TSTF", 40);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBSIX_TSTFv, UINT32, "FIBSIX_TSTFv", "FIBSIX_TSTF");         \
                                                                                                              \
    } while (0);

#endif /*__GUARD_H101_FIBSIX_EXT_H101_FIBSIX_H__*/

/*******************************************************/
