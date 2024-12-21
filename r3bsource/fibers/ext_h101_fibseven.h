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

#ifndef __GUARD_H101_FIBSEVEN_EXT_H101_FIBSEVEN_H__
#define __GUARD_H101_FIBSEVEN_EXT_H101_FIBSEVEN_H__

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

typedef struct EXT_STR_h101_FIBSEVEN_t
{
    /* RAW */
    uint32_t FIBSEVEN_TMLCM /* [1,512] */;
    uint32_t FIBSEVEN_TMLCMI[512 EXT_STRUCT_CTRL(FIBSEVEN_TMLCM)] /* [1,512] */;
    uint32_t FIBSEVEN_TMLCME[512 EXT_STRUCT_CTRL(FIBSEVEN_TMLCM)] /* [1,10240] */;
    uint32_t FIBSEVEN_TMLC /* [0,10240] */;
    uint32_t FIBSEVEN_TMLCv[10240 EXT_STRUCT_CTRL(FIBSEVEN_TMLC)] /* [0,65535] */;
    uint32_t FIBSEVEN_TMLFM /* [1,512] */;
    uint32_t FIBSEVEN_TMLFMI[512 EXT_STRUCT_CTRL(FIBSEVEN_TMLFM)] /* [1,512] */;
    uint32_t FIBSEVEN_TMLFME[512 EXT_STRUCT_CTRL(FIBSEVEN_TMLFM)] /* [1,10240] */;
    uint32_t FIBSEVEN_TMLF /* [0,10240] */;
    uint32_t FIBSEVEN_TMLFv[10240 EXT_STRUCT_CTRL(FIBSEVEN_TMLF)] /* [0,65535] */;
    uint32_t FIBSEVEN_TMTCM /* [1,512] */;
    uint32_t FIBSEVEN_TMTCMI[512 EXT_STRUCT_CTRL(FIBSEVEN_TMTCM)] /* [1,512] */;
    uint32_t FIBSEVEN_TMTCME[512 EXT_STRUCT_CTRL(FIBSEVEN_TMTCM)] /* [1,10240] */;
    uint32_t FIBSEVEN_TMTC /* [0,10240] */;
    uint32_t FIBSEVEN_TMTCv[10240 EXT_STRUCT_CTRL(FIBSEVEN_TMTC)] /* [0,65535] */;
    uint32_t FIBSEVEN_TMTFM /* [1,512] */;
    uint32_t FIBSEVEN_TMTFMI[512 EXT_STRUCT_CTRL(FIBSEVEN_TMTFM)] /* [1,512] */;
    uint32_t FIBSEVEN_TMTFME[512 EXT_STRUCT_CTRL(FIBSEVEN_TMTFM)] /* [1,10240] */;
    uint32_t FIBSEVEN_TMTF /* [0,10240] */;
    uint32_t FIBSEVEN_TMTFv[10240 EXT_STRUCT_CTRL(FIBSEVEN_TMTF)] /* [0,65535] */;
    uint32_t FIBSEVEN_TSLCM /* [1,2] */;
    uint32_t FIBSEVEN_TSLCMI[2 EXT_STRUCT_CTRL(FIBSEVEN_TSLCM)] /* [1,2] */;
    uint32_t FIBSEVEN_TSLCME[2 EXT_STRUCT_CTRL(FIBSEVEN_TSLCM)] /* [1,40] */;
    uint32_t FIBSEVEN_TSLC /* [0,40] */;
    uint32_t FIBSEVEN_TSLCv[40 EXT_STRUCT_CTRL(FIBSEVEN_TSLC)] /* [0,65535] */;
    uint32_t FIBSEVEN_TSLFM /* [1,2] */;
    uint32_t FIBSEVEN_TSLFMI[2 EXT_STRUCT_CTRL(FIBSEVEN_TSLFM)] /* [1,2] */;
    uint32_t FIBSEVEN_TSLFME[2 EXT_STRUCT_CTRL(FIBSEVEN_TSLFM)] /* [1,40] */;
    uint32_t FIBSEVEN_TSLF /* [0,40] */;
    uint32_t FIBSEVEN_TSLFv[40 EXT_STRUCT_CTRL(FIBSEVEN_TSLF)] /* [0,65535] */;
    uint32_t FIBSEVEN_TSTCM /* [1,2] */;
    uint32_t FIBSEVEN_TSTCMI[2 EXT_STRUCT_CTRL(FIBSEVEN_TSTCM)] /* [1,2] */;
    uint32_t FIBSEVEN_TSTCME[2 EXT_STRUCT_CTRL(FIBSEVEN_TSTCM)] /* [1,40] */;
    uint32_t FIBSEVEN_TSTC /* [0,40] */;
    uint32_t FIBSEVEN_TSTCv[40 EXT_STRUCT_CTRL(FIBSEVEN_TSTC)] /* [0,65535] */;
    uint32_t FIBSEVEN_TSTFM /* [1,2] */;
    uint32_t FIBSEVEN_TSTFMI[2 EXT_STRUCT_CTRL(FIBSEVEN_TSTFM)] /* [1,2] */;
    uint32_t FIBSEVEN_TSTFME[2 EXT_STRUCT_CTRL(FIBSEVEN_TSTFM)] /* [1,40] */;
    uint32_t FIBSEVEN_TSTF /* [0,40] */;
    uint32_t FIBSEVEN_TSTFv[40 EXT_STRUCT_CTRL(FIBSEVEN_TSTF)] /* [0,65535] */;

} EXT_STR_h101_FIBSEVEN;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBSEVEN_onion_t
{
    /* RAW */
    uint32_t FIBSEVEN_TMLCM;
    uint32_t FIBSEVEN_TMLCMI[512 /* FIBSEVEN_TMLCM */];
    uint32_t FIBSEVEN_TMLCME[512 /* FIBSEVEN_TMLCM */];
    uint32_t FIBSEVEN_TMLC;
    uint32_t FIBSEVEN_TMLCv[10240 /* FIBSEVEN_TMLC */];
    uint32_t FIBSEVEN_TMLFM;
    uint32_t FIBSEVEN_TMLFMI[512 /* FIBSEVEN_TMLFM */];
    uint32_t FIBSEVEN_TMLFME[512 /* FIBSEVEN_TMLFM */];
    uint32_t FIBSEVEN_TMLF;
    uint32_t FIBSEVEN_TMLFv[10240 /* FIBSEVEN_TMLF */];
    uint32_t FIBSEVEN_TMTCM;
    uint32_t FIBSEVEN_TMTCMI[512 /* FIBSEVEN_TMTCM */];
    uint32_t FIBSEVEN_TMTCME[512 /* FIBSEVEN_TMTCM */];
    uint32_t FIBSEVEN_TMTC;
    uint32_t FIBSEVEN_TMTCv[10240 /* FIBSEVEN_TMTC */];
    uint32_t FIBSEVEN_TMTFM;
    uint32_t FIBSEVEN_TMTFMI[512 /* FIBSEVEN_TMTFM */];
    uint32_t FIBSEVEN_TMTFME[512 /* FIBSEVEN_TMTFM */];
    uint32_t FIBSEVEN_TMTF;
    uint32_t FIBSEVEN_TMTFv[10240 /* FIBSEVEN_TMTF */];
    uint32_t FIBSEVEN_TSLCM;
    uint32_t FIBSEVEN_TSLCMI[2 /* FIBSEVEN_TSLCM */];
    uint32_t FIBSEVEN_TSLCME[2 /* FIBSEVEN_TSLCM */];
    uint32_t FIBSEVEN_TSLC;
    uint32_t FIBSEVEN_TSLCv[40 /* FIBSEVEN_TSLC */];
    uint32_t FIBSEVEN_TSLFM;
    uint32_t FIBSEVEN_TSLFMI[2 /* FIBSEVEN_TSLFM */];
    uint32_t FIBSEVEN_TSLFME[2 /* FIBSEVEN_TSLFM */];
    uint32_t FIBSEVEN_TSLF;
    uint32_t FIBSEVEN_TSLFv[40 /* FIBSEVEN_TSLF */];
    uint32_t FIBSEVEN_TSTCM;
    uint32_t FIBSEVEN_TSTCMI[2 /* FIBSEVEN_TSTCM */];
    uint32_t FIBSEVEN_TSTCME[2 /* FIBSEVEN_TSTCM */];
    uint32_t FIBSEVEN_TSTC;
    uint32_t FIBSEVEN_TSTCv[40 /* FIBSEVEN_TSTC */];
    uint32_t FIBSEVEN_TSTFM;
    uint32_t FIBSEVEN_TSTFMI[2 /* FIBSEVEN_TSTFM */];
    uint32_t FIBSEVEN_TSTFME[2 /* FIBSEVEN_TSTFM */];
    uint32_t FIBSEVEN_TSTF;
    uint32_t FIBSEVEN_TSTFv[40 /* FIBSEVEN_TSTF */];

} EXT_STR_h101_FIBSEVEN_onion;

/*******************************************************/

#define EXT_STR_h101_FIBSEVEN_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                      \
    do                                                                                                            \
    {                                                                                                             \
        ok = 1;                                                                                                   \
        /* RAW */                                                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TMLCM, UINT32, "FIBSEVEN_TMLCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMLCMI, UINT32, "FIBSEVEN_TMLCMI", "FIBSEVEN_TMLCM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMLCME, UINT32, "FIBSEVEN_TMLCME", "FIBSEVEN_TMLCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TMLC, UINT32, "FIBSEVEN_TMLC", 10240); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMLCv, UINT32, "FIBSEVEN_TMLCv", "FIBSEVEN_TMLC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TMLFM, UINT32, "FIBSEVEN_TMLFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMLFMI, UINT32, "FIBSEVEN_TMLFMI", "FIBSEVEN_TMLFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMLFME, UINT32, "FIBSEVEN_TMLFME", "FIBSEVEN_TMLFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TMLF, UINT32, "FIBSEVEN_TMLF", 10240); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMLFv, UINT32, "FIBSEVEN_TMLFv", "FIBSEVEN_TMLF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TMTCM, UINT32, "FIBSEVEN_TMTCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMTCMI, UINT32, "FIBSEVEN_TMTCMI", "FIBSEVEN_TMTCM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMTCME, UINT32, "FIBSEVEN_TMTCME", "FIBSEVEN_TMTCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TMTC, UINT32, "FIBSEVEN_TMTC", 10240); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMTCv, UINT32, "FIBSEVEN_TMTCv", "FIBSEVEN_TMTC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TMTFM, UINT32, "FIBSEVEN_TMTFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMTFMI, UINT32, "FIBSEVEN_TMTFMI", "FIBSEVEN_TMTFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMTFME, UINT32, "FIBSEVEN_TMTFME", "FIBSEVEN_TMTFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TMTF, UINT32, "FIBSEVEN_TMTF", 10240); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TMTFv, UINT32, "FIBSEVEN_TMTFv", "FIBSEVEN_TMTF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TSLCM, UINT32, "FIBSEVEN_TSLCM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSLCMI, UINT32, "FIBSEVEN_TSLCMI", "FIBSEVEN_TSLCM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSLCME, UINT32, "FIBSEVEN_TSLCME", "FIBSEVEN_TSLCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TSLC, UINT32, "FIBSEVEN_TSLC", 40);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSLCv, UINT32, "FIBSEVEN_TSLCv", "FIBSEVEN_TSLC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TSLFM, UINT32, "FIBSEVEN_TSLFM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSLFMI, UINT32, "FIBSEVEN_TSLFMI", "FIBSEVEN_TSLFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSLFME, UINT32, "FIBSEVEN_TSLFME", "FIBSEVEN_TSLFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TSLF, UINT32, "FIBSEVEN_TSLF", 40);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSLFv, UINT32, "FIBSEVEN_TSLFv", "FIBSEVEN_TSLF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TSTCM, UINT32, "FIBSEVEN_TSTCM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSTCMI, UINT32, "FIBSEVEN_TSTCMI", "FIBSEVEN_TSTCM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSTCME, UINT32, "FIBSEVEN_TSTCME", "FIBSEVEN_TSTCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TSTC, UINT32, "FIBSEVEN_TSTC", 40);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSTCv, UINT32, "FIBSEVEN_TSTCv", "FIBSEVEN_TSTC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TSTFM, UINT32, "FIBSEVEN_TSTFM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSTFMI, UINT32, "FIBSEVEN_TSTFMI", "FIBSEVEN_TSTFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSTFME, UINT32, "FIBSEVEN_TSTFME", "FIBSEVEN_TSTFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBSEVEN_TSTF, UINT32, "FIBSEVEN_TSTF", 40);    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBSEVEN_TSTFv, UINT32, "FIBSEVEN_TSTFv", "FIBSEVEN_TSTF");       \
                                                                                                                  \
    } while (0);

#endif /*__GUARD_H101_FIBSEVEN_EXT_H101_FIBSEVEN_H__*/

/*******************************************************/
