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

#ifndef __GUARD_H101_FIBEIGHT_EXT_H101_FIBEIGHT_H__
#define __GUARD_H101_FIBEIGHT_EXT_H101_FIBEIGHT_H__

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

typedef struct EXT_STR_h101_FIBEIGHT_t
{
    /* RAW */
    uint32_t FIBEIGHT_TMLCM /* [1,512] */;
    uint32_t FIBEIGHT_TMLCMI[512 EXT_STRUCT_CTRL(FIBEIGHT_TMLCM)] /* [1,512] */;
    uint32_t FIBEIGHT_TMLCME[512 EXT_STRUCT_CTRL(FIBEIGHT_TMLCM)] /* [1,32768] */;
    uint32_t FIBEIGHT_TMLC /* [0,32768] */;
    uint32_t FIBEIGHT_TMLCv[32768 EXT_STRUCT_CTRL(FIBEIGHT_TMLC)] /* [0,65535] */;
    uint32_t FIBEIGHT_TMLFM /* [1,512] */;
    uint32_t FIBEIGHT_TMLFMI[512 EXT_STRUCT_CTRL(FIBEIGHT_TMLFM)] /* [1,512] */;
    uint32_t FIBEIGHT_TMLFME[512 EXT_STRUCT_CTRL(FIBEIGHT_TMLFM)] /* [1,32768] */;
    uint32_t FIBEIGHT_TMLF /* [0,32768] */;
    uint32_t FIBEIGHT_TMLFv[32768 EXT_STRUCT_CTRL(FIBEIGHT_TMLF)] /* [0,65535] */;
    uint32_t FIBEIGHT_TMTCM /* [1,512] */;
    uint32_t FIBEIGHT_TMTCMI[512 EXT_STRUCT_CTRL(FIBEIGHT_TMTCM)] /* [1,512] */;
    uint32_t FIBEIGHT_TMTCME[512 EXT_STRUCT_CTRL(FIBEIGHT_TMTCM)] /* [1,32768] */;
    uint32_t FIBEIGHT_TMTC /* [0,32768] */;
    uint32_t FIBEIGHT_TMTCv[32768 EXT_STRUCT_CTRL(FIBEIGHT_TMTC)] /* [0,65535] */;
    uint32_t FIBEIGHT_TMTFM /* [1,512] */;
    uint32_t FIBEIGHT_TMTFMI[512 EXT_STRUCT_CTRL(FIBEIGHT_TMTFM)] /* [1,512] */;
    uint32_t FIBEIGHT_TMTFME[512 EXT_STRUCT_CTRL(FIBEIGHT_TMTFM)] /* [1,32768] */;
    uint32_t FIBEIGHT_TMTF /* [0,32768] */;
    uint32_t FIBEIGHT_TMTFv[32768 EXT_STRUCT_CTRL(FIBEIGHT_TMTF)] /* [0,65535] */;
    uint32_t FIBEIGHT_TSLCM /* [1,2] */;
    uint32_t FIBEIGHT_TSLCMI[2 EXT_STRUCT_CTRL(FIBEIGHT_TSLCM)] /* [1,2] */;
    uint32_t FIBEIGHT_TSLCME[2 EXT_STRUCT_CTRL(FIBEIGHT_TSLCM)] /* [1,128] */;
    uint32_t FIBEIGHT_TSLC /* [0,128] */;
    uint32_t FIBEIGHT_TSLCv[128 EXT_STRUCT_CTRL(FIBEIGHT_TSLC)] /* [0,65535] */;
    uint32_t FIBEIGHT_TSLFM /* [1,2] */;
    uint32_t FIBEIGHT_TSLFMI[2 EXT_STRUCT_CTRL(FIBEIGHT_TSLFM)] /* [1,2] */;
    uint32_t FIBEIGHT_TSLFME[2 EXT_STRUCT_CTRL(FIBEIGHT_TSLFM)] /* [1,128] */;
    uint32_t FIBEIGHT_TSLF /* [0,128] */;
    uint32_t FIBEIGHT_TSLFv[128 EXT_STRUCT_CTRL(FIBEIGHT_TSLF)] /* [0,65535] */;
    uint32_t FIBEIGHT_TSTCM /* [1,2] */;
    uint32_t FIBEIGHT_TSTCMI[2 EXT_STRUCT_CTRL(FIBEIGHT_TSTCM)] /* [1,2] */;
    uint32_t FIBEIGHT_TSTCME[2 EXT_STRUCT_CTRL(FIBEIGHT_TSTCM)] /* [1,128] */;
    uint32_t FIBEIGHT_TSTC /* [0,128] */;
    uint32_t FIBEIGHT_TSTCv[128 EXT_STRUCT_CTRL(FIBEIGHT_TSTC)] /* [0,65535] */;
    uint32_t FIBEIGHT_TSTFM /* [1,2] */;
    uint32_t FIBEIGHT_TSTFMI[2 EXT_STRUCT_CTRL(FIBEIGHT_TSTFM)] /* [1,2] */;
    uint32_t FIBEIGHT_TSTFME[2 EXT_STRUCT_CTRL(FIBEIGHT_TSTFM)] /* [1,128] */;
    uint32_t FIBEIGHT_TSTF /* [0,128] */;
    uint32_t FIBEIGHT_TSTFv[128 EXT_STRUCT_CTRL(FIBEIGHT_TSTF)] /* [0,65535] */;

} EXT_STR_h101_FIBEIGHT;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBEIGHT_onion_t
{
    /* RAW */
    uint32_t FIBEIGHT_TMLCM;
    uint32_t FIBEIGHT_TMLCMI[512 /* FIBEIGHT_TMLCM */];
    uint32_t FIBEIGHT_TMLCME[512 /* FIBEIGHT_TMLCM */];
    uint32_t FIBEIGHT_TMLC;
    uint32_t FIBEIGHT_TMLCv[32768 /* FIBEIGHT_TMLC */];
    uint32_t FIBEIGHT_TMLFM;
    uint32_t FIBEIGHT_TMLFMI[512 /* FIBEIGHT_TMLFM */];
    uint32_t FIBEIGHT_TMLFME[512 /* FIBEIGHT_TMLFM */];
    uint32_t FIBEIGHT_TMLF;
    uint32_t FIBEIGHT_TMLFv[32768 /* FIBEIGHT_TMLF */];
    uint32_t FIBEIGHT_TMTCM;
    uint32_t FIBEIGHT_TMTCMI[512 /* FIBEIGHT_TMTCM */];
    uint32_t FIBEIGHT_TMTCME[512 /* FIBEIGHT_TMTCM */];
    uint32_t FIBEIGHT_TMTC;
    uint32_t FIBEIGHT_TMTCv[32768 /* FIBEIGHT_TMTC */];
    uint32_t FIBEIGHT_TMTFM;
    uint32_t FIBEIGHT_TMTFMI[512 /* FIBEIGHT_TMTFM */];
    uint32_t FIBEIGHT_TMTFME[512 /* FIBEIGHT_TMTFM */];
    uint32_t FIBEIGHT_TMTF;
    uint32_t FIBEIGHT_TMTFv[32768 /* FIBEIGHT_TMTF */];
    uint32_t FIBEIGHT_TSLCM;
    uint32_t FIBEIGHT_TSLCMI[2 /* FIBEIGHT_TSLCM */];
    uint32_t FIBEIGHT_TSLCME[2 /* FIBEIGHT_TSLCM */];
    uint32_t FIBEIGHT_TSLC;
    uint32_t FIBEIGHT_TSLCv[128 /* FIBEIGHT_TSLC */];
    uint32_t FIBEIGHT_TSLFM;
    uint32_t FIBEIGHT_TSLFMI[2 /* FIBEIGHT_TSLFM */];
    uint32_t FIBEIGHT_TSLFME[2 /* FIBEIGHT_TSLFM */];
    uint32_t FIBEIGHT_TSLF;
    uint32_t FIBEIGHT_TSLFv[128 /* FIBEIGHT_TSLF */];
    uint32_t FIBEIGHT_TSTCM;
    uint32_t FIBEIGHT_TSTCMI[2 /* FIBEIGHT_TSTCM */];
    uint32_t FIBEIGHT_TSTCME[2 /* FIBEIGHT_TSTCM */];
    uint32_t FIBEIGHT_TSTC;
    uint32_t FIBEIGHT_TSTCv[128 /* FIBEIGHT_TSTC */];
    uint32_t FIBEIGHT_TSTFM;
    uint32_t FIBEIGHT_TSTFMI[2 /* FIBEIGHT_TSTFM */];
    uint32_t FIBEIGHT_TSTFME[2 /* FIBEIGHT_TSTFM */];
    uint32_t FIBEIGHT_TSTF;
    uint32_t FIBEIGHT_TSTFv[128 /* FIBEIGHT_TSTF */];

} EXT_STR_h101_FIBEIGHT_onion;

/*******************************************************/

#define EXT_STR_h101_FIBEIGHT_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                      \
    do                                                                                                            \
    {                                                                                                             \
        ok = 1;                                                                                                   \
        /* RAW */                                                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TMLCM, UINT32, "FIBEIGHT_TMLCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMLCMI, UINT32, "FIBEIGHT_TMLCMI", "FIBEIGHT_TMLCM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMLCME, UINT32, "FIBEIGHT_TMLCME", "FIBEIGHT_TMLCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TMLC, UINT32, "FIBEIGHT_TMLC", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMLCv, UINT32, "FIBEIGHT_TMLCv", "FIBEIGHT_TMLC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TMLFM, UINT32, "FIBEIGHT_TMLFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMLFMI, UINT32, "FIBEIGHT_TMLFMI", "FIBEIGHT_TMLFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMLFME, UINT32, "FIBEIGHT_TMLFME", "FIBEIGHT_TMLFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TMLF, UINT32, "FIBEIGHT_TMLF", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMLFv, UINT32, "FIBEIGHT_TMLFv", "FIBEIGHT_TMLF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TMTCM, UINT32, "FIBEIGHT_TMTCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMTCMI, UINT32, "FIBEIGHT_TMTCMI", "FIBEIGHT_TMTCM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMTCME, UINT32, "FIBEIGHT_TMTCME", "FIBEIGHT_TMTCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TMTC, UINT32, "FIBEIGHT_TMTC", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMTCv, UINT32, "FIBEIGHT_TMTCv", "FIBEIGHT_TMTC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TMTFM, UINT32, "FIBEIGHT_TMTFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMTFMI, UINT32, "FIBEIGHT_TMTFMI", "FIBEIGHT_TMTFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMTFME, UINT32, "FIBEIGHT_TMTFME", "FIBEIGHT_TMTFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TMTF, UINT32, "FIBEIGHT_TMTF", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TMTFv, UINT32, "FIBEIGHT_TMTFv", "FIBEIGHT_TMTF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TSLCM, UINT32, "FIBEIGHT_TSLCM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSLCMI, UINT32, "FIBEIGHT_TSLCMI", "FIBEIGHT_TSLCM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSLCME, UINT32, "FIBEIGHT_TSLCME", "FIBEIGHT_TSLCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TSLC, UINT32, "FIBEIGHT_TSLC", 128);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSLCv, UINT32, "FIBEIGHT_TSLCv", "FIBEIGHT_TSLC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TSLFM, UINT32, "FIBEIGHT_TSLFM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSLFMI, UINT32, "FIBEIGHT_TSLFMI", "FIBEIGHT_TSLFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSLFME, UINT32, "FIBEIGHT_TSLFME", "FIBEIGHT_TSLFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TSLF, UINT32, "FIBEIGHT_TSLF", 128);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSLFv, UINT32, "FIBEIGHT_TSLFv", "FIBEIGHT_TSLF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TSTCM, UINT32, "FIBEIGHT_TSTCM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSTCMI, UINT32, "FIBEIGHT_TSTCMI", "FIBEIGHT_TSTCM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSTCME, UINT32, "FIBEIGHT_TSTCME", "FIBEIGHT_TSTCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TSTC, UINT32, "FIBEIGHT_TSTC", 128);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSTCv, UINT32, "FIBEIGHT_TSTCv", "FIBEIGHT_TSTC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TSTFM, UINT32, "FIBEIGHT_TSTFM", 2);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSTFMI, UINT32, "FIBEIGHT_TSTFMI", "FIBEIGHT_TSTFM");    \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSTFME, UINT32, "FIBEIGHT_TSTFME", "FIBEIGHT_TSTFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEIGHT_TSTF, UINT32, "FIBEIGHT_TSTF", 128);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBEIGHT_TSTFv, UINT32, "FIBEIGHT_TSTFv", "FIBEIGHT_TSTF");       \
                                                                                                                  \
    } while (0);

#endif /*__GUARD_H101_FIBEIGHT_EXT_H101_FIBEIGHT_H__*/

/*******************************************************/
