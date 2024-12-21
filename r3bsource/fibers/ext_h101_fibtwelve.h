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

#ifndef __GUARD_H101_FIBTWELVE_EXT_H101_FIBTWELVE_H__
#define __GUARD_H101_FIBTWELVE_EXT_H101_FIBTWELVE_H__

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

typedef struct EXT_STR_h101_FIBTWELVE_t
{
    /* RAW */
    uint32_t FIBTWELVE_TMLCM /* [1,512] */;
    uint32_t FIBTWELVE_TMLCMI[512 EXT_STRUCT_CTRL(FIBTWELVE_TMLCM)] /* [1,512] */;
    uint32_t FIBTWELVE_TMLCME[512 EXT_STRUCT_CTRL(FIBTWELVE_TMLCM)] /* [1,512] */;
    uint32_t FIBTWELVE_TMLC /* [0,32768] */;
    uint32_t FIBTWELVE_TMLCv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TMLC)] /* [0,32768] */;
    uint32_t FIBTWELVE_TMTCM /* [1,512] */;
    uint32_t FIBTWELVE_TMTCMI[512 EXT_STRUCT_CTRL(FIBTWELVE_TMTCM)] /* [1,512] */;
    uint32_t FIBTWELVE_TMTCME[512 EXT_STRUCT_CTRL(FIBTWELVE_TMTCM)] /* [1,512] */;
    uint32_t FIBTWELVE_TMTC /* [0,32768] */;
    uint32_t FIBTWELVE_TMTCv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TMTC)] /* [0,32768] */;
    uint32_t FIBTWELVE_TMLFM /* [1,512] */;
    uint32_t FIBTWELVE_TMLFMI[512 EXT_STRUCT_CTRL(FIBTWELVE_TMLFM)] /* [1,512] */;
    uint32_t FIBTWELVE_TMLFME[512 EXT_STRUCT_CTRL(FIBTWELVE_TMLFM)] /* [1,512] */;
    uint32_t FIBTWELVE_TMLF /* [0,32768] */;
    uint32_t FIBTWELVE_TMLFv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TMLF)] /* [0,32768] */;
    uint32_t FIBTWELVE_TMTFM /* [1,512] */;
    uint32_t FIBTWELVE_TMTFMI[512 EXT_STRUCT_CTRL(FIBTWELVE_TMTFM)] /* [1,512] */;
    uint32_t FIBTWELVE_TMTFME[512 EXT_STRUCT_CTRL(FIBTWELVE_TMTFM)] /* [1,512] */;
    uint32_t FIBTWELVE_TMTF /* [0,32768] */;
    uint32_t FIBTWELVE_TMTFv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TMTF)] /* [0,32768] */;
    uint32_t FIBTWELVE_TRIGMLCM /* [1,4] */;
    uint32_t FIBTWELVE_TRIGMLCMI[4 EXT_STRUCT_CTRL(FIBTWELVE_TRIGMLCM)] /* [1,4] */;
    uint32_t FIBTWELVE_TRIGMLCME[4 EXT_STRUCT_CTRL(FIBTWELVE_TRIGMLCM)] /* [1,4] */;
    uint32_t FIBTWELVE_TRIGMLC /* [0,32768] */;
    uint32_t FIBTWELVE_TRIGMLCv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TRIGMLC)] /* [0,32768] */;
    uint32_t FIBTWELVE_TRIGMLFM /* [1,4] */;
    uint32_t FIBTWELVE_TRIGMLFMI[4 EXT_STRUCT_CTRL(FIBTWELVE_TRIGMLFM)] /* [1,4] */;
    uint32_t FIBTWELVE_TRIGMLFME[4 EXT_STRUCT_CTRL(FIBTWELVE_TRIGMLFM)] /* [1,4] */;
    uint32_t FIBTWELVE_TRIGMLF /* [0,32768] */;
    uint32_t FIBTWELVE_TRIGMLFv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TRIGMLF)] /* [0,32768] */;
    uint32_t FIBTWELVE_TSLCM /* [1,4] */;
    uint32_t FIBTWELVE_TSLCMI[4 EXT_STRUCT_CTRL(FIBTWELVE_TSLCM)] /* [1,4] */;
    uint32_t FIBTWELVE_TSLCME[4 EXT_STRUCT_CTRL(FIBTWELVE_TSLCM)] /* [1,4] */;
    uint32_t FIBTWELVE_TSLC /* [0,32768] */;
    uint32_t FIBTWELVE_TSLCv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TSLC)] /* [0,32768] */;
    uint32_t FIBTWELVE_TSTCM /* [1,4] */;
    uint32_t FIBTWELVE_TSTCMI[4 EXT_STRUCT_CTRL(FIBTWELVE_TSTCM)] /* [1,4] */;
    uint32_t FIBTWELVE_TSTCME[4 EXT_STRUCT_CTRL(FIBTWELVE_TSTCM)] /* [1,4] */;
    uint32_t FIBTWELVE_TSTC /* [0,32768] */;
    uint32_t FIBTWELVE_TSTCv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TSTC)] /* [0,65535] */;
    uint32_t FIBTWELVE_TSLFM /* [1,4] */;
    uint32_t FIBTWELVE_TSLFMI[4 EXT_STRUCT_CTRL(FIBTWELVE_TSLFM)] /* [1,4] */;
    uint32_t FIBTWELVE_TSLFME[4 EXT_STRUCT_CTRL(FIBTWELVE_TSLFM)] /* [1,4] */;
    uint32_t FIBTWELVE_TSLF /* [0,32768] */;
    uint32_t FIBTWELVE_TSLFv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TSLF)] /* [0,65535] */;
    uint32_t FIBTWELVE_TSTFM /* [1,4] */;
    uint32_t FIBTWELVE_TSTFMI[4 EXT_STRUCT_CTRL(FIBTWELVE_TSTFM)] /* [1,4] */;
    uint32_t FIBTWELVE_TSTFME[4 EXT_STRUCT_CTRL(FIBTWELVE_TSTFM)] /* [1,4] */;
    uint32_t FIBTWELVE_TSTF /* [0,32768] */;
    uint32_t FIBTWELVE_TSTFv[32768 EXT_STRUCT_CTRL(FIBTWELVE_TSTF)] /* [0,32768] */;

} EXT_STR_h101_FIBTWELVE;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBTWELVE_onion_t
{
    /* RAW */
    uint32_t FIBTWELVE_TMLCM;
    uint32_t FIBTWELVE_TMLCMI[512 /* FIBTWELVE_TMLCM */];
    uint32_t FIBTWELVE_TMLCME[512 /* FIBTWELVE_TMLCM */];
    uint32_t FIBTWELVE_TMLC;
    uint32_t FIBTWELVE_TMLCv[32768 /* FIBTWELVE_TMLC */];
    uint32_t FIBTWELVE_TMTCM;
    uint32_t FIBTWELVE_TMTCMI[512 /* FIBTWELVE_TMTCM */];
    uint32_t FIBTWELVE_TMTCME[512 /* FIBTWELVE_TMTCM */];
    uint32_t FIBTWELVE_TMTC;
    uint32_t FIBTWELVE_TMTCv[32768 /* FIBTWELVE_TMTC */];
    uint32_t FIBTWELVE_TMLFM;
    uint32_t FIBTWELVE_TMLFMI[512 /* FIBTWELVE_TMLFM */];
    uint32_t FIBTWELVE_TMLFME[512 /* FIBTWELVE_TMLFM */];
    uint32_t FIBTWELVE_TMLF;
    uint32_t FIBTWELVE_TMLFv[32768 /* FIBTWELVE_TMLF */];
    uint32_t FIBTWELVE_TMTFM;
    uint32_t FIBTWELVE_TMTFMI[512 /* FIBTWELVE_TMTFM */];
    uint32_t FIBTWELVE_TMTFME[512 /* FIBTWELVE_TMTFM */];
    uint32_t FIBTWELVE_TMTF;
    uint32_t FIBTWELVE_TMTFv[32768 /* FIBTWELVE_TMTF */];
    uint32_t FIBTWELVE_TRIGMLCM;
    uint32_t FIBTWELVE_TRIGMLCMI[4 /* FIBTWELVE_TRIGMLCM */];
    uint32_t FIBTWELVE_TRIGMLCME[4 /* FIBTWELVE_TRIGMLCM */];
    uint32_t FIBTWELVE_TRIGMLC;
    uint32_t FIBTWELVE_TRIGMLCv[32768 /* FIBTWELVE_TRIGMLC */];
    uint32_t FIBTWELVE_TRIGMLFM;
    uint32_t FIBTWELVE_TRIGMLFMI[4 /* FIBTWELVE_TRIGMLFM */];
    uint32_t FIBTWELVE_TRIGMLFME[4 /* FIBTWELVE_TRIGMLFM */];
    uint32_t FIBTWELVE_TRIGMLF;
    uint32_t FIBTWELVE_TRIGMLFv[32768 /* FIBTWELVE_TRIGMLF */];
    uint32_t FIBTWELVE_TSLCM;
    uint32_t FIBTWELVE_TSLCMI[4 /* FIBTWELVE_TSLCM */];
    uint32_t FIBTWELVE_TSLCME[4 /* FIBTWELVE_TSLCM */];
    uint32_t FIBTWELVE_TSLC;
    uint32_t FIBTWELVE_TSLCv[32768 /* FIBTWELVE_TSLC */];
    uint32_t FIBTWELVE_TSTCM;
    uint32_t FIBTWELVE_TSTCMI[4 /* FIBTWELVE_TSTCM */];
    uint32_t FIBTWELVE_TSTCME[4 /* FIBTWELVE_TSTCM */];
    uint32_t FIBTWELVE_TSTC;
    uint32_t FIBTWELVE_TSTCv[32768 /* FIBTWELVE_TSTC */];
    uint32_t FIBTWELVE_TSLFM;
    uint32_t FIBTWELVE_TSLFMI[4 /* FIBTWELVE_TSLFM */];
    uint32_t FIBTWELVE_TSLFME[4 /* FIBTWELVE_TSLFM */];
    uint32_t FIBTWELVE_TSLF;
    uint32_t FIBTWELVE_TSLFv[32768 /* FIBTWELVE_TSLF */];
    uint32_t FIBTWELVE_TSTFM;
    uint32_t FIBTWELVE_TSTFMI[4 /* FIBTWELVE_TSTFM */];
    uint32_t FIBTWELVE_TSTFME[4 /* FIBTWELVE_TSTFM */];
    uint32_t FIBTWELVE_TSTF;
    uint32_t FIBTWELVE_TSTFv[32768 /* FIBTWELVE_TSTF */];

} EXT_STR_h101_FIBTWELVE_onion;

/*******************************************************/

#define EXT_STR_h101_FIBTWELVE_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                       \
    do                                                                                                              \
    {                                                                                                               \
        ok = 1;                                                                                                     \
        /* RAW */                                                                                                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TMLCM, UINT32, "FIBTWELVE_TMLCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMLCMI, UINT32, "FIBTWELVE_TMLCMI", "FIBTWELVE_TMLCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMLCME, UINT32, "FIBTWELVE_TMLCME", "FIBTWELVE_TMLCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TMLC, UINT32, "FIBTWELVE_TMLC", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMLCv, UINT32, "FIBTWELVE_TMLCv", "FIBTWELVE_TMLC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TMTCM, UINT32, "FIBTWELVE_TMTCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMTCMI, UINT32, "FIBTWELVE_TMTCMI", "FIBTWELVE_TMTCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMTCME, UINT32, "FIBTWELVE_TMTCME", "FIBTWELVE_TMTCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TMTC, UINT32, "FIBTWELVE_TMTC", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMTCv, UINT32, "FIBTWELVE_TMTCv", "FIBTWELVE_TMTC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TMLFM, UINT32, "FIBTWELVE_TMLFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMLFMI, UINT32, "FIBTWELVE_TMLFMI", "FIBTWELVE_TMLFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMLFME, UINT32, "FIBTWELVE_TMLFME", "FIBTWELVE_TMLFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TMLF, UINT32, "FIBTWELVE_TMLF", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMLFv, UINT32, "FIBTWELVE_TMLFv", "FIBTWELVE_TMLF");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TMTFM, UINT32, "FIBTWELVE_TMTFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMTFMI, UINT32, "FIBTWELVE_TMTFMI", "FIBTWELVE_TMTFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMTFME, UINT32, "FIBTWELVE_TMTFME", "FIBTWELVE_TMTFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TMTF, UINT32, "FIBTWELVE_TMTF", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TMTFv, UINT32, "FIBTWELVE_TMTFv", "FIBTWELVE_TMTF");      \
        EXT_STR_ITEM_INFO_LIM(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TRIGMLCM, UINT32, "FIBTWELVE_TRIGMLCM", 4);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBTWELVE_TRIGMLCMI,                                                                  \
                              UINT32,                                                                               \
                              "FIBTWELVE_TRIGMLCMI",                                                                \
                              "FIBTWELVE_TRIGMLCM");                                                                \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBTWELVE_TRIGMLCME,                                                                  \
                              UINT32,                                                                               \
                              "FIBTWELVE_TRIGMLCME",                                                                \
                              "FIBTWELVE_TRIGMLCM");                                                                \
        EXT_STR_ITEM_INFO_LIM(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TRIGMLC, UINT32, "FIBTWELVE_TRIGMLC", 256);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBTWELVE_TRIGMLCv,                                                                   \
                              UINT32,                                                                               \
                              "FIBTWELVE_TRIGMLCv",                                                                 \
                              "FIBTWELVE_TRIGMLC");                                                                 \
        EXT_STR_ITEM_INFO_LIM(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TRIGMLFM, UINT32, "FIBTWELVE_TRIGMLFM", 4);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBTWELVE_TRIGMLFMI,                                                                  \
                              UINT32,                                                                               \
                              "FIBTWELVE_TRIGMLFMI",                                                                \
                              "FIBTWELVE_TRIGMLFM");                                                                \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBTWELVE_TRIGMLFME,                                                                  \
                              UINT32,                                                                               \
                              "FIBTWELVE_TRIGMLFME",                                                                \
                              "FIBTWELVE_TRIGMLFM");                                                                \
        EXT_STR_ITEM_INFO_LIM(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TRIGMLF, UINT32, "FIBTWELVE_TRIGMLF", 256);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBTWELVE_TRIGMLFv,                                                                   \
                              UINT32,                                                                               \
                              "FIBTWELVE_TRIGMLFv",                                                                 \
                              "FIBTWELVE_TRIGMLF");                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TSLCM, UINT32, "FIBTWELVE_TSLCM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSLCMI, UINT32, "FIBTWELVE_TSLCMI", "FIBTWELVE_TSLCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSLCME, UINT32, "FIBTWELVE_TSLCME", "FIBTWELVE_TSLCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TSLC, UINT32, "FIBTWELVE_TSLC", 256);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSLCv, UINT32, "FIBTWELVE_TSLCv", "FIBTWELVE_TSLC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TSTCM, UINT32, "FIBTWELVE_TSTCM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSTCMI, UINT32, "FIBTWELVE_TSTCMI", "FIBTWELVE_TSTCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSTCME, UINT32, "FIBTWELVE_TSTCME", "FIBTWELVE_TSTCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TSTC, UINT32, "FIBTWELVE_TSTC", 256);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSTCv, UINT32, "FIBTWELVE_TSTCv", "FIBTWELVE_TSTC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TSLFM, UINT32, "FIBTWELVE_TSLFM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSLFMI, UINT32, "FIBTWELVE_TSLFMI", "FIBTWELVE_TSLFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSLFME, UINT32, "FIBTWELVE_TSLFME", "FIBTWELVE_TSLFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TSLF, UINT32, "FIBTWELVE_TSLF", 256);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSLFv, UINT32, "FIBTWELVE_TSLFv", "FIBTWELVE_TSLF");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TSTFM, UINT32, "FIBTWELVE_TSTFM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSTFMI, UINT32, "FIBTWELVE_TSTFMI", "FIBTWELVE_TSTFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSTFME, UINT32, "FIBTWELVE_TSTFME", "FIBTWELVE_TSTFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTWELVE_TSTF, UINT32, "FIBTWELVE_TSTF", 256);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBTWELVE_TSTFv, UINT32, "FIBTWELVE_TSTFv", "FIBTWELVE_TSTF");      \
                                                                                                                    \
    } while (0);

#endif /*__GUARD_H101_FIBTWELVE_EXT_H101_FIBTWELVE_H__*/

/*******************************************************/
