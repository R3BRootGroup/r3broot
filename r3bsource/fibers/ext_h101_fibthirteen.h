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

#ifndef __GUARD_H101_FIBTHIRTEEN_EXT_H101_FIBTHIRTEEN_H__
#define __GUARD_H101_FIBTHIRTEEN_EXT_H101_FIBTHIRTEEN_H__

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

typedef struct EXT_STR_h101_FIBTHIRTEEN_t
{
    /* RAW */
    uint32_t FIBTHIRTEEN_TMLCM /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMLCMI[512 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMLCM)] /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMLCME[512 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMLCM)] /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMLC /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TMLCv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMLC)] /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TMTCM /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMTCMI[512 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMTCM)] /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMTCME[512 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMTCM)] /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMTC /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TMTCv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMTC)] /* [0,65535] */;
    uint32_t FIBTHIRTEEN_TMLFM /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMLFMI[512 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMLFM)] /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMLFME[512 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMLFM)] /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMLF /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TMLFv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMLF)] /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TMTFM /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMTFMI[512 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMTFM)] /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMTFME[512 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMTFM)] /* [1,512] */;
    uint32_t FIBTHIRTEEN_TMTF /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TMTFv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TMTF)] /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TRIGMLCM /* [1,4] */;
    uint32_t FIBTHIRTEEN_TRIGMLCMI[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TRIGMLCM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TRIGMLCME[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TRIGMLCM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TRIGMLC /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TRIGMLCv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TRIGMLC)] /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TRIGMLFM /* [1,4] */;
    uint32_t FIBTHIRTEEN_TRIGMLFMI[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TRIGMLFM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TRIGMLFME[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TRIGMLFM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TRIGMLF /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TRIGMLFv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TRIGMLF)] /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TSLCM /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSLCMI[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSLCM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSLCME[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSLCM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSLC /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TSLCv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSLC)] /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TSTCM /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSTCMI[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSTCM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSTCME[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSTCM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSTC /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TSTCv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSTC)] /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TSLFM /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSLFMI[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSLFM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSLFME[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSLFM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSLF /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TSLFv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSLF)] /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TSTFM /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSTFMI[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSTFM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSTFME[4 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSTFM)] /* [1,4] */;
    uint32_t FIBTHIRTEEN_TSTF /* [0,32768] */;
    uint32_t FIBTHIRTEEN_TSTFv[32768 EXT_STRUCT_CTRL(FIBTHIRTEEN_TSTF)] /* [0,32768] */;

} EXT_STR_h101_FIBTHIRTEEN;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBTHIRTEEN_onion_t
{
    /* RAW */
    uint32_t FIBTHIRTEEN_TMLCM;
    uint32_t FIBTHIRTEEN_TMLCMI[512 /* FIBTHIRTEEN_TMLCM */];
    uint32_t FIBTHIRTEEN_TMLCME[512 /* FIBTHIRTEEN_TMLCM */];
    uint32_t FIBTHIRTEEN_TMLC;
    uint32_t FIBTHIRTEEN_TMLCv[32768 /* FIBTHIRTEEN_TMLC */];
    uint32_t FIBTHIRTEEN_TMTCM;
    uint32_t FIBTHIRTEEN_TMTCMI[512 /* FIBTHIRTEEN_TMTCM */];
    uint32_t FIBTHIRTEEN_TMTCME[512 /* FIBTHIRTEEN_TMTCM */];
    uint32_t FIBTHIRTEEN_TMTC;
    uint32_t FIBTHIRTEEN_TMTCv[32768 /* FIBTHIRTEEN_TMTC */];
    uint32_t FIBTHIRTEEN_TMLFM;
    uint32_t FIBTHIRTEEN_TMLFMI[512 /* FIBTHIRTEEN_TMLFM */];
    uint32_t FIBTHIRTEEN_TMLFME[512 /* FIBTHIRTEEN_TMLFM */];
    uint32_t FIBTHIRTEEN_TMLF;
    uint32_t FIBTHIRTEEN_TMLFv[32768 /* FIBTHIRTEEN_TMLF */];
    uint32_t FIBTHIRTEEN_TMTFM;
    uint32_t FIBTHIRTEEN_TMTFMI[512 /* FIBTHIRTEEN_TMTFM */];
    uint32_t FIBTHIRTEEN_TMTFME[512 /* FIBTHIRTEEN_TMTFM */];
    uint32_t FIBTHIRTEEN_TMTF;
    uint32_t FIBTHIRTEEN_TMTFv[32768 /* FIBTHIRTEEN_TMTF */];
    uint32_t FIBTHIRTEEN_TRIGMLCM;
    uint32_t FIBTHIRTEEN_TRIGMLCMI[4 /* FIBTHIRTEEN_TRIGMLCM */];
    uint32_t FIBTHIRTEEN_TRIGMLCME[4 /* FIBTHIRTEEN_TRIGMLCM */];
    uint32_t FIBTHIRTEEN_TRIGMLC;
    uint32_t FIBTHIRTEEN_TRIGMLCv[32768 /* FIBTHIRTEEN_TRIGMLC */];
    uint32_t FIBTHIRTEEN_TRIGMLFM;
    uint32_t FIBTHIRTEEN_TRIGMLFMI[4 /* FIBTHIRTEEN_TRIGMLFM */];
    uint32_t FIBTHIRTEEN_TRIGMLFME[4 /* FIBTHIRTEEN_TRIGMLFM */];
    uint32_t FIBTHIRTEEN_TRIGMLF;
    uint32_t FIBTHIRTEEN_TRIGMLFv[32768 /* FIBTHIRTEEN_TRIGMLF */];
    uint32_t FIBTHIRTEEN_TSLCM;
    uint32_t FIBTHIRTEEN_TSLCMI[4 /* FIBTHIRTEEN_TSLCM */];
    uint32_t FIBTHIRTEEN_TSLCME[4 /* FIBTHIRTEEN_TSLCM */];
    uint32_t FIBTHIRTEEN_TSLC;
    uint32_t FIBTHIRTEEN_TSLCv[32768 /* FIBTHIRTEEN_TSLC */];
    uint32_t FIBTHIRTEEN_TSTCM;
    uint32_t FIBTHIRTEEN_TSTCMI[4 /* FIBTHIRTEEN_TSTCM */];
    uint32_t FIBTHIRTEEN_TSTCME[4 /* FIBTHIRTEEN_TSTCM */];
    uint32_t FIBTHIRTEEN_TSTC;
    uint32_t FIBTHIRTEEN_TSTCv[32768 /* FIBTHIRTEEN_TSTC */];
    uint32_t FIBTHIRTEEN_TSLFM;
    uint32_t FIBTHIRTEEN_TSLFMI[4 /* FIBTHIRTEEN_TSLFM */];
    uint32_t FIBTHIRTEEN_TSLFME[4 /* FIBTHIRTEEN_TSLFM */];
    uint32_t FIBTHIRTEEN_TSLF;
    uint32_t FIBTHIRTEEN_TSLFv[32768 /* FIBTHIRTEEN_TSLF */];
    uint32_t FIBTHIRTEEN_TSTFM;
    uint32_t FIBTHIRTEEN_TSTFMI[4 /* FIBTHIRTEEN_TSTFM */];
    uint32_t FIBTHIRTEEN_TSTFME[4 /* FIBTHIRTEEN_TSTFM */];
    uint32_t FIBTHIRTEEN_TSTF;
    uint32_t FIBTHIRTEEN_TSTFv[32768 /* FIBTHIRTEEN_TSTF */];

} EXT_STR_h101_FIBTHIRTEEN_onion;

/*******************************************************/

#define EXT_STR_h101_FIBTHIRTEEN_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                       \
    do                                                                                                                \
    {                                                                                                                 \
        ok = 1;                                                                                                       \
        /* RAW */                                                                                                     \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMLCM, UINT32, "FIBTHIRTEEN_TMLCM", 512);                 \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TMLCMI,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TMLCMI",                                                                   \
                              "FIBTHIRTEEN_TMLCM");                                                                   \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TMLCME,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TMLCME",                                                                   \
                              "FIBTHIRTEEN_TMLCM");                                                                   \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMLC, UINT32, "FIBTHIRTEEN_TMLC", 32768);                 \
        EXT_STR_ITEM_INFO_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMLCv, UINT32, "FIBTHIRTEEN_TMLCv", "FIBTHIRTEEN_TMLC");  \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMTCM, UINT32, "FIBTHIRTEEN_TMTCM", 512);                 \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TMTCMI,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TMTCMI",                                                                   \
                              "FIBTHIRTEEN_TMTCM");                                                                   \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TMTCME,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TMTCME",                                                                   \
                              "FIBTHIRTEEN_TMTCM");                                                                   \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMTC, UINT32, "FIBTHIRTEEN_TMTC", 32768);                 \
        EXT_STR_ITEM_INFO_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMTCv, UINT32, "FIBTHIRTEEN_TMTCv", "FIBTHIRTEEN_TMTC");  \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMLFM, UINT32, "FIBTHIRTEEN_TMLFM", 512);                 \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TMLFMI,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TMLFMI",                                                                   \
                              "FIBTHIRTEEN_TMLFM");                                                                   \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TMLFME,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TMLFME",                                                                   \
                              "FIBTHIRTEEN_TMLFM");                                                                   \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMLF, UINT32, "FIBTHIRTEEN_TMLF", 32768);                 \
        EXT_STR_ITEM_INFO_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMLFv, UINT32, "FIBTHIRTEEN_TMLFv", "FIBTHIRTEEN_TMLF");  \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMTFM, UINT32, "FIBTHIRTEEN_TMTFM", 512);                 \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TMTFMI,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TMTFMI",                                                                   \
                              "FIBTHIRTEEN_TMTFM");                                                                   \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TMTFME,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TMTFME",                                                                   \
                              "FIBTHIRTEEN_TMTFM");                                                                   \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMTF, UINT32, "FIBTHIRTEEN_TMTF", 32768);                 \
        EXT_STR_ITEM_INFO_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TMTFv, UINT32, "FIBTHIRTEEN_TMTFv", "FIBTHIRTEEN_TMTF");  \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TRIGMLCM, UINT32, "FIBTHIRTEEN_TRIGMLCM", 4);             \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TRIGMLCMI,                                                                  \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TRIGMLCMI",                                                                \
                              "FIBTHIRTEEN_TRIGMLCM");                                                                \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TRIGMLCME,                                                                  \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TRIGMLCME",                                                                \
                              "FIBTHIRTEEN_TRIGMLCM");                                                                \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TRIGMLC, UINT32, "FIBTHIRTEEN_TRIGMLC", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TRIGMLCv,                                                                   \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TRIGMLCv",                                                                 \
                              "FIBTHIRTEEN_TRIGMLC");                                                                 \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TRIGMLFM, UINT32, "FIBTHIRTEEN_TRIGMLFM", 4);             \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TRIGMLFMI,                                                                  \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TRIGMLFMI",                                                                \
                              "FIBTHIRTEEN_TRIGMLFM");                                                                \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TRIGMLFME,                                                                  \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TRIGMLFME",                                                                \
                              "FIBTHIRTEEN_TRIGMLFM");                                                                \
        EXT_STR_ITEM_INFO_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TRIGMLF, UINT32, "FIBTHIRTEEN_TRIGMLF", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TRIGMLFv,                                                                   \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TRIGMLFv",                                                                 \
                              "FIBTHIRTEEN_TRIGMLF");                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSLCM, UINT32, "FIBTHIRTEEN_TSLCM", 4); \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TSLCMI,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TSLCMI",                                                                   \
                              "FIBTHIRTEEN_TSLCM");                                                                   \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TSLCME,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TSLCME",                                                                   \
                              "FIBTHIRTEEN_TSLCM");                                                                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSLC, UINT32, "FIBTHIRTEEN_TSLC", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSLCv, UINT32, "FIBTHIRTEEN_TSLCv", "FIBTHIRTEEN_TSLC");  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSTCM, UINT32, "FIBTHIRTEEN_TSTCM", 4); \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TSTCMI,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TSTCMI",                                                                   \
                              "FIBTHIRTEEN_TSTCM");                                                                   \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TSTCME,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TSTCME",                                                                   \
                              "FIBTHIRTEEN_TSTCM");                                                                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSTC, UINT32, "FIBTHIRTEEN_TSTC", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSTCv, UINT32, "FIBTHIRTEEN_TSTCv", "FIBTHIRTEEN_TSTC");  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSLFM, UINT32, "FIBTHIRTEEN_TSLFM", 4); \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TSLFMI,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TSLFMI",                                                                   \
                              "FIBTHIRTEEN_TSLFM");                                                                   \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TSLFME,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TSLFME",                                                                   \
                              "FIBTHIRTEEN_TSLFM");                                                                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSLF, UINT32, "FIBTHIRTEEN_TSLF", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSLFv, UINT32, "FIBTHIRTEEN_TSLFv", "FIBTHIRTEEN_TSLF");  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSTFM, UINT32, "FIBTHIRTEEN_TSTFM", 4); \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TSTFMI,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TSTFMI",                                                                   \
                              "FIBTHIRTEEN_TSTFM");                                                                   \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                     \
                              si,                                                                                     \
                              offset,                                                                                 \
                              struct_t,                                                                               \
                              printerr,                                                                               \
                              FIBTHIRTEEN_TSTFME,                                                                     \
                              UINT32,                                                                                 \
                              "FIBTHIRTEEN_TSTFME",                                                                   \
                              "FIBTHIRTEEN_TSTFM");                                                                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSTF, UINT32, "FIBTHIRTEEN_TSTF", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, FIBTHIRTEEN_TSTFv, UINT32, "FIBTHIRTEEN_TSTFv", "FIBTHIRTEEN_TSTF");  \
                                                                                                                      \
    } while (0);

#endif /*__GUARD_H101_FIBTHIRTEEN_EXT_H101_FIBTHIRTEEN_H__*/

/*******************************************************/
