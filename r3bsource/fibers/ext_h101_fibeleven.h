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

#ifndef __GUARD_H101_FIBELEVEN_EXT_H101_FIBELEVEN_H__
#define __GUARD_H101_FIBELEVEN_EXT_H101_FIBELEVEN_H__

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

typedef struct EXT_STR_h101_FIBELEVEN_t
{
    /* RAW */
    uint32_t FIBELEVEN_TMLCM /* [1,512] */;
    uint32_t FIBELEVEN_TMLCMI[512 EXT_STRUCT_CTRL(FIBELEVEN_TMLCM)] /* [1,512] */;
    uint32_t FIBELEVEN_TMLCME[512 EXT_STRUCT_CTRL(FIBELEVEN_TMLCM)] /* [1,512] */;
    uint32_t FIBELEVEN_TMLC /* [0,32768] */;
    uint32_t FIBELEVEN_TMLCv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TMLC)] /* [0,32768] */;
    uint32_t FIBELEVEN_TMTCM /* [1,512] */;
    uint32_t FIBELEVEN_TMTCMI[512 EXT_STRUCT_CTRL(FIBELEVEN_TMTCM)] /* [1,512] */;
    uint32_t FIBELEVEN_TMTCME[512 EXT_STRUCT_CTRL(FIBELEVEN_TMTCM)] /* [1,512] */;
    uint32_t FIBELEVEN_TMTC /* [0,32768] */;
    uint32_t FIBELEVEN_TMTCv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TMTC)] /* [0,32768] */;
    uint32_t FIBELEVEN_TMLFM /* [1,512] */;
    uint32_t FIBELEVEN_TMLFMI[512 EXT_STRUCT_CTRL(FIBELEVEN_TMLFM)] /* [1,512] */;
    uint32_t FIBELEVEN_TMLFME[512 EXT_STRUCT_CTRL(FIBELEVEN_TMLFM)] /* [1,512] */;
    uint32_t FIBELEVEN_TMLF /* [0,32768] */;
    uint32_t FIBELEVEN_TMLFv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TMLF)] /* [0,32768] */;
    uint32_t FIBELEVEN_TMTFM /* [1,512] */;
    uint32_t FIBELEVEN_TMTFMI[512 EXT_STRUCT_CTRL(FIBELEVEN_TMTFM)] /* [1,512] */;
    uint32_t FIBELEVEN_TMTFME[512 EXT_STRUCT_CTRL(FIBELEVEN_TMTFM)] /* [1,512] */;
    uint32_t FIBELEVEN_TMTF /* [0,32768] */;
    uint32_t FIBELEVEN_TMTFv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TMTF)] /* [0,32768] */;
    uint32_t FIBELEVEN_TRIGMLCM /* [1,4] */;
    uint32_t FIBELEVEN_TRIGMLCMI[4 EXT_STRUCT_CTRL(FIBELEVEN_TRIGMLCM)] /* [1,4] */;
    uint32_t FIBELEVEN_TRIGMLCME[4 EXT_STRUCT_CTRL(FIBELEVEN_TRIGMLCM)] /* [1,4] */;
    uint32_t FIBELEVEN_TRIGMLC /* [0,32768] */;
    uint32_t FIBELEVEN_TRIGMLCv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TRIGMLC)] /* [0,32768] */;
    uint32_t FIBELEVEN_TRIGMLFM /* [1,4] */;
    uint32_t FIBELEVEN_TRIGMLFMI[4 EXT_STRUCT_CTRL(FIBELEVEN_TRIGMLFM)] /* [1,4] */;
    uint32_t FIBELEVEN_TRIGMLFME[4 EXT_STRUCT_CTRL(FIBELEVEN_TRIGMLFM)] /* [1,4] */;
    uint32_t FIBELEVEN_TRIGMLF /* [0,32768] */;
    uint32_t FIBELEVEN_TRIGMLFv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TRIGMLF)] /* [0,32768] */;
    uint32_t FIBELEVEN_TSLCM /* [1,4] */;
    uint32_t FIBELEVEN_TSLCMI[4 EXT_STRUCT_CTRL(FIBELEVEN_TSLCM)] /* [1,4] */;
    uint32_t FIBELEVEN_TSLCME[4 EXT_STRUCT_CTRL(FIBELEVEN_TSLCM)] /* [1,4] */;
    uint32_t FIBELEVEN_TSLC /* [0,32768] */;
    uint32_t FIBELEVEN_TSLCv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TSLC)] /* [0,32768] */;
    uint32_t FIBELEVEN_TSTCM /* [1,4] */;
    uint32_t FIBELEVEN_TSTCMI[4 EXT_STRUCT_CTRL(FIBELEVEN_TSTCM)] /* [1,4] */;
    uint32_t FIBELEVEN_TSTCME[4 EXT_STRUCT_CTRL(FIBELEVEN_TSTCM)] /* [1,4] */;
    uint32_t FIBELEVEN_TSTC /* [0,32768] */;
    uint32_t FIBELEVEN_TSTCv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TSTC)] /* [0,32768] */;
    uint32_t FIBELEVEN_TSLFM /* [1,4] */;
    uint32_t FIBELEVEN_TSLFMI[4 EXT_STRUCT_CTRL(FIBELEVEN_TSLFM)] /* [1,4] */;
    uint32_t FIBELEVEN_TSLFME[4 EXT_STRUCT_CTRL(FIBELEVEN_TSLFM)] /* [1,4] */;
    uint32_t FIBELEVEN_TSLF /* [0,32768] */;
    uint32_t FIBELEVEN_TSLFv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TSLF)] /* [0,32768] */;
    uint32_t FIBELEVEN_TSTFM /* [1,4] */;
    uint32_t FIBELEVEN_TSTFMI[4 EXT_STRUCT_CTRL(FIBELEVEN_TSTFM)] /* [1,4] */;
    uint32_t FIBELEVEN_TSTFME[4 EXT_STRUCT_CTRL(FIBELEVEN_TSTFM)] /* [1,4] */;
    uint32_t FIBELEVEN_TSTF /* [0,32768] */;
    uint32_t FIBELEVEN_TSTFv[32768 EXT_STRUCT_CTRL(FIBELEVEN_TSTF)] /* [0,32768] */;

} EXT_STR_h101_FIBELEVEN;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBELEVEN_onion_t
{
    /* RAW */
    uint32_t FIBELEVEN_TMLCM;
    uint32_t FIBELEVEN_TMLCMI[512 /* FIBELEVEN_TMLCM */];
    uint32_t FIBELEVEN_TMLCME[512 /* FIBELEVEN_TMLCM */];
    uint32_t FIBELEVEN_TMLC;
    uint32_t FIBELEVEN_TMLCv[32768 /* FIBELEVEN_TMLC */];
    uint32_t FIBELEVEN_TMTCM;
    uint32_t FIBELEVEN_TMTCMI[512 /* FIBELEVEN_TMTCM */];
    uint32_t FIBELEVEN_TMTCME[512 /* FIBELEVEN_TMTCM */];
    uint32_t FIBELEVEN_TMTC;
    uint32_t FIBELEVEN_TMTCv[32768 /* FIBELEVEN_TMTC */];
    uint32_t FIBELEVEN_TMLFM;
    uint32_t FIBELEVEN_TMLFMI[512 /* FIBELEVEN_TMLFM */];
    uint32_t FIBELEVEN_TMLFME[512 /* FIBELEVEN_TMLFM */];
    uint32_t FIBELEVEN_TMLF;
    uint32_t FIBELEVEN_TMLFv[32768 /* FIBELEVEN_TMLF */];
    uint32_t FIBELEVEN_TMTFM;
    uint32_t FIBELEVEN_TMTFMI[512 /* FIBELEVEN_TMTFM */];
    uint32_t FIBELEVEN_TMTFME[512 /* FIBELEVEN_TMTFM */];
    uint32_t FIBELEVEN_TMTF;
    uint32_t FIBELEVEN_TMTFv[32768 /* FIBELEVEN_TMTF */];
    uint32_t FIBELEVEN_TRIGMLCM;
    uint32_t FIBELEVEN_TRIGMLCMI[4 /* FIBELEVEN_TRIGMLCM */];
    uint32_t FIBELEVEN_TRIGMLCME[4 /* FIBELEVEN_TRIGMLCM */];
    uint32_t FIBELEVEN_TRIGMLC;
    uint32_t FIBELEVEN_TRIGMLCv[32768 /* FIBELEVEN_TRIGMLC */];
    uint32_t FIBELEVEN_TRIGMLFM;
    uint32_t FIBELEVEN_TRIGMLFMI[4 /* FIBELEVEN_TRIGMLFM */];
    uint32_t FIBELEVEN_TRIGMLFME[4 /* FIBELEVEN_TRIGMLFM */];
    uint32_t FIBELEVEN_TRIGMLF;
    uint32_t FIBELEVEN_TRIGMLFv[32768 /* FIBELEVEN_TRIGMLF */];
    uint32_t FIBELEVEN_TSLCM;
    uint32_t FIBELEVEN_TSLCMI[4 /* FIBELEVEN_TSLCM */];
    uint32_t FIBELEVEN_TSLCME[4 /* FIBELEVEN_TSLCM */];
    uint32_t FIBELEVEN_TSLC;
    uint32_t FIBELEVEN_TSLCv[32768 /* FIBELEVEN_TSLC */];
    uint32_t FIBELEVEN_TSTCM;
    uint32_t FIBELEVEN_TSTCMI[4 /* FIBELEVEN_TSTCM */];
    uint32_t FIBELEVEN_TSTCME[4 /* FIBELEVEN_TSTCM */];
    uint32_t FIBELEVEN_TSTC;
    uint32_t FIBELEVEN_TSTCv[32768 /* FIBELEVEN_TSTC */];
    uint32_t FIBELEVEN_TSLFM;
    uint32_t FIBELEVEN_TSLFMI[4 /* FIBELEVEN_TSLFM */];
    uint32_t FIBELEVEN_TSLFME[4 /* FIBELEVEN_TSLFM */];
    uint32_t FIBELEVEN_TSLF;
    uint32_t FIBELEVEN_TSLFv[32768 /* FIBELEVEN_TSLF */];
    uint32_t FIBELEVEN_TSTFM;
    uint32_t FIBELEVEN_TSTFMI[4 /* FIBELEVEN_TSTFM */];
    uint32_t FIBELEVEN_TSTFME[4 /* FIBELEVEN_TSTFM */];
    uint32_t FIBELEVEN_TSTF;
    uint32_t FIBELEVEN_TSTFv[32768 /* FIBELEVEN_TSTF */];

} EXT_STR_h101_FIBELEVEN_onion;

/*******************************************************/

#define EXT_STR_h101_FIBELEVEN_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                       \
    do                                                                                                              \
    {                                                                                                               \
        ok = 1;                                                                                                     \
        /* RAW */                                                                                                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TMLCM, UINT32, "FIBELEVEN_TMLCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMLCMI, UINT32, "FIBELEVEN_TMLCMI", "FIBELEVEN_TMLCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMLCME, UINT32, "FIBELEVEN_TMLCME", "FIBELEVEN_TMLCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TMLC, UINT32, "FIBELEVEN_TMLC", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMLCv, UINT32, "FIBELEVEN_TMLCv", "FIBELEVEN_TMLC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TMTCM, UINT32, "FIBELEVEN_TMTCM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMTCMI, UINT32, "FIBELEVEN_TMTCMI", "FIBELEVEN_TMTCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMTCME, UINT32, "FIBELEVEN_TMTCME", "FIBELEVEN_TMTCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TMTC, UINT32, "FIBELEVEN_TMTC", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMTCv, UINT32, "FIBELEVEN_TMTCv", "FIBELEVEN_TMTC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TMLFM, UINT32, "FIBELEVEN_TMLFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMLFMI, UINT32, "FIBELEVEN_TMLFMI", "FIBELEVEN_TMLFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMLFME, UINT32, "FIBELEVEN_TMLFME", "FIBELEVEN_TMLFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TMLF, UINT32, "FIBELEVEN_TMLF", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMLFv, UINT32, "FIBELEVEN_TMLFv", "FIBELEVEN_TMLF");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TMTFM, UINT32, "FIBELEVEN_TMTFM", 512); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMTFMI, UINT32, "FIBELEVEN_TMTFMI", "FIBELEVEN_TMTFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMTFME, UINT32, "FIBELEVEN_TMTFME", "FIBELEVEN_TMTFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TMTF, UINT32, "FIBELEVEN_TMTF", 32768); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TMTFv, UINT32, "FIBELEVEN_TMTFv", "FIBELEVEN_TMTF");      \
        EXT_STR_ITEM_INFO_LIM(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TRIGMLCM, UINT32, "FIBELEVEN_TRIGMLCM", 4);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBELEVEN_TRIGMLCMI,                                                                  \
                              UINT32,                                                                               \
                              "FIBELEVEN_TRIGMLCMI",                                                                \
                              "FIBELEVEN_TRIGMLCM");                                                                \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBELEVEN_TRIGMLCME,                                                                  \
                              UINT32,                                                                               \
                              "FIBELEVEN_TRIGMLCME",                                                                \
                              "FIBELEVEN_TRIGMLCM");                                                                \
        EXT_STR_ITEM_INFO_LIM(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TRIGMLC, UINT32, "FIBELEVEN_TRIGMLC", 256);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBELEVEN_TRIGMLCv,                                                                   \
                              UINT32,                                                                               \
                              "FIBELEVEN_TRIGMLCv",                                                                 \
                              "FIBELEVEN_TRIGMLC");                                                                 \
        EXT_STR_ITEM_INFO_LIM(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TRIGMLFM, UINT32, "FIBELEVEN_TRIGMLFM", 4);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBELEVEN_TRIGMLFMI,                                                                  \
                              UINT32,                                                                               \
                              "FIBELEVEN_TRIGMLFMI",                                                                \
                              "FIBELEVEN_TRIGMLFM");                                                                \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBELEVEN_TRIGMLFME,                                                                  \
                              UINT32,                                                                               \
                              "FIBELEVEN_TRIGMLFME",                                                                \
                              "FIBELEVEN_TRIGMLFM");                                                                \
        EXT_STR_ITEM_INFO_LIM(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TRIGMLF, UINT32, "FIBELEVEN_TRIGMLF", 256);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                   \
                              si,                                                                                   \
                              offset,                                                                               \
                              struct_t,                                                                             \
                              printerr,                                                                             \
                              FIBELEVEN_TRIGMLFv,                                                                   \
                              UINT32,                                                                               \
                              "FIBELEVEN_TRIGMLFv",                                                                 \
                              "FIBELEVEN_TRIGMLF");                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TSLCM, UINT32, "FIBELEVEN_TSLCM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSLCMI, UINT32, "FIBELEVEN_TSLCMI", "FIBELEVEN_TSLCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSLCME, UINT32, "FIBELEVEN_TSLCME", "FIBELEVEN_TSLCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TSLC, UINT32, "FIBELEVEN_TSLC", 256);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSLCv, UINT32, "FIBELEVEN_TSLCv", "FIBELEVEN_TSLC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TSTCM, UINT32, "FIBELEVEN_TSTCM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSTCMI, UINT32, "FIBELEVEN_TSTCMI", "FIBELEVEN_TSTCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSTCME, UINT32, "FIBELEVEN_TSTCME", "FIBELEVEN_TSTCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TSTC, UINT32, "FIBELEVEN_TSTC", 256);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSTCv, UINT32, "FIBELEVEN_TSTCv", "FIBELEVEN_TSTC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TSLFM, UINT32, "FIBELEVEN_TSLFM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSLFMI, UINT32, "FIBELEVEN_TSLFMI", "FIBELEVEN_TSLFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSLFME, UINT32, "FIBELEVEN_TSLFME", "FIBELEVEN_TSLFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TSLF, UINT32, "FIBELEVEN_TSLF", 256);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSLFv, UINT32, "FIBELEVEN_TSLFv", "FIBELEVEN_TSLF");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TSTFM, UINT32, "FIBELEVEN_TSTFM", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSTFMI, UINT32, "FIBELEVEN_TSTFMI", "FIBELEVEN_TSTFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSTFME, UINT32, "FIBELEVEN_TSTFME", "FIBELEVEN_TSTFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBELEVEN_TSTF, UINT32, "FIBELEVEN_TSTF", 256);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, FIBELEVEN_TSTFv, UINT32, "FIBELEVEN_TSTFv", "FIBELEVEN_TSTF");      \
                                                                                                                    \
    } while (0);

#endif /*__GUARD_H101_FIBELEVEN_EXT_H101_FIBELEVEN_H__*/

/*******************************************************/
