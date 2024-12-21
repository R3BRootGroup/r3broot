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

#ifndef __GUARD_H101_FIBTEN_EXT_H101_FIBTEN_H__
#define __GUARD_H101_FIBTEN_EXT_H101_FIBTEN_H__

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

typedef struct EXT_STR_h101_FIBTEN_t
{
    /* RAW */
    uint32_t FIBTEN_TMLCM /* [1,512] */;
    uint32_t FIBTEN_TMLCMI[512 EXT_STRUCT_CTRL(FIBTEN_TMLCM)] /* [1,512] */;
    uint32_t FIBTEN_TMLCME[512 EXT_STRUCT_CTRL(FIBTEN_TMLCM)] /* [1,512] */;
    uint32_t FIBTEN_TMLC /* [0,32768] */;
    uint32_t FIBTEN_TMLCv[32768 EXT_STRUCT_CTRL(FIBTEN_TMLC)] /* [0,32768] */;
    uint32_t FIBTEN_TMTCM /* [1,512] */;
    uint32_t FIBTEN_TMTCMI[512 EXT_STRUCT_CTRL(FIBTEN_TMTCM)] /* [1,512] */;
    uint32_t FIBTEN_TMTCME[512 EXT_STRUCT_CTRL(FIBTEN_TMTCM)] /* [1,512] */;
    uint32_t FIBTEN_TMTC /* [0,32768] */;
    uint32_t FIBTEN_TMTCv[32768 EXT_STRUCT_CTRL(FIBTEN_TMTC)] /* [0,32768] */;
    uint32_t FIBTEN_TMLFM /* [1,512] */;
    uint32_t FIBTEN_TMLFMI[512 EXT_STRUCT_CTRL(FIBTEN_TMLFM)] /* [1,512] */;
    uint32_t FIBTEN_TMLFME[512 EXT_STRUCT_CTRL(FIBTEN_TMLFM)] /* [1,512] */;
    uint32_t FIBTEN_TMLF /* [0,32768] */;
    uint32_t FIBTEN_TMLFv[32768 EXT_STRUCT_CTRL(FIBTEN_TMLF)] /* [0,32768] */;
    uint32_t FIBTEN_TMTFM /* [1,512] */;
    uint32_t FIBTEN_TMTFMI[512 EXT_STRUCT_CTRL(FIBTEN_TMTFM)] /* [1,512] */;
    uint32_t FIBTEN_TMTFME[512 EXT_STRUCT_CTRL(FIBTEN_TMTFM)] /* [1,512] */;
    uint32_t FIBTEN_TMTF /* [0,32768] */;
    uint32_t FIBTEN_TMTFv[32768 EXT_STRUCT_CTRL(FIBTEN_TMTF)] /* [0,32768] */;
    uint32_t FIBTEN_TRIGMLCM /* [1,4] */;
    uint32_t FIBTEN_TRIGMLCMI[4 EXT_STRUCT_CTRL(FIBTEN_TRIGMLCM)] /* [1,4] */;
    uint32_t FIBTEN_TRIGMLCME[4 EXT_STRUCT_CTRL(FIBTEN_TRIGMLCM)] /* [1,4] */;
    uint32_t FIBTEN_TRIGMLC /* [0,32768] */;
    uint32_t FIBTEN_TRIGMLCv[32768 EXT_STRUCT_CTRL(FIBTEN_TRIGMLC)] /* [0,32767] */;
    uint32_t FIBTEN_TRIGMLFM /* [1,4] */;
    uint32_t FIBTEN_TRIGMLFMI[4 EXT_STRUCT_CTRL(FIBTEN_TRIGMLFM)] /* [1,4] */;
    uint32_t FIBTEN_TRIGMLFME[4 EXT_STRUCT_CTRL(FIBTEN_TRIGMLFM)] /* [1,4] */;
    uint32_t FIBTEN_TRIGMLF /* [0,32768] */;
    uint32_t FIBTEN_TRIGMLFv[32768 EXT_STRUCT_CTRL(FIBTEN_TRIGMLF)] /* [0,32767] */;
    uint32_t FIBTEN_TSLCM /* [1,4] */;
    uint32_t FIBTEN_TSLCMI[4 EXT_STRUCT_CTRL(FIBTEN_TSLCM)] /* [1,4] */;
    uint32_t FIBTEN_TSLCME[4 EXT_STRUCT_CTRL(FIBTEN_TSLCM)] /* [1,4] */;
    uint32_t FIBTEN_TSLC /* [0,32768] */;
    uint32_t FIBTEN_TSLCv[32768 EXT_STRUCT_CTRL(FIBTEN_TSLC)] /* [0,32768] */;
    uint32_t FIBTEN_TSTCM /* [1,4] */;
    uint32_t FIBTEN_TSTCMI[4 EXT_STRUCT_CTRL(FIBTEN_TSTCM)] /* [1,4] */;
    uint32_t FIBTEN_TSTCME[4 EXT_STRUCT_CTRL(FIBTEN_TSTCM)] /* [1,4] */;
    uint32_t FIBTEN_TSTC /* [0,32768] */;
    uint32_t FIBTEN_TSTCv[32768 EXT_STRUCT_CTRL(FIBTEN_TSTC)] /* [0,32768] */;
    uint32_t FIBTEN_TSLFM /* [1,4] */;
    uint32_t FIBTEN_TSLFMI[4 EXT_STRUCT_CTRL(FIBTEN_TSLFM)] /* [1,4] */;
    uint32_t FIBTEN_TSLFME[4 EXT_STRUCT_CTRL(FIBTEN_TSLFM)] /* [1,4] */;
    uint32_t FIBTEN_TSLF /* [0,32768] */;
    uint32_t FIBTEN_TSLFv[32768 EXT_STRUCT_CTRL(FIBTEN_TSLF)] /* [0,32768] */;
    uint32_t FIBTEN_TSTFM /* [1,4] */;
    uint32_t FIBTEN_TSTFMI[4 EXT_STRUCT_CTRL(FIBTEN_TSTFM)] /* [1,4] */;
    uint32_t FIBTEN_TSTFME[4 EXT_STRUCT_CTRL(FIBTEN_TSTFM)] /* [1,4] */;
    uint32_t FIBTEN_TSTF /* [0,32768] */;
    uint32_t FIBTEN_TSTFv[32768 EXT_STRUCT_CTRL(FIBTEN_TSTF)] /* [0,32768] */;

} EXT_STR_h101_FIBTEN;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBTEN_onion_t
{
    /* RAW */
    uint32_t FIBTEN_TMLCM;
    uint32_t FIBTEN_TMLCMI[512 /* FIBTEN_TMLCM */];
    uint32_t FIBTEN_TMLCME[512 /* FIBTEN_TMLCM */];
    uint32_t FIBTEN_TMLC;
    uint32_t FIBTEN_TMLCv[32768 /* FIBTEN_TMLC */];
    uint32_t FIBTEN_TMTCM;
    uint32_t FIBTEN_TMTCMI[512 /* FIBTEN_TMTCM */];
    uint32_t FIBTEN_TMTCME[512 /* FIBTEN_TMTCM */];
    uint32_t FIBTEN_TMTC;
    uint32_t FIBTEN_TMTCv[32768 /* FIBTEN_TMTC */];
    uint32_t FIBTEN_TMLFM;
    uint32_t FIBTEN_TMLFMI[512 /* FIBTEN_TMLFM */];
    uint32_t FIBTEN_TMLFME[512 /* FIBTEN_TMLFM */];
    uint32_t FIBTEN_TMLF;
    uint32_t FIBTEN_TMLFv[32768 /* FIBTEN_TMLF */];
    uint32_t FIBTEN_TMTFM;
    uint32_t FIBTEN_TMTFMI[512 /* FIBTEN_TMTFM */];
    uint32_t FIBTEN_TMTFME[512 /* FIBTEN_TMTFM */];
    uint32_t FIBTEN_TMTF;
    uint32_t FIBTEN_TMTFv[32768 /* FIBTEN_TMTF */];
    uint32_t FIBTEN_TRIGMLCM;
    uint32_t FIBTEN_TRIGMLCMI[4 /* FIBTEN_TRIGMLCM */];
    uint32_t FIBTEN_TRIGMLCME[4 /* FIBTEN_TRIGMLCM */];
    uint32_t FIBTEN_TRIGMLC;
    uint32_t FIBTEN_TRIGMLCv[32768 /* FIBTEN_TRIGMLC */];
    uint32_t FIBTEN_TRIGMLFM;
    uint32_t FIBTEN_TRIGMLFMI[4 /* FIBTEN_TRIGMLFM */];
    uint32_t FIBTEN_TRIGMLFME[4 /* FIBTEN_TRIGMLFM */];
    uint32_t FIBTEN_TRIGMLF;
    uint32_t FIBTEN_TRIGMLFv[32768 /* FIBTEN_TRIGMLF */];
    uint32_t FIBTEN_TSLCM;
    uint32_t FIBTEN_TSLCMI[4 /* FIBTEN_TSLCM */];
    uint32_t FIBTEN_TSLCME[4 /* FIBTEN_TSLCM */];
    uint32_t FIBTEN_TSLC;
    uint32_t FIBTEN_TSLCv[32768 /* FIBTEN_TSLC */];
    uint32_t FIBTEN_TSTCM;
    uint32_t FIBTEN_TSTCMI[4 /* FIBTEN_TSTCM */];
    uint32_t FIBTEN_TSTCME[4 /* FIBTEN_TSTCM */];
    uint32_t FIBTEN_TSTC;
    uint32_t FIBTEN_TSTCv[32768 /* FIBTEN_TSTC */];
    uint32_t FIBTEN_TSLFM;
    uint32_t FIBTEN_TSLFMI[4 /* FIBTEN_TSLFM */];
    uint32_t FIBTEN_TSLFME[4 /* FIBTEN_TSLFM */];
    uint32_t FIBTEN_TSLF;
    uint32_t FIBTEN_TSLFv[32768 /* FIBTEN_TSLF */];
    uint32_t FIBTEN_TSTFM;
    uint32_t FIBTEN_TSTFMI[4 /* FIBTEN_TSTFM */];
    uint32_t FIBTEN_TSTFME[4 /* FIBTEN_TSTFM */];
    uint32_t FIBTEN_TSTF;
    uint32_t FIBTEN_TSTFv[32768 /* FIBTEN_TSTF */];

} EXT_STR_h101_FIBTEN_onion;

/*******************************************************/

#define EXT_STR_h101_FIBTEN_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                        \
    do                                                                                                            \
    {                                                                                                             \
        ok = 1;                                                                                                   \
        /* RAW */                                                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TMLCM, UINT32, "FIBTEN_TMLCM", 512);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMLCMI, UINT32, "FIBTEN_TMLCMI", "FIBTEN_TMLCM");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMLCME, UINT32, "FIBTEN_TMLCME", "FIBTEN_TMLCM");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TMLC, UINT32, "FIBTEN_TMLC", 32768);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMLCv, UINT32, "FIBTEN_TMLCv", "FIBTEN_TMLC");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TMTCM, UINT32, "FIBTEN_TMTCM", 512);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMTCMI, UINT32, "FIBTEN_TMTCMI", "FIBTEN_TMTCM");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMTCME, UINT32, "FIBTEN_TMTCME", "FIBTEN_TMTCM");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TMTC, UINT32, "FIBTEN_TMTC", 32768);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMTCv, UINT32, "FIBTEN_TMTCv", "FIBTEN_TMTC");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TMLFM, UINT32, "FIBTEN_TMLFM", 512);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMLFMI, UINT32, "FIBTEN_TMLFMI", "FIBTEN_TMLFM");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMLFME, UINT32, "FIBTEN_TMLFME", "FIBTEN_TMLFM");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TMLF, UINT32, "FIBTEN_TMLF", 32768);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMLFv, UINT32, "FIBTEN_TMLFv", "FIBTEN_TMLF");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TMTFM, UINT32, "FIBTEN_TMTFM", 512);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMTFMI, UINT32, "FIBTEN_TMTFMI", "FIBTEN_TMTFM");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMTFME, UINT32, "FIBTEN_TMTFME", "FIBTEN_TMTFM");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TMTF, UINT32, "FIBTEN_TMTF", 32768);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TMTFv, UINT32, "FIBTEN_TMTFv", "FIBTEN_TMTF");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLCM, UINT32, "FIBTEN_TRIGMLCM", 4); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLCMI, UINT32, "FIBTEN_TRIGMLCMI", "FIBTEN_TRIGMLCM"); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLCME, UINT32, "FIBTEN_TRIGMLCME", "FIBTEN_TRIGMLCM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLC, UINT32, "FIBTEN_TRIGMLC", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLCv, UINT32, "FIBTEN_TRIGMLCv", "FIBTEN_TRIGMLC");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLFM, UINT32, "FIBTEN_TRIGMLFM", 4); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLFMI, UINT32, "FIBTEN_TRIGMLFMI", "FIBTEN_TRIGMLFM"); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLFME, UINT32, "FIBTEN_TRIGMLFME", "FIBTEN_TRIGMLFM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLF, UINT32, "FIBTEN_TRIGMLF", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TRIGMLFv, UINT32, "FIBTEN_TRIGMLFv", "FIBTEN_TRIGMLF");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TSLCM, UINT32, "FIBTEN_TSLCM", 4);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSLCMI, UINT32, "FIBTEN_TSLCMI", "FIBTEN_TSLCM");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSLCME, UINT32, "FIBTEN_TSLCME", "FIBTEN_TSLCM");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TSLC, UINT32, "FIBTEN_TSLC", 256);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSLCv, UINT32, "FIBTEN_TSLCv", "FIBTEN_TSLC");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TSTCM, UINT32, "FIBTEN_TSTCM", 4);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSTCMI, UINT32, "FIBTEN_TSTCMI", "FIBTEN_TSTCM");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSTCME, UINT32, "FIBTEN_TSTCME", "FIBTEN_TSTCM");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TSTC, UINT32, "FIBTEN_TSTC", 256);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSTCv, UINT32, "FIBTEN_TSTCv", "FIBTEN_TSTC");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TSLFM, UINT32, "FIBTEN_TSLFM", 4);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSLFMI, UINT32, "FIBTEN_TSLFMI", "FIBTEN_TSLFM");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSLFME, UINT32, "FIBTEN_TSLFME", "FIBTEN_TSLFM");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TSLF, UINT32, "FIBTEN_TSLF", 256);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSLFv, UINT32, "FIBTEN_TSLFv", "FIBTEN_TSLF");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TSTFM, UINT32, "FIBTEN_TSTFM", 4);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSTFMI, UINT32, "FIBTEN_TSTFMI", "FIBTEN_TSTFM");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSTFME, UINT32, "FIBTEN_TSTFME", "FIBTEN_TSTFM");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBTEN_TSTF, UINT32, "FIBTEN_TSTF", 256);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FIBTEN_TSTFv, UINT32, "FIBTEN_TSTFv", "FIBTEN_TSTF");             \
                                                                                                                  \
    } while (0);

#endif /*__GUARD_H101_FIBTEN_EXT_H101_FIBTEN_H__*/

/*******************************************************/
