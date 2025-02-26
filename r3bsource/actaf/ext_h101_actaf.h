/******************************************************************************
 *   Copyright (C) 2025 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2025 Members of R3B Collaboration                          *
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

#pragma once

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

typedef struct EXT_STR_h101_ACTAF_t
{
    /* RAW */
    uint32_t ACTAF1TSLO /* [0,20] */;
    uint32_t ACTAF1TSLOv[20 EXT_STRUCT_CTRL(ACTAF1TSLO)] /* [-1,-1] */;
    uint32_t ACTAF1TSHI /* [0,20] */;
    uint32_t ACTAF1TSHIv[20 EXT_STRUCT_CTRL(ACTAF1TSHI)] /* [-1,-1] */;
    uint32_t ACTAF1ENERGY /* [0,20] */;
    uint32_t ACTAF1ENERGYv[20 EXT_STRUCT_CTRL(ACTAF1ENERGY)] /* [-1,-1] */;
    uint32_t ACTAF1BASELINE /* [0,20] */;
    uint32_t ACTAF1BASELINEv[20 EXT_STRUCT_CTRL(ACTAF1BASELINE)] /* [-1,-1] */;
    uint32_t ACTAF1MAW1 /* [0,20] */;
    uint32_t ACTAF1MAW1v[20 EXT_STRUCT_CTRL(ACTAF1MAW1)] /* [-1,-1] */;
    uint32_t ACTAF1MAW2 /* [0,20] */;
    uint32_t ACTAF1MAW2v[20 EXT_STRUCT_CTRL(ACTAF1MAW2)] /* [-1,-1] */;
    uint32_t ACTAF1MAW3 /* [0,20] */;
    uint32_t ACTAF1MAW3v[20 EXT_STRUCT_CTRL(ACTAF1MAW3)] /* [-1,-1] */;
    uint32_t ACTAF1CH /* [0,20] */;
    uint32_t ACTAF1CHv[20 EXT_STRUCT_CTRL(ACTAF1CH)] /* [0,255] */;
    uint32_t ACTAF1TRACECH /* [0,60000] */;
    uint32_t ACTAF1TRACECHv[60000 EXT_STRUCT_CTRL(ACTAF1TRACECH)] /* [0,255] */;
    uint32_t ACTAF1TRACERAW /* [0,60000] */;
    uint32_t ACTAF1TRACERAWv[60000 EXT_STRUCT_CTRL(ACTAF1TRACERAW)] /* [0,65535] */;
    uint32_t ACTAF1TRACECHMAW /* [0,20440] */;
    uint32_t ACTAF1TRACECHMAWv[20440 EXT_STRUCT_CTRL(ACTAF1TRACECHMAW)] /* [0,255] */;
    uint32_t ACTAF1TRACEMAW /* [0,20440] */;
    uint32_t ACTAF1TRACEMAWv[20440 EXT_STRUCT_CTRL(ACTAF1TRACEMAW)] /* [-1,-1] */;
    uint32_t ACTAF2TSLO /* [0,20] */;
    uint32_t ACTAF2TSLOv[20 EXT_STRUCT_CTRL(ACTAF2TSLO)] /* [-1,-1] */;
    uint32_t ACTAF2TSHI /* [0,20] */;
    uint32_t ACTAF2TSHIv[20 EXT_STRUCT_CTRL(ACTAF2TSHI)] /* [-1,-1] */;
    uint32_t ACTAF2ENERGY /* [0,20] */;
    uint32_t ACTAF2ENERGYv[20 EXT_STRUCT_CTRL(ACTAF2ENERGY)] /* [-1,-1] */;
    uint32_t ACTAF2BASELINE /* [0,20] */;
    uint32_t ACTAF2BASELINEv[20 EXT_STRUCT_CTRL(ACTAF2BASELINE)] /* [-1,-1] */;
    uint32_t ACTAF2MAW1 /* [0,20] */;
    uint32_t ACTAF2MAW1v[20 EXT_STRUCT_CTRL(ACTAF2MAW1)] /* [-1,-1] */;
    uint32_t ACTAF2MAW2 /* [0,20] */;
    uint32_t ACTAF2MAW2v[20 EXT_STRUCT_CTRL(ACTAF2MAW2)] /* [-1,-1] */;
    uint32_t ACTAF2MAW3 /* [0,20] */;
    uint32_t ACTAF2MAW3v[20 EXT_STRUCT_CTRL(ACTAF2MAW3)] /* [-1,-1] */;
    uint32_t ACTAF2CH /* [0,20] */;
    uint32_t ACTAF2CHv[20 EXT_STRUCT_CTRL(ACTAF2CH)] /* [0,255] */;
    uint32_t ACTAF2TRACECH /* [0,60000] */;
    uint32_t ACTAF2TRACECHv[60000 EXT_STRUCT_CTRL(ACTAF2TRACECH)] /* [0,255] */;
    uint32_t ACTAF2TRACERAW /* [0,60000] */;
    uint32_t ACTAF2TRACERAWv[60000 EXT_STRUCT_CTRL(ACTAF2TRACERAW)] /* [0,65535] */;
    uint32_t ACTAF2TRACECHMAW /* [0,20440] */;
    uint32_t ACTAF2TRACECHMAWv[20440 EXT_STRUCT_CTRL(ACTAF2TRACECHMAW)] /* [0,255] */;
    uint32_t ACTAF2TRACEMAW /* [0,20440] */;
    uint32_t ACTAF2TRACEMAWv[20440 EXT_STRUCT_CTRL(ACTAF2TRACEMAW)] /* [-1,-1] */;
    uint32_t ACTAF3TSLO /* [0,20] */;
    uint32_t ACTAF3TSLOv[20 EXT_STRUCT_CTRL(ACTAF3TSLO)] /* [-1,-1] */;
    uint32_t ACTAF3TSHI /* [0,20] */;
    uint32_t ACTAF3TSHIv[20 EXT_STRUCT_CTRL(ACTAF3TSHI)] /* [-1,-1] */;
    uint32_t ACTAF3ENERGY /* [0,20] */;
    uint32_t ACTAF3ENERGYv[20 EXT_STRUCT_CTRL(ACTAF3ENERGY)] /* [-1,-1] */;
    uint32_t ACTAF3BASELINE /* [0,20] */;
    uint32_t ACTAF3BASELINEv[20 EXT_STRUCT_CTRL(ACTAF3BASELINE)] /* [-1,-1] */;
    uint32_t ACTAF3MAW1 /* [0,20] */;
    uint32_t ACTAF3MAW1v[20 EXT_STRUCT_CTRL(ACTAF3MAW1)] /* [-1,-1] */;
    uint32_t ACTAF3MAW2 /* [0,20] */;
    uint32_t ACTAF3MAW2v[20 EXT_STRUCT_CTRL(ACTAF3MAW2)] /* [-1,-1] */;
    uint32_t ACTAF3MAW3 /* [0,20] */;
    uint32_t ACTAF3MAW3v[20 EXT_STRUCT_CTRL(ACTAF3MAW3)] /* [-1,-1] */;
    uint32_t ACTAF3CH /* [0,20] */;
    uint32_t ACTAF3CHv[20 EXT_STRUCT_CTRL(ACTAF3CH)] /* [0,255] */;
    uint32_t ACTAF3TRACECH /* [0,60000] */;
    uint32_t ACTAF3TRACECHv[60000 EXT_STRUCT_CTRL(ACTAF3TRACECH)] /* [0,255] */;
    uint32_t ACTAF3TRACERAW /* [0,60000] */;
    uint32_t ACTAF3TRACERAWv[60000 EXT_STRUCT_CTRL(ACTAF3TRACERAW)] /* [0,65535] */;
    uint32_t ACTAF3TRACECHMAW /* [0,20440] */;
    uint32_t ACTAF3TRACECHMAWv[20440 EXT_STRUCT_CTRL(ACTAF3TRACECHMAW)] /* [0,255] */;
    uint32_t ACTAF3TRACEMAW /* [0,20440] */;
    uint32_t ACTAF3TRACEMAWv[20440 EXT_STRUCT_CTRL(ACTAF3TRACEMAW)] /* [-1,-1] */;
    uint32_t ACTAF4TSLO /* [0,20] */;
    uint32_t ACTAF4TSLOv[20 EXT_STRUCT_CTRL(ACTAF4TSLO)] /* [-1,-1] */;
    uint32_t ACTAF4TSHI /* [0,20] */;
    uint32_t ACTAF4TSHIv[20 EXT_STRUCT_CTRL(ACTAF4TSHI)] /* [-1,-1] */;
    uint32_t ACTAF4ENERGY /* [0,20] */;
    uint32_t ACTAF4ENERGYv[20 EXT_STRUCT_CTRL(ACTAF4ENERGY)] /* [-1,-1] */;
    uint32_t ACTAF4BASELINE /* [0,20] */;
    uint32_t ACTAF4BASELINEv[20 EXT_STRUCT_CTRL(ACTAF4BASELINE)] /* [-1,-1] */;
    uint32_t ACTAF4MAW1 /* [0,20] */;
    uint32_t ACTAF4MAW1v[20 EXT_STRUCT_CTRL(ACTAF4MAW1)] /* [-1,-1] */;
    uint32_t ACTAF4MAW2 /* [0,20] */;
    uint32_t ACTAF4MAW2v[20 EXT_STRUCT_CTRL(ACTAF4MAW2)] /* [-1,-1] */;
    uint32_t ACTAF4MAW3 /* [0,20] */;
    uint32_t ACTAF4MAW3v[20 EXT_STRUCT_CTRL(ACTAF4MAW3)] /* [-1,-1] */;
    uint32_t ACTAF4CH /* [0,20] */;
    uint32_t ACTAF4CHv[20 EXT_STRUCT_CTRL(ACTAF4CH)] /* [0,255] */;
    uint32_t ACTAF4TRACECH /* [0,60000] */;
    uint32_t ACTAF4TRACECHv[60000 EXT_STRUCT_CTRL(ACTAF4TRACECH)] /* [0,255] */;
    uint32_t ACTAF4TRACERAW /* [0,60000] */;
    uint32_t ACTAF4TRACERAWv[60000 EXT_STRUCT_CTRL(ACTAF4TRACERAW)] /* [0,65535] */;
    uint32_t ACTAF4TRACECHMAW /* [0,20440] */;
    uint32_t ACTAF4TRACECHMAWv[20440 EXT_STRUCT_CTRL(ACTAF4TRACECHMAW)] /* [0,255] */;
    uint32_t ACTAF4TRACEMAW /* [0,20440] */;
    uint32_t ACTAF4TRACEMAWv[20440 EXT_STRUCT_CTRL(ACTAF4TRACEMAW)] /* [-1,-1] */;

} EXT_STR_h101_ACTAF;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_ACTAF_onion_t
{
    /* RAW */
    struct
    {
        uint32_t TSLO;
        uint32_t TSLOv[20 /* TSLO */];
        uint32_t TSHI;
        uint32_t TSHIv[20 /* TSHI */];
        uint32_t ENERGY;
        uint32_t ENERGYv[20 /* ENERGY */];
        uint32_t BASELINE;
        uint32_t BASELINEv[20 /* BASELINE */];
        struct
        {
            uint32_t _;
            uint32_t v[20 /*  */];
        } MAW[3];
        uint32_t CH;
        uint32_t CHv[20 /* CH */];
        uint32_t TRACECH;
        uint32_t TRACECHv[60000 /* TRACECH */];
        uint32_t TRACERAW;
        uint32_t TRACERAWv[60000 /* TRACERAW */];
        uint32_t TRACECHMAW;
        uint32_t TRACECHMAWv[20440 /* TRACECHMAW */];
        uint32_t TRACEMAW;
        uint32_t TRACEMAWv[20440 /* TRACEMAW */];
    } ACTAF[4];

} EXT_STR_h101_ACTAF_onion;

/*******************************************************/

#define EXT_STR_h101_ACTAF_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                              \
    do                                                                                                                 \
    {                                                                                                                  \
        ok = 1;                                                                                                        \
        /* RAW */                                                                                                      \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF1TSLO, UINT32, "ACTAF1TSLO", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1TSLOv, UINT32, "ACTAF1TSLOv", "ACTAF1TSLO", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF1TSHI, UINT32, "ACTAF1TSHI", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1TSHIv, UINT32, "ACTAF1TSHIv", "ACTAF1TSHI", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1ENERGY, UINT32, "ACTAF1ENERGY", 20, 0 /*flags*/);                \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1ENERGYv, UINT32, "ACTAF1ENERGYv", "ACTAF1ENERGY", 0 /*flags*/);  \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1BASELINE, UINT32, "ACTAF1BASELINE", 20, 0 /*flags*/);            \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF1BASELINEv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF1BASELINEv",                                                                      \
                               "ACTAF1BASELINE",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF1MAW1, UINT32, "ACTAF1MAW1", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1MAW1v, UINT32, "ACTAF1MAW1v", "ACTAF1MAW1", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF1MAW2, UINT32, "ACTAF1MAW2", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1MAW2v, UINT32, "ACTAF1MAW2v", "ACTAF1MAW2", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF1MAW3, UINT32, "ACTAF1MAW3", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1MAW3v, UINT32, "ACTAF1MAW3v", "ACTAF1MAW3", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF1CH, UINT32, "ACTAF1CH", 20, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1CHv, UINT32, "ACTAF1CHv", "ACTAF1CH", 0 /*flags*/);              \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1TRACECH, UINT32, "ACTAF1TRACECH", 60000, 0 /*flags*/);           \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF1TRACECHv,                                                                         \
                               UINT32,                                                                                 \
                               "ACTAF1TRACECHv",                                                                       \
                               "ACTAF1TRACECH",                                                                        \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1TRACERAW, UINT32, "ACTAF1TRACERAW", 60000, 0 /*flags*/);         \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF1TRACERAWv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF1TRACERAWv",                                                                      \
                               "ACTAF1TRACERAW",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1TRACECHMAW, UINT32, "ACTAF1TRACECHMAW", 20440, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF1TRACECHMAWv,                                                                      \
                               UINT32,                                                                                 \
                               "ACTAF1TRACECHMAWv",                                                                    \
                               "ACTAF1TRACECHMAW",                                                                     \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF1TRACEMAW, UINT32, "ACTAF1TRACEMAW", 20440, 0 /*flags*/);         \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF1TRACEMAWv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF1TRACEMAWv",                                                                      \
                               "ACTAF1TRACEMAW",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF2TSLO, UINT32, "ACTAF2TSLO", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2TSLOv, UINT32, "ACTAF2TSLOv", "ACTAF2TSLO", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF2TSHI, UINT32, "ACTAF2TSHI", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2TSHIv, UINT32, "ACTAF2TSHIv", "ACTAF2TSHI", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2ENERGY, UINT32, "ACTAF2ENERGY", 20, 0 /*flags*/);                \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2ENERGYv, UINT32, "ACTAF2ENERGYv", "ACTAF2ENERGY", 0 /*flags*/);  \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2BASELINE, UINT32, "ACTAF2BASELINE", 20, 0 /*flags*/);            \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF2BASELINEv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF2BASELINEv",                                                                      \
                               "ACTAF2BASELINE",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF2MAW1, UINT32, "ACTAF2MAW1", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2MAW1v, UINT32, "ACTAF2MAW1v", "ACTAF2MAW1", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF2MAW2, UINT32, "ACTAF2MAW2", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2MAW2v, UINT32, "ACTAF2MAW2v", "ACTAF2MAW2", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF2MAW3, UINT32, "ACTAF2MAW3", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2MAW3v, UINT32, "ACTAF2MAW3v", "ACTAF2MAW3", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF2CH, UINT32, "ACTAF2CH", 20, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2CHv, UINT32, "ACTAF2CHv", "ACTAF2CH", 0 /*flags*/);              \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2TRACECH, UINT32, "ACTAF2TRACECH", 60000, 0 /*flags*/);           \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF2TRACECHv,                                                                         \
                               UINT32,                                                                                 \
                               "ACTAF2TRACECHv",                                                                       \
                               "ACTAF2TRACECH",                                                                        \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2TRACERAW, UINT32, "ACTAF2TRACERAW", 60000, 0 /*flags*/);         \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF2TRACERAWv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF2TRACERAWv",                                                                      \
                               "ACTAF2TRACERAW",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2TRACECHMAW, UINT32, "ACTAF2TRACECHMAW", 20440, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF2TRACECHMAWv,                                                                      \
                               UINT32,                                                                                 \
                               "ACTAF2TRACECHMAWv",                                                                    \
                               "ACTAF2TRACECHMAW",                                                                     \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF2TRACEMAW, UINT32, "ACTAF2TRACEMAW", 20440, 0 /*flags*/);         \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF2TRACEMAWv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF2TRACEMAWv",                                                                      \
                               "ACTAF2TRACEMAW",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF3TSLO, UINT32, "ACTAF3TSLO", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3TSLOv, UINT32, "ACTAF3TSLOv", "ACTAF3TSLO", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF3TSHI, UINT32, "ACTAF3TSHI", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3TSHIv, UINT32, "ACTAF3TSHIv", "ACTAF3TSHI", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3ENERGY, UINT32, "ACTAF3ENERGY", 20, 0 /*flags*/);                \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3ENERGYv, UINT32, "ACTAF3ENERGYv", "ACTAF3ENERGY", 0 /*flags*/);  \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3BASELINE, UINT32, "ACTAF3BASELINE", 20, 0 /*flags*/);            \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF3BASELINEv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF3BASELINEv",                                                                      \
                               "ACTAF3BASELINE",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF3MAW1, UINT32, "ACTAF3MAW1", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3MAW1v, UINT32, "ACTAF3MAW1v", "ACTAF3MAW1", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF3MAW2, UINT32, "ACTAF3MAW2", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3MAW2v, UINT32, "ACTAF3MAW2v", "ACTAF3MAW2", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF3MAW3, UINT32, "ACTAF3MAW3", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3MAW3v, UINT32, "ACTAF3MAW3v", "ACTAF3MAW3", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF3CH, UINT32, "ACTAF3CH", 20, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3CHv, UINT32, "ACTAF3CHv", "ACTAF3CH", 0 /*flags*/);              \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3TRACECH, UINT32, "ACTAF3TRACECH", 60000, 0 /*flags*/);           \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF3TRACECHv,                                                                         \
                               UINT32,                                                                                 \
                               "ACTAF3TRACECHv",                                                                       \
                               "ACTAF3TRACECH",                                                                        \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3TRACERAW, UINT32, "ACTAF3TRACERAW", 60000, 0 /*flags*/);         \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF3TRACERAWv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF3TRACERAWv",                                                                      \
                               "ACTAF3TRACERAW",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3TRACECHMAW, UINT32, "ACTAF3TRACECHMAW", 20440, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF3TRACECHMAWv,                                                                      \
                               UINT32,                                                                                 \
                               "ACTAF3TRACECHMAWv",                                                                    \
                               "ACTAF3TRACECHMAW",                                                                     \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF3TRACEMAW, UINT32, "ACTAF3TRACEMAW", 20440, 0 /*flags*/);         \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF3TRACEMAWv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF3TRACEMAWv",                                                                      \
                               "ACTAF3TRACEMAW",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF4TSLO, UINT32, "ACTAF4TSLO", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4TSLOv, UINT32, "ACTAF4TSLOv", "ACTAF4TSLO", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF4TSHI, UINT32, "ACTAF4TSHI", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4TSHIv, UINT32, "ACTAF4TSHIv", "ACTAF4TSHI", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4ENERGY, UINT32, "ACTAF4ENERGY", 20, 0 /*flags*/);                \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4ENERGYv, UINT32, "ACTAF4ENERGYv", "ACTAF4ENERGY", 0 /*flags*/);  \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4BASELINE, UINT32, "ACTAF4BASELINE", 20, 0 /*flags*/);            \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF4BASELINEv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF4BASELINEv",                                                                      \
                               "ACTAF4BASELINE",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF4MAW1, UINT32, "ACTAF4MAW1", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4MAW1v, UINT32, "ACTAF4MAW1v", "ACTAF4MAW1", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF4MAW2, UINT32, "ACTAF4MAW2", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4MAW2v, UINT32, "ACTAF4MAW2v", "ACTAF4MAW2", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF4MAW3, UINT32, "ACTAF4MAW3", 20, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4MAW3v, UINT32, "ACTAF4MAW3v", "ACTAF4MAW3", 0 /*flags*/);        \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ACTAF4CH, UINT32, "ACTAF4CH", 20, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4CHv, UINT32, "ACTAF4CHv", "ACTAF4CH", 0 /*flags*/);              \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4TRACECH, UINT32, "ACTAF4TRACECH", 60000, 0 /*flags*/);           \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF4TRACECHv,                                                                         \
                               UINT32,                                                                                 \
                               "ACTAF4TRACECHv",                                                                       \
                               "ACTAF4TRACECH",                                                                        \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4TRACERAW, UINT32, "ACTAF4TRACERAW", 60000, 0 /*flags*/);         \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF4TRACERAWv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF4TRACERAWv",                                                                      \
                               "ACTAF4TRACERAW",                                                                       \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4TRACECHMAW, UINT32, "ACTAF4TRACECHMAW", 20440, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF4TRACECHMAWv,                                                                      \
                               UINT32,                                                                                 \
                               "ACTAF4TRACECHMAWv",                                                                    \
                               "ACTAF4TRACECHMAW",                                                                     \
                               0 /*flags*/);                                                                           \
        EXT_STR_ITEM_INFO2_LIM(                                                                                        \
            ok, si, offset, struct_t, printerr, ACTAF4TRACEMAW, UINT32, "ACTAF4TRACEMAW", 20440, 0 /*flags*/);         \
        EXT_STR_ITEM_INFO2_ZZP(ok,                                                                                     \
                               si,                                                                                     \
                               offset,                                                                                 \
                               struct_t,                                                                               \
                               printerr,                                                                               \
                               ACTAF4TRACEMAWv,                                                                        \
                               UINT32,                                                                                 \
                               "ACTAF4TRACEMAWv",                                                                      \
                               "ACTAF4TRACEMAW",                                                                       \
                               0 /*flags*/);                                                                           \
                                                                                                                       \
    } while (0);

/*******************************************************/
