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

typedef struct EXT_STR_h101_ROLU_NOTRIG_t
{
    /* RAW */
    uint32_t ROLU1TTFLM /* [1,4] */;
    uint32_t ROLU1TTFLMI[4 EXT_STRUCT_CTRL(ROLU1TTFLM)] /* [1,4] */;
    uint32_t ROLU1TTFLME[4 EXT_STRUCT_CTRL(ROLU1TTFLM)] /* [1,128] */;
    uint32_t ROLU1TTFL /* [0,128] */;
    uint32_t ROLU1TTFLv[128 EXT_STRUCT_CTRL(ROLU1TTFL)] /* [0,65535] */;
    uint32_t ROLU1TTFTM /* [1,4] */;
    uint32_t ROLU1TTFTMI[4 EXT_STRUCT_CTRL(ROLU1TTFTM)] /* [1,4] */;
    uint32_t ROLU1TTFTME[4 EXT_STRUCT_CTRL(ROLU1TTFTM)] /* [1,128] */;
    uint32_t ROLU1TTFT /* [0,128] */;
    uint32_t ROLU1TTFTv[128 EXT_STRUCT_CTRL(ROLU1TTFT)] /* [0,65535] */;
    uint32_t ROLU1TTCLM /* [1,4] */;
    uint32_t ROLU1TTCLMI[4 EXT_STRUCT_CTRL(ROLU1TTCLM)] /* [1,4] */;
    uint32_t ROLU1TTCLME[4 EXT_STRUCT_CTRL(ROLU1TTCLM)] /* [1,128] */;
    uint32_t ROLU1TTCL /* [0,128] */;
    uint32_t ROLU1TTCLv[128 EXT_STRUCT_CTRL(ROLU1TTCL)] /* [0,65535] */;
    uint32_t ROLU1TTCTM /* [1,4] */;
    uint32_t ROLU1TTCTMI[4 EXT_STRUCT_CTRL(ROLU1TTCTM)] /* [1,4] */;
    uint32_t ROLU1TTCTME[4 EXT_STRUCT_CTRL(ROLU1TTCTM)] /* [1,128] */;
    uint32_t ROLU1TTCT /* [0,128] */;
    uint32_t ROLU1TTCTv[128 EXT_STRUCT_CTRL(ROLU1TTCT)] /* [0,65535] */;
    uint32_t ROLU2TTFLM /* [1,4] */;
    uint32_t ROLU2TTFLMI[4 EXT_STRUCT_CTRL(ROLU2TTFLM)] /* [1,4] */;
    uint32_t ROLU2TTFLME[4 EXT_STRUCT_CTRL(ROLU2TTFLM)] /* [1,128] */;
    uint32_t ROLU2TTFL /* [0,128] */;
    uint32_t ROLU2TTFLv[128 EXT_STRUCT_CTRL(ROLU2TTFL)] /* [0,65535] */;
    uint32_t ROLU2TTFTM /* [1,4] */;
    uint32_t ROLU2TTFTMI[4 EXT_STRUCT_CTRL(ROLU2TTFTM)] /* [1,4] */;
    uint32_t ROLU2TTFTME[4 EXT_STRUCT_CTRL(ROLU2TTFTM)] /* [1,128] */;
    uint32_t ROLU2TTFT /* [0,128] */;
    uint32_t ROLU2TTFTv[128 EXT_STRUCT_CTRL(ROLU2TTFT)] /* [0,65535] */;
    uint32_t ROLU2TTCLM /* [1,4] */;
    uint32_t ROLU2TTCLMI[4 EXT_STRUCT_CTRL(ROLU2TTCLM)] /* [1,4] */;
    uint32_t ROLU2TTCLME[4 EXT_STRUCT_CTRL(ROLU2TTCLM)] /* [1,128] */;
    uint32_t ROLU2TTCL /* [0,128] */;
    uint32_t ROLU2TTCLv[128 EXT_STRUCT_CTRL(ROLU2TTCL)] /* [0,65535] */;
    uint32_t ROLU2TTCTM /* [1,4] */;
    uint32_t ROLU2TTCTMI[4 EXT_STRUCT_CTRL(ROLU2TTCTM)] /* [1,4] */;
    uint32_t ROLU2TTCTME[4 EXT_STRUCT_CTRL(ROLU2TTCTM)] /* [1,128] */;
    uint32_t ROLU2TTCT /* [0,128] */;
    uint32_t ROLU2TTCTv[128 EXT_STRUCT_CTRL(ROLU2TTCT)] /* [0,65535] */;

} EXT_STR_h101_ROLU_NOTRIG;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_ROLU_NOTRIG_onion_t
{
    /* RAW */
    struct
    {
        uint32_t TTFLM;
        uint32_t TTFLMI[4 /* TTFLM */];
        uint32_t TTFLME[4 /* TTFLM */];
        uint32_t TTFL;
        uint32_t TTFLv[128 /* TTFL */];
        uint32_t TTFTM;
        uint32_t TTFTMI[4 /* TTFTM */];
        uint32_t TTFTME[4 /* TTFTM */];
        uint32_t TTFT;
        uint32_t TTFTv[128 /* TTFT */];
        uint32_t TTCLM;
        uint32_t TTCLMI[4 /* TTCLM */];
        uint32_t TTCLME[4 /* TTCLM */];
        uint32_t TTCL;
        uint32_t TTCLv[128 /* TTCL */];
        uint32_t TTCTM;
        uint32_t TTCTMI[4 /* TTCTM */];
        uint32_t TTCTME[4 /* TTCTM */];
        uint32_t TTCT;
        uint32_t TTCTv[128 /* TTCT */];
    } ROLU[2];

} EXT_STR_h101_ROLU_NOTRIG_onion;

/*******************************************************/

#define EXT_STR_h101_ROLU_NOTRIG_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                       \
    do                                                                                                                \
    {                                                                                                                 \
        ok = 1;                                                                                                       \
        /* RAW */                                                                                                     \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFLM, UINT32, "ROLU1TTFLM", 4, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTFLMI, UINT32, "ROLU1TTFLMI", "ROLU1TTFLM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTFLME, UINT32, "ROLU1TTFLME", "ROLU1TTFLM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFL, UINT32, "ROLU1TTFL", 128, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTFLv, UINT32, "ROLU1TTFLv", "ROLU1TTFL", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFTM, UINT32, "ROLU1TTFTM", 4, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTFTMI, UINT32, "ROLU1TTFTMI", "ROLU1TTFTM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTFTME, UINT32, "ROLU1TTFTME", "ROLU1TTFTM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFT, UINT32, "ROLU1TTFT", 128, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTFTv, UINT32, "ROLU1TTFTv", "ROLU1TTFT", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCLM, UINT32, "ROLU1TTCLM", 4, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTCLMI, UINT32, "ROLU1TTCLMI", "ROLU1TTCLM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTCLME, UINT32, "ROLU1TTCLME", "ROLU1TTCLM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCL, UINT32, "ROLU1TTCL", 128, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTCLv, UINT32, "ROLU1TTCLv", "ROLU1TTCL", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCTM, UINT32, "ROLU1TTCTM", 4, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTCTMI, UINT32, "ROLU1TTCTMI", "ROLU1TTCTM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTCTME, UINT32, "ROLU1TTCTME", "ROLU1TTCTM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCT, UINT32, "ROLU1TTCT", 128, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TTCTv, UINT32, "ROLU1TTCTv", "ROLU1TTCT", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU2TTFLM, UINT32, "ROLU2TTFLM", 4, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTFLMI, UINT32, "ROLU2TTFLMI", "ROLU2TTFLM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTFLME, UINT32, "ROLU2TTFLME", "ROLU2TTFLM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU2TTFL, UINT32, "ROLU2TTFL", 128, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTFLv, UINT32, "ROLU2TTFLv", "ROLU2TTFL", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU2TTFTM, UINT32, "ROLU2TTFTM", 4, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTFTMI, UINT32, "ROLU2TTFTMI", "ROLU2TTFTM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTFTME, UINT32, "ROLU2TTFTME", "ROLU2TTFTM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU2TTFT, UINT32, "ROLU2TTFT", 128, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTFTv, UINT32, "ROLU2TTFTv", "ROLU2TTFT", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU2TTCLM, UINT32, "ROLU2TTCLM", 4, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTCLMI, UINT32, "ROLU2TTCLMI", "ROLU2TTCLM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTCLME, UINT32, "ROLU2TTCLME", "ROLU2TTCLM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU2TTCL, UINT32, "ROLU2TTCL", 128, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTCLv, UINT32, "ROLU2TTCLv", "ROLU2TTCL", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU2TTCTM, UINT32, "ROLU2TTCTM", 4, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTCTMI, UINT32, "ROLU2TTCTMI", "ROLU2TTCTM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTCTME, UINT32, "ROLU2TTCTME", "ROLU2TTCTM", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, ROLU2TTCT, UINT32, "ROLU2TTCT", 128, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TTCTv, UINT32, "ROLU2TTCTv", "ROLU2TTCT", 0 /*flags*/);          \
                                                                                                                      \
    } while (0);

/*******************************************************/
