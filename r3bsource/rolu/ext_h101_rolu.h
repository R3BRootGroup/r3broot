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

#ifndef __GUARD_H101_ROLU_EXT_H101_ROLU_H__
#define __GUARD_H101_ROLU_EXT_H101_ROLU_H__

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

typedef struct EXT_STR_h101_ROLU_t
{
    /* RAW */
    uint32_t ROLU1TTFLM /* [1,4] */;
    uint32_t ROLU1TTFLMI[4 EXT_STRUCT_CTRL(ROLU1TTFLM)] /* [1,4] */;
    uint32_t ROLU1TTFLME[4 EXT_STRUCT_CTRL(ROLU1TTFLM)] /* [1,256] */;
    uint32_t ROLU1TTFL /* [0,256] */;
    uint32_t ROLU1TTFLv[256 EXT_STRUCT_CTRL(ROLU1TTFL)] /* [0,65535] */;
    uint32_t ROLU1TTFTM /* [1,4] */;
    uint32_t ROLU1TTFTMI[4 EXT_STRUCT_CTRL(ROLU1TTFTM)] /* [1,4] */;
    uint32_t ROLU1TTFTME[4 EXT_STRUCT_CTRL(ROLU1TTFTM)] /* [1,256] */;
    uint32_t ROLU1TTFT /* [0,256] */;
    uint32_t ROLU1TTFTv[256 EXT_STRUCT_CTRL(ROLU1TTFT)] /* [0,65535] */;
    uint32_t ROLU1TTCLM /* [1,4] */;
    uint32_t ROLU1TTCLMI[4 EXT_STRUCT_CTRL(ROLU1TTCLM)] /* [1,4] */;
    uint32_t ROLU1TTCLME[4 EXT_STRUCT_CTRL(ROLU1TTCLM)] /* [1,256] */;
    uint32_t ROLU1TTCL /* [0,256] */;
    uint32_t ROLU1TTCLv[256 EXT_STRUCT_CTRL(ROLU1TTCL)] /* [0,65535] */;
    uint32_t ROLU1TTCTM /* [1,4] */;
    uint32_t ROLU1TTCTMI[4 EXT_STRUCT_CTRL(ROLU1TTCTM)] /* [1,4] */;
    uint32_t ROLU1TTCTME[4 EXT_STRUCT_CTRL(ROLU1TTCTM)] /* [1,256] */;
    uint32_t ROLU1TTCT /* [0,256] */;
    uint32_t ROLU1TTCTv[256 EXT_STRUCT_CTRL(ROLU1TTCT)] /* [0,65535] */;
    uint32_t ROLU1TRIGCL /* [0,1] */;
    uint32_t ROLU1TRIGCLI[1 EXT_STRUCT_CTRL(ROLU1TRIGCL)] /* [1,1] */;
    uint32_t ROLU1TRIGCLv[1 EXT_STRUCT_CTRL(ROLU1TRIGCL)] /* [0,65535] */;
    uint32_t ROLU1TRIGFL /* [0,1] */;
    uint32_t ROLU1TRIGFLI[1 EXT_STRUCT_CTRL(ROLU1TRIGFL)] /* [1,1] */;
    uint32_t ROLU1TRIGFLv[1 EXT_STRUCT_CTRL(ROLU1TRIGFL)] /* [0,65535] */;
    uint32_t ROLU2TTFLM /* [1,4] */;
    uint32_t ROLU2TTFLMI[4 EXT_STRUCT_CTRL(ROLU2TTFLM)] /* [1,4] */;
    uint32_t ROLU2TTFLME[4 EXT_STRUCT_CTRL(ROLU2TTFLM)] /* [1,256] */;
    uint32_t ROLU2TTFL /* [0,256] */;
    uint32_t ROLU2TTFLv[256 EXT_STRUCT_CTRL(ROLU2TTFL)] /* [0,65535] */;
    uint32_t ROLU2TTFTM /* [1,4] */;
    uint32_t ROLU2TTFTMI[4 EXT_STRUCT_CTRL(ROLU2TTFTM)] /* [1,4] */;
    uint32_t ROLU2TTFTME[4 EXT_STRUCT_CTRL(ROLU2TTFTM)] /* [1,256] */;
    uint32_t ROLU2TTFT /* [0,256] */;
    uint32_t ROLU2TTFTv[256 EXT_STRUCT_CTRL(ROLU2TTFT)] /* [0,65535] */;
    uint32_t ROLU2TTCLM /* [1,4] */;
    uint32_t ROLU2TTCLMI[4 EXT_STRUCT_CTRL(ROLU2TTCLM)] /* [1,4] */;
    uint32_t ROLU2TTCLME[4 EXT_STRUCT_CTRL(ROLU2TTCLM)] /* [1,256] */;
    uint32_t ROLU2TTCL /* [0,256] */;
    uint32_t ROLU2TTCLv[256 EXT_STRUCT_CTRL(ROLU2TTCL)] /* [0,65535] */;
    uint32_t ROLU2TTCTM /* [1,4] */;
    uint32_t ROLU2TTCTMI[4 EXT_STRUCT_CTRL(ROLU2TTCTM)] /* [1,4] */;
    uint32_t ROLU2TTCTME[4 EXT_STRUCT_CTRL(ROLU2TTCTM)] /* [1,256] */;
    uint32_t ROLU2TTCT /* [0,256] */;
    uint32_t ROLU2TTCTv[256 EXT_STRUCT_CTRL(ROLU2TTCT)] /* [0,65535] */;
    uint32_t ROLU2TRIGCL /* [0,1] */;
    uint32_t ROLU2TRIGCLI[1 EXT_STRUCT_CTRL(ROLU2TRIGCL)] /* [1,1] */;
    uint32_t ROLU2TRIGCLv[1 EXT_STRUCT_CTRL(ROLU2TRIGCL)] /* [0,65535] */;
    uint32_t ROLU2TRIGFL /* [0,1] */;
    uint32_t ROLU2TRIGFLI[1 EXT_STRUCT_CTRL(ROLU2TRIGFL)] /* [1,1] */;
    uint32_t ROLU2TRIGFLv[1 EXT_STRUCT_CTRL(ROLU2TRIGFL)] /* [0,65535] */;

} EXT_STR_h101_ROLU;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_ROLU_onion_t
{
    /* RAW */
    struct
    {
        uint32_t TTFLM;
        uint32_t TTFLMI[4 /* TTFLM */];
        uint32_t TTFLME[4 /* TTFLM */];
        uint32_t TTFL;
        uint32_t TTFLv[256 /* TTFL */];
        uint32_t TTFTM;
        uint32_t TTFTMI[4 /* TTFTM */];
        uint32_t TTFTME[4 /* TTFTM */];
        uint32_t TTFT;
        uint32_t TTFTv[256 /* TTFT */];
        uint32_t TTCLM;
        uint32_t TTCLMI[4 /* TTCLM */];
        uint32_t TTCLME[4 /* TTCLM */];
        uint32_t TTCL;
        uint32_t TTCLv[256 /* TTCL */];
        uint32_t TTCTM;
        uint32_t TTCTMI[4 /* TTCTM */];
        uint32_t TTCTME[4 /* TTCTM */];
        uint32_t TTCT;
        uint32_t TTCTv[256 /* TTCT */];
        uint32_t TRIGCL;
        uint32_t TRIGCLI[1 /* TRIGCL */];
        uint32_t TRIGCLv[1 /* TRIGCL */];
        uint32_t TRIGFL;
        uint32_t TRIGFLI[1 /* TRIGFL */];
        uint32_t TRIGFLv[1 /* TRIGFL */];
    } ROLU[2];

} EXT_STR_h101_ROLU_onion;

/*******************************************************/

#define EXT_STR_h101_ROLU_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                             \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFLM, UINT32, "ROLU1TTFLM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFLMI, UINT32, "ROLU1TTFLMI", "ROLU1TTFLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFLME, UINT32, "ROLU1TTFLME", "ROLU1TTFLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFL, UINT32, "ROLU1TTFL", 256);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFLv, UINT32, "ROLU1TTFLv", "ROLU1TTFL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFTM, UINT32, "ROLU1TTFTM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFTMI, UINT32, "ROLU1TTFTMI", "ROLU1TTFTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFTME, UINT32, "ROLU1TTFTME", "ROLU1TTFTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFT, UINT32, "ROLU1TTFT", 256);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFTv, UINT32, "ROLU1TTFTv", "ROLU1TTFT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCLM, UINT32, "ROLU1TTCLM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCLMI, UINT32, "ROLU1TTCLMI", "ROLU1TTCLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCLME, UINT32, "ROLU1TTCLME", "ROLU1TTCLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCL, UINT32, "ROLU1TTCL", 256);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCLv, UINT32, "ROLU1TTCLv", "ROLU1TTCL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCTM, UINT32, "ROLU1TTCTM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCTMI, UINT32, "ROLU1TTCTMI", "ROLU1TTCTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCTME, UINT32, "ROLU1TTCTME", "ROLU1TTCTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCT, UINT32, "ROLU1TTCT", 256);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCTv, UINT32, "ROLU1TTCTv", "ROLU1TTCT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TRIGCL, UINT32, "ROLU1TRIGCL", 1);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TRIGCLI, UINT32, "ROLU1TRIGCLI", "ROLU1TRIGCL");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TRIGCLv, UINT32, "ROLU1TRIGCLv", "ROLU1TRIGCL");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TRIGFL, UINT32, "ROLU1TRIGFL", 1);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TRIGFLI, UINT32, "ROLU1TRIGFLI", "ROLU1TRIGFL");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU1TRIGFLv, UINT32, "ROLU1TRIGFLv", "ROLU1TRIGFL");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TTFLM, UINT32, "ROLU2TTFLM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTFLMI, UINT32, "ROLU2TTFLMI", "ROLU2TTFLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTFLME, UINT32, "ROLU2TTFLME", "ROLU2TTFLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TTFL, UINT32, "ROLU2TTFL", 256);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTFLv, UINT32, "ROLU2TTFLv", "ROLU2TTFL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TTFTM, UINT32, "ROLU2TTFTM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTFTMI, UINT32, "ROLU2TTFTMI", "ROLU2TTFTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTFTME, UINT32, "ROLU2TTFTME", "ROLU2TTFTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TTFT, UINT32, "ROLU2TTFT", 256);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTFTv, UINT32, "ROLU2TTFTv", "ROLU2TTFT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TTCLM, UINT32, "ROLU2TTCLM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTCLMI, UINT32, "ROLU2TTCLMI", "ROLU2TTCLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTCLME, UINT32, "ROLU2TTCLME", "ROLU2TTCLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TTCL, UINT32, "ROLU2TTCL", 256);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTCLv, UINT32, "ROLU2TTCLv", "ROLU2TTCL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TTCTM, UINT32, "ROLU2TTCTM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTCTMI, UINT32, "ROLU2TTCTMI", "ROLU2TTCTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTCTME, UINT32, "ROLU2TTCTME", "ROLU2TTCTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TTCT, UINT32, "ROLU2TTCT", 256);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU2TTCTv, UINT32, "ROLU2TTCTv", "ROLU2TTCT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TRIGCL, UINT32, "ROLU2TRIGCL", 1);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TRIGCLI, UINT32, "ROLU2TRIGCLI", "ROLU2TRIGCL");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TRIGCLv, UINT32, "ROLU2TRIGCLv", "ROLU2TRIGCL");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU2TRIGFL, UINT32, "ROLU2TRIGFL", 1);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TRIGFLI, UINT32, "ROLU2TRIGFLI", "ROLU2TRIGFL");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ROLU2TRIGFLv, UINT32, "ROLU2TRIGFLv", "ROLU2TRIGFL");                \
                                                                                                                     \
    } while (0);
#endif /*__GUARD_H101_ROLU_EXT_H101_ROLU_H__*/

/*******************************************************/
