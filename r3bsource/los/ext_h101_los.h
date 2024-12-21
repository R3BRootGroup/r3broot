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

typedef struct EXT_STR_h101_LOS_t
{
    /* RAW */
    uint32_t LOS1VTFM /* [1,8] */;
    uint32_t LOS1VTFMI[8 EXT_STRUCT_CTRL(LOS1VTFM)] /* [1,8] */;
    uint32_t LOS1VTFME[8 EXT_STRUCT_CTRL(LOS1VTFM)] /* [1,256] */;
    uint32_t LOS1VTF /* [0,256] */;
    uint32_t LOS1VTFv[256 EXT_STRUCT_CTRL(LOS1VTF)] /* [0,65535] */;
    uint32_t LOS1VTCM /* [1,8] */;
    uint32_t LOS1VTCMI[8 EXT_STRUCT_CTRL(LOS1VTCM)] /* [1,8] */;
    uint32_t LOS1VTCME[8 EXT_STRUCT_CTRL(LOS1VTCM)] /* [1,256] */;
    uint32_t LOS1VTC /* [0,256] */;
    uint32_t LOS1VTCv[256 EXT_STRUCT_CTRL(LOS1VTC)] /* [0,65535] */;
    uint32_t LOS1VTRIGF /* [0,1] */;
    uint32_t LOS1VTRIGFI[1 EXT_STRUCT_CTRL(LOS1VTRIGF)] /* [1,1] */;
    uint32_t LOS1VTRIGFv[1 EXT_STRUCT_CTRL(LOS1VTRIGF)] /* [0,65535] */;
    uint32_t LOS1VTRIGC /* [0,1] */;
    uint32_t LOS1VTRIGCI[1 EXT_STRUCT_CTRL(LOS1VTRIGC)] /* [1,1] */;
    uint32_t LOS1VTRIGCv[1 EXT_STRUCT_CTRL(LOS1VTRIGC)] /* [0,65535] */;
    uint32_t LOS1TTFLM /* [1,8] */;
    uint32_t LOS1TTFLMI[8 EXT_STRUCT_CTRL(LOS1TTFLM)] /* [1,8] */;
    uint32_t LOS1TTFLME[8 EXT_STRUCT_CTRL(LOS1TTFLM)] /* [1,256] */;
    uint32_t LOS1TTFL /* [0,256] */;
    uint32_t LOS1TTFLv[256 EXT_STRUCT_CTRL(LOS1TTFL)] /* [0,65535] */;
    uint32_t LOS1TTFTM /* [1,8] */;
    uint32_t LOS1TTFTMI[8 EXT_STRUCT_CTRL(LOS1TTFTM)] /* [1,8] */;
    uint32_t LOS1TTFTME[8 EXT_STRUCT_CTRL(LOS1TTFTM)] /* [1,256] */;
    uint32_t LOS1TTFT /* [0,256] */;
    uint32_t LOS1TTFTv[256 EXT_STRUCT_CTRL(LOS1TTFT)] /* [0,65535] */;
    uint32_t LOS1TTCLM /* [1,8] */;
    uint32_t LOS1TTCLMI[8 EXT_STRUCT_CTRL(LOS1TTCLM)] /* [1,8] */;
    uint32_t LOS1TTCLME[8 EXT_STRUCT_CTRL(LOS1TTCLM)] /* [1,256] */;
    uint32_t LOS1TTCL /* [0,256] */;
    uint32_t LOS1TTCLv[256 EXT_STRUCT_CTRL(LOS1TTCL)] /* [0,65535] */;
    uint32_t LOS1TTCTM /* [1,8] */;
    uint32_t LOS1TTCTMI[8 EXT_STRUCT_CTRL(LOS1TTCTM)] /* [1,8] */;
    uint32_t LOS1TTCTME[8 EXT_STRUCT_CTRL(LOS1TTCTM)] /* [1,256] */;
    uint32_t LOS1TTCT /* [0,256] */;
    uint32_t LOS1TTCTv[256 EXT_STRUCT_CTRL(LOS1TTCT)] /* [0,65535] */;
    uint32_t LOS1TTRIGFL /* [0,1] */;
    uint32_t LOS1TTRIGFLI[1 EXT_STRUCT_CTRL(LOS1TTRIGFL)] /* [1,1] */;
    uint32_t LOS1TTRIGFLv[1 EXT_STRUCT_CTRL(LOS1TTRIGFL)] /* [0,65535] */;
    uint32_t LOS1TTRIGCL /* [0,1] */;
    uint32_t LOS1TTRIGCLI[1 EXT_STRUCT_CTRL(LOS1TTRIGCL)] /* [1,1] */;
    uint32_t LOS1TTRIGCLv[1 EXT_STRUCT_CTRL(LOS1TTRIGCL)] /* [0,65535] */;
    uint32_t LOS1TTRIGFT /* [0,1] */;
    uint32_t LOS1TTRIGFTI[1 EXT_STRUCT_CTRL(LOS1TTRIGFT)] /* [1,1] */;
    uint32_t LOS1TTRIGFTv[1 EXT_STRUCT_CTRL(LOS1TTRIGFT)] /* [0,65535] */;
    uint32_t LOS1TTRIGCT /* [0,1] */;
    uint32_t LOS1TTRIGCTI[1 EXT_STRUCT_CTRL(LOS1TTRIGCT)] /* [1,1] */;
    uint32_t LOS1TTRIGCTv[1 EXT_STRUCT_CTRL(LOS1TTRIGCT)] /* [0,65535] */;

} EXT_STR_h101_LOS;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_LOS_onion_t
{
    /* RAW */
    struct
    {
        uint32_t VTFM;
        uint32_t VTFMI[8 /* VTFM */];
        uint32_t VTFME[8 /* VTFM */];
        uint32_t VTF;
        uint32_t VTFv[256 /* VTF */];
        uint32_t VTCM;
        uint32_t VTCMI[8 /* VTCM */];
        uint32_t VTCME[8 /* VTCM */];
        uint32_t VTC;
        uint32_t VTCv[256 /* VTC */];
        uint32_t VTRIGF;
        uint32_t VTRIGFI[1 /* VTRIGF */];
        uint32_t VTRIGFv[1 /* VTRIGF */];
        uint32_t VTRIGC;
        uint32_t VTRIGCI[1 /* VTRIGC */];
        uint32_t VTRIGCv[1 /* VTRIGC */];
        uint32_t TTFLM;
        uint32_t TTFLMI[8 /* TTFLM */];
        uint32_t TTFLME[8 /* TTFLM */];
        uint32_t TTFL;
        uint32_t TTFLv[256 /* TTFL */];
        uint32_t TTFTM;
        uint32_t TTFTMI[8 /* TTFTM */];
        uint32_t TTFTME[8 /* TTFTM */];
        uint32_t TTFT;
        uint32_t TTFTv[256 /* TTFT */];
        uint32_t TTCLM;
        uint32_t TTCLMI[8 /* TTCLM */];
        uint32_t TTCLME[8 /* TTCLM */];
        uint32_t TTCL;
        uint32_t TTCLv[256 /* TTCL */];
        uint32_t TTCTM;
        uint32_t TTCTMI[8 /* TTCTM */];
        uint32_t TTCTME[8 /* TTCTM */];
        uint32_t TTCT;
        uint32_t TTCTv[256 /* TTCT */];
        uint32_t TTRIGFL;
        uint32_t TTRIGFLI[1 /* TTRIGFL */];
        uint32_t TTRIGFLv[1 /* TTRIGFL */];
        uint32_t TTRIGCL;
        uint32_t TTRIGCLI[1 /* TTRIGCL */];
        uint32_t TTRIGCLv[1 /* TTRIGCL */];
        uint32_t TTRIGFT;
        uint32_t TTRIGFTI[1 /* TTRIGFT */];
        uint32_t TTRIGFTv[1 /* TTRIGFT */];
        uint32_t TTRIGCT;
        uint32_t TTRIGCTI[1 /* TTRIGCT */];
        uint32_t TTRIGCTv[1 /* TTRIGCT */];
    } LOS[1];

} EXT_STR_h101_LOS_onion;

/*******************************************************/

#define EXT_STR_h101_LOS_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                               \
    do                                                                                                                \
    {                                                                                                                 \
        ok = 1;                                                                                                       \
        /* RAW */                                                                                                     \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1VTFM, UINT32, "LOS1VTFM", 8, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTFMI, UINT32, "LOS1VTFMI", "LOS1VTFM", 0 /*flags*/);             \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTFME, UINT32, "LOS1VTFME", "LOS1VTFM", 0 /*flags*/);             \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1VTF, UINT32, "LOS1VTF", 256, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTFv, UINT32, "LOS1VTFv", "LOS1VTF", 0 /*flags*/);                \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1VTCM, UINT32, "LOS1VTCM", 8, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTCMI, UINT32, "LOS1VTCMI", "LOS1VTCM", 0 /*flags*/);             \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTCME, UINT32, "LOS1VTCME", "LOS1VTCM", 0 /*flags*/);             \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1VTC, UINT32, "LOS1VTC", 256, 0 /*flags*/);     \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTCv, UINT32, "LOS1VTCv", "LOS1VTC", 0 /*flags*/);                \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1VTRIGF, UINT32, "LOS1VTRIGF", 1, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTRIGFI, UINT32, "LOS1VTRIGFI", "LOS1VTRIGF", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTRIGFv, UINT32, "LOS1VTRIGFv", "LOS1VTRIGF", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1VTRIGC, UINT32, "LOS1VTRIGC", 1, 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTRIGCI, UINT32, "LOS1VTRIGCI", "LOS1VTRIGC", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1VTRIGCv, UINT32, "LOS1VTRIGCv", "LOS1VTRIGC", 0 /*flags*/);       \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1TTFLM, UINT32, "LOS1TTFLM", 8, 0 /*flags*/);   \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTFLMI, UINT32, "LOS1TTFLMI", "LOS1TTFLM", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTFLME, UINT32, "LOS1TTFLME", "LOS1TTFLM", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1TTFL, UINT32, "LOS1TTFL", 256, 0 /*flags*/);   \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTFLv, UINT32, "LOS1TTFLv", "LOS1TTFL", 0 /*flags*/);             \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1TTFTM, UINT32, "LOS1TTFTM", 8, 0 /*flags*/);   \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTFTMI, UINT32, "LOS1TTFTMI", "LOS1TTFTM", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTFTME, UINT32, "LOS1TTFTME", "LOS1TTFTM", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1TTFT, UINT32, "LOS1TTFT", 256, 0 /*flags*/);   \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTFTv, UINT32, "LOS1TTFTv", "LOS1TTFT", 0 /*flags*/);             \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1TTCLM, UINT32, "LOS1TTCLM", 8, 0 /*flags*/);   \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTCLMI, UINT32, "LOS1TTCLMI", "LOS1TTCLM", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTCLME, UINT32, "LOS1TTCLME", "LOS1TTCLM", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1TTCL, UINT32, "LOS1TTCL", 256, 0 /*flags*/);   \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTCLv, UINT32, "LOS1TTCLv", "LOS1TTCL", 0 /*flags*/);             \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1TTCTM, UINT32, "LOS1TTCTM", 8, 0 /*flags*/);   \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTCTMI, UINT32, "LOS1TTCTMI", "LOS1TTCTM", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTCTME, UINT32, "LOS1TTCTME", "LOS1TTCTM", 0 /*flags*/);          \
        EXT_STR_ITEM_INFO2_LIM(ok, si, offset, struct_t, printerr, LOS1TTCT, UINT32, "LOS1TTCT", 256, 0 /*flags*/);   \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTCTv, UINT32, "LOS1TTCTv", "LOS1TTCT", 0 /*flags*/);             \
        EXT_STR_ITEM_INFO2_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGFL, UINT32, "LOS1TTRIGFL", 1, 0 /*flags*/);                  \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGFLI, UINT32, "LOS1TTRIGFLI", "LOS1TTRIGFL", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGFLv, UINT32, "LOS1TTRIGFLv", "LOS1TTRIGFL", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGCL, UINT32, "LOS1TTRIGCL", 1, 0 /*flags*/);                  \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGCLI, UINT32, "LOS1TTRIGCLI", "LOS1TTRIGCL", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGCLv, UINT32, "LOS1TTRIGCLv", "LOS1TTRIGCL", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGFT, UINT32, "LOS1TTRIGFT", 1, 0 /*flags*/);                  \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGFTI, UINT32, "LOS1TTRIGFTI", "LOS1TTRIGFT", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGFTv, UINT32, "LOS1TTRIGFTv", "LOS1TTRIGFT", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGCT, UINT32, "LOS1TTRIGCT", 1, 0 /*flags*/);                  \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGCTI, UINT32, "LOS1TTRIGCTI", "LOS1TTRIGCT", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, LOS1TTRIGCTv, UINT32, "LOS1TTRIGCTv", "LOS1TTRIGCT", 0 /*flags*/);    \
                                                                                                                      \
    } while (0);

/*******************************************************/
