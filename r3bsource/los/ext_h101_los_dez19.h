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

#ifndef __GUARD_H101_LOS_EXT_H101_LOS_DEZ19_H__
#define __GUARD_H101_LOS_EXT_H101_LOS_DEZ19_H__

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
    uint32_t LOS1MTM /* [1,8] */;
    uint32_t LOS1MTMI[8 EXT_STRUCT_CTRL(LOS1MTM)] /* [1,8] */;
    uint32_t LOS1MTME[8 EXT_STRUCT_CTRL(LOS1MTM)] /* [1,256] */;
    uint32_t LOS1MT /* [0,256] */;
    uint32_t LOS1MTv[256 EXT_STRUCT_CTRL(LOS1MT)] /* [0,65535] */;
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
    uint32_t LOS2VTFM /* [1,8] */;
    uint32_t LOS2VTFMI[8 EXT_STRUCT_CTRL(LOS2VTFM)] /* [1,8] */;
    uint32_t LOS2VTFME[8 EXT_STRUCT_CTRL(LOS2VTFM)] /* [1,256] */;
    uint32_t LOS2VTF /* [0,256] */;
    uint32_t LOS2VTFv[256 EXT_STRUCT_CTRL(LOS2VTF)] /* [0,65535] */;
    uint32_t LOS2VTCM /* [1,8] */;
    uint32_t LOS2VTCMI[8 EXT_STRUCT_CTRL(LOS2VTCM)] /* [1,8] */;
    uint32_t LOS2VTCME[8 EXT_STRUCT_CTRL(LOS2VTCM)] /* [1,256] */;
    uint32_t LOS2VTC /* [0,256] */;
    uint32_t LOS2VTCv[256 EXT_STRUCT_CTRL(LOS2VTC)] /* [0,65535] */;
    uint32_t LOS2MTM /* [1,8] */;
    uint32_t LOS2MTMI[8 EXT_STRUCT_CTRL(LOS2MTM)] /* [1,8] */;
    uint32_t LOS2MTME[8 EXT_STRUCT_CTRL(LOS2MTM)] /* [1,256] */;
    uint32_t LOS2MT /* [0,256] */;
    uint32_t LOS2MTv[256 EXT_STRUCT_CTRL(LOS2MT)] /* [0,65535] */;
    uint32_t LOS2TTFLM /* [1,8] */;
    uint32_t LOS2TTFLMI[8 EXT_STRUCT_CTRL(LOS2TTFLM)] /* [1,8] */;
    uint32_t LOS2TTFLME[8 EXT_STRUCT_CTRL(LOS2TTFLM)] /* [1,256] */;
    uint32_t LOS2TTFL /* [0,256] */;
    uint32_t LOS2TTFLv[256 EXT_STRUCT_CTRL(LOS2TTFL)] /* [0,65535] */;
    uint32_t LOS2TTFTM /* [1,8] */;
    uint32_t LOS2TTFTMI[8 EXT_STRUCT_CTRL(LOS2TTFTM)] /* [1,8] */;
    uint32_t LOS2TTFTME[8 EXT_STRUCT_CTRL(LOS2TTFTM)] /* [1,256] */;
    uint32_t LOS2TTFT /* [0,256] */;
    uint32_t LOS2TTFTv[256 EXT_STRUCT_CTRL(LOS2TTFT)] /* [0,65535] */;
    uint32_t LOS2TTCLM /* [1,8] */;
    uint32_t LOS2TTCLMI[8 EXT_STRUCT_CTRL(LOS2TTCLM)] /* [1,8] */;
    uint32_t LOS2TTCLME[8 EXT_STRUCT_CTRL(LOS2TTCLM)] /* [1,256] */;
    uint32_t LOS2TTCL /* [0,256] */;
    uint32_t LOS2TTCLv[256 EXT_STRUCT_CTRL(LOS2TTCL)] /* [0,65535] */;
    uint32_t LOS2TTCTM /* [1,8] */;
    uint32_t LOS2TTCTMI[8 EXT_STRUCT_CTRL(LOS2TTCTM)] /* [1,8] */;
    uint32_t LOS2TTCTME[8 EXT_STRUCT_CTRL(LOS2TTCTM)] /* [1,256] */;
    uint32_t LOS2TTCT /* [0,256] */;
    uint32_t LOS2TTCTv[256 EXT_STRUCT_CTRL(LOS2TTCT)] /* [0,65535] */;

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
        uint32_t MTM;
        uint32_t MTMI[8 /* MTM */];
        uint32_t MTME[8 /* MTM */];
        uint32_t MT;
        uint32_t MTv[256 /* MT */];
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
    } LOS[2];

} EXT_STR_h101_LOS_onion;

/*******************************************************/

#define EXT_STR_h101_LOS_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                           \
    do                                                                                                            \
    {                                                                                                             \
        ok = 1;                                                                                                   \
        /* RAW */                                                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1VTFM, UINT32, "LOS1VTFM", 8);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTFMI, UINT32, "LOS1VTFMI", "LOS1VTFM");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTFME, UINT32, "LOS1VTFME", "LOS1VTFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1VTF, UINT32, "LOS1VTF", 256);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTFv, UINT32, "LOS1VTFv", "LOS1VTF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1VTCM, UINT32, "LOS1VTCM", 8);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTCMI, UINT32, "LOS1VTCMI", "LOS1VTCM");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTCME, UINT32, "LOS1VTCME", "LOS1VTCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1VTC, UINT32, "LOS1VTC", 256);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTCv, UINT32, "LOS1VTCv", "LOS1VTC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1MTM, UINT32, "LOS1MTM", 8);                 \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1MTMI, UINT32, "LOS1MTMI", "LOS1MTM");       \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1MTME, UINT32, "LOS1MTME", "LOS1MTM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1MT, UINT32, "LOS1MT", 256);                 \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1MTv, UINT32, "LOS1MTv", "LOS1MT");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTFLM, UINT32, "LOS1TTFLM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFLMI, UINT32, "LOS1TTFLMI", "LOS1TTFLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFLME, UINT32, "LOS1TTFLME", "LOS1TTFLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTFL, UINT32, "LOS1TTFL", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFLv, UINT32, "LOS1TTFLv", "LOS1TTFL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTFTM, UINT32, "LOS1TTFTM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFTMI, UINT32, "LOS1TTFTMI", "LOS1TTFTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFTME, UINT32, "LOS1TTFTME", "LOS1TTFTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTFT, UINT32, "LOS1TTFT", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFTv, UINT32, "LOS1TTFTv", "LOS1TTFT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTCLM, UINT32, "LOS1TTCLM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCLMI, UINT32, "LOS1TTCLMI", "LOS1TTCLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCLME, UINT32, "LOS1TTCLME", "LOS1TTCLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTCL, UINT32, "LOS1TTCL", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCLv, UINT32, "LOS1TTCLv", "LOS1TTCL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTCTM, UINT32, "LOS1TTCTM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCTMI, UINT32, "LOS1TTCTMI", "LOS1TTCTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCTME, UINT32, "LOS1TTCTME", "LOS1TTCTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTCT, UINT32, "LOS1TTCT", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCTv, UINT32, "LOS1TTCTv", "LOS1TTCT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2VTFM, UINT32, "LOS2VTFM", 8);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2VTFMI, UINT32, "LOS2VTFMI", "LOS2VTFM");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2VTFME, UINT32, "LOS2VTFME", "LOS2VTFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2VTF, UINT32, "LOS2VTF", 256);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2VTFv, UINT32, "LOS2VTFv", "LOS2VTF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2VTCM, UINT32, "LOS2VTCM", 8);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2VTCMI, UINT32, "LOS2VTCMI", "LOS2VTCM");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2VTCME, UINT32, "LOS2VTCME", "LOS2VTCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2VTC, UINT32, "LOS2VTC", 256);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2VTCv, UINT32, "LOS2VTCv", "LOS2VTC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2MTM, UINT32, "LOS2MTM", 8);                 \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2MTMI, UINT32, "LOS2MTMI", "LOS2MTM");       \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2MTME, UINT32, "LOS2MTME", "LOS2MTM");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2MT, UINT32, "LOS2MT", 256);                 \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2MTv, UINT32, "LOS2MTv", "LOS2MT");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2TTFLM, UINT32, "LOS2TTFLM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTFLMI, UINT32, "LOS2TTFLMI", "LOS2TTFLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTFLME, UINT32, "LOS2TTFLME", "LOS2TTFLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2TTFL, UINT32, "LOS2TTFL", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTFLv, UINT32, "LOS2TTFLv", "LOS2TTFL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2TTFTM, UINT32, "LOS2TTFTM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTFTMI, UINT32, "LOS2TTFTMI", "LOS2TTFTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTFTME, UINT32, "LOS2TTFTME", "LOS2TTFTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2TTFT, UINT32, "LOS2TTFT", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTFTv, UINT32, "LOS2TTFTv", "LOS2TTFT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2TTCLM, UINT32, "LOS2TTCLM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTCLMI, UINT32, "LOS2TTCLMI", "LOS2TTCLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTCLME, UINT32, "LOS2TTCLME", "LOS2TTCLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2TTCL, UINT32, "LOS2TTCL", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTCLv, UINT32, "LOS2TTCLv", "LOS2TTCL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2TTCTM, UINT32, "LOS2TTCTM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTCTMI, UINT32, "LOS2TTCTMI", "LOS2TTCTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTCTME, UINT32, "LOS2TTCTME", "LOS2TTCTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS2TTCT, UINT32, "LOS2TTCT", 256);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS2TTCTv, UINT32, "LOS2TTCTv", "LOS2TTCT");    \
                                                                                                                  \
    } while (0);

#endif /*__GUARD_H101_LOS_EXT_H101_LOS_DEZ19_H__*/

/*******************************************************/
