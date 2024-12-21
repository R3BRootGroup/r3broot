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

#ifndef __GUARD_H101_EXT_H101_LOS_TAMEX_H__
#define __GUARD_H101_EXT_H101_LOS_TAMEX_H__

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

typedef struct EXT_STR_h101_LOS_TAMEX_t
{
    /* RAW */
    uint32_t LOS1VTFM /* [1,16] */;
    uint32_t LOS1VTFMI[16 EXT_STRUCT_CTRL(LOS1VTFM)] /* [1,16] */;
    uint32_t LOS1VTFME[16 EXT_STRUCT_CTRL(LOS1VTFM)] /* [1,160] */;
    uint32_t LOS1VTF /* [0,160] */;
    uint32_t LOS1VTFv[160 EXT_STRUCT_CTRL(LOS1VTF)] /* [0,65535] */;
    uint32_t LOS1VTCM /* [1,16] */;
    uint32_t LOS1VTCMI[16 EXT_STRUCT_CTRL(LOS1VTCM)] /* [1,16] */;
    uint32_t LOS1VTCME[16 EXT_STRUCT_CTRL(LOS1VTCM)] /* [1,160] */;
    uint32_t LOS1VTC /* [0,160] */;
    uint32_t LOS1VTCv[160 EXT_STRUCT_CTRL(LOS1VTC)] /* [0,65535] */;
    uint32_t LOS1TTFLM /* [1,16] */;
    uint32_t LOS1TTFLMI[16 EXT_STRUCT_CTRL(LOS1TTFLM)] /* [1,16] */;
    uint32_t LOS1TTFLME[16 EXT_STRUCT_CTRL(LOS1TTFLM)] /* [1,160] */;
    uint32_t LOS1TTFL /* [0,160] */;
    uint32_t LOS1TTFLv[160 EXT_STRUCT_CTRL(LOS1TTFL)] /* [0,65535] */;
    uint32_t LOS1TTFTM /* [1,16] */;
    uint32_t LOS1TTFTMI[16 EXT_STRUCT_CTRL(LOS1TTFTM)] /* [1,16] */;
    uint32_t LOS1TTFTME[16 EXT_STRUCT_CTRL(LOS1TTFTM)] /* [1,160] */;
    uint32_t LOS1TTFT /* [0,160] */;
    uint32_t LOS1TTFTv[160 EXT_STRUCT_CTRL(LOS1TTFT)] /* [0,65535] */;
    uint32_t LOS1TTCLM /* [1,16] */;
    uint32_t LOS1TTCLMI[16 EXT_STRUCT_CTRL(LOS1TTCLM)] /* [1,16] */;
    uint32_t LOS1TTCLME[16 EXT_STRUCT_CTRL(LOS1TTCLM)] /* [1,160] */;
    uint32_t LOS1TTCL /* [0,160] */;
    uint32_t LOS1TTCLv[160 EXT_STRUCT_CTRL(LOS1TTCL)] /* [0,65535] */;
    uint32_t LOS1TTCTM /* [1,16] */;
    uint32_t LOS1TTCTMI[16 EXT_STRUCT_CTRL(LOS1TTCTM)] /* [1,16] */;
    uint32_t LOS1TTCTME[16 EXT_STRUCT_CTRL(LOS1TTCTM)] /* [1,160] */;
    uint32_t LOS1TTCT /* [0,160] */;
    uint32_t LOS1TTCTv[160 EXT_STRUCT_CTRL(LOS1TTCT)] /* [0,65535] */;

} EXT_STR_h101_LOS_TAMEX;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_LOS_TAMEX_onion_t
{
    /* RAW */
    struct
    {
        uint32_t VTFM;
        uint32_t VTFMI[16 /* VTFM */];
        uint32_t VTFME[16 /* VTFM */];
        uint32_t VTF;
        uint32_t VTFv[160 /* VTF */];
        uint32_t VTCM;
        uint32_t VTCMI[16 /* VTCM */];
        uint32_t VTCME[16 /* VTCM */];
        uint32_t VTC;
        uint32_t VTCv[160 /* VTC */];
        uint32_t TTFLM;
        uint32_t TTFLMI[16 /* TTFLM */];
        uint32_t TTFLME[16 /* TTFLM */];
        uint32_t TTFL;
        uint32_t TTFLv[160 /* TTFL */];
        uint32_t TTFTM;
        uint32_t TTFTMI[16 /* TTFTM */];
        uint32_t TTFTME[16 /* TTFTM */];
        uint32_t TTFT;
        uint32_t TTFTv[160 /* TTFT */];
        uint32_t TTCLM;
        uint32_t TTCLMI[16 /* TTCLM */];
        uint32_t TTCLME[16 /* TTCLM */];
        uint32_t TTCL;
        uint32_t TTCLv[160 /* TTCL */];
        uint32_t TTCTM;
        uint32_t TTCTMI[16 /* TTCTM */];
        uint32_t TTCTME[16 /* TTCTM */];
        uint32_t TTCT;
        uint32_t TTCTv[160 /* TTCT */];
    } LOS[1];

} EXT_STR_h101_LOS_TAMEX_onion;

/*******************************************************/

#define EXT_STR_h101_LOS_TAMEX_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                     \
    do                                                                                                            \
    {                                                                                                             \
        ok = 1;                                                                                                   \
        /* RAW */                                                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1VTFM, UINT32, "LOS1VTFM", 16);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTFMI, UINT32, "LOS1VTFMI", "LOS1VTFM");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTFME, UINT32, "LOS1VTFME", "LOS1VTFM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1VTF, UINT32, "LOS1VTF", 160);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTFv, UINT32, "LOS1VTFv", "LOS1VTF");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1VTCM, UINT32, "LOS1VTCM", 16);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTCMI, UINT32, "LOS1VTCMI", "LOS1VTCM");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTCME, UINT32, "LOS1VTCME", "LOS1VTCM");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1VTC, UINT32, "LOS1VTC", 160);               \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1VTCv, UINT32, "LOS1VTCv", "LOS1VTC");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTFLM, UINT32, "LOS1TTFLM", 16);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFLMI, UINT32, "LOS1TTFLMI", "LOS1TTFLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFLME, UINT32, "LOS1TTFLME", "LOS1TTFLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTFL, UINT32, "LOS1TTFL", 160);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFLv, UINT32, "LOS1TTFLv", "LOS1TTFL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTFTM, UINT32, "LOS1TTFTM", 16);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFTMI, UINT32, "LOS1TTFTMI", "LOS1TTFTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFTME, UINT32, "LOS1TTFTME", "LOS1TTFTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTFT, UINT32, "LOS1TTFT", 160);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTFTv, UINT32, "LOS1TTFTv", "LOS1TTFT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTCLM, UINT32, "LOS1TTCLM", 16);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCLMI, UINT32, "LOS1TTCLMI", "LOS1TTCLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCLME, UINT32, "LOS1TTCLME", "LOS1TTCLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTCL, UINT32, "LOS1TTCL", 160);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCLv, UINT32, "LOS1TTCLv", "LOS1TTCL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTCTM, UINT32, "LOS1TTCTM", 16);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCTMI, UINT32, "LOS1TTCTMI", "LOS1TTCTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCTME, UINT32, "LOS1TTCTME", "LOS1TTCTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, LOS1TTCT, UINT32, "LOS1TTCT", 160);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, LOS1TTCTv, UINT32, "LOS1TTCTv", "LOS1TTCT");    \
                                                                                                                  \
    } while (0);

#endif /*__GUARD_H101_EXT_H101_LOS_TAMEX_H__*/
