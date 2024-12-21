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

#ifndef __GUARD_H101_SCI8_EXT_H101_SCI8_H__
#define __GUARD_H101_SCI8_EXT_H101_SCI8_H__

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

typedef struct EXT_STR_h101_SCI8_t
{
    /* RAW */
    uint32_t SCIEIGHT_VTFM /* [1,8] */;
    uint32_t SCIEIGHT_VTFMI[8 EXT_STRUCT_CTRL(SCIEIGHT_VTFM)] /* [1,8] */;
    uint32_t SCIEIGHT_VTFME[8 EXT_STRUCT_CTRL(SCIEIGHT_VTFM)] /* [1,256] */;
    uint32_t SCIEIGHT_VTF /* [0,256] */;
    uint32_t SCIEIGHT_VTFv[256 EXT_STRUCT_CTRL(SCIEIGHT_VTF)] /* [0,65535] */;
    uint32_t SCIEIGHT_VTCM /* [1,8] */;
    uint32_t SCIEIGHT_VTCMI[8 EXT_STRUCT_CTRL(SCIEIGHT_VTCM)] /* [1,8] */;
    uint32_t SCIEIGHT_VTCME[8 EXT_STRUCT_CTRL(SCIEIGHT_VTCM)] /* [1,256] */;
    uint32_t SCIEIGHT_VTC /* [0,256] */;
    uint32_t SCIEIGHT_VTCv[256 EXT_STRUCT_CTRL(SCIEIGHT_VTC)] /* [0,65535] */;
} EXT_STR_h101_SCI8;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_SCI8_onion_t
{
    /* RAW */
    uint32_t SCIEIGHT_VTFM;
    uint32_t SCIEIGHT_VTFMI[8 /* SCIEIGHT_VTFM */];
    uint32_t SCIEIGHT_VTFME[8 /* SCIEIGHT_VTFM */];
    uint32_t SCIEIGHT_VTF;
    uint32_t SCIEIGHT_VTFv[256 /* SCIEIGHT_VTF */];
    uint32_t SCIEIGHT_VTCM;
    uint32_t SCIEIGHT_VTCMI[8 /* SCIEIGHT_VTCM */];
    uint32_t SCIEIGHT_VTCME[8 /* SCIEIGHT_VTCM */];
    uint32_t SCIEIGHT_VTC;
    uint32_t SCIEIGHT_VTCv[256 /* SCIEIGHT_VTC */];
} EXT_STR_h101_SCI8_onion;

/*******************************************************/

#define EXT_STR_h101_SCI8_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                      \
    do                                                                                                        \
    {                                                                                                         \
        ok = 1;                                                                                               \
        /* RAW */                                                                                             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SCIEIGHT_VTFM, UINT32, "SCIEIGHT_VTFM", 8); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, SCIEIGHT_VTFMI, UINT32, "SCIEIGHT_VTFMI", "SCIEIGHT_VTFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, SCIEIGHT_VTFME, UINT32, "SCIEIGHT_VTFME", "SCIEIGHT_VTFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SCIEIGHT_VTF, UINT32, "SCIEIGHT_VTF", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, SCIEIGHT_VTFv, UINT32, "SCIEIGHT_VTFv", "SCIEIGHT_VTF");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SCIEIGHT_VTCM, UINT32, "SCIEIGHT_VTCM", 8); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, SCIEIGHT_VTCMI, UINT32, "SCIEIGHT_VTCMI", "SCIEIGHT_VTCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, SCIEIGHT_VTCME, UINT32, "SCIEIGHT_VTCME", "SCIEIGHT_VTCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SCIEIGHT_VTC, UINT32, "SCIEIGHT_VTC", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, SCIEIGHT_VTCv, UINT32, "SCIEIGHT_VTCv", "SCIEIGHT_VTC");      \
                                                                                                              \
    } while (0);

#endif /*__GUARD_H101_SCI8_EXT_H101_SCI8_H__*/

/*******************************************************/
