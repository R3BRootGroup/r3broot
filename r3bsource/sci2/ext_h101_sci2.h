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

#ifndef __GUARD_H101_SCI2_EXT_H101_SCI2_H__
#define __GUARD_H101_SCI2_EXT_H101_SCI2_H__

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

typedef struct EXT_STR_h101_SCI2_t
{
    /* RAW */
    uint32_t SCITWO_VTFM /* [1,8] */;
    uint32_t SCITWO_VTFMI[8 EXT_STRUCT_CTRL(SCITWO_VTFM)] /* [1,8] */;
    uint32_t SCITWO_VTFME[8 EXT_STRUCT_CTRL(SCITWO_VTFM)] /* [1,1024] */;
    uint32_t SCITWO_VTF /* [0,1024] */;
    uint32_t SCITWO_VTFv[1024 EXT_STRUCT_CTRL(SCITWO_VTF)] /* [0,65535] */;
    uint32_t SCITWO_VTCM /* [1,8] */;
    uint32_t SCITWO_VTCMI[8 EXT_STRUCT_CTRL(SCITWO_VTCM)] /* [1,8] */;
    uint32_t SCITWO_VTCME[8 EXT_STRUCT_CTRL(SCITWO_VTCM)] /* [1,1024] */;
    uint32_t SCITWO_VTC /* [0,1024] */;
    uint32_t SCITWO_VTCv[1024 EXT_STRUCT_CTRL(SCITWO_VTC)] /* [0,65535] */;

} EXT_STR_h101_SCI2;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_SCI2_onion_t
{
    /* RAW */
    uint32_t SCITWO_VTFM;
    uint32_t SCITWO_VTFMI[8 /* SCITWO_VTFM */];
    uint32_t SCITWO_VTFME[8 /* SCITWO_VTFM */];
    uint32_t SCITWO_VTF;
    uint32_t SCITWO_VTFv[1024 /* SCITWO_VTF */];
    uint32_t SCITWO_VTCM;
    uint32_t SCITWO_VTCMI[8 /* SCITWO_VTCM */];
    uint32_t SCITWO_VTCME[8 /* SCITWO_VTCM */];
    uint32_t SCITWO_VTC;
    uint32_t SCITWO_VTCv[1024 /* SCITWO_VTC */];

} EXT_STR_h101_SCI2_onion;

/*******************************************************/

#define EXT_STR_h101_SCI2_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                             \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SCITWO_VTFM, UINT32, "SCITWO_VTFM", 8);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SCITWO_VTFMI, UINT32, "SCITWO_VTFMI", "SCITWO_VTFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SCITWO_VTFME, UINT32, "SCITWO_VTFME", "SCITWO_VTFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SCITWO_VTF, UINT32, "SCITWO_VTF", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SCITWO_VTFv, UINT32, "SCITWO_VTFv", "SCITWO_VTF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SCITWO_VTCM, UINT32, "SCITWO_VTCM", 8);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SCITWO_VTCMI, UINT32, "SCITWO_VTCMI", "SCITWO_VTCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, SCITWO_VTCME, UINT32, "SCITWO_VTCME", "SCITWO_VTCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SCITWO_VTC, UINT32, "SCITWO_VTC", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SCITWO_VTCv, UINT32, "SCITWO_VTCv", "SCITWO_VTC"); \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_SCI2_EXT_H101_SCI2_H__*/

/*******************************************************/
