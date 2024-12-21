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

#ifndef __GUARD_H101_AMS_EXT_H101_AMS_H__
#define __GUARD_H101_AMS_EXT_H101_AMS_H__

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

typedef struct EXT_STR_h101_AMS_t
{
    /* RAW */
    uint32_t SST1 /* [0,1024] */;
    uint32_t SST1I[1024 EXT_STRUCT_CTRL(SST1)] /* [1,1024] */;
    uint32_t SST1E[1024 EXT_STRUCT_CTRL(SST1)] /* [0,65535] */;
    uint32_t SST2 /* [0,1024] */;
    uint32_t SST2I[1024 EXT_STRUCT_CTRL(SST2)] /* [1,1024] */;
    uint32_t SST2E[1024 EXT_STRUCT_CTRL(SST2)] /* [0,65535] */;
    uint32_t SST3 /* [0,1024] */;
    uint32_t SST3I[1024 EXT_STRUCT_CTRL(SST3)] /* [1,1024] */;
    uint32_t SST3E[1024 EXT_STRUCT_CTRL(SST3)] /* [0,65535] */;
    uint32_t SST4 /* [0,1024] */;
    uint32_t SST4I[1024 EXT_STRUCT_CTRL(SST4)] /* [1,1024] */;
    uint32_t SST4E[1024 EXT_STRUCT_CTRL(SST4)] /* [0,65535] */;
    uint32_t SST5 /* [0,1024] */;
    uint32_t SST5I[1024 EXT_STRUCT_CTRL(SST5)] /* [1,1024] */;
    uint32_t SST5E[1024 EXT_STRUCT_CTRL(SST5)] /* [0,65535] */;
    uint32_t SST6 /* [0,1024] */;
    uint32_t SST6I[1024 EXT_STRUCT_CTRL(SST6)] /* [1,1024] */;
    uint32_t SST6E[1024 EXT_STRUCT_CTRL(SST6)] /* [0,65535] */;

} EXT_STR_h101_AMS;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_AMS_onion_t
{
    /* RAW */
    struct
    {
        uint32_t _;
        uint32_t I[1024 /*  */];
        uint32_t E[1024 /*  */];
    } SST[6];

} EXT_STR_h101_AMS_onion;

/*******************************************************/

#define EXT_STR_h101_AMS_ITEMS_INFO(ok, si, offset, struct_t, printerr)                            \
    do                                                                                             \
    {                                                                                              \
        ok = 1;                                                                                    \
        /* RAW */                                                                                  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SST1, UINT32, "SST1", 1024);     \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST1I, UINT32, "SST1I", "SST1"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST1E, UINT32, "SST1E", "SST1"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SST2, UINT32, "SST2", 1024);     \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST2I, UINT32, "SST2I", "SST2"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST2E, UINT32, "SST2E", "SST2"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SST3, UINT32, "SST3", 1024);     \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST3I, UINT32, "SST3I", "SST3"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST3E, UINT32, "SST3E", "SST3"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SST4, UINT32, "SST4", 1024);     \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST4I, UINT32, "SST4I", "SST4"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST4E, UINT32, "SST4E", "SST4"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SST5, UINT32, "SST5", 1024);     \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST5I, UINT32, "SST5I", "SST5"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST5E, UINT32, "SST5E", "SST5"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SST6, UINT32, "SST6", 1024);     \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST6I, UINT32, "SST6I", "SST6"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SST6E, UINT32, "SST6E", "SST6"); \
                                                                                                   \
    } while (0);

#endif /*__GUARD_H101_AMS_EXT_H101_AMS_H__*/

/*******************************************************/
