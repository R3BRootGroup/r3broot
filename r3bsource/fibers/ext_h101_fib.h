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

#ifndef __GUARD_H101_FIB_EXT_H101_FIB_H__
#define __GUARD_H101_FIB_EXT_H101_FIB_H__

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

typedef struct EXT_STR_h101_FIB_t
{
    /* RAW */
    uint32_t FIB_TRIGSLC /* [0,3] */;
    uint32_t FIB_TRIGSLCI[3 EXT_STRUCT_CTRL(FIB_TRIGSLC)] /* [1,3] */;
    uint32_t FIB_TRIGSLCv[3 EXT_STRUCT_CTRL(FIB_TRIGSLC)] /* [0,65535] */;
    uint32_t FIB_TRIGSLF /* [0,3] */;
    uint32_t FIB_TRIGSLFI[3 EXT_STRUCT_CTRL(FIB_TRIGSLF)] /* [1,3] */;
    uint32_t FIB_TRIGSLFv[3 EXT_STRUCT_CTRL(FIB_TRIGSLF)] /* [0,65535] */;

} EXT_STR_h101_FIB;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIB_onion_t
{
    /* RAW */
    uint32_t FIB_TRIGSLC;
    uint32_t FIB_TRIGSLCI[3 /* FIB_TRIGSLC */];
    uint32_t FIB_TRIGSLCv[3 /* FIB_TRIGSLC */];
    uint32_t FIB_TRIGSLF;
    uint32_t FIB_TRIGSLFI[3 /* FIB_TRIGSLF */];
    uint32_t FIB_TRIGSLFv[3 /* FIB_TRIGSLF */];

} EXT_STR_h101_FIB_onion;

/*******************************************************/

#define EXT_STR_h101_FIB_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                   \
    do                                                                                                    \
    {                                                                                                     \
        ok = 1;                                                                                           \
        /* RAW */                                                                                         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIB_TRIGSLC, UINT32, "FIB_TRIGSLC", 3); \
        EXT_STR_ITEM_INFO_ZZP(                                                                            \
            ok, si, offset, struct_t, printerr, FIB_TRIGSLCI, UINT32, "FIB_TRIGSLCI", "FIB_TRIGSLC");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                            \
            ok, si, offset, struct_t, printerr, FIB_TRIGSLCv, UINT32, "FIB_TRIGSLCv", "FIB_TRIGSLC");     \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIB_TRIGSLF, UINT32, "FIB_TRIGSLF", 3); \
        EXT_STR_ITEM_INFO_ZZP(                                                                            \
            ok, si, offset, struct_t, printerr, FIB_TRIGSLFI, UINT32, "FIB_TRIGSLFI", "FIB_TRIGSLF");     \
        EXT_STR_ITEM_INFO_ZZP(                                                                            \
            ok, si, offset, struct_t, printerr, FIB_TRIGSLFv, UINT32, "FIB_TRIGSLFv", "FIB_TRIGSLF");     \
                                                                                                          \
    } while (0);

#endif /*__GUARD_H101_FIB_EXT_H101_FIB_H__*/

/*******************************************************/
