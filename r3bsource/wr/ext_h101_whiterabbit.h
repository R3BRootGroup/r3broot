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

#ifndef __GUARD_H101_EXT_H101_WHITERABBIT_H__
#define __GUARD_H101_EXT_H101_WHITERABBIT_H__

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

typedef struct EXT_STR_h101_whiterabbit_t
{
    uint32_t TIMESTAMP_MASTER_ID;
    uint32_t TIMESTAMP_MASTER_WR_T1;
    uint32_t TIMESTAMP_MASTER_WR_T2;
    uint32_t TIMESTAMP_MASTER_WR_T3;
    uint32_t TIMESTAMP_MASTER_WR_T4;
} EXT_STR_h101_whiterabbit;

typedef struct EXT_STR_h101_whiterabbit_onion_t
{
    uint32_t TIMESTAMP_MASTER_ID;
    uint32_t TIMESTAMP_MASTER_WR_T[4];
} EXT_STR_h101_whiterabbit_onion;

#define EXT_STR_h101_whiterabbit_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                   \
    do                                                                                                            \
    {                                                                                                             \
        ok = 1;                                                                                                   \
        /* Whiterabbit timestamp */                                                                               \
        EXT_STR_ITEM_INFO_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, TIMESTAMP_MASTER_ID, UINT32, "TIMESTAMP_MASTER_ID", 65535);       \
        EXT_STR_ITEM_INFO_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, TIMESTAMP_MASTER_WR_T1, UINT32, "TIMESTAMP_MASTER_WR_T1", 65535); \
        EXT_STR_ITEM_INFO_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, TIMESTAMP_MASTER_WR_T2, UINT32, "TIMESTAMP_MASTER_WR_T2", 65535); \
        EXT_STR_ITEM_INFO_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, TIMESTAMP_MASTER_WR_T3, UINT32, "TIMESTAMP_MASTER_WR_T3", 65535); \
        EXT_STR_ITEM_INFO_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, TIMESTAMP_MASTER_WR_T4, UINT32, "TIMESTAMP_MASTER_WR_T4", 65535); \
    } while (0);

#endif
