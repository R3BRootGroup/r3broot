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

#ifndef __GUARD_H101_SOFMWPC_EXT_H101_H__
#define __GUARD_H101_SOFMWPC_EXT_H101_H__

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

typedef struct EXT_STR_h101_SOFMWPC_t
{
    /* RAW */
    uint32_t SOFMWPC1Plane1Q /* [0,320] */;
    uint32_t SOFMWPC1Plane1QI[320 EXT_STRUCT_CTRL(SOFMWPC1Plane1Q)] /* [1,320] */;
    uint32_t SOFMWPC1Plane1Qv[320 EXT_STRUCT_CTRL(SOFMWPC1Plane1Q)] /* [0,65535] */;
    uint32_t SOFMWPC1Plane2Q /* [0,320] */;
    uint32_t SOFMWPC1Plane2QI[320 EXT_STRUCT_CTRL(SOFMWPC1Plane2Q)] /* [1,320] */;
    uint32_t SOFMWPC1Plane2Qv[320 EXT_STRUCT_CTRL(SOFMWPC1Plane2Q)] /* [0,65535] */;
    uint32_t SOFMWPC1Plane3Q /* [0,320] */;
    uint32_t SOFMWPC1Plane3QI[320 EXT_STRUCT_CTRL(SOFMWPC1Plane3Q)] /* [1,320] */;
    uint32_t SOFMWPC1Plane3Qv[320 EXT_STRUCT_CTRL(SOFMWPC1Plane3Q)] /* [0,65535] */;
    uint32_t SOFMWPC2Plane1Q /* [0,320] */;
    uint32_t SOFMWPC2Plane1QI[320 EXT_STRUCT_CTRL(SOFMWPC2Plane1Q)] /* [1,320] */;
    uint32_t SOFMWPC2Plane1Qv[320 EXT_STRUCT_CTRL(SOFMWPC2Plane1Q)] /* [0,65535] */;
    uint32_t SOFMWPC2Plane2Q /* [0,320] */;
    uint32_t SOFMWPC2Plane2QI[320 EXT_STRUCT_CTRL(SOFMWPC2Plane2Q)] /* [1,320] */;
    uint32_t SOFMWPC2Plane2Qv[320 EXT_STRUCT_CTRL(SOFMWPC2Plane2Q)] /* [0,65535] */;
    uint32_t SOFMWPC2Plane3Q /* [0,320] */;
    uint32_t SOFMWPC2Plane3QI[320 EXT_STRUCT_CTRL(SOFMWPC2Plane3Q)] /* [1,320] */;
    uint32_t SOFMWPC2Plane3Qv[320 EXT_STRUCT_CTRL(SOFMWPC2Plane3Q)] /* [0,65535] */;
    uint32_t SOFMWPC3Plane1Q /* [0,320] */;
    uint32_t SOFMWPC3Plane1QI[320 EXT_STRUCT_CTRL(SOFMWPC3Plane1Q)] /* [1,320] */;
    uint32_t SOFMWPC3Plane1Qv[320 EXT_STRUCT_CTRL(SOFMWPC3Plane1Q)] /* [0,65535] */;
    uint32_t SOFMWPC3Plane2Q /* [0,320] */;
    uint32_t SOFMWPC3Plane2QI[320 EXT_STRUCT_CTRL(SOFMWPC3Plane2Q)] /* [1,320] */;
    uint32_t SOFMWPC3Plane2Qv[320 EXT_STRUCT_CTRL(SOFMWPC3Plane2Q)] /* [0,65535] */;
    uint32_t SOFMWPC3Plane3Q /* [0,320] */;
    uint32_t SOFMWPC3Plane3QI[320 EXT_STRUCT_CTRL(SOFMWPC3Plane3Q)] /* [1,320] */;
    uint32_t SOFMWPC3Plane3Qv[320 EXT_STRUCT_CTRL(SOFMWPC3Plane3Q)] /* [0,65535] */;
    uint32_t SOFMWPC4Plane1Q /* [0,320] */;
    uint32_t SOFMWPC4Plane1QI[320 EXT_STRUCT_CTRL(SOFMWPC4Plane1Q)] /* [1,320] */;
    uint32_t SOFMWPC4Plane1Qv[320 EXT_STRUCT_CTRL(SOFMWPC4Plane1Q)] /* [0,65535] */;
    uint32_t SOFMWPC4Plane2Q /* [0,320] */;
    uint32_t SOFMWPC4Plane2QI[320 EXT_STRUCT_CTRL(SOFMWPC4Plane2Q)] /* [1,320] */;
    uint32_t SOFMWPC4Plane2Qv[320 EXT_STRUCT_CTRL(SOFMWPC4Plane2Q)] /* [0,65535] */;
    uint32_t SOFMWPC4Plane3Q /* [0,320] */;
    uint32_t SOFMWPC4Plane3QI[320 EXT_STRUCT_CTRL(SOFMWPC4Plane3Q)] /* [1,320] */;
    uint32_t SOFMWPC4Plane3Qv[320 EXT_STRUCT_CTRL(SOFMWPC4Plane3Q)] /* [0,65535] */;

} EXT_STR_h101_SOFMWPC;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_SOFMWPC_onion_t
{
    /* RAW */
    struct
    {
        struct
        {
            uint32_t Q;
            uint32_t QI[320 /* Q */];
            uint32_t Qv[320 /* Q */];
        } Plane[3];
    } SOFMWPC[4];

} EXT_STR_h101_SOFMWPC_onion;

/*******************************************************/

#define EXT_STR_h101_SOFMWPC_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                         \
    do                                                                                                              \
    {                                                                                                               \
        ok = 1;                                                                                                     \
        /* RAW */                                                                                                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC1Plane1Q, UINT32, "SOFMWPC1Plane1Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC1Plane1QI, UINT32, "SOFMWPC1Plane1QI", "SOFMWPC1Plane1Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC1Plane1Qv, UINT32, "SOFMWPC1Plane1Qv", "SOFMWPC1Plane1Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC1Plane2Q, UINT32, "SOFMWPC1Plane2Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC1Plane2QI, UINT32, "SOFMWPC1Plane2QI", "SOFMWPC1Plane2Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC1Plane2Qv, UINT32, "SOFMWPC1Plane2Qv", "SOFMWPC1Plane2Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC1Plane3Q, UINT32, "SOFMWPC1Plane3Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC1Plane3QI, UINT32, "SOFMWPC1Plane3QI", "SOFMWPC1Plane3Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC1Plane3Qv, UINT32, "SOFMWPC1Plane3Qv", "SOFMWPC1Plane3Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC2Plane1Q, UINT32, "SOFMWPC2Plane1Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC2Plane1QI, UINT32, "SOFMWPC2Plane1QI", "SOFMWPC2Plane1Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC2Plane1Qv, UINT32, "SOFMWPC2Plane1Qv", "SOFMWPC2Plane1Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC2Plane2Q, UINT32, "SOFMWPC2Plane2Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC2Plane2QI, UINT32, "SOFMWPC2Plane2QI", "SOFMWPC2Plane2Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC2Plane2Qv, UINT32, "SOFMWPC2Plane2Qv", "SOFMWPC2Plane2Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC2Plane3Q, UINT32, "SOFMWPC2Plane3Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC2Plane3QI, UINT32, "SOFMWPC2Plane3QI", "SOFMWPC2Plane3Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC2Plane3Qv, UINT32, "SOFMWPC2Plane3Qv", "SOFMWPC2Plane3Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC3Plane1Q, UINT32, "SOFMWPC3Plane1Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC3Plane1QI, UINT32, "SOFMWPC3Plane1QI", "SOFMWPC3Plane1Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC3Plane1Qv, UINT32, "SOFMWPC3Plane1Qv", "SOFMWPC3Plane1Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC3Plane2Q, UINT32, "SOFMWPC3Plane2Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC3Plane2QI, UINT32, "SOFMWPC3Plane2QI", "SOFMWPC3Plane2Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC3Plane2Qv, UINT32, "SOFMWPC3Plane2Qv", "SOFMWPC3Plane2Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC3Plane3Q, UINT32, "SOFMWPC3Plane3Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC3Plane3QI, UINT32, "SOFMWPC3Plane3QI", "SOFMWPC3Plane3Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC3Plane3Qv, UINT32, "SOFMWPC3Plane3Qv", "SOFMWPC3Plane3Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC4Plane1Q, UINT32, "SOFMWPC4Plane1Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC4Plane1QI, UINT32, "SOFMWPC4Plane1QI", "SOFMWPC4Plane1Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC4Plane1Qv, UINT32, "SOFMWPC4Plane1Qv", "SOFMWPC4Plane1Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC4Plane2Q, UINT32, "SOFMWPC4Plane2Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC4Plane2QI, UINT32, "SOFMWPC4Plane2QI", "SOFMWPC4Plane2Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC4Plane2Qv, UINT32, "SOFMWPC4Plane2Qv", "SOFMWPC4Plane2Q");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SOFMWPC4Plane3Q, UINT32, "SOFMWPC4Plane3Q", 320); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC4Plane3QI, UINT32, "SOFMWPC4Plane3QI", "SOFMWPC4Plane3Q");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                      \
            ok, si, offset, struct_t, printerr, SOFMWPC4Plane3Qv, UINT32, "SOFMWPC4Plane3Qv", "SOFMWPC4Plane3Q");   \
                                                                                                                    \
    } while (0);

/********************************************************
 *
 * For internal use by the network data reader:
 * (version checks, etc)
 */

#endif /*__GUARD_H101_SOFMWPC_EXT_H101_SOFMWPC_H__*/

/*******************************************************/
