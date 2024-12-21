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

#ifndef __GUARD_H101_EXT_H101_H__
#define __GUARD_H101_EXT_H101_H__

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

typedef struct EXT_STR_h101_t
{
    /* UNPACK */
    uint32_t TRIGGER /* [0,15] */;
    uint32_t EVENTNO /* [-1,-1] */;
    /* RAW */
    uint32_t TOFD_P1B1T1tfl /* [0,65535] */;
    uint32_t TOFD_P1B1T1tft /* [0,65535] */;
    uint32_t TOFD_P1B1T1tcl /* [0,65535] */;
    uint32_t TOFD_P1B1T1tct /* [0,65535] */;
    uint32_t TOFD_P1B1T2tfl /* [0,65535] */;
    uint32_t TOFD_P1B1T2tft /* [0,65535] */;
    uint32_t TOFD_P1B1T2tcl /* [0,65535] */;
    uint32_t TOFD_P1B1T2tct /* [0,65535] */;
    uint32_t TOFD_P1B2T1tfl /* [0,65535] */;
    uint32_t TOFD_P1B2T1tft /* [0,65535] */;
    uint32_t TOFD_P1B2T1tcl /* [0,65535] */;
    uint32_t TOFD_P1B2T1tct /* [0,65535] */;
    uint32_t TOFD_P1B2T2tfl /* [0,65535] */;
    uint32_t TOFD_P1B2T2tft /* [0,65535] */;
    uint32_t TOFD_P1B2T2tcl /* [0,65535] */;
    uint32_t TOFD_P1B2T2tct /* [0,65535] */;
    uint32_t TOFD_P1B3T1tfl /* [0,65535] */;
    uint32_t TOFD_P1B3T1tft /* [0,65535] */;
    uint32_t TOFD_P1B3T1tcl /* [0,65535] */;
    uint32_t TOFD_P1B3T1tct /* [0,65535] */;
    uint32_t TOFD_P1B3T2tfl /* [0,65535] */;
    uint32_t TOFD_P1B3T2tft /* [0,65535] */;
    uint32_t TOFD_P1B3T2tcl /* [0,65535] */;
    uint32_t TOFD_P1B3T2tct /* [0,65535] */;
    uint32_t TOFD_P1B4T1tfl /* [0,65535] */;
    uint32_t TOFD_P1B4T1tft /* [0,65535] */;
    uint32_t TOFD_P1B4T1tcl /* [0,65535] */;
    uint32_t TOFD_P1B4T1tct /* [0,65535] */;
    uint32_t TOFD_P1B4T2tfl /* [0,65535] */;
    uint32_t TOFD_P1B4T2tft /* [0,65535] */;
    uint32_t TOFD_P1B4T2tcl /* [0,65535] */;
    uint32_t TOFD_P1B4T2tct /* [0,65535] */;

} EXT_STR_h101;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_onion_t
{
    /* UNPACK */
    uint32_t TRIGGER;
    uint32_t EVENTNO;
    /* RAW */
    struct
    {
        struct
        {
            struct
            {
                uint32_t tfl;
                uint32_t tft;
                uint32_t tcl;
                uint32_t tct;
            } T[2];
        } B[4];
    } TOFD_P[1];

} EXT_STR_h101_onion;

/*******************************************************/

#define EXT_STR_h101_ITEMS_INFO(ok, si, struct_t, printerr)                                                 \
    do                                                                                                      \
    {                                                                                                       \
        ok = 1;                                                                                             \
        /* UNPACK */                                                                                        \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TRIGGER, UINT32, "TRIGGER", 15);                  \
        EXT_STR_ITEM_INFO(ok, si, struct_t, printerr, EVENTNO, UINT32, "EVENTNO");                          \
        /* RAW */                                                                                           \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B1T1tfl, UINT32, "TOFD_P1B1T1tfl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B1T1tft, UINT32, "TOFD_P1B1T1tft", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B1T1tcl, UINT32, "TOFD_P1B1T1tcl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B1T1tct, UINT32, "TOFD_P1B1T1tct", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B1T2tfl, UINT32, "TOFD_P1B1T2tfl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B1T2tft, UINT32, "TOFD_P1B1T2tft", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B1T2tcl, UINT32, "TOFD_P1B1T2tcl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B1T2tct, UINT32, "TOFD_P1B1T2tct", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B2T1tfl, UINT32, "TOFD_P1B2T1tfl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B2T1tft, UINT32, "TOFD_P1B2T1tft", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B2T1tcl, UINT32, "TOFD_P1B2T1tcl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B2T1tct, UINT32, "TOFD_P1B2T1tct", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B2T2tfl, UINT32, "TOFD_P1B2T2tfl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B2T2tft, UINT32, "TOFD_P1B2T2tft", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B2T2tcl, UINT32, "TOFD_P1B2T2tcl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B2T2tct, UINT32, "TOFD_P1B2T2tct", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B3T1tfl, UINT32, "TOFD_P1B3T1tfl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B3T1tft, UINT32, "TOFD_P1B3T1tft", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B3T1tcl, UINT32, "TOFD_P1B3T1tcl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B3T1tct, UINT32, "TOFD_P1B3T1tct", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B3T2tfl, UINT32, "TOFD_P1B3T2tfl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B3T2tft, UINT32, "TOFD_P1B3T2tft", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B3T2tcl, UINT32, "TOFD_P1B3T2tcl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B3T2tct, UINT32, "TOFD_P1B3T2tct", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B4T1tfl, UINT32, "TOFD_P1B4T1tfl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B4T1tft, UINT32, "TOFD_P1B4T1tft", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B4T1tcl, UINT32, "TOFD_P1B4T1tcl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B4T1tct, UINT32, "TOFD_P1B4T1tct", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B4T2tfl, UINT32, "TOFD_P1B4T2tfl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B4T2tft, UINT32, "TOFD_P1B4T2tft", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B4T2tcl, UINT32, "TOFD_P1B4T2tcl", 65535); \
        EXT_STR_ITEM_INFO_LIM(ok, si, struct_t, printerr, TOFD_P1B4T2tct, UINT32, "TOFD_P1B4T2tct", 65535); \
                                                                                                            \
    } while (0);

#endif /*__GUARD_H101_EXT_H101_H__*/

/*******************************************************/
