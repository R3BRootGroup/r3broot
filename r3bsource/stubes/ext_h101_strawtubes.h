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

#ifndef __GUARD_H101_EXT_H101_STRAWTUBES_H__
#define __GUARD_H101_EXT_H101_STRAWTUBES_H__

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

typedef struct EXT_STR_h101_STRAWTUBES_t
{
    uint32_t STRAW1S1TFM /* [1,32] */;
    uint32_t STRAW1S1TFMI[32 EXT_STRUCT_CTRL(STRAW1S1TFM)] /* [1,32] */;
    uint32_t STRAW1S1TFME[32 EXT_STRUCT_CTRL(STRAW1S1TFM)] /* [1,6400] */;
    uint32_t STRAW1S1TF /* [0,6400] */;
    uint32_t STRAW1S1TFv[6400 EXT_STRUCT_CTRL(STRAW1S1TF)] /* [0,65535] */;
    uint32_t STRAW1S1TCM /* [1,32] */;
    uint32_t STRAW1S1TCMI[32 EXT_STRUCT_CTRL(STRAW1S1TCM)] /* [1,32] */;
    uint32_t STRAW1S1TCME[32 EXT_STRUCT_CTRL(STRAW1S1TCM)] /* [1,6400] */;
    uint32_t STRAW1S1TC /* [0,6400] */;
    uint32_t STRAW1S1TCv[6400 EXT_STRUCT_CTRL(STRAW1S1TC)] /* [0,65535] */;
    uint32_t STRAW1S1TRIG /* [0,65535] */;

} EXT_STR_h101_STRAWTUBES;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_STRAWTUBES_onion_t
{
    struct
    {
        struct
        {
            uint32_t TFM;
            uint32_t TFMI[32 /* TFM */];
            uint32_t TFME[32 /* TFM */];
            uint32_t TF;
            uint32_t TFv[6400 /* TF */];
            uint32_t TCM;
            uint32_t TCMI[32 /* TCM */];
            uint32_t TCME[32 /* TCM */];
            uint32_t TC;
            uint32_t TCv[6400 /* TC */];
            uint32_t TRIG;
        } S[1];
    } STRAW[1];

} EXT_STR_h101_STRAWTUBES_onion;

/*******************************************************/

#define EXT_STR_h101_STRAWTUBES_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                       \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, STRAW1S1TFM, UINT32, "STRAW1S1TFM", 32);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, STRAW1S1TFMI, UINT32, "STRAW1S1TFMI", "STRAW1S1TFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, STRAW1S1TFME, UINT32, "STRAW1S1TFME", "STRAW1S1TFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, STRAW1S1TF, UINT32, "STRAW1S1TF", 6400);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, STRAW1S1TFv, UINT32, "STRAW1S1TFv", "STRAW1S1TF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, STRAW1S1TCM, UINT32, "STRAW1S1TCM", 32);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, STRAW1S1TCMI, UINT32, "STRAW1S1TCMI", "STRAW1S1TCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, STRAW1S1TCME, UINT32, "STRAW1S1TCME", "STRAW1S1TCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, STRAW1S1TC, UINT32, "STRAW1S1TC", 6400);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, STRAW1S1TCv, UINT32, "STRAW1S1TCv", "STRAW1S1TC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, STRAW1S1TRIG, UINT32, "STRAW1S1TRIG", 65535);      \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_EXT_H101_STRAWTUBES_H__*/

/*******************************************************/
