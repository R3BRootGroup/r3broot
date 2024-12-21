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

#ifndef __GUARD_H101_PTOF_EXT_H101_PTOF_H__
#define __GUARD_H101_PTOF_EXT_H101_PTOF_H__

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

typedef struct EXT_STR_h101_PTOF_t
{
    /* RAW */
    uint32_t PTOF_TFLM /* [1,8] */;
    uint32_t PTOF_TFLMI[8 EXT_STRUCT_CTRL(PTOF_TFLM)] /* [1,8] */;
    uint32_t PTOF_TFLME[8 EXT_STRUCT_CTRL(PTOF_TFLM)] /* [1,80] */;
    uint32_t PTOF_TFL /* [0,80] */;
    uint32_t PTOF_TFLv[80 EXT_STRUCT_CTRL(PTOF_TFL)] /* [0,65535] */;
    uint32_t PTOF_TFTM /* [1,8] */;
    uint32_t PTOF_TFTMI[8 EXT_STRUCT_CTRL(PTOF_TFTM)] /* [1,8] */;
    uint32_t PTOF_TFTME[8 EXT_STRUCT_CTRL(PTOF_TFTM)] /* [1,80] */;
    uint32_t PTOF_TFT /* [0,80] */;
    uint32_t PTOF_TFTv[80 EXT_STRUCT_CTRL(PTOF_TFT)] /* [0,65535] */;
    uint32_t PTOF_TCLM /* [1,8] */;
    uint32_t PTOF_TCLMI[8 EXT_STRUCT_CTRL(PTOF_TCLM)] /* [1,8] */;
    uint32_t PTOF_TCLME[8 EXT_STRUCT_CTRL(PTOF_TCLM)] /* [1,80] */;
    uint32_t PTOF_TCL /* [0,80] */;
    uint32_t PTOF_TCLv[80 EXT_STRUCT_CTRL(PTOF_TCL)] /* [0,65535] */;
    uint32_t PTOF_TCTM /* [1,8] */;
    uint32_t PTOF_TCTMI[8 EXT_STRUCT_CTRL(PTOF_TCTM)] /* [1,8] */;
    uint32_t PTOF_TCTME[8 EXT_STRUCT_CTRL(PTOF_TCTM)] /* [1,80] */;
    uint32_t PTOF_TCT /* [0,80] */;
    uint32_t PTOF_TCTv[80 EXT_STRUCT_CTRL(PTOF_TCT)] /* [0,65535] */;

} EXT_STR_h101_PTOF;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_PTOF_onion_t
{
    /* RAW */
    uint32_t PTOF_TFLM;
    uint32_t PTOF_TFLMI[8 /* PTOF_TFLM */];
    uint32_t PTOF_TFLME[8 /* PTOF_TFLM */];
    uint32_t PTOF_TFL;
    uint32_t PTOF_TFLv[80 /* PTOF_TFL */];
    uint32_t PTOF_TFTM;
    uint32_t PTOF_TFTMI[8 /* PTOF_TFTM */];
    uint32_t PTOF_TFTME[8 /* PTOF_TFTM */];
    uint32_t PTOF_TFT;
    uint32_t PTOF_TFTv[80 /* PTOF_TFT */];
    uint32_t PTOF_TCLM;
    uint32_t PTOF_TCLMI[8 /* PTOF_TCLM */];
    uint32_t PTOF_TCLME[8 /* PTOF_TCLM */];
    uint32_t PTOF_TCL;
    uint32_t PTOF_TCLv[80 /* PTOF_TCL */];
    uint32_t PTOF_TCTM;
    uint32_t PTOF_TCTMI[8 /* PTOF_TCTM */];
    uint32_t PTOF_TCTME[8 /* PTOF_TCTM */];
    uint32_t PTOF_TCT;
    uint32_t PTOF_TCTv[80 /* PTOF_TCT */];

} EXT_STR_h101_PTOF_onion;

/*******************************************************/

#define EXT_STR_h101_PTOF_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                          \
    do                                                                                                            \
    {                                                                                                             \
        ok = 1;                                                                                                   \
        /* RAW */                                                                                                 \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, PTOF_TFLM, UINT32, "PTOF_TFLM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TFLMI, UINT32, "PTOF_TFLMI", "PTOF_TFLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TFLME, UINT32, "PTOF_TFLME", "PTOF_TFLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, PTOF_TFL, UINT32, "PTOF_TFL", 80);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TFLv, UINT32, "PTOF_TFLv", "PTOF_TFL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, PTOF_TFTM, UINT32, "PTOF_TFTM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TFTMI, UINT32, "PTOF_TFTMI", "PTOF_TFTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TFTME, UINT32, "PTOF_TFTME", "PTOF_TFTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, PTOF_TFT, UINT32, "PTOF_TFT", 80);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TFTv, UINT32, "PTOF_TFTv", "PTOF_TFT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, PTOF_TCLM, UINT32, "PTOF_TCLM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TCLMI, UINT32, "PTOF_TCLMI", "PTOF_TCLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TCLME, UINT32, "PTOF_TCLME", "PTOF_TCLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, PTOF_TCL, UINT32, "PTOF_TCL", 80);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TCLv, UINT32, "PTOF_TCLv", "PTOF_TCL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, PTOF_TCTM, UINT32, "PTOF_TCTM", 8);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TCTMI, UINT32, "PTOF_TCTMI", "PTOF_TCTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TCTME, UINT32, "PTOF_TCTME", "PTOF_TCTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, PTOF_TCT, UINT32, "PTOF_TCT", 80);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, PTOF_TCTv, UINT32, "PTOF_TCTv", "PTOF_TCT");    \
                                                                                                                  \
    } while (0);

#endif /*__GUARD_H101_PTOF_EXT_H101_PTOF_H__*/

/*******************************************************/
