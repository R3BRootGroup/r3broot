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

#ifndef __GUARD_H101_TOFI_EXT_H101_TOFI_H__
#define __GUARD_H101_TOFI_EXT_H101_TOFI_H__

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

typedef struct EXT_STR_h101_TOFI_t
{
    /* RAW */
    uint32_t TOFI_TRIGCL /* [0,3] */;
    uint32_t TOFI_TRIGCLI[3 EXT_STRUCT_CTRL(TOFI_TRIGCL)] /* [1,3] */;
    uint32_t TOFI_TRIGCLv[3 EXT_STRUCT_CTRL(TOFI_TRIGCL)] /* [0,65535] */;
    uint32_t TOFI_TRIGFL /* [0,3] */;
    uint32_t TOFI_TRIGFLI[3 EXT_STRUCT_CTRL(TOFI_TRIGFL)] /* [1,3] */;
    uint32_t TOFI_TRIGFLv[3 EXT_STRUCT_CTRL(TOFI_TRIGFL)] /* [0,65535] */;
    uint32_t TOFI_T1TCLM /* [1,24] */;
    uint32_t TOFI_T1TCLMI[24 EXT_STRUCT_CTRL(TOFI_T1TCLM)] /* [1,24] */;
    uint32_t TOFI_T1TCLME[24 EXT_STRUCT_CTRL(TOFI_T1TCLM)] /* [1,1536] */;
    uint32_t TOFI_T1TCL /* [0,1536] */;
    uint32_t TOFI_T1TCLv[1536 EXT_STRUCT_CTRL(TOFI_T1TCL)] /* [0,65535] */;
    uint32_t TOFI_T1TFLM /* [1,24] */;
    uint32_t TOFI_T1TFLMI[24 EXT_STRUCT_CTRL(TOFI_T1TFLM)] /* [1,24] */;
    uint32_t TOFI_T1TFLME[24 EXT_STRUCT_CTRL(TOFI_T1TFLM)] /* [1,1536] */;
    uint32_t TOFI_T1TFL /* [0,1536] */;
    uint32_t TOFI_T1TFLv[1536 EXT_STRUCT_CTRL(TOFI_T1TFL)] /* [0,65535] */;
    uint32_t TOFI_T1TCTM /* [1,24] */;
    uint32_t TOFI_T1TCTMI[24 EXT_STRUCT_CTRL(TOFI_T1TCTM)] /* [1,24] */;
    uint32_t TOFI_T1TCTME[24 EXT_STRUCT_CTRL(TOFI_T1TCTM)] /* [1,1536] */;
    uint32_t TOFI_T1TCT /* [0,1536] */;
    uint32_t TOFI_T1TCTv[1536 EXT_STRUCT_CTRL(TOFI_T1TCT)] /* [0,65535] */;
    uint32_t TOFI_T1TFTM /* [1,24] */;
    uint32_t TOFI_T1TFTMI[24 EXT_STRUCT_CTRL(TOFI_T1TFTM)] /* [1,24] */;
    uint32_t TOFI_T1TFTME[24 EXT_STRUCT_CTRL(TOFI_T1TFTM)] /* [1,1536] */;
    uint32_t TOFI_T1TFT /* [0,1536] */;
    uint32_t TOFI_T1TFTv[1536 EXT_STRUCT_CTRL(TOFI_T1TFT)] /* [0,65535] */;
    uint32_t TOFI_T2TCLM /* [1,24] */;
    uint32_t TOFI_T2TCLMI[24 EXT_STRUCT_CTRL(TOFI_T2TCLM)] /* [1,24] */;
    uint32_t TOFI_T2TCLME[24 EXT_STRUCT_CTRL(TOFI_T2TCLM)] /* [1,1536] */;
    uint32_t TOFI_T2TCL /* [0,1536] */;
    uint32_t TOFI_T2TCLv[1536 EXT_STRUCT_CTRL(TOFI_T2TCL)] /* [0,65535] */;
    uint32_t TOFI_T2TFLM /* [1,24] */;
    uint32_t TOFI_T2TFLMI[24 EXT_STRUCT_CTRL(TOFI_T2TFLM)] /* [1,24] */;
    uint32_t TOFI_T2TFLME[24 EXT_STRUCT_CTRL(TOFI_T2TFLM)] /* [1,1536] */;
    uint32_t TOFI_T2TFL /* [0,1536] */;
    uint32_t TOFI_T2TFLv[1536 EXT_STRUCT_CTRL(TOFI_T2TFL)] /* [0,65535] */;
    uint32_t TOFI_T2TCTM /* [1,24] */;
    uint32_t TOFI_T2TCTMI[24 EXT_STRUCT_CTRL(TOFI_T2TCTM)] /* [1,24] */;
    uint32_t TOFI_T2TCTME[24 EXT_STRUCT_CTRL(TOFI_T2TCTM)] /* [1,1536] */;
    uint32_t TOFI_T2TCT /* [0,1536] */;
    uint32_t TOFI_T2TCTv[1536 EXT_STRUCT_CTRL(TOFI_T2TCT)] /* [0,65535] */;
    uint32_t TOFI_T2TFTM /* [1,24] */;
    uint32_t TOFI_T2TFTMI[24 EXT_STRUCT_CTRL(TOFI_T2TFTM)] /* [1,24] */;
    uint32_t TOFI_T2TFTME[24 EXT_STRUCT_CTRL(TOFI_T2TFTM)] /* [1,1536] */;
    uint32_t TOFI_T2TFT /* [0,1536] */;
    uint32_t TOFI_T2TFTv[1536 EXT_STRUCT_CTRL(TOFI_T2TFT)] /* [0,65535] */;

} EXT_STR_h101_TOFI;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_TOFI_onion_t
{
    /* RAW */
    uint32_t TOFI_TRIGCL;
    uint32_t TOFI_TRIGCLI[3 /* TOFI_TRIGCL */];
    uint32_t TOFI_TRIGCLv[3 /* TOFI_TRIGCL */];
    uint32_t TOFI_TRIGFL;
    uint32_t TOFI_TRIGFLI[3 /* TOFI_TRIGFL */];
    uint32_t TOFI_TRIGFLv[3 /* TOFI_TRIGFL */];
    struct
    {
        uint32_t TCLM;
        uint32_t TCLMI[24 /* TCLM */];
        uint32_t TCLME[24 /* TCLM */];
        uint32_t TCL;
        uint32_t TCLv[1536 /* TCL */];
        uint32_t TFLM;
        uint32_t TFLMI[24 /* TFLM */];
        uint32_t TFLME[24 /* TFLM */];
        uint32_t TFL;
        uint32_t TFLv[1536 /* TFL */];
        uint32_t TCTM;
        uint32_t TCTMI[24 /* TCTM */];
        uint32_t TCTME[24 /* TCTM */];
        uint32_t TCT;
        uint32_t TCTv[1536 /* TCT */];
        uint32_t TFTM;
        uint32_t TFTMI[24 /* TFTM */];
        uint32_t TFTME[24 /* TFTM */];
        uint32_t TFT;
        uint32_t TFTv[1536 /* TFT */];
    } TOFI_T[2];

} EXT_STR_h101_TOFI_onion;

/*******************************************************/

#define EXT_STR_h101_TOFI_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                             \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_TRIGCL, UINT32, "TOFI_TRIGCL", 3);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_TRIGCLI, UINT32, "TOFI_TRIGCLI", "TOFI_TRIGCL");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_TRIGCLv, UINT32, "TOFI_TRIGCLv", "TOFI_TRIGCL");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_TRIGFL, UINT32, "TOFI_TRIGFL", 3);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_TRIGFLI, UINT32, "TOFI_TRIGFLI", "TOFI_TRIGFL");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_TRIGFLv, UINT32, "TOFI_TRIGFLv", "TOFI_TRIGFL");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T1TCLM, UINT32, "TOFI_T1TCLM", 24);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T1TCLMI, UINT32, "TOFI_T1TCLMI", "TOFI_T1TCLM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T1TCLME, UINT32, "TOFI_T1TCLME", "TOFI_T1TCLM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T1TCL, UINT32, "TOFI_T1TCL", 1536);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, TOFI_T1TCLv, UINT32, "TOFI_T1TCLv", "TOFI_T1TCL"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T1TFLM, UINT32, "TOFI_T1TFLM", 24);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T1TFLMI, UINT32, "TOFI_T1TFLMI", "TOFI_T1TFLM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T1TFLME, UINT32, "TOFI_T1TFLME", "TOFI_T1TFLM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T1TFL, UINT32, "TOFI_T1TFL", 1536);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, TOFI_T1TFLv, UINT32, "TOFI_T1TFLv", "TOFI_T1TFL"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T1TCTM, UINT32, "TOFI_T1TCTM", 24);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T1TCTMI, UINT32, "TOFI_T1TCTMI", "TOFI_T1TCTM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T1TCTME, UINT32, "TOFI_T1TCTME", "TOFI_T1TCTM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T1TCT, UINT32, "TOFI_T1TCT", 1536);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, TOFI_T1TCTv, UINT32, "TOFI_T1TCTv", "TOFI_T1TCT"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T1TFTM, UINT32, "TOFI_T1TFTM", 24);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T1TFTMI, UINT32, "TOFI_T1TFTMI", "TOFI_T1TFTM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T1TFTME, UINT32, "TOFI_T1TFTME", "TOFI_T1TFTM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T1TFT, UINT32, "TOFI_T1TFT", 1536);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, TOFI_T1TFTv, UINT32, "TOFI_T1TFTv", "TOFI_T1TFT"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T2TCLM, UINT32, "TOFI_T2TCLM", 24);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T2TCLMI, UINT32, "TOFI_T2TCLMI", "TOFI_T2TCLM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T2TCLME, UINT32, "TOFI_T2TCLME", "TOFI_T2TCLM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T2TCL, UINT32, "TOFI_T2TCL", 1536);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, TOFI_T2TCLv, UINT32, "TOFI_T2TCLv", "TOFI_T2TCL"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T2TFLM, UINT32, "TOFI_T2TFLM", 24);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T2TFLMI, UINT32, "TOFI_T2TFLMI", "TOFI_T2TFLM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T2TFLME, UINT32, "TOFI_T2TFLME", "TOFI_T2TFLM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T2TFL, UINT32, "TOFI_T2TFL", 1536);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, TOFI_T2TFLv, UINT32, "TOFI_T2TFLv", "TOFI_T2TFL"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T2TCTM, UINT32, "TOFI_T2TCTM", 24);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T2TCTMI, UINT32, "TOFI_T2TCTMI", "TOFI_T2TCTM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T2TCTME, UINT32, "TOFI_T2TCTME", "TOFI_T2TCTM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T2TCT, UINT32, "TOFI_T2TCT", 1536);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, TOFI_T2TCTv, UINT32, "TOFI_T2TCTv", "TOFI_T2TCT"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T2TFTM, UINT32, "TOFI_T2TFTM", 24);           \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T2TFTMI, UINT32, "TOFI_T2TFTMI", "TOFI_T2TFTM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, TOFI_T2TFTME, UINT32, "TOFI_T2TFTME", "TOFI_T2TFTM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, TOFI_T2TFT, UINT32, "TOFI_T2TFT", 1536);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, TOFI_T2TFTv, UINT32, "TOFI_T2TFTv", "TOFI_T2TFT"); \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_TOFI_EXT_H101_TOFI_H__*/

/*******************************************************/
