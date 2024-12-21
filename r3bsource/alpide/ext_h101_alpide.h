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

#ifndef __GUARD_H101_ALPIDE_EXT_H101_ALPIDE_H__
#define __GUARD_H101_ALPIDE_EXT_H101_ALPIDE_H__

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

typedef struct EXT_STR_h101_ALPIDE_t
{
    /* RAW */
    uint32_t ALPIDE1T_HI /* [-1,-1] */;
    uint32_t ALPIDE1T_LO /* [-1,-1] */;
    uint32_t ALPIDE1PRECHIP /* [0,2000] */;
    uint32_t ALPIDE1PRECHIPv[2000 EXT_STRUCT_CTRL(ALPIDE1PRECHIP)] /* [0,255] */;
    uint32_t ALPIDE1REGION /* [0,2000] */;
    uint32_t ALPIDE1REGIONv[2000 EXT_STRUCT_CTRL(ALPIDE1REGION)] /* [0,255] */;
    uint32_t ALPIDE1ENCODER_ID /* [0,2000] */;
    uint32_t ALPIDE1ENCODER_IDv[2000 EXT_STRUCT_CTRL(ALPIDE1ENCODER_ID)] /* [0,255] */;
    uint32_t ALPIDE1ADDRESS /* [0,2000] */;
    uint32_t ALPIDE1ADDRESSv[2000 EXT_STRUCT_CTRL(ALPIDE1ADDRESS)] /* [0,65535] */;
    uint32_t ALPIDE1HITMAP /* [0,2000] */;
    uint32_t ALPIDE1HITMAPv[2000 EXT_STRUCT_CTRL(ALPIDE1HITMAP)] /* [0,255] */;
    uint32_t ALPIDE1CHIP /* [0,2000] */;
    uint32_t ALPIDE1CHIPv[2000 EXT_STRUCT_CTRL(ALPIDE1CHIP)] /* [0,255] */;
    uint32_t ALPIDE1ROW /* [0,2000] */;
    uint32_t ALPIDE1ROWv[2000 EXT_STRUCT_CTRL(ALPIDE1ROW)] /* [0,65535] */;
    uint32_t ALPIDE1COL /* [0,2000] */;
    uint32_t ALPIDE1COLv[2000 EXT_STRUCT_CTRL(ALPIDE1COL)] /* [0,65535] */;
    uint32_t ALPIDE2T_HI /* [-1,-1] */;
    uint32_t ALPIDE2T_LO /* [-1,-1] */;
    uint32_t ALPIDE2PRECHIP /* [0,2000] */;
    uint32_t ALPIDE2PRECHIPv[2000 EXT_STRUCT_CTRL(ALPIDE2PRECHIP)] /* [0,255] */;
    uint32_t ALPIDE2REGION /* [0,2000] */;
    uint32_t ALPIDE2REGIONv[2000 EXT_STRUCT_CTRL(ALPIDE2REGION)] /* [0,255] */;
    uint32_t ALPIDE2ENCODER_ID /* [0,2000] */;
    uint32_t ALPIDE2ENCODER_IDv[2000 EXT_STRUCT_CTRL(ALPIDE2ENCODER_ID)] /* [0,255] */;
    uint32_t ALPIDE2ADDRESS /* [0,2000] */;
    uint32_t ALPIDE2ADDRESSv[2000 EXT_STRUCT_CTRL(ALPIDE2ADDRESS)] /* [0,65535] */;
    uint32_t ALPIDE2HITMAP /* [0,2000] */;
    uint32_t ALPIDE2HITMAPv[2000 EXT_STRUCT_CTRL(ALPIDE2HITMAP)] /* [0,255] */;
    uint32_t ALPIDE2CHIP /* [0,2000] */;
    uint32_t ALPIDE2CHIPv[2000 EXT_STRUCT_CTRL(ALPIDE2CHIP)] /* [0,255] */;
    uint32_t ALPIDE2ROW /* [0,2000] */;
    uint32_t ALPIDE2ROWv[2000 EXT_STRUCT_CTRL(ALPIDE2ROW)] /* [0,65535] */;
    uint32_t ALPIDE2COL /* [0,2000] */;
    uint32_t ALPIDE2COLv[2000 EXT_STRUCT_CTRL(ALPIDE2COL)] /* [0,65535] */;
    uint32_t ALPIDE3T_HI /* [-1,-1] */;
    uint32_t ALPIDE3T_LO /* [-1,-1] */;
    uint32_t ALPIDE3PRECHIP /* [0,2000] */;
    uint32_t ALPIDE3PRECHIPv[2000 EXT_STRUCT_CTRL(ALPIDE3PRECHIP)] /* [0,255] */;
    uint32_t ALPIDE3REGION /* [0,2000] */;
    uint32_t ALPIDE3REGIONv[2000 EXT_STRUCT_CTRL(ALPIDE3REGION)] /* [0,255] */;
    uint32_t ALPIDE3ENCODER_ID /* [0,2000] */;
    uint32_t ALPIDE3ENCODER_IDv[2000 EXT_STRUCT_CTRL(ALPIDE3ENCODER_ID)] /* [0,255] */;
    uint32_t ALPIDE3ADDRESS /* [0,2000] */;
    uint32_t ALPIDE3ADDRESSv[2000 EXT_STRUCT_CTRL(ALPIDE3ADDRESS)] /* [0,65535] */;
    uint32_t ALPIDE3HITMAP /* [0,2000] */;
    uint32_t ALPIDE3HITMAPv[2000 EXT_STRUCT_CTRL(ALPIDE3HITMAP)] /* [0,255] */;
    uint32_t ALPIDE3CHIP /* [0,2000] */;
    uint32_t ALPIDE3CHIPv[2000 EXT_STRUCT_CTRL(ALPIDE3CHIP)] /* [0,255] */;
    uint32_t ALPIDE3ROW /* [0,2000] */;
    uint32_t ALPIDE3ROWv[2000 EXT_STRUCT_CTRL(ALPIDE3ROW)] /* [0,65535] */;
    uint32_t ALPIDE3COL /* [0,2000] */;
    uint32_t ALPIDE3COLv[2000 EXT_STRUCT_CTRL(ALPIDE3COL)] /* [0,65535] */;
    uint32_t ALPIDE4T_HI /* [-1,-1] */;
    uint32_t ALPIDE4T_LO /* [-1,-1] */;
    uint32_t ALPIDE4PRECHIP /* [0,2000] */;
    uint32_t ALPIDE4PRECHIPv[2000 EXT_STRUCT_CTRL(ALPIDE4PRECHIP)] /* [0,255] */;
    uint32_t ALPIDE4REGION /* [0,2000] */;
    uint32_t ALPIDE4REGIONv[2000 EXT_STRUCT_CTRL(ALPIDE4REGION)] /* [0,255] */;
    uint32_t ALPIDE4ENCODER_ID /* [0,2000] */;
    uint32_t ALPIDE4ENCODER_IDv[2000 EXT_STRUCT_CTRL(ALPIDE4ENCODER_ID)] /* [0,255] */;
    uint32_t ALPIDE4ADDRESS /* [0,2000] */;
    uint32_t ALPIDE4ADDRESSv[2000 EXT_STRUCT_CTRL(ALPIDE4ADDRESS)] /* [0,65535] */;
    uint32_t ALPIDE4HITMAP /* [0,2000] */;
    uint32_t ALPIDE4HITMAPv[2000 EXT_STRUCT_CTRL(ALPIDE4HITMAP)] /* [0,255] */;
    uint32_t ALPIDE4CHIP /* [0,2000] */;
    uint32_t ALPIDE4CHIPv[2000 EXT_STRUCT_CTRL(ALPIDE4CHIP)] /* [0,255] */;
    uint32_t ALPIDE4ROW /* [0,2000] */;
    uint32_t ALPIDE4ROWv[2000 EXT_STRUCT_CTRL(ALPIDE4ROW)] /* [0,65535] */;
    uint32_t ALPIDE4COL /* [0,2000] */;
    uint32_t ALPIDE4COLv[2000 EXT_STRUCT_CTRL(ALPIDE4COL)] /* [0,65535] */;
    uint32_t ALPIDE5T_HI /* [-1,-1] */;
    uint32_t ALPIDE5T_LO /* [-1,-1] */;
    uint32_t ALPIDE5PRECHIP /* [0,2000] */;
    uint32_t ALPIDE5PRECHIPv[2000 EXT_STRUCT_CTRL(ALPIDE5PRECHIP)] /* [0,255] */;
    uint32_t ALPIDE5REGION /* [0,2000] */;
    uint32_t ALPIDE5REGIONv[2000 EXT_STRUCT_CTRL(ALPIDE5REGION)] /* [0,255] */;
    uint32_t ALPIDE5ENCODER_ID /* [0,2000] */;
    uint32_t ALPIDE5ENCODER_IDv[2000 EXT_STRUCT_CTRL(ALPIDE5ENCODER_ID)] /* [0,255] */;
    uint32_t ALPIDE5ADDRESS /* [0,2000] */;
    uint32_t ALPIDE5ADDRESSv[2000 EXT_STRUCT_CTRL(ALPIDE5ADDRESS)] /* [0,65535] */;
    uint32_t ALPIDE5HITMAP /* [0,2000] */;
    uint32_t ALPIDE5HITMAPv[2000 EXT_STRUCT_CTRL(ALPIDE5HITMAP)] /* [0,255] */;
    uint32_t ALPIDE5CHIP /* [0,2000] */;
    uint32_t ALPIDE5CHIPv[2000 EXT_STRUCT_CTRL(ALPIDE5CHIP)] /* [0,255] */;
    uint32_t ALPIDE5ROW /* [0,2000] */;
    uint32_t ALPIDE5ROWv[2000 EXT_STRUCT_CTRL(ALPIDE5ROW)] /* [0,65535] */;
    uint32_t ALPIDE5COL /* [0,2000] */;
    uint32_t ALPIDE5COLv[2000 EXT_STRUCT_CTRL(ALPIDE5COL)] /* [0,65535] */;
    uint32_t ALPIDE6T_HI /* [-1,-1] */;
    uint32_t ALPIDE6T_LO /* [-1,-1] */;
    uint32_t ALPIDE6PRECHIP /* [0,2000] */;
    uint32_t ALPIDE6PRECHIPv[2000 EXT_STRUCT_CTRL(ALPIDE6PRECHIP)] /* [0,255] */;
    uint32_t ALPIDE6REGION /* [0,2000] */;
    uint32_t ALPIDE6REGIONv[2000 EXT_STRUCT_CTRL(ALPIDE6REGION)] /* [0,255] */;
    uint32_t ALPIDE6ENCODER_ID /* [0,2000] */;
    uint32_t ALPIDE6ENCODER_IDv[2000 EXT_STRUCT_CTRL(ALPIDE6ENCODER_ID)] /* [0,255] */;
    uint32_t ALPIDE6ADDRESS /* [0,2000] */;
    uint32_t ALPIDE6ADDRESSv[2000 EXT_STRUCT_CTRL(ALPIDE6ADDRESS)] /* [0,65535] */;
    uint32_t ALPIDE6HITMAP /* [0,2000] */;
    uint32_t ALPIDE6HITMAPv[2000 EXT_STRUCT_CTRL(ALPIDE6HITMAP)] /* [0,255] */;
    uint32_t ALPIDE6CHIP /* [0,2000] */;
    uint32_t ALPIDE6CHIPv[2000 EXT_STRUCT_CTRL(ALPIDE6CHIP)] /* [0,255] */;
    uint32_t ALPIDE6ROW /* [0,2000] */;
    uint32_t ALPIDE6ROWv[2000 EXT_STRUCT_CTRL(ALPIDE6ROW)] /* [0,65535] */;
    uint32_t ALPIDE6COL /* [0,2000] */;
    uint32_t ALPIDE6COLv[2000 EXT_STRUCT_CTRL(ALPIDE6COL)] /* [0,65535] */;

} EXT_STR_h101_ALPIDE;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_ALPIDE_onion_t
{
    /* RAW */
    struct
    {
        uint32_t T_HI;
        uint32_t T_LO;
        uint32_t PRECHIP;
        uint32_t PRECHIPv[2000 /* PRECHIP */];
        uint32_t REGION;
        uint32_t REGIONv[2000 /* REGION */];
        uint32_t ENCODER_ID;
        uint32_t ENCODER_IDv[2000 /* ENCODER_ID */];
        uint32_t ADDRESS;
        uint32_t ADDRESSv[2000 /* ADDRESS */];
        uint32_t HITMAP;
        uint32_t HITMAPv[2000 /* HITMAP */];
        uint32_t CHIP;
        uint32_t CHIPv[2000 /* CHIP */];
        uint32_t ROW;
        uint32_t ROWv[2000 /* ROW */];
        uint32_t COL;
        uint32_t COLv[2000 /* COL */];
    } ALPIDE[6];

} EXT_STR_h101_ALPIDE_onion;

/*******************************************************/

#define EXT_STR_h101_ALPIDE_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                           \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE1T_HI, UINT32, "ALPIDE1T_HI");                   \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE1T_LO, UINT32, "ALPIDE1T_LO");                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE1PRECHIP, UINT32, "ALPIDE1PRECHIP", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE1PRECHIPv, UINT32, "ALPIDE1PRECHIPv", "ALPIDE1PRECHIP");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE1REGION, UINT32, "ALPIDE1REGION", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE1REGIONv, UINT32, "ALPIDE1REGIONv", "ALPIDE1REGION");          \
        EXT_STR_ITEM_INFO_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE1ENCODER_ID, UINT32, "ALPIDE1ENCODER_ID", 2000);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                    \
                              si,                                                                                    \
                              offset,                                                                                \
                              struct_t,                                                                              \
                              printerr,                                                                              \
                              ALPIDE1ENCODER_IDv,                                                                    \
                              UINT32,                                                                                \
                              "ALPIDE1ENCODER_IDv",                                                                  \
                              "ALPIDE1ENCODER_ID");                                                                  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE1ADDRESS, UINT32, "ALPIDE1ADDRESS", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE1ADDRESSv, UINT32, "ALPIDE1ADDRESSv", "ALPIDE1ADDRESS");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE1HITMAP, UINT32, "ALPIDE1HITMAP", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE1HITMAPv, UINT32, "ALPIDE1HITMAPv", "ALPIDE1HITMAP");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE1CHIP, UINT32, "ALPIDE1CHIP", 2000);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE1CHIPv, UINT32, "ALPIDE1CHIPv", "ALPIDE1CHIP");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE1ROW, UINT32, "ALPIDE1ROW", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE1ROWv, UINT32, "ALPIDE1ROWv", "ALPIDE1ROW"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE1COL, UINT32, "ALPIDE1COL", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE1COLv, UINT32, "ALPIDE1COLv", "ALPIDE1COL"); \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE2T_HI, UINT32, "ALPIDE2T_HI");                   \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE2T_LO, UINT32, "ALPIDE2T_LO");                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE2PRECHIP, UINT32, "ALPIDE2PRECHIP", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE2PRECHIPv, UINT32, "ALPIDE2PRECHIPv", "ALPIDE2PRECHIP");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE2REGION, UINT32, "ALPIDE2REGION", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE2REGIONv, UINT32, "ALPIDE2REGIONv", "ALPIDE2REGION");          \
        EXT_STR_ITEM_INFO_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE2ENCODER_ID, UINT32, "ALPIDE2ENCODER_ID", 2000);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                    \
                              si,                                                                                    \
                              offset,                                                                                \
                              struct_t,                                                                              \
                              printerr,                                                                              \
                              ALPIDE2ENCODER_IDv,                                                                    \
                              UINT32,                                                                                \
                              "ALPIDE2ENCODER_IDv",                                                                  \
                              "ALPIDE2ENCODER_ID");                                                                  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE2ADDRESS, UINT32, "ALPIDE2ADDRESS", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE2ADDRESSv, UINT32, "ALPIDE2ADDRESSv", "ALPIDE2ADDRESS");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE2HITMAP, UINT32, "ALPIDE2HITMAP", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE2HITMAPv, UINT32, "ALPIDE2HITMAPv", "ALPIDE2HITMAP");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE2CHIP, UINT32, "ALPIDE2CHIP", 2000);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE2CHIPv, UINT32, "ALPIDE2CHIPv", "ALPIDE2CHIP");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE2ROW, UINT32, "ALPIDE2ROW", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE2ROWv, UINT32, "ALPIDE2ROWv", "ALPIDE2ROW"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE2COL, UINT32, "ALPIDE2COL", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE2COLv, UINT32, "ALPIDE2COLv", "ALPIDE2COL"); \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE3T_HI, UINT32, "ALPIDE3T_HI");                   \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE3T_LO, UINT32, "ALPIDE3T_LO");                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE3PRECHIP, UINT32, "ALPIDE3PRECHIP", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE3PRECHIPv, UINT32, "ALPIDE3PRECHIPv", "ALPIDE3PRECHIP");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE3REGION, UINT32, "ALPIDE3REGION", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE3REGIONv, UINT32, "ALPIDE3REGIONv", "ALPIDE3REGION");          \
        EXT_STR_ITEM_INFO_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE3ENCODER_ID, UINT32, "ALPIDE3ENCODER_ID", 2000);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                    \
                              si,                                                                                    \
                              offset,                                                                                \
                              struct_t,                                                                              \
                              printerr,                                                                              \
                              ALPIDE3ENCODER_IDv,                                                                    \
                              UINT32,                                                                                \
                              "ALPIDE3ENCODER_IDv",                                                                  \
                              "ALPIDE3ENCODER_ID");                                                                  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE3ADDRESS, UINT32, "ALPIDE3ADDRESS", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE3ADDRESSv, UINT32, "ALPIDE3ADDRESSv", "ALPIDE3ADDRESS");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE3HITMAP, UINT32, "ALPIDE3HITMAP", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE3HITMAPv, UINT32, "ALPIDE3HITMAPv", "ALPIDE3HITMAP");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE3CHIP, UINT32, "ALPIDE3CHIP", 2000);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE3CHIPv, UINT32, "ALPIDE3CHIPv", "ALPIDE3CHIP");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE3ROW, UINT32, "ALPIDE3ROW", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE3ROWv, UINT32, "ALPIDE3ROWv", "ALPIDE3ROW"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE3COL, UINT32, "ALPIDE3COL", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE3COLv, UINT32, "ALPIDE3COLv", "ALPIDE3COL"); \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE4T_HI, UINT32, "ALPIDE4T_HI");                   \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE4T_LO, UINT32, "ALPIDE4T_LO");                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE4PRECHIP, UINT32, "ALPIDE4PRECHIP", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE4PRECHIPv, UINT32, "ALPIDE4PRECHIPv", "ALPIDE4PRECHIP");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE4REGION, UINT32, "ALPIDE4REGION", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE4REGIONv, UINT32, "ALPIDE4REGIONv", "ALPIDE4REGION");          \
        EXT_STR_ITEM_INFO_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE4ENCODER_ID, UINT32, "ALPIDE4ENCODER_ID", 2000);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                    \
                              si,                                                                                    \
                              offset,                                                                                \
                              struct_t,                                                                              \
                              printerr,                                                                              \
                              ALPIDE4ENCODER_IDv,                                                                    \
                              UINT32,                                                                                \
                              "ALPIDE4ENCODER_IDv",                                                                  \
                              "ALPIDE4ENCODER_ID");                                                                  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE4ADDRESS, UINT32, "ALPIDE4ADDRESS", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE4ADDRESSv, UINT32, "ALPIDE4ADDRESSv", "ALPIDE4ADDRESS");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE4HITMAP, UINT32, "ALPIDE4HITMAP", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE4HITMAPv, UINT32, "ALPIDE4HITMAPv", "ALPIDE4HITMAP");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE4CHIP, UINT32, "ALPIDE4CHIP", 2000);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE4CHIPv, UINT32, "ALPIDE4CHIPv", "ALPIDE4CHIP");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE4ROW, UINT32, "ALPIDE4ROW", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE4ROWv, UINT32, "ALPIDE4ROWv", "ALPIDE4ROW"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE4COL, UINT32, "ALPIDE4COL", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE4COLv, UINT32, "ALPIDE4COLv", "ALPIDE4COL"); \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE5T_HI, UINT32, "ALPIDE5T_HI");                   \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE5T_LO, UINT32, "ALPIDE5T_LO");                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE5PRECHIP, UINT32, "ALPIDE5PRECHIP", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE5PRECHIPv, UINT32, "ALPIDE5PRECHIPv", "ALPIDE5PRECHIP");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE5REGION, UINT32, "ALPIDE5REGION", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE5REGIONv, UINT32, "ALPIDE5REGIONv", "ALPIDE5REGION");          \
        EXT_STR_ITEM_INFO_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE5ENCODER_ID, UINT32, "ALPIDE5ENCODER_ID", 2000);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                    \
                              si,                                                                                    \
                              offset,                                                                                \
                              struct_t,                                                                              \
                              printerr,                                                                              \
                              ALPIDE5ENCODER_IDv,                                                                    \
                              UINT32,                                                                                \
                              "ALPIDE5ENCODER_IDv",                                                                  \
                              "ALPIDE5ENCODER_ID");                                                                  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE5ADDRESS, UINT32, "ALPIDE5ADDRESS", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE5ADDRESSv, UINT32, "ALPIDE5ADDRESSv", "ALPIDE5ADDRESS");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE5HITMAP, UINT32, "ALPIDE5HITMAP", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE5HITMAPv, UINT32, "ALPIDE5HITMAPv", "ALPIDE5HITMAP");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE5CHIP, UINT32, "ALPIDE5CHIP", 2000);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE5CHIPv, UINT32, "ALPIDE5CHIPv", "ALPIDE5CHIP");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE5ROW, UINT32, "ALPIDE5ROW", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE5ROWv, UINT32, "ALPIDE5ROWv", "ALPIDE5ROW"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE5COL, UINT32, "ALPIDE5COL", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE5COLv, UINT32, "ALPIDE5COLv", "ALPIDE5COL"); \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE6T_HI, UINT32, "ALPIDE6T_HI");                   \
        EXT_STR_ITEM_INFO(ok, si, offset, struct_t, printerr, ALPIDE6T_LO, UINT32, "ALPIDE6T_LO");                   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE6PRECHIP, UINT32, "ALPIDE6PRECHIP", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE6PRECHIPv, UINT32, "ALPIDE6PRECHIPv", "ALPIDE6PRECHIP");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE6REGION, UINT32, "ALPIDE6REGION", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE6REGIONv, UINT32, "ALPIDE6REGIONv", "ALPIDE6REGION");          \
        EXT_STR_ITEM_INFO_LIM(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE6ENCODER_ID, UINT32, "ALPIDE6ENCODER_ID", 2000);               \
        EXT_STR_ITEM_INFO_ZZP(ok,                                                                                    \
                              si,                                                                                    \
                              offset,                                                                                \
                              struct_t,                                                                              \
                              printerr,                                                                              \
                              ALPIDE6ENCODER_IDv,                                                                    \
                              UINT32,                                                                                \
                              "ALPIDE6ENCODER_IDv",                                                                  \
                              "ALPIDE6ENCODER_ID");                                                                  \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE6ADDRESS, UINT32, "ALPIDE6ADDRESS", 2000);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE6ADDRESSv, UINT32, "ALPIDE6ADDRESSv", "ALPIDE6ADDRESS");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE6HITMAP, UINT32, "ALPIDE6HITMAP", 2000);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE6HITMAPv, UINT32, "ALPIDE6HITMAPv", "ALPIDE6HITMAP");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE6CHIP, UINT32, "ALPIDE6CHIP", 2000);         \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, ALPIDE6CHIPv, UINT32, "ALPIDE6CHIPv", "ALPIDE6CHIP");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE6ROW, UINT32, "ALPIDE6ROW", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE6ROWv, UINT32, "ALPIDE6ROWv", "ALPIDE6ROW"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ALPIDE6COL, UINT32, "ALPIDE6COL", 2000);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ALPIDE6COLv, UINT32, "ALPIDE6COLv", "ALPIDE6COL"); \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_ALPIDE_EXT_H101_ALPIDE_H__*/

/*******************************************************/
