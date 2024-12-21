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

#ifndef __GUARD_H101_EXT_H101_FRSSCI_H__
#define __GUARD_H101_EXT_H101_FRSSCI_H__

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

typedef struct EXT_STR_h101_FRSSCI_t
{
    /* RAW */
    uint32_t FRSSCI1VTFM /* [1,3] */;
    uint32_t FRSSCI1VTFMI[3 EXT_STRUCT_CTRL(FRSSCI1VTFM)] /* [1,3] */;
    uint32_t FRSSCI1VTFME[3 EXT_STRUCT_CTRL(FRSSCI1VTFM)] /* [1,384] */;
    uint32_t FRSSCI1VTF /* [0,384] */;
    uint32_t FRSSCI1VTFv[384 EXT_STRUCT_CTRL(FRSSCI1VTF)] /* [0,65535] */;
    uint32_t FRSSCI1VTCM /* [1,3] */;
    uint32_t FRSSCI1VTCMI[3 EXT_STRUCT_CTRL(FRSSCI1VTCM)] /* [1,3] */;
    uint32_t FRSSCI1VTCME[3 EXT_STRUCT_CTRL(FRSSCI1VTCM)] /* [1,384] */;
    uint32_t FRSSCI1VTC /* [0,384] */;
    uint32_t FRSSCI1VTCv[384 EXT_STRUCT_CTRL(FRSSCI1VTC)] /* [0,65535] */;
    uint32_t FRSSCI2VTFM /* [1,3] */;
    uint32_t FRSSCI2VTFMI[3 EXT_STRUCT_CTRL(FRSSCI2VTFM)] /* [1,3] */;
    uint32_t FRSSCI2VTFME[3 EXT_STRUCT_CTRL(FRSSCI2VTFM)] /* [1,384] */;
    uint32_t FRSSCI2VTF /* [0,384] */;
    uint32_t FRSSCI2VTFv[384 EXT_STRUCT_CTRL(FRSSCI2VTF)] /* [0,65535] */;
    uint32_t FRSSCI2VTCM /* [1,3] */;
    uint32_t FRSSCI2VTCMI[3 EXT_STRUCT_CTRL(FRSSCI2VTCM)] /* [1,3] */;
    uint32_t FRSSCI2VTCME[3 EXT_STRUCT_CTRL(FRSSCI2VTCM)] /* [1,384] */;
    uint32_t FRSSCI2VTC /* [0,384] */;
    uint32_t FRSSCI2VTCv[384 EXT_STRUCT_CTRL(FRSSCI2VTC)] /* [0,65535] */;
    uint32_t FRSSCI3VTFM /* [1,3] */;
    uint32_t FRSSCI3VTFMI[3 EXT_STRUCT_CTRL(FRSSCI3VTFM)] /* [1,3] */;
    uint32_t FRSSCI3VTFME[3 EXT_STRUCT_CTRL(FRSSCI3VTFM)] /* [1,384] */;
    uint32_t FRSSCI3VTF /* [0,384] */;
    uint32_t FRSSCI3VTFv[384 EXT_STRUCT_CTRL(FRSSCI3VTF)] /* [0,65535] */;
    uint32_t FRSSCI3VTCM /* [1,3] */;
    uint32_t FRSSCI3VTCMI[3 EXT_STRUCT_CTRL(FRSSCI3VTCM)] /* [1,3] */;
    uint32_t FRSSCI3VTCME[3 EXT_STRUCT_CTRL(FRSSCI3VTCM)] /* [1,384] */;
    uint32_t FRSSCI3VTC /* [0,384] */;
    uint32_t FRSSCI3VTCv[384 EXT_STRUCT_CTRL(FRSSCI3VTC)] /* [0,65535] */;

} EXT_STR_h101_FRSSCI;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FRSSCI_onion_t
{
    /* RAW */
    struct
    {
        uint32_t VTFM;
        uint32_t VTFMI[3 /* VTFM */];
        uint32_t VTFME[3 /* VTFM */];
        uint32_t VTF;
        uint32_t VTFv[384 /* VTF */];
        uint32_t VTCM;
        uint32_t VTCMI[3 /* VTCM */];
        uint32_t VTCME[3 /* VTCM */];
        uint32_t VTC;
        uint32_t VTCv[384 /* VTC */];
    } FRSSCI[3];

} EXT_STR_h101_FRSSCI_onion;

/*******************************************************/

#define EXT_STR_h101_FRSSCI_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                         \
    do                                                                                                             \
    {                                                                                                              \
        ok = 1;                                                                                                    \
        /* RAW */                                                                                                  \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTFM, UINT32, "FRSSCI1VTFM", 3, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTFMI, UINT32, "FRSSCI1VTFMI", "FRSSCI1VTFM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTFME, UINT32, "FRSSCI1VTFME", "FRSSCI1VTFM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTF, UINT32, "FRSSCI1VTF", 384, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTFv, UINT32, "FRSSCI1VTFv", "FRSSCI1VTF", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTCM, UINT32, "FRSSCI1VTCM", 3, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTCMI, UINT32, "FRSSCI1VTCMI", "FRSSCI1VTCM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTCME, UINT32, "FRSSCI1VTCME", "FRSSCI1VTCM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTC, UINT32, "FRSSCI1VTC", 384, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI1VTCv, UINT32, "FRSSCI1VTCv", "FRSSCI1VTC", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTFM, UINT32, "FRSSCI2VTFM", 3, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTFMI, UINT32, "FRSSCI2VTFMI", "FRSSCI2VTFM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTFME, UINT32, "FRSSCI2VTFME", "FRSSCI2VTFM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTF, UINT32, "FRSSCI2VTF", 384, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTFv, UINT32, "FRSSCI2VTFv", "FRSSCI2VTF", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTCM, UINT32, "FRSSCI2VTCM", 3, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTCMI, UINT32, "FRSSCI2VTCMI", "FRSSCI2VTCM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTCME, UINT32, "FRSSCI2VTCME", "FRSSCI2VTCM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTC, UINT32, "FRSSCI2VTC", 384, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI2VTCv, UINT32, "FRSSCI2VTCv", "FRSSCI2VTC", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTFM, UINT32, "FRSSCI3VTFM", 3, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTFMI, UINT32, "FRSSCI3VTFMI", "FRSSCI3VTFM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTFME, UINT32, "FRSSCI3VTFME", "FRSSCI3VTFM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTF, UINT32, "FRSSCI3VTF", 384, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTFv, UINT32, "FRSSCI3VTFv", "FRSSCI3VTF", 0 /*flags*/);    \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTCM, UINT32, "FRSSCI3VTCM", 3, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTCMI, UINT32, "FRSSCI3VTCMI", "FRSSCI3VTCM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTCME, UINT32, "FRSSCI3VTCME", "FRSSCI3VTCM", 0 /*flags*/); \
        EXT_STR_ITEM_INFO2_LIM(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTC, UINT32, "FRSSCI3VTC", 384, 0 /*flags*/);               \
        EXT_STR_ITEM_INFO2_ZZP(                                                                                    \
            ok, si, offset, struct_t, printerr, FRSSCI3VTCv, UINT32, "FRSSCI3VTCv", "FRSSCI3VTC", 0 /*flags*/);    \
                                                                                                                   \
    } while (0);

#endif /*__GUARD_H101_EXT_H101_FRSSCI_H__*/

/*******************************************************/
