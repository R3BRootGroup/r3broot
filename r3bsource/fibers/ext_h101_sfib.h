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

#ifndef __GUARD_H101_SFIB_EXT_H101_SFIB_H__
#define __GUARD_H101_SFIB_EXT_H101_SFIB_H__

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

typedef struct EXT_STR_h101_sfib_t
{
    /* RAW */
    uint32_t SFIB_TTLCM /* [1,256] */;
    uint32_t SFIB_TTLCMI[256 EXT_STRUCT_CTRL(SFIB_TTLCM)] /* [1,256] */;
    uint32_t SFIB_TTLCME[256 EXT_STRUCT_CTRL(SFIB_TTLCM)] /* [1,16384] */;
    uint32_t SFIB_TTLC /* [0,16384] */;
    uint32_t SFIB_TTLCv[16384 EXT_STRUCT_CTRL(SFIB_TTLC)] /* [0,65535] */;
    uint32_t SFIB_BTLCM /* [1,256] */;
    uint32_t SFIB_BTLCMI[256 EXT_STRUCT_CTRL(SFIB_BTLCM)] /* [1,256] */;
    uint32_t SFIB_BTLCME[256 EXT_STRUCT_CTRL(SFIB_BTLCM)] /* [1,16384] */;
    uint32_t SFIB_BTLC /* [0,16384] */;
    uint32_t SFIB_BTLCv[16384 EXT_STRUCT_CTRL(SFIB_BTLC)] /* [0,65535] */;
    uint32_t SFIB_TTLFM /* [1,256] */;
    uint32_t SFIB_TTLFMI[256 EXT_STRUCT_CTRL(SFIB_TTLFM)] /* [1,256] */;
    uint32_t SFIB_TTLFME[256 EXT_STRUCT_CTRL(SFIB_TTLFM)] /* [1,16384] */;
    uint32_t SFIB_TTLF /* [0,16384] */;
    uint32_t SFIB_TTLFv[16384 EXT_STRUCT_CTRL(SFIB_TTLF)] /* [0,65535] */;
    uint32_t SFIB_BTLFM /* [1,256] */;
    uint32_t SFIB_BTLFMI[256 EXT_STRUCT_CTRL(SFIB_BTLFM)] /* [1,256] */;
    uint32_t SFIB_BTLFME[256 EXT_STRUCT_CTRL(SFIB_BTLFM)] /* [1,16384] */;
    uint32_t SFIB_BTLF /* [0,16384] */;
    uint32_t SFIB_BTLFv[16384 EXT_STRUCT_CTRL(SFIB_BTLF)] /* [0,65535] */;
    uint32_t SFIB_TTTCM /* [1,256] */;
    uint32_t SFIB_TTTCMI[256 EXT_STRUCT_CTRL(SFIB_TTTCM)] /* [1,256] */;
    uint32_t SFIB_TTTCME[256 EXT_STRUCT_CTRL(SFIB_TTTCM)] /* [1,16384] */;
    uint32_t SFIB_TTTC /* [0,16384] */;
    uint32_t SFIB_TTTCv[16384 EXT_STRUCT_CTRL(SFIB_TTTC)] /* [0,65535] */;
    uint32_t SFIB_BTTCM /* [1,256] */;
    uint32_t SFIB_BTTCMI[256 EXT_STRUCT_CTRL(SFIB_BTTCM)] /* [1,256] */;
    uint32_t SFIB_BTTCME[256 EXT_STRUCT_CTRL(SFIB_BTTCM)] /* [1,16384] */;
    uint32_t SFIB_BTTC /* [0,16384] */;
    uint32_t SFIB_BTTCv[16384 EXT_STRUCT_CTRL(SFIB_BTTC)] /* [0,65535] */;
    uint32_t SFIB_TTTFM /* [1,256] */;
    uint32_t SFIB_TTTFMI[256 EXT_STRUCT_CTRL(SFIB_TTTFM)] /* [1,256] */;
    uint32_t SFIB_TTTFME[256 EXT_STRUCT_CTRL(SFIB_TTTFM)] /* [1,16384] */;
    uint32_t SFIB_TTTF /* [0,16384] */;
    uint32_t SFIB_TTTFv[16384 EXT_STRUCT_CTRL(SFIB_TTTF)] /* [0,65535] */;
    uint32_t SFIB_BTTFM /* [1,256] */;
    uint32_t SFIB_BTTFMI[256 EXT_STRUCT_CTRL(SFIB_BTTFM)] /* [1,256] */;
    uint32_t SFIB_BTTFME[256 EXT_STRUCT_CTRL(SFIB_BTTFM)] /* [1,16384] */;
    uint32_t SFIB_BTTF /* [0,16384] */;
    uint32_t SFIB_BTTFv[16384 EXT_STRUCT_CTRL(SFIB_BTTF)] /* [0,65535] */;

} EXT_STR_h101_sfib;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_sfib_onion_t
{
    /* RAW */
    uint32_t SFIB_TTLCM;
    uint32_t SFIB_TTLCMI[256 /* SFIB_TTLCM */];
    uint32_t SFIB_TTLCME[256 /* SFIB_TTLCM */];
    uint32_t SFIB_TTLC;
    uint32_t SFIB_TTLCv[16384 /* SFIB_TTLC */];
    uint32_t SFIB_BTLCM;
    uint32_t SFIB_BTLCMI[256 /* SFIB_BTLCM */];
    uint32_t SFIB_BTLCME[256 /* SFIB_BTLCM */];
    uint32_t SFIB_BTLC;
    uint32_t SFIB_BTLCv[16384 /* SFIB_BTLC */];
    uint32_t SFIB_TTLFM;
    uint32_t SFIB_TTLFMI[256 /* SFIB_TTLFM */];
    uint32_t SFIB_TTLFME[256 /* SFIB_TTLFM */];
    uint32_t SFIB_TTLF;
    uint32_t SFIB_TTLFv[16384 /* SFIB_TTLF */];
    uint32_t SFIB_BTLFM;
    uint32_t SFIB_BTLFMI[256 /* SFIB_BTLFM */];
    uint32_t SFIB_BTLFME[256 /* SFIB_BTLFM */];
    uint32_t SFIB_BTLF;
    uint32_t SFIB_BTLFv[16384 /* SFIB_BTLF */];
    uint32_t SFIB_TTTCM;
    uint32_t SFIB_TTTCMI[256 /* SFIB_TTTCM */];
    uint32_t SFIB_TTTCME[256 /* SFIB_TTTCM */];
    uint32_t SFIB_TTTC;
    uint32_t SFIB_TTTCv[16384 /* SFIB_TTTC */];
    uint32_t SFIB_BTTCM;
    uint32_t SFIB_BTTCMI[256 /* SFIB_BTTCM */];
    uint32_t SFIB_BTTCME[256 /* SFIB_BTTCM */];
    uint32_t SFIB_BTTC;
    uint32_t SFIB_BTTCv[16384 /* SFIB_BTTC */];
    uint32_t SFIB_TTTFM;
    uint32_t SFIB_TTTFMI[256 /* SFIB_TTTFM */];
    uint32_t SFIB_TTTFME[256 /* SFIB_TTTFM */];
    uint32_t SFIB_TTTF;
    uint32_t SFIB_TTTFv[16384 /* SFIB_TTTF */];
    uint32_t SFIB_BTTFM;
    uint32_t SFIB_BTTFMI[256 /* SFIB_BTTFM */];
    uint32_t SFIB_BTTFME[256 /* SFIB_BTTFM */];
    uint32_t SFIB_BTTF;
    uint32_t SFIB_BTTFv[16384 /* SFIB_BTTF */];

} EXT_STR_h101_sfib_onion;

/*******************************************************/

#define EXT_STR_h101_sfib_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                             \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_TTLCM, UINT32, "SFIB_TTLCM", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTLCMI, UINT32, "SFIB_TTLCMI", "SFIB_TTLCM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTLCME, UINT32, "SFIB_TTLCME", "SFIB_TTLCM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_TTLC, UINT32, "SFIB_TTLC", 16384);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTLCv, UINT32, "SFIB_TTLCv", "SFIB_TTLC");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_BTLCM, UINT32, "SFIB_BTLCM", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTLCMI, UINT32, "SFIB_BTLCMI", "SFIB_BTLCM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTLCME, UINT32, "SFIB_BTLCME", "SFIB_BTLCM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_BTLC, UINT32, "SFIB_BTLC", 16384);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTLCv, UINT32, "SFIB_BTLCv", "SFIB_BTLC");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_TTLFM, UINT32, "SFIB_TTLFM", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTLFMI, UINT32, "SFIB_TTLFMI", "SFIB_TTLFM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTLFME, UINT32, "SFIB_TTLFME", "SFIB_TTLFM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_TTLF, UINT32, "SFIB_TTLF", 16384);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTLFv, UINT32, "SFIB_TTLFv", "SFIB_TTLF");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_BTLFM, UINT32, "SFIB_BTLFM", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTLFMI, UINT32, "SFIB_BTLFMI", "SFIB_BTLFM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTLFME, UINT32, "SFIB_BTLFME", "SFIB_BTLFM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_BTLF, UINT32, "SFIB_BTLF", 16384);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTLFv, UINT32, "SFIB_BTLFv", "SFIB_BTLF");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_TTTCM, UINT32, "SFIB_TTTCM", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTTCMI, UINT32, "SFIB_TTTCMI", "SFIB_TTTCM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTTCME, UINT32, "SFIB_TTTCME", "SFIB_TTTCM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_TTTC, UINT32, "SFIB_TTTC", 16384);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTTCv, UINT32, "SFIB_TTTCv", "SFIB_TTTC");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_BTTCM, UINT32, "SFIB_BTTCM", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTTCMI, UINT32, "SFIB_BTTCMI", "SFIB_BTTCM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTTCME, UINT32, "SFIB_BTTCME", "SFIB_BTTCM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_BTTC, UINT32, "SFIB_BTTC", 16384);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTTCv, UINT32, "SFIB_BTTCv", "SFIB_BTTC");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_TTTFM, UINT32, "SFIB_TTTFM", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTTFMI, UINT32, "SFIB_TTTFMI", "SFIB_TTTFM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTTFME, UINT32, "SFIB_TTTFME", "SFIB_TTTFM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_TTTF, UINT32, "SFIB_TTTF", 16384);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_TTTFv, UINT32, "SFIB_TTTFv", "SFIB_TTTF");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_BTTFM, UINT32, "SFIB_BTTFM", 256);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTTFMI, UINT32, "SFIB_BTTFMI", "SFIB_BTTFM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTTFME, UINT32, "SFIB_BTTFME", "SFIB_BTTFM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, SFIB_BTTF, UINT32, "SFIB_BTTF", 16384);            \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, SFIB_BTTFv, UINT32, "SFIB_BTTFv", "SFIB_BTTF");    \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_SFIB_EXT_H101_SFIB_H__*/

/*******************************************************/
