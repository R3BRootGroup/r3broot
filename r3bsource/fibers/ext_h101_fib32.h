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

#ifndef __GUARD_H101_FIBEZ_EXT_H101_FIB32_H__
#define __GUARD_H101_FIBEZ_EXT_H101_FIB32_H__

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

typedef struct EXT_STR_h101_FIBEZ_t
{
    /* RAW */
    uint32_t FIBEZ_TBLCM /* [1,512] */;
    uint32_t FIBEZ_TBLCMI[512 EXT_STRUCT_CTRL(FIBEZ_TBLCM)] /* [1,512] */;
    uint32_t FIBEZ_TBLCME[512 EXT_STRUCT_CTRL(FIBEZ_TBLCM)] /* [1,32768] */;
    uint32_t FIBEZ_TBLC /* [0,32768] */;
    uint32_t FIBEZ_TBLCv[32768 EXT_STRUCT_CTRL(FIBEZ_TBLC)] /* [0,65535] */;
    uint32_t FIBEZ_TBLFM /* [1,512] */;
    uint32_t FIBEZ_TBLFMI[512 EXT_STRUCT_CTRL(FIBEZ_TBLFM)] /* [1,512] */;
    uint32_t FIBEZ_TBLFME[512 EXT_STRUCT_CTRL(FIBEZ_TBLFM)] /* [1,32768] */;
    uint32_t FIBEZ_TBLF /* [0,32768] */;
    uint32_t FIBEZ_TBLFv[32768 EXT_STRUCT_CTRL(FIBEZ_TBLF)] /* [0,65535] */;
    uint32_t FIBEZ_TBTCM /* [1,512] */;
    uint32_t FIBEZ_TBTCMI[512 EXT_STRUCT_CTRL(FIBEZ_TBTCM)] /* [1,512] */;
    uint32_t FIBEZ_TBTCME[512 EXT_STRUCT_CTRL(FIBEZ_TBTCM)] /* [1,32768] */;
    uint32_t FIBEZ_TBTC /* [0,32768] */;
    uint32_t FIBEZ_TBTCv[32768 EXT_STRUCT_CTRL(FIBEZ_TBTC)] /* [0,65535] */;
    uint32_t FIBEZ_TBTFM /* [1,512] */;
    uint32_t FIBEZ_TBTFMI[512 EXT_STRUCT_CTRL(FIBEZ_TBTFM)] /* [1,512] */;
    uint32_t FIBEZ_TBTFME[512 EXT_STRUCT_CTRL(FIBEZ_TBTFM)] /* [1,32768] */;
    uint32_t FIBEZ_TBTF /* [0,32768] */;
    uint32_t FIBEZ_TBTFv[32768 EXT_STRUCT_CTRL(FIBEZ_TBTF)] /* [0,65535] */;
    uint32_t FIBEZ_TTLCM /* [1,512] */;
    uint32_t FIBEZ_TTLCMI[512 EXT_STRUCT_CTRL(FIBEZ_TTLCM)] /* [1,512] */;
    uint32_t FIBEZ_TTLCME[512 EXT_STRUCT_CTRL(FIBEZ_TTLCM)] /* [1,32768] */;
    uint32_t FIBEZ_TTLC /* [0,32768] */;
    uint32_t FIBEZ_TTLCv[32768 EXT_STRUCT_CTRL(FIBEZ_TTLC)] /* [0,65535] */;
    uint32_t FIBEZ_TTLFM /* [1,512] */;
    uint32_t FIBEZ_TTLFMI[512 EXT_STRUCT_CTRL(FIBEZ_TTLFM)] /* [1,512] */;
    uint32_t FIBEZ_TTLFME[512 EXT_STRUCT_CTRL(FIBEZ_TTLFM)] /* [1,32768] */;
    uint32_t FIBEZ_TTLF /* [0,32768] */;
    uint32_t FIBEZ_TTLFv[32768 EXT_STRUCT_CTRL(FIBEZ_TTLF)] /* [0,65535] */;
    uint32_t FIBEZ_TTTCM /* [1,512] */;
    uint32_t FIBEZ_TTTCMI[512 EXT_STRUCT_CTRL(FIBEZ_TTTCM)] /* [1,512] */;
    uint32_t FIBEZ_TTTCME[512 EXT_STRUCT_CTRL(FIBEZ_TTTCM)] /* [1,32768] */;
    uint32_t FIBEZ_TTTC /* [0,32768] */;
    uint32_t FIBEZ_TTTCv[32768 EXT_STRUCT_CTRL(FIBEZ_TTTC)] /* [0,65535] */;
    uint32_t FIBEZ_TTTFM /* [1,512] */;
    uint32_t FIBEZ_TTTFMI[512 EXT_STRUCT_CTRL(FIBEZ_TTTFM)] /* [1,512] */;
    uint32_t FIBEZ_TTTFME[512 EXT_STRUCT_CTRL(FIBEZ_TTTFM)] /* [1,32768] */;
    uint32_t FIBEZ_TTTF /* [0,32768] */;
    uint32_t FIBEZ_TTTFv[32768 EXT_STRUCT_CTRL(FIBEZ_TTTF)] /* [0,65535] */;
    uint32_t FIBEZ_TRIGCM /* [1,8] */;
    uint32_t FIBEZ_TRIGCMI[8 EXT_STRUCT_CTRL(FIBEZ_TRIGCM)] /* [1,8] */;
    uint32_t FIBEZ_TRIGCME[8 EXT_STRUCT_CTRL(FIBEZ_TRIGCM)] /* [1,8] */;
    uint32_t FIBEZ_TRIGC /* [0,8] */;
    uint32_t FIBEZ_TRIGCv[8 EXT_STRUCT_CTRL(FIBEZ_TRIGC)] /* [0,65535] */;
    uint32_t FIBEZ_TRIGFM /* [1,8] */;
    uint32_t FIBEZ_TRIGFMI[8 EXT_STRUCT_CTRL(FIBEZ_TRIGFM)] /* [1,8] */;
    uint32_t FIBEZ_TRIGFME[8 EXT_STRUCT_CTRL(FIBEZ_TRIGFM)] /* [1,8] */;
    uint32_t FIBEZ_TRIGF /* [0,8] */;
    uint32_t FIBEZ_TRIGFv[8 EXT_STRUCT_CTRL(FIBEZ_TRIGF)] /* [0,65535] */;

} EXT_STR_h101_FIBEZ;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBEZ_onion_t
{
    /* RAW */
    uint32_t FIBEZ_TBLCM;
    uint32_t FIBEZ_TBLCMI[512 /* FIBEZ_TBLCM */];
    uint32_t FIBEZ_TBLCME[512 /* FIBEZ_TBLCM */];
    uint32_t FIBEZ_TBLC;
    uint32_t FIBEZ_TBLCv[32768 /* FIBEZ_TBLC */];
    uint32_t FIBEZ_TBLFM;
    uint32_t FIBEZ_TBLFMI[512 /* FIBEZ_TBLFM */];
    uint32_t FIBEZ_TBLFME[512 /* FIBEZ_TBLFM */];
    uint32_t FIBEZ_TBLF;
    uint32_t FIBEZ_TBLFv[32768 /* FIBEZ_TBLF */];
    uint32_t FIBEZ_TBTCM;
    uint32_t FIBEZ_TBTCMI[512 /* FIBEZ_TBTCM */];
    uint32_t FIBEZ_TBTCME[512 /* FIBEZ_TBTCM */];
    uint32_t FIBEZ_TBTC;
    uint32_t FIBEZ_TBTCv[32768 /* FIBEZ_TBTC */];
    uint32_t FIBEZ_TBTFM;
    uint32_t FIBEZ_TBTFMI[512 /* FIBEZ_TBTFM */];
    uint32_t FIBEZ_TBTFME[512 /* FIBEZ_TBTFM */];
    uint32_t FIBEZ_TBTF;
    uint32_t FIBEZ_TBTFv[32768 /* FIBEZ_TBTF */];
    uint32_t FIBEZ_TTLCM;
    uint32_t FIBEZ_TTLCMI[512 /* FIBEZ_TTLCM */];
    uint32_t FIBEZ_TTLCME[512 /* FIBEZ_TTLCM */];
    uint32_t FIBEZ_TTLC;
    uint32_t FIBEZ_TTLCv[32768 /* FIBEZ_TTLC */];
    uint32_t FIBEZ_TTLFM;
    uint32_t FIBEZ_TTLFMI[512 /* FIBEZ_TTLFM */];
    uint32_t FIBEZ_TTLFME[512 /* FIBEZ_TTLFM */];
    uint32_t FIBEZ_TTLF;
    uint32_t FIBEZ_TTLFv[32768 /* FIBEZ_TTLF */];
    uint32_t FIBEZ_TTTCM;
    uint32_t FIBEZ_TTTCMI[512 /* FIBEZ_TTTCM */];
    uint32_t FIBEZ_TTTCME[512 /* FIBEZ_TTTCM */];
    uint32_t FIBEZ_TTTC;
    uint32_t FIBEZ_TTTCv[32768 /* FIBEZ_TTTC */];
    uint32_t FIBEZ_TTTFM;
    uint32_t FIBEZ_TTTFMI[512 /* FIBEZ_TTTFM */];
    uint32_t FIBEZ_TTTFME[512 /* FIBEZ_TTTFM */];
    uint32_t FIBEZ_TTTF;
    uint32_t FIBEZ_TTTFv[32768 /* FIBEZ_TTTF */];
    uint32_t FIBEZ_TRIGCM;
    uint32_t FIBEZ_TRIGCMI[8 /* FIBEZ_TRIGCM */];
    uint32_t FIBEZ_TRIGCME[8 /* FIBEZ_TRIGCM */];
    uint32_t FIBEZ_TRIGC;
    uint32_t FIBEZ_TRIGCv[8 /* FIBEZ_TRIGC */];
    uint32_t FIBEZ_TRIGFM;
    uint32_t FIBEZ_TRIGFMI[8 /* FIBEZ_TRIGFM */];
    uint32_t FIBEZ_TRIGFME[8 /* FIBEZ_TRIGFM */];
    uint32_t FIBEZ_TRIGF;
    uint32_t FIBEZ_TRIGFv[8 /* FIBEZ_TRIGF */];

} EXT_STR_h101_FIBEZ_onion;

/*******************************************************/

#define EXT_STR_h101_FIBEZ_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                            \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TBLCM, UINT32, "FIBEZ_TBLCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TBLCMI, UINT32, "FIBEZ_TBLCMI", "FIBEZ_TBLCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TBLCME, UINT32, "FIBEZ_TBLCME", "FIBEZ_TBLCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TBLC, UINT32, "FIBEZ_TBLC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEZ_TBLCv, UINT32, "FIBEZ_TBLCv", "FIBEZ_TBLC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TBLFM, UINT32, "FIBEZ_TBLFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TBLFMI, UINT32, "FIBEZ_TBLFMI", "FIBEZ_TBLFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TBLFME, UINT32, "FIBEZ_TBLFME", "FIBEZ_TBLFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TBLF, UINT32, "FIBEZ_TBLF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEZ_TBLFv, UINT32, "FIBEZ_TBLFv", "FIBEZ_TBLF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TBTCM, UINT32, "FIBEZ_TBTCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TBTCMI, UINT32, "FIBEZ_TBTCMI", "FIBEZ_TBTCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TBTCME, UINT32, "FIBEZ_TBTCME", "FIBEZ_TBTCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TBTC, UINT32, "FIBEZ_TBTC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEZ_TBTCv, UINT32, "FIBEZ_TBTCv", "FIBEZ_TBTC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TBTFM, UINT32, "FIBEZ_TBTFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TBTFMI, UINT32, "FIBEZ_TBTFMI", "FIBEZ_TBTFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TBTFME, UINT32, "FIBEZ_TBTFME", "FIBEZ_TBTFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TBTF, UINT32, "FIBEZ_TBTF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEZ_TBTFv, UINT32, "FIBEZ_TBTFv", "FIBEZ_TBTF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TTLCM, UINT32, "FIBEZ_TTLCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TTLCMI, UINT32, "FIBEZ_TTLCMI", "FIBEZ_TTLCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TTLCME, UINT32, "FIBEZ_TTLCME", "FIBEZ_TTLCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TTLC, UINT32, "FIBEZ_TTLC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEZ_TTLCv, UINT32, "FIBEZ_TTLCv", "FIBEZ_TTLC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TTLFM, UINT32, "FIBEZ_TTLFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TTLFMI, UINT32, "FIBEZ_TTLFMI", "FIBEZ_TTLFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TTLFME, UINT32, "FIBEZ_TTLFME", "FIBEZ_TTLFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TTLF, UINT32, "FIBEZ_TTLF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEZ_TTLFv, UINT32, "FIBEZ_TTLFv", "FIBEZ_TTLF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TTTCM, UINT32, "FIBEZ_TTTCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TTTCMI, UINT32, "FIBEZ_TTTCMI", "FIBEZ_TTTCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TTTCME, UINT32, "FIBEZ_TTTCME", "FIBEZ_TTTCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TTTC, UINT32, "FIBEZ_TTTC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEZ_TTTCv, UINT32, "FIBEZ_TTTCv", "FIBEZ_TTTC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TTTFM, UINT32, "FIBEZ_TTTFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TTTFMI, UINT32, "FIBEZ_TTTFMI", "FIBEZ_TTTFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TTTFME, UINT32, "FIBEZ_TTTFME", "FIBEZ_TTTFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TTTF, UINT32, "FIBEZ_TTTF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEZ_TTTFv, UINT32, "FIBEZ_TTTFv", "FIBEZ_TTTF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TRIGCM, UINT32, "FIBEZ_TRIGCM", 8);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TRIGCMI, UINT32, "FIBEZ_TRIGCMI", "FIBEZ_TRIGCM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TRIGCME, UINT32, "FIBEZ_TRIGCME", "FIBEZ_TRIGCM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TRIGC, UINT32, "FIBEZ_TRIGC", 8);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TRIGCv, UINT32, "FIBEZ_TRIGCv", "FIBEZ_TRIGC");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TRIGFM, UINT32, "FIBEZ_TRIGFM", 8);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TRIGFMI, UINT32, "FIBEZ_TRIGFMI", "FIBEZ_TRIGFM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TRIGFME, UINT32, "FIBEZ_TRIGFME", "FIBEZ_TRIGFM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEZ_TRIGF, UINT32, "FIBEZ_TRIGF", 8);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEZ_TRIGFv, UINT32, "FIBEZ_TRIGFv", "FIBEZ_TRIGF");                \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_FIBEZ_EXT_H101_FIB32_H__*/

/*******************************************************/
