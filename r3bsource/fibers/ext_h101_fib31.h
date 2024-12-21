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

#ifndef __GUARD_H101_FIB31_EXT_H101_FIB31_H__
#define __GUARD_H101_FIB31_EXT_H101_FIB31_H__

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

typedef struct EXT_STR_h101_FIBEI_t
{
    /* RAW */
    uint32_t FIBEI_TBLCM /* [1,512] */;
    uint32_t FIBEI_TBLCMI[512 EXT_STRUCT_CTRL(FIBEI_TBLCM)] /* [1,512] */;
    uint32_t FIBEI_TBLCME[512 EXT_STRUCT_CTRL(FIBEI_TBLCM)] /* [1,32768] */;
    uint32_t FIBEI_TBLC /* [0,32768] */;
    uint32_t FIBEI_TBLCv[32768 EXT_STRUCT_CTRL(FIBEI_TBLC)] /* [0,65535] */;
    uint32_t FIBEI_TBLFM /* [1,512] */;
    uint32_t FIBEI_TBLFMI[512 EXT_STRUCT_CTRL(FIBEI_TBLFM)] /* [1,512] */;
    uint32_t FIBEI_TBLFME[512 EXT_STRUCT_CTRL(FIBEI_TBLFM)] /* [1,32768] */;
    uint32_t FIBEI_TBLF /* [0,32768] */;
    uint32_t FIBEI_TBLFv[32768 EXT_STRUCT_CTRL(FIBEI_TBLF)] /* [0,65535] */;
    uint32_t FIBEI_TBTCM /* [1,512] */;
    uint32_t FIBEI_TBTCMI[512 EXT_STRUCT_CTRL(FIBEI_TBTCM)] /* [1,512] */;
    uint32_t FIBEI_TBTCME[512 EXT_STRUCT_CTRL(FIBEI_TBTCM)] /* [1,32768] */;
    uint32_t FIBEI_TBTC /* [0,32768] */;
    uint32_t FIBEI_TBTCv[32768 EXT_STRUCT_CTRL(FIBEI_TBTC)] /* [0,65535] */;
    uint32_t FIBEI_TBTFM /* [1,512] */;
    uint32_t FIBEI_TBTFMI[512 EXT_STRUCT_CTRL(FIBEI_TBTFM)] /* [1,512] */;
    uint32_t FIBEI_TBTFME[512 EXT_STRUCT_CTRL(FIBEI_TBTFM)] /* [1,32768] */;
    uint32_t FIBEI_TBTF /* [0,32768] */;
    uint32_t FIBEI_TBTFv[32768 EXT_STRUCT_CTRL(FIBEI_TBTF)] /* [0,65535] */;
    uint32_t FIBEI_TTLCM /* [1,512] */;
    uint32_t FIBEI_TTLCMI[512 EXT_STRUCT_CTRL(FIBEI_TTLCM)] /* [1,512] */;
    uint32_t FIBEI_TTLCME[512 EXT_STRUCT_CTRL(FIBEI_TTLCM)] /* [1,32768] */;
    uint32_t FIBEI_TTLC /* [0,32768] */;
    uint32_t FIBEI_TTLCv[32768 EXT_STRUCT_CTRL(FIBEI_TTLC)] /* [0,65535] */;
    uint32_t FIBEI_TTLFM /* [1,512] */;
    uint32_t FIBEI_TTLFMI[512 EXT_STRUCT_CTRL(FIBEI_TTLFM)] /* [1,512] */;
    uint32_t FIBEI_TTLFME[512 EXT_STRUCT_CTRL(FIBEI_TTLFM)] /* [1,32768] */;
    uint32_t FIBEI_TTLF /* [0,32768] */;
    uint32_t FIBEI_TTLFv[32768 EXT_STRUCT_CTRL(FIBEI_TTLF)] /* [0,65535] */;
    uint32_t FIBEI_TTTCM /* [1,512] */;
    uint32_t FIBEI_TTTCMI[512 EXT_STRUCT_CTRL(FIBEI_TTTCM)] /* [1,512] */;
    uint32_t FIBEI_TTTCME[512 EXT_STRUCT_CTRL(FIBEI_TTTCM)] /* [1,32768] */;
    uint32_t FIBEI_TTTC /* [0,32768] */;
    uint32_t FIBEI_TTTCv[32768 EXT_STRUCT_CTRL(FIBEI_TTTC)] /* [0,65535] */;
    uint32_t FIBEI_TTTFM /* [1,512] */;
    uint32_t FIBEI_TTTFMI[512 EXT_STRUCT_CTRL(FIBEI_TTTFM)] /* [1,512] */;
    uint32_t FIBEI_TTTFME[512 EXT_STRUCT_CTRL(FIBEI_TTTFM)] /* [1,32768] */;
    uint32_t FIBEI_TTTF /* [0,32768] */;
    uint32_t FIBEI_TTTFv[32768 EXT_STRUCT_CTRL(FIBEI_TTTF)] /* [0,65535] */;
    uint32_t FIBEI_TRIGCM /* [1,8] */;
    uint32_t FIBEI_TRIGCMI[8 EXT_STRUCT_CTRL(FIBEI_TRIGCM)] /* [1,8] */;
    uint32_t FIBEI_TRIGCME[8 EXT_STRUCT_CTRL(FIBEI_TRIGCM)] /* [1,8] */;
    uint32_t FIBEI_TRIGC /* [0,8] */;
    uint32_t FIBEI_TRIGCv[8 EXT_STRUCT_CTRL(FIBEI_TRIGC)] /* [0,65535] */;
    uint32_t FIBEI_TRIGFM /* [1,8] */;
    uint32_t FIBEI_TRIGFMI[8 EXT_STRUCT_CTRL(FIBEI_TRIGFM)] /* [1,8] */;
    uint32_t FIBEI_TRIGFME[8 EXT_STRUCT_CTRL(FIBEI_TRIGFM)] /* [1,8] */;
    uint32_t FIBEI_TRIGF /* [0,8] */;
    uint32_t FIBEI_TRIGFv[8 EXT_STRUCT_CTRL(FIBEI_TRIGF)] /* [0,65535] */;

} EXT_STR_h101_FIBEI;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBEI_onion_t
{
    /* RAW */
    uint32_t FIBEI_TBLCM;
    uint32_t FIBEI_TBLCMI[512 /* FIBEI_TBLCM */];
    uint32_t FIBEI_TBLCME[512 /* FIBEI_TBLCM */];
    uint32_t FIBEI_TBLC;
    uint32_t FIBEI_TBLCv[32768 /* FIBEI_TBLC */];
    uint32_t FIBEI_TBLFM;
    uint32_t FIBEI_TBLFMI[512 /* FIBEI_TBLFM */];
    uint32_t FIBEI_TBLFME[512 /* FIBEI_TBLFM */];
    uint32_t FIBEI_TBLF;
    uint32_t FIBEI_TBLFv[32768 /* FIBEI_TBLF */];
    uint32_t FIBEI_TBTCM;
    uint32_t FIBEI_TBTCMI[512 /* FIBEI_TBTCM */];
    uint32_t FIBEI_TBTCME[512 /* FIBEI_TBTCM */];
    uint32_t FIBEI_TBTC;
    uint32_t FIBEI_TBTCv[32768 /* FIBEI_TBTC */];
    uint32_t FIBEI_TBTFM;
    uint32_t FIBEI_TBTFMI[512 /* FIBEI_TBTFM */];
    uint32_t FIBEI_TBTFME[512 /* FIBEI_TBTFM */];
    uint32_t FIBEI_TBTF;
    uint32_t FIBEI_TBTFv[32768 /* FIBEI_TBTF */];
    uint32_t FIBEI_TTLCM;
    uint32_t FIBEI_TTLCMI[512 /* FIBEI_TTLCM */];
    uint32_t FIBEI_TTLCME[512 /* FIBEI_TTLCM */];
    uint32_t FIBEI_TTLC;
    uint32_t FIBEI_TTLCv[32768 /* FIBEI_TTLC */];
    uint32_t FIBEI_TTLFM;
    uint32_t FIBEI_TTLFMI[512 /* FIBEI_TTLFM */];
    uint32_t FIBEI_TTLFME[512 /* FIBEI_TTLFM */];
    uint32_t FIBEI_TTLF;
    uint32_t FIBEI_TTLFv[32768 /* FIBEI_TTLF */];
    uint32_t FIBEI_TTTCM;
    uint32_t FIBEI_TTTCMI[512 /* FIBEI_TTTCM */];
    uint32_t FIBEI_TTTCME[512 /* FIBEI_TTTCM */];
    uint32_t FIBEI_TTTC;
    uint32_t FIBEI_TTTCv[32768 /* FIBEI_TTTC */];
    uint32_t FIBEI_TTTFM;
    uint32_t FIBEI_TTTFMI[512 /* FIBEI_TTTFM */];
    uint32_t FIBEI_TTTFME[512 /* FIBEI_TTTFM */];
    uint32_t FIBEI_TTTF;
    uint32_t FIBEI_TTTFv[32768 /* FIBEI_TTTF */];
    uint32_t FIBEI_TRIGCM;
    uint32_t FIBEI_TRIGCMI[8 /* FIBEI_TRIGCM */];
    uint32_t FIBEI_TRIGCME[8 /* FIBEI_TRIGCM */];
    uint32_t FIBEI_TRIGC;
    uint32_t FIBEI_TRIGCv[8 /* FIBEI_TRIGC */];
    uint32_t FIBEI_TRIGFM;
    uint32_t FIBEI_TRIGFMI[8 /* FIBEI_TRIGFM */];
    uint32_t FIBEI_TRIGFME[8 /* FIBEI_TRIGFM */];
    uint32_t FIBEI_TRIGF;
    uint32_t FIBEI_TRIGFv[8 /* FIBEI_TRIGF */];

} EXT_STR_h101_FIBEI_onion;

/*******************************************************/

#define EXT_STR_h101_FIBEI_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                            \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TBLCM, UINT32, "FIBEI_TBLCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TBLCMI, UINT32, "FIBEI_TBLCMI", "FIBEI_TBLCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TBLCME, UINT32, "FIBEI_TBLCME", "FIBEI_TBLCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TBLC, UINT32, "FIBEI_TBLC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEI_TBLCv, UINT32, "FIBEI_TBLCv", "FIBEI_TBLC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TBLFM, UINT32, "FIBEI_TBLFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TBLFMI, UINT32, "FIBEI_TBLFMI", "FIBEI_TBLFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TBLFME, UINT32, "FIBEI_TBLFME", "FIBEI_TBLFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TBLF, UINT32, "FIBEI_TBLF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEI_TBLFv, UINT32, "FIBEI_TBLFv", "FIBEI_TBLF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TBTCM, UINT32, "FIBEI_TBTCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TBTCMI, UINT32, "FIBEI_TBTCMI", "FIBEI_TBTCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TBTCME, UINT32, "FIBEI_TBTCME", "FIBEI_TBTCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TBTC, UINT32, "FIBEI_TBTC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEI_TBTCv, UINT32, "FIBEI_TBTCv", "FIBEI_TBTC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TBTFM, UINT32, "FIBEI_TBTFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TBTFMI, UINT32, "FIBEI_TBTFMI", "FIBEI_TBTFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TBTFME, UINT32, "FIBEI_TBTFME", "FIBEI_TBTFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TBTF, UINT32, "FIBEI_TBTF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEI_TBTFv, UINT32, "FIBEI_TBTFv", "FIBEI_TBTF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TTLCM, UINT32, "FIBEI_TTLCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TTLCMI, UINT32, "FIBEI_TTLCMI", "FIBEI_TTLCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TTLCME, UINT32, "FIBEI_TTLCME", "FIBEI_TTLCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TTLC, UINT32, "FIBEI_TTLC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEI_TTLCv, UINT32, "FIBEI_TTLCv", "FIBEI_TTLC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TTLFM, UINT32, "FIBEI_TTLFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TTLFMI, UINT32, "FIBEI_TTLFMI", "FIBEI_TTLFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TTLFME, UINT32, "FIBEI_TTLFME", "FIBEI_TTLFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TTLF, UINT32, "FIBEI_TTLF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEI_TTLFv, UINT32, "FIBEI_TTLFv", "FIBEI_TTLF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TTTCM, UINT32, "FIBEI_TTTCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TTTCMI, UINT32, "FIBEI_TTTCMI", "FIBEI_TTTCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TTTCME, UINT32, "FIBEI_TTTCME", "FIBEI_TTTCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TTTC, UINT32, "FIBEI_TTTC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEI_TTTCv, UINT32, "FIBEI_TTTCv", "FIBEI_TTTC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TTTFM, UINT32, "FIBEI_TTTFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TTTFMI, UINT32, "FIBEI_TTTFMI", "FIBEI_TTTFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TTTFME, UINT32, "FIBEI_TTTFME", "FIBEI_TTTFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TTTF, UINT32, "FIBEI_TTTF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEI_TTTFv, UINT32, "FIBEI_TTTFv", "FIBEI_TTTF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TRIGCM, UINT32, "FIBEI_TRIGCM", 8);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TRIGCMI, UINT32, "FIBEI_TRIGCMI", "FIBEI_TRIGCM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TRIGCME, UINT32, "FIBEI_TRIGCME", "FIBEI_TRIGCM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TRIGC, UINT32, "FIBEI_TRIGC", 8);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TRIGCv, UINT32, "FIBEI_TRIGCv", "FIBEI_TRIGC");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TRIGFM, UINT32, "FIBEI_TRIGFM", 8);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TRIGFMI, UINT32, "FIBEI_TRIGFMI", "FIBEI_TRIGFM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TRIGFME, UINT32, "FIBEI_TRIGFME", "FIBEI_TRIGFM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEI_TRIGF, UINT32, "FIBEI_TRIGF", 8);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEI_TRIGFv, UINT32, "FIBEI_TRIGFv", "FIBEI_TRIGF");                \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_FIB31_EXT_H101_FIB31_H__*/

/*******************************************************/
