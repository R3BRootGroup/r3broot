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

#ifndef __GUARD_H101_FIB30_EXT_H101_FIB30_H__
#define __GUARD_H101_FIB30_EXT_H101_FIB30_H__

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

typedef struct EXT_STR_h101_FIBEO_t
{
    /* RAW */
    uint32_t FIBEO_TBLCM /* [1,512] */;
    uint32_t FIBEO_TBLCMI[512 EXT_STRUCT_CTRL(FIBEO_TBLCM)] /* [1,512] */;
    uint32_t FIBEO_TBLCME[512 EXT_STRUCT_CTRL(FIBEO_TBLCM)] /* [1,32768] */;
    uint32_t FIBEO_TBLC /* [0,32768] */;
    uint32_t FIBEO_TBLCv[32768 EXT_STRUCT_CTRL(FIBEO_TBLC)] /* [0,65535] */;
    uint32_t FIBEO_TBLFM /* [1,512] */;
    uint32_t FIBEO_TBLFMI[512 EXT_STRUCT_CTRL(FIBEO_TBLFM)] /* [1,512] */;
    uint32_t FIBEO_TBLFME[512 EXT_STRUCT_CTRL(FIBEO_TBLFM)] /* [1,32768] */;
    uint32_t FIBEO_TBLF /* [0,32768] */;
    uint32_t FIBEO_TBLFv[32768 EXT_STRUCT_CTRL(FIBEO_TBLF)] /* [0,65535] */;
    uint32_t FIBEO_TBTCM /* [1,512] */;
    uint32_t FIBEO_TBTCMI[512 EXT_STRUCT_CTRL(FIBEO_TBTCM)] /* [1,512] */;
    uint32_t FIBEO_TBTCME[512 EXT_STRUCT_CTRL(FIBEO_TBTCM)] /* [1,32768] */;
    uint32_t FIBEO_TBTC /* [0,32768] */;
    uint32_t FIBEO_TBTCv[32768 EXT_STRUCT_CTRL(FIBEO_TBTC)] /* [0,65535] */;
    uint32_t FIBEO_TBTFM /* [1,512] */;
    uint32_t FIBEO_TBTFMI[512 EXT_STRUCT_CTRL(FIBEO_TBTFM)] /* [1,512] */;
    uint32_t FIBEO_TBTFME[512 EXT_STRUCT_CTRL(FIBEO_TBTFM)] /* [1,32768] */;
    uint32_t FIBEO_TBTF /* [0,32768] */;
    uint32_t FIBEO_TBTFv[32768 EXT_STRUCT_CTRL(FIBEO_TBTF)] /* [0,65535] */;
    uint32_t FIBEO_TTLCM /* [1,512] */;
    uint32_t FIBEO_TTLCMI[512 EXT_STRUCT_CTRL(FIBEO_TTLCM)] /* [1,512] */;
    uint32_t FIBEO_TTLCME[512 EXT_STRUCT_CTRL(FIBEO_TTLCM)] /* [1,32768] */;
    uint32_t FIBEO_TTLC /* [0,32768] */;
    uint32_t FIBEO_TTLCv[32768 EXT_STRUCT_CTRL(FIBEO_TTLC)] /* [0,65535] */;
    uint32_t FIBEO_TTLFM /* [1,512] */;
    uint32_t FIBEO_TTLFMI[512 EXT_STRUCT_CTRL(FIBEO_TTLFM)] /* [1,512] */;
    uint32_t FIBEO_TTLFME[512 EXT_STRUCT_CTRL(FIBEO_TTLFM)] /* [1,32768] */;
    uint32_t FIBEO_TTLF /* [0,32768] */;
    uint32_t FIBEO_TTLFv[32768 EXT_STRUCT_CTRL(FIBEO_TTLF)] /* [0,65535] */;
    uint32_t FIBEO_TTTCM /* [1,512] */;
    uint32_t FIBEO_TTTCMI[512 EXT_STRUCT_CTRL(FIBEO_TTTCM)] /* [1,512] */;
    uint32_t FIBEO_TTTCME[512 EXT_STRUCT_CTRL(FIBEO_TTTCM)] /* [1,32768] */;
    uint32_t FIBEO_TTTC /* [0,32768] */;
    uint32_t FIBEO_TTTCv[32768 EXT_STRUCT_CTRL(FIBEO_TTTC)] /* [0,65535] */;
    uint32_t FIBEO_TTTFM /* [1,512] */;
    uint32_t FIBEO_TTTFMI[512 EXT_STRUCT_CTRL(FIBEO_TTTFM)] /* [1,512] */;
    uint32_t FIBEO_TTTFME[512 EXT_STRUCT_CTRL(FIBEO_TTTFM)] /* [1,32768] */;
    uint32_t FIBEO_TTTF /* [0,32768] */;
    uint32_t FIBEO_TTTFv[32768 EXT_STRUCT_CTRL(FIBEO_TTTF)] /* [0,65535] */;
    uint32_t FIBEO_TRIGCM /* [1,8] */;
    uint32_t FIBEO_TRIGCMI[8 EXT_STRUCT_CTRL(FIBEO_TRIGCM)] /* [1,8] */;
    uint32_t FIBEO_TRIGCME[8 EXT_STRUCT_CTRL(FIBEO_TRIGCM)] /* [1,8] */;
    uint32_t FIBEO_TRIGC /* [0,8] */;
    uint32_t FIBEO_TRIGCv[8 EXT_STRUCT_CTRL(FIBEO_TRIGC)] /* [0,65535] */;
    uint32_t FIBEO_TRIGFM /* [1,8] */;
    uint32_t FIBEO_TRIGFMI[8 EXT_STRUCT_CTRL(FIBEO_TRIGFM)] /* [1,8] */;
    uint32_t FIBEO_TRIGFME[8 EXT_STRUCT_CTRL(FIBEO_TRIGFM)] /* [1,8] */;
    uint32_t FIBEO_TRIGF /* [0,8] */;
    uint32_t FIBEO_TRIGFv[8 EXT_STRUCT_CTRL(FIBEO_TRIGF)] /* [0,65535] */;

} EXT_STR_h101_FIBEO;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBEO_onion_t
{
    /* RAW */
    uint32_t FIBEO_TBLCM;
    uint32_t FIBEO_TBLCMI[512 /* FIBEO_TBLCM */];
    uint32_t FIBEO_TBLCME[512 /* FIBEO_TBLCM */];
    uint32_t FIBEO_TBLC;
    uint32_t FIBEO_TBLCv[32768 /* FIBEO_TBLC */];
    uint32_t FIBEO_TBLFM;
    uint32_t FIBEO_TBLFMI[512 /* FIBEO_TBLFM */];
    uint32_t FIBEO_TBLFME[512 /* FIBEO_TBLFM */];
    uint32_t FIBEO_TBLF;
    uint32_t FIBEO_TBLFv[32768 /* FIBEO_TBLF */];
    uint32_t FIBEO_TBTCM;
    uint32_t FIBEO_TBTCMI[512 /* FIBEO_TBTCM */];
    uint32_t FIBEO_TBTCME[512 /* FIBEO_TBTCM */];
    uint32_t FIBEO_TBTC;
    uint32_t FIBEO_TBTCv[32768 /* FIBEO_TBTC */];
    uint32_t FIBEO_TBTFM;
    uint32_t FIBEO_TBTFMI[512 /* FIBEO_TBTFM */];
    uint32_t FIBEO_TBTFME[512 /* FIBEO_TBTFM */];
    uint32_t FIBEO_TBTF;
    uint32_t FIBEO_TBTFv[32768 /* FIBEO_TBTF */];
    uint32_t FIBEO_TTLCM;
    uint32_t FIBEO_TTLCMI[512 /* FIBEO_TTLCM */];
    uint32_t FIBEO_TTLCME[512 /* FIBEO_TTLCM */];
    uint32_t FIBEO_TTLC;
    uint32_t FIBEO_TTLCv[32768 /* FIBEO_TTLC */];
    uint32_t FIBEO_TTLFM;
    uint32_t FIBEO_TTLFMI[512 /* FIBEO_TTLFM */];
    uint32_t FIBEO_TTLFME[512 /* FIBEO_TTLFM */];
    uint32_t FIBEO_TTLF;
    uint32_t FIBEO_TTLFv[32768 /* FIBEO_TTLF */];
    uint32_t FIBEO_TTTCM;
    uint32_t FIBEO_TTTCMI[512 /* FIBEO_TTTCM */];
    uint32_t FIBEO_TTTCME[512 /* FIBEO_TTTCM */];
    uint32_t FIBEO_TTTC;
    uint32_t FIBEO_TTTCv[32768 /* FIBEO_TTTC */];
    uint32_t FIBEO_TTTFM;
    uint32_t FIBEO_TTTFMI[512 /* FIBEO_TTTFM */];
    uint32_t FIBEO_TTTFME[512 /* FIBEO_TTTFM */];
    uint32_t FIBEO_TTTF;
    uint32_t FIBEO_TTTFv[32768 /* FIBEO_TTTF */];
    uint32_t FIBEO_TRIGCM;
    uint32_t FIBEO_TRIGCMI[8 /* FIBEO_TRIGCM */];
    uint32_t FIBEO_TRIGCME[8 /* FIBEO_TRIGCM */];
    uint32_t FIBEO_TRIGC;
    uint32_t FIBEO_TRIGCv[8 /* FIBEO_TRIGC */];
    uint32_t FIBEO_TRIGFM;
    uint32_t FIBEO_TRIGFMI[8 /* FIBEO_TRIGFM */];
    uint32_t FIBEO_TRIGFME[8 /* FIBEO_TRIGFM */];
    uint32_t FIBEO_TRIGF;
    uint32_t FIBEO_TRIGFv[8 /* FIBEO_TRIGF */];

} EXT_STR_h101_FIBEO_onion;

/*******************************************************/

#define EXT_STR_h101_FIBEO_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                            \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TBLCM, UINT32, "FIBEO_TBLCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TBLCMI, UINT32, "FIBEO_TBLCMI", "FIBEO_TBLCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TBLCME, UINT32, "FIBEO_TBLCME", "FIBEO_TBLCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TBLC, UINT32, "FIBEO_TBLC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEO_TBLCv, UINT32, "FIBEO_TBLCv", "FIBEO_TBLC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TBLFM, UINT32, "FIBEO_TBLFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TBLFMI, UINT32, "FIBEO_TBLFMI", "FIBEO_TBLFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TBLFME, UINT32, "FIBEO_TBLFME", "FIBEO_TBLFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TBLF, UINT32, "FIBEO_TBLF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEO_TBLFv, UINT32, "FIBEO_TBLFv", "FIBEO_TBLF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TBTCM, UINT32, "FIBEO_TBTCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TBTCMI, UINT32, "FIBEO_TBTCMI", "FIBEO_TBTCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TBTCME, UINT32, "FIBEO_TBTCME", "FIBEO_TBTCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TBTC, UINT32, "FIBEO_TBTC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEO_TBTCv, UINT32, "FIBEO_TBTCv", "FIBEO_TBTC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TBTFM, UINT32, "FIBEO_TBTFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TBTFMI, UINT32, "FIBEO_TBTFMI", "FIBEO_TBTFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TBTFME, UINT32, "FIBEO_TBTFME", "FIBEO_TBTFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TBTF, UINT32, "FIBEO_TBTF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEO_TBTFv, UINT32, "FIBEO_TBTFv", "FIBEO_TBTF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TTLCM, UINT32, "FIBEO_TTLCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TTLCMI, UINT32, "FIBEO_TTLCMI", "FIBEO_TTLCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TTLCME, UINT32, "FIBEO_TTLCME", "FIBEO_TTLCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TTLC, UINT32, "FIBEO_TTLC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEO_TTLCv, UINT32, "FIBEO_TTLCv", "FIBEO_TTLC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TTLFM, UINT32, "FIBEO_TTLFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TTLFMI, UINT32, "FIBEO_TTLFMI", "FIBEO_TTLFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TTLFME, UINT32, "FIBEO_TTLFME", "FIBEO_TTLFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TTLF, UINT32, "FIBEO_TTLF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEO_TTLFv, UINT32, "FIBEO_TTLFv", "FIBEO_TTLF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TTTCM, UINT32, "FIBEO_TTTCM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TTTCMI, UINT32, "FIBEO_TTTCMI", "FIBEO_TTTCM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TTTCME, UINT32, "FIBEO_TTTCME", "FIBEO_TTTCM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TTTC, UINT32, "FIBEO_TTTC", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEO_TTTCv, UINT32, "FIBEO_TTTCv", "FIBEO_TTTC"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TTTFM, UINT32, "FIBEO_TTTFM", 512);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TTTFMI, UINT32, "FIBEO_TTTFMI", "FIBEO_TTTFM");                \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TTTFME, UINT32, "FIBEO_TTTFME", "FIBEO_TTTFM");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TTTF, UINT32, "FIBEO_TTTF", 32768);          \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, FIBEO_TTTFv, UINT32, "FIBEO_TTTFv", "FIBEO_TTTF"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TRIGCM, UINT32, "FIBEO_TRIGCM", 8);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TRIGCMI, UINT32, "FIBEO_TRIGCMI", "FIBEO_TRIGCM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TRIGCME, UINT32, "FIBEO_TRIGCME", "FIBEO_TRIGCM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TRIGC, UINT32, "FIBEO_TRIGC", 8);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TRIGCv, UINT32, "FIBEO_TRIGCv", "FIBEO_TRIGC");                \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TRIGFM, UINT32, "FIBEO_TRIGFM", 8);          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TRIGFMI, UINT32, "FIBEO_TRIGFMI", "FIBEO_TRIGFM");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TRIGFME, UINT32, "FIBEO_TRIGFME", "FIBEO_TRIGFM");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBEO_TRIGF, UINT32, "FIBEO_TRIGF", 8);            \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, FIBEO_TRIGFv, UINT32, "FIBEO_TRIGFv", "FIBEO_TRIGF");                \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_FIB30_EXT_H101_FIB30_H__*/

/*******************************************************/
