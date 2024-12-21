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

#ifndef __GUARD_H101_FIBZEA_EXT_H101_FIB23A_H__
#define __GUARD_H101_FIBZEA_EXT_H101_FIB23A_H__

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

typedef struct EXT_STR_h101_FIBZEA_t
{
    /* RAW */
    uint32_t FIBZEA_TBLCM /* [1,256] */;
    uint32_t FIBZEA_TBLCMI[256 EXT_STRUCT_CTRL(FIBZEA_TBLCM)] /* [1,256] */;
    uint32_t FIBZEA_TBLCME[256 EXT_STRUCT_CTRL(FIBZEA_TBLCM)] /* [1,16384] */;
    uint32_t FIBZEA_TBLC /* [0,16384] */;
    uint32_t FIBZEA_TBLCv[16384 EXT_STRUCT_CTRL(FIBZEA_TBLC)] /* [0,65535] */;
    uint32_t FIBZEA_TBLFM /* [1,256] */;
    uint32_t FIBZEA_TBLFMI[256 EXT_STRUCT_CTRL(FIBZEA_TBLFM)] /* [1,256] */;
    uint32_t FIBZEA_TBLFME[256 EXT_STRUCT_CTRL(FIBZEA_TBLFM)] /* [1,16384] */;
    uint32_t FIBZEA_TBLF /* [0,16384] */;
    uint32_t FIBZEA_TBLFv[16384 EXT_STRUCT_CTRL(FIBZEA_TBLF)] /* [0,65535] */;
    uint32_t FIBZEA_TBTCM /* [1,256] */;
    uint32_t FIBZEA_TBTCMI[256 EXT_STRUCT_CTRL(FIBZEA_TBTCM)] /* [1,256] */;
    uint32_t FIBZEA_TBTCME[256 EXT_STRUCT_CTRL(FIBZEA_TBTCM)] /* [1,16384] */;
    uint32_t FIBZEA_TBTC /* [0,16384] */;
    uint32_t FIBZEA_TBTCv[16384 EXT_STRUCT_CTRL(FIBZEA_TBTC)] /* [0,65535] */;
    uint32_t FIBZEA_TBTFM /* [1,256] */;
    uint32_t FIBZEA_TBTFMI[256 EXT_STRUCT_CTRL(FIBZEA_TBTFM)] /* [1,256] */;
    uint32_t FIBZEA_TBTFME[256 EXT_STRUCT_CTRL(FIBZEA_TBTFM)] /* [1,16384] */;
    uint32_t FIBZEA_TBTF /* [0,16384] */;
    uint32_t FIBZEA_TBTFv[16384 EXT_STRUCT_CTRL(FIBZEA_TBTF)] /* [0,65535] */;
    uint32_t FIBZEA_TTLCM /* [1,256] */;
    uint32_t FIBZEA_TTLCMI[256 EXT_STRUCT_CTRL(FIBZEA_TTLCM)] /* [1,256] */;
    uint32_t FIBZEA_TTLCME[256 EXT_STRUCT_CTRL(FIBZEA_TTLCM)] /* [1,16384] */;
    uint32_t FIBZEA_TTLC /* [0,16384] */;
    uint32_t FIBZEA_TTLCv[16384 EXT_STRUCT_CTRL(FIBZEA_TTLC)] /* [0,65535] */;
    uint32_t FIBZEA_TTLFM /* [1,256] */;
    uint32_t FIBZEA_TTLFMI[256 EXT_STRUCT_CTRL(FIBZEA_TTLFM)] /* [1,256] */;
    uint32_t FIBZEA_TTLFME[256 EXT_STRUCT_CTRL(FIBZEA_TTLFM)] /* [1,16384] */;
    uint32_t FIBZEA_TTLF /* [0,16384] */;
    uint32_t FIBZEA_TTLFv[16384 EXT_STRUCT_CTRL(FIBZEA_TTLF)] /* [0,65535] */;
    uint32_t FIBZEA_TTTCM /* [1,256] */;
    uint32_t FIBZEA_TTTCMI[256 EXT_STRUCT_CTRL(FIBZEA_TTTCM)] /* [1,256] */;
    uint32_t FIBZEA_TTTCME[256 EXT_STRUCT_CTRL(FIBZEA_TTTCM)] /* [1,16384] */;
    uint32_t FIBZEA_TTTC /* [0,16384] */;
    uint32_t FIBZEA_TTTCv[16384 EXT_STRUCT_CTRL(FIBZEA_TTTC)] /* [0,65535] */;
    uint32_t FIBZEA_TTTFM /* [1,256] */;
    uint32_t FIBZEA_TTTFMI[256 EXT_STRUCT_CTRL(FIBZEA_TTTFM)] /* [1,256] */;
    uint32_t FIBZEA_TTTFME[256 EXT_STRUCT_CTRL(FIBZEA_TTTFM)] /* [1,16384] */;
    uint32_t FIBZEA_TTTF /* [0,16384] */;
    uint32_t FIBZEA_TTTFv[16384 EXT_STRUCT_CTRL(FIBZEA_TTTF)] /* [0,65535] */;
    uint32_t FIBZEA_TRIGCM /* [1,4] */;
    uint32_t FIBZEA_TRIGCMI[4 EXT_STRUCT_CTRL(FIBZEA_TRIGCM)] /* [1,4] */;
    uint32_t FIBZEA_TRIGCME[4 EXT_STRUCT_CTRL(FIBZEA_TRIGCM)] /* [1,4] */;
    uint32_t FIBZEA_TRIGC /* [0,4] */;
    uint32_t FIBZEA_TRIGCv[4 EXT_STRUCT_CTRL(FIBZEA_TRIGC)] /* [0,65535] */;
    uint32_t FIBZEA_TRIGFM /* [1,4] */;
    uint32_t FIBZEA_TRIGFMI[4 EXT_STRUCT_CTRL(FIBZEA_TRIGFM)] /* [1,4] */;
    uint32_t FIBZEA_TRIGFME[4 EXT_STRUCT_CTRL(FIBZEA_TRIGFM)] /* [1,4] */;
    uint32_t FIBZEA_TRIGF /* [0,4] */;
    uint32_t FIBZEA_TRIGFv[4 EXT_STRUCT_CTRL(FIBZEA_TRIGF)] /* [0,65535] */;

} EXT_STR_h101_FIBZEA;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBZEA_onion_t
{
    /* RAW */
    uint32_t FIBZEA_TBLCM;
    uint32_t FIBZEA_TBLCMI[256 /* FIBZEA_TBLCM */];
    uint32_t FIBZEA_TBLCME[256 /* FIBZEA_TBLCM */];
    uint32_t FIBZEA_TBLC;
    uint32_t FIBZEA_TBLCv[16384 /* FIBZEA_TBLC */];
    uint32_t FIBZEA_TBLFM;
    uint32_t FIBZEA_TBLFMI[256 /* FIBZEA_TBLFM */];
    uint32_t FIBZEA_TBLFME[256 /* FIBZEA_TBLFM */];
    uint32_t FIBZEA_TBLF;
    uint32_t FIBZEA_TBLFv[16384 /* FIBZEA_TBLF */];
    uint32_t FIBZEA_TBTCM;
    uint32_t FIBZEA_TBTCMI[256 /* FIBZEA_TBTCM */];
    uint32_t FIBZEA_TBTCME[256 /* FIBZEA_TBTCM */];
    uint32_t FIBZEA_TBTC;
    uint32_t FIBZEA_TBTCv[16384 /* FIBZEA_TBTC */];
    uint32_t FIBZEA_TBTFM;
    uint32_t FIBZEA_TBTFMI[256 /* FIBZEA_TBTFM */];
    uint32_t FIBZEA_TBTFME[256 /* FIBZEA_TBTFM */];
    uint32_t FIBZEA_TBTF;
    uint32_t FIBZEA_TBTFv[16384 /* FIBZEA_TBTF */];
    uint32_t FIBZEA_TTLCM;
    uint32_t FIBZEA_TTLCMI[256 /* FIBZEA_TTLCM */];
    uint32_t FIBZEA_TTLCME[256 /* FIBZEA_TTLCM */];
    uint32_t FIBZEA_TTLC;
    uint32_t FIBZEA_TTLCv[16384 /* FIBZEA_TTLC */];
    uint32_t FIBZEA_TTLFM;
    uint32_t FIBZEA_TTLFMI[256 /* FIBZEA_TTLFM */];
    uint32_t FIBZEA_TTLFME[256 /* FIBZEA_TTLFM */];
    uint32_t FIBZEA_TTLF;
    uint32_t FIBZEA_TTLFv[16384 /* FIBZEA_TTLF */];
    uint32_t FIBZEA_TTTCM;
    uint32_t FIBZEA_TTTCMI[256 /* FIBZEA_TTTCM */];
    uint32_t FIBZEA_TTTCME[256 /* FIBZEA_TTTCM */];
    uint32_t FIBZEA_TTTC;
    uint32_t FIBZEA_TTTCv[16384 /* FIBZEA_TTTC */];
    uint32_t FIBZEA_TTTFM;
    uint32_t FIBZEA_TTTFMI[256 /* FIBZEA_TTTFM */];
    uint32_t FIBZEA_TTTFME[256 /* FIBZEA_TTTFM */];
    uint32_t FIBZEA_TTTF;
    uint32_t FIBZEA_TTTFv[16384 /* FIBZEA_TTTF */];
    uint32_t FIBZEA_TRIGCM;
    uint32_t FIBZEA_TRIGCMI[4 /* FIBZEA_TRIGCM */];
    uint32_t FIBZEA_TRIGCME[4 /* FIBZEA_TRIGCM */];
    uint32_t FIBZEA_TRIGC;
    uint32_t FIBZEA_TRIGCv[4 /* FIBZEA_TRIGC */];
    uint32_t FIBZEA_TRIGFM;
    uint32_t FIBZEA_TRIGFMI[4 /* FIBZEA_TRIGFM */];
    uint32_t FIBZEA_TRIGFME[4 /* FIBZEA_TRIGFM */];
    uint32_t FIBZEA_TRIGF;
    uint32_t FIBZEA_TRIGFv[4 /* FIBZEA_TRIGF */];

} EXT_STR_h101_FIBZEA_onion;

/*******************************************************/

#define EXT_STR_h101_FIBZEA_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                    \
    do                                                                                                        \
    {                                                                                                         \
        ok = 1;                                                                                               \
        /* RAW */                                                                                             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TBLCM, UINT32, "FIBZEA_TBLCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBLCMI, UINT32, "FIBZEA_TBLCMI", "FIBZEA_TBLCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBLCME, UINT32, "FIBZEA_TBLCME", "FIBZEA_TBLCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TBLC, UINT32, "FIBZEA_TBLC", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBLCv, UINT32, "FIBZEA_TBLCv", "FIBZEA_TBLC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TBLFM, UINT32, "FIBZEA_TBLFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBLFMI, UINT32, "FIBZEA_TBLFMI", "FIBZEA_TBLFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBLFME, UINT32, "FIBZEA_TBLFME", "FIBZEA_TBLFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TBLF, UINT32, "FIBZEA_TBLF", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBLFv, UINT32, "FIBZEA_TBLFv", "FIBZEA_TBLF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TBTCM, UINT32, "FIBZEA_TBTCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBTCMI, UINT32, "FIBZEA_TBTCMI", "FIBZEA_TBTCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBTCME, UINT32, "FIBZEA_TBTCME", "FIBZEA_TBTCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TBTC, UINT32, "FIBZEA_TBTC", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBTCv, UINT32, "FIBZEA_TBTCv", "FIBZEA_TBTC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TBTFM, UINT32, "FIBZEA_TBTFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBTFMI, UINT32, "FIBZEA_TBTFMI", "FIBZEA_TBTFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBTFME, UINT32, "FIBZEA_TBTFME", "FIBZEA_TBTFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TBTF, UINT32, "FIBZEA_TBTF", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TBTFv, UINT32, "FIBZEA_TBTFv", "FIBZEA_TBTF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TTLCM, UINT32, "FIBZEA_TTLCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTLCMI, UINT32, "FIBZEA_TTLCMI", "FIBZEA_TTLCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTLCME, UINT32, "FIBZEA_TTLCME", "FIBZEA_TTLCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TTLC, UINT32, "FIBZEA_TTLC", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTLCv, UINT32, "FIBZEA_TTLCv", "FIBZEA_TTLC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TTLFM, UINT32, "FIBZEA_TTLFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTLFMI, UINT32, "FIBZEA_TTLFMI", "FIBZEA_TTLFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTLFME, UINT32, "FIBZEA_TTLFME", "FIBZEA_TTLFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TTLF, UINT32, "FIBZEA_TTLF", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTLFv, UINT32, "FIBZEA_TTLFv", "FIBZEA_TTLF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TTTCM, UINT32, "FIBZEA_TTTCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTTCMI, UINT32, "FIBZEA_TTTCMI", "FIBZEA_TTTCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTTCME, UINT32, "FIBZEA_TTTCME", "FIBZEA_TTTCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TTTC, UINT32, "FIBZEA_TTTC", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTTCv, UINT32, "FIBZEA_TTTCv", "FIBZEA_TTTC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TTTFM, UINT32, "FIBZEA_TTTFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTTFMI, UINT32, "FIBZEA_TTTFMI", "FIBZEA_TTTFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTTFME, UINT32, "FIBZEA_TTTFME", "FIBZEA_TTTFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TTTF, UINT32, "FIBZEA_TTTF", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TTTFv, UINT32, "FIBZEA_TTTFv", "FIBZEA_TTTF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TRIGCM, UINT32, "FIBZEA_TRIGCM", 4); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TRIGCMI, UINT32, "FIBZEA_TRIGCMI", "FIBZEA_TRIGCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TRIGCME, UINT32, "FIBZEA_TRIGCME", "FIBZEA_TRIGCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TRIGC, UINT32, "FIBZEA_TRIGC", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TRIGCv, UINT32, "FIBZEA_TRIGCv", "FIBZEA_TRIGC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TRIGFM, UINT32, "FIBZEA_TRIGFM", 4); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TRIGFMI, UINT32, "FIBZEA_TRIGFMI", "FIBZEA_TRIGFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TRIGFME, UINT32, "FIBZEA_TRIGFME", "FIBZEA_TRIGFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEA_TRIGF, UINT32, "FIBZEA_TRIGF", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEA_TRIGFv, UINT32, "FIBZEA_TRIGFv", "FIBZEA_TRIGF");      \
                                                                                                              \
    } while (0);

#endif /*__GUARD_H101_FIBZEA_EXT_H101_FIB23A_H__*/

/*******************************************************/
