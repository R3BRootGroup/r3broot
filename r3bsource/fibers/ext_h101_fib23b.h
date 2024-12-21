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

#ifndef __GUARD_H101_FIBZEB_EXT_H101_FIB23B_H__
#define __GUARD_H101_FIBZEB_EXT_H101_FIB23B_H__

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

typedef struct EXT_STR_h101_FIBZEB_t
{
    /* RAW */
    uint32_t FIBZEB_TBLCM /* [1,256] */;
    uint32_t FIBZEB_TBLCMI[256 EXT_STRUCT_CTRL(FIBZEB_TBLCM)] /* [1,256] */;
    uint32_t FIBZEB_TBLCME[256 EXT_STRUCT_CTRL(FIBZEB_TBLCM)] /* [1,16384] */;
    uint32_t FIBZEB_TBLC /* [0,16384] */;
    uint32_t FIBZEB_TBLCv[16384 EXT_STRUCT_CTRL(FIBZEB_TBLC)] /* [0,65535] */;
    uint32_t FIBZEB_TBLFM /* [1,256] */;
    uint32_t FIBZEB_TBLFMI[256 EXT_STRUCT_CTRL(FIBZEB_TBLFM)] /* [1,256] */;
    uint32_t FIBZEB_TBLFME[256 EXT_STRUCT_CTRL(FIBZEB_TBLFM)] /* [1,16384] */;
    uint32_t FIBZEB_TBLF /* [0,16384] */;
    uint32_t FIBZEB_TBLFv[16384 EXT_STRUCT_CTRL(FIBZEB_TBLF)] /* [0,65535] */;
    uint32_t FIBZEB_TBTCM /* [1,256] */;
    uint32_t FIBZEB_TBTCMI[256 EXT_STRUCT_CTRL(FIBZEB_TBTCM)] /* [1,256] */;
    uint32_t FIBZEB_TBTCME[256 EXT_STRUCT_CTRL(FIBZEB_TBTCM)] /* [1,16384] */;
    uint32_t FIBZEB_TBTC /* [0,16384] */;
    uint32_t FIBZEB_TBTCv[16384 EXT_STRUCT_CTRL(FIBZEB_TBTC)] /* [0,65535] */;
    uint32_t FIBZEB_TBTFM /* [1,256] */;
    uint32_t FIBZEB_TBTFMI[256 EXT_STRUCT_CTRL(FIBZEB_TBTFM)] /* [1,256] */;
    uint32_t FIBZEB_TBTFME[256 EXT_STRUCT_CTRL(FIBZEB_TBTFM)] /* [1,16384] */;
    uint32_t FIBZEB_TBTF /* [0,16384] */;
    uint32_t FIBZEB_TBTFv[16384 EXT_STRUCT_CTRL(FIBZEB_TBTF)] /* [0,65535] */;
    uint32_t FIBZEB_TTLCM /* [1,256] */;
    uint32_t FIBZEB_TTLCMI[256 EXT_STRUCT_CTRL(FIBZEB_TTLCM)] /* [1,256] */;
    uint32_t FIBZEB_TTLCME[256 EXT_STRUCT_CTRL(FIBZEB_TTLCM)] /* [1,16384] */;
    uint32_t FIBZEB_TTLC /* [0,16384] */;
    uint32_t FIBZEB_TTLCv[16384 EXT_STRUCT_CTRL(FIBZEB_TTLC)] /* [0,65535] */;
    uint32_t FIBZEB_TTLFM /* [1,256] */;
    uint32_t FIBZEB_TTLFMI[256 EXT_STRUCT_CTRL(FIBZEB_TTLFM)] /* [1,256] */;
    uint32_t FIBZEB_TTLFME[256 EXT_STRUCT_CTRL(FIBZEB_TTLFM)] /* [1,16384] */;
    uint32_t FIBZEB_TTLF /* [0,16384] */;
    uint32_t FIBZEB_TTLFv[16384 EXT_STRUCT_CTRL(FIBZEB_TTLF)] /* [0,65535] */;
    uint32_t FIBZEB_TTTCM /* [1,256] */;
    uint32_t FIBZEB_TTTCMI[256 EXT_STRUCT_CTRL(FIBZEB_TTTCM)] /* [1,256] */;
    uint32_t FIBZEB_TTTCME[256 EXT_STRUCT_CTRL(FIBZEB_TTTCM)] /* [1,16384] */;
    uint32_t FIBZEB_TTTC /* [0,16384] */;
    uint32_t FIBZEB_TTTCv[16384 EXT_STRUCT_CTRL(FIBZEB_TTTC)] /* [0,65535] */;
    uint32_t FIBZEB_TTTFM /* [1,256] */;
    uint32_t FIBZEB_TTTFMI[256 EXT_STRUCT_CTRL(FIBZEB_TTTFM)] /* [1,256] */;
    uint32_t FIBZEB_TTTFME[256 EXT_STRUCT_CTRL(FIBZEB_TTTFM)] /* [1,16384] */;
    uint32_t FIBZEB_TTTF /* [0,16384] */;
    uint32_t FIBZEB_TTTFv[16384 EXT_STRUCT_CTRL(FIBZEB_TTTF)] /* [0,65535] */;
    uint32_t FIBZEB_TRIGCM /* [1,4] */;
    uint32_t FIBZEB_TRIGCMI[4 EXT_STRUCT_CTRL(FIBZEB_TRIGCM)] /* [1,4] */;
    uint32_t FIBZEB_TRIGCME[4 EXT_STRUCT_CTRL(FIBZEB_TRIGCM)] /* [1,4] */;
    uint32_t FIBZEB_TRIGC /* [0,4] */;
    uint32_t FIBZEB_TRIGCv[4 EXT_STRUCT_CTRL(FIBZEB_TRIGC)] /* [0,65535] */;
    uint32_t FIBZEB_TRIGFM /* [1,4] */;
    uint32_t FIBZEB_TRIGFMI[4 EXT_STRUCT_CTRL(FIBZEB_TRIGFM)] /* [1,4] */;
    uint32_t FIBZEB_TRIGFME[4 EXT_STRUCT_CTRL(FIBZEB_TRIGFM)] /* [1,4] */;
    uint32_t FIBZEB_TRIGF /* [0,4] */;
    uint32_t FIBZEB_TRIGFv[4 EXT_STRUCT_CTRL(FIBZEB_TRIGF)] /* [0,65535] */;

} EXT_STR_h101_FIBZEB;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_FIBZEB_onion_t
{
    /* RAW */
    uint32_t FIBZEB_TBLCM;
    uint32_t FIBZEB_TBLCMI[256 /* FIBZEB_TBLCM */];
    uint32_t FIBZEB_TBLCME[256 /* FIBZEB_TBLCM */];
    uint32_t FIBZEB_TBLC;
    uint32_t FIBZEB_TBLCv[16384 /* FIBZEB_TBLC */];
    uint32_t FIBZEB_TBLFM;
    uint32_t FIBZEB_TBLFMI[256 /* FIBZEB_TBLFM */];
    uint32_t FIBZEB_TBLFME[256 /* FIBZEB_TBLFM */];
    uint32_t FIBZEB_TBLF;
    uint32_t FIBZEB_TBLFv[16384 /* FIBZEB_TBLF */];
    uint32_t FIBZEB_TBTCM;
    uint32_t FIBZEB_TBTCMI[256 /* FIBZEB_TBTCM */];
    uint32_t FIBZEB_TBTCME[256 /* FIBZEB_TBTCM */];
    uint32_t FIBZEB_TBTC;
    uint32_t FIBZEB_TBTCv[16384 /* FIBZEB_TBTC */];
    uint32_t FIBZEB_TBTFM;
    uint32_t FIBZEB_TBTFMI[256 /* FIBZEB_TBTFM */];
    uint32_t FIBZEB_TBTFME[256 /* FIBZEB_TBTFM */];
    uint32_t FIBZEB_TBTF;
    uint32_t FIBZEB_TBTFv[16384 /* FIBZEB_TBTF */];
    uint32_t FIBZEB_TTLCM;
    uint32_t FIBZEB_TTLCMI[256 /* FIBZEB_TTLCM */];
    uint32_t FIBZEB_TTLCME[256 /* FIBZEB_TTLCM */];
    uint32_t FIBZEB_TTLC;
    uint32_t FIBZEB_TTLCv[16384 /* FIBZEB_TTLC */];
    uint32_t FIBZEB_TTLFM;
    uint32_t FIBZEB_TTLFMI[256 /* FIBZEB_TTLFM */];
    uint32_t FIBZEB_TTLFME[256 /* FIBZEB_TTLFM */];
    uint32_t FIBZEB_TTLF;
    uint32_t FIBZEB_TTLFv[16384 /* FIBZEB_TTLF */];
    uint32_t FIBZEB_TTTCM;
    uint32_t FIBZEB_TTTCMI[256 /* FIBZEB_TTTCM */];
    uint32_t FIBZEB_TTTCME[256 /* FIBZEB_TTTCM */];
    uint32_t FIBZEB_TTTC;
    uint32_t FIBZEB_TTTCv[16384 /* FIBZEB_TTTC */];
    uint32_t FIBZEB_TTTFM;
    uint32_t FIBZEB_TTTFMI[256 /* FIBZEB_TTTFM */];
    uint32_t FIBZEB_TTTFME[256 /* FIBZEB_TTTFM */];
    uint32_t FIBZEB_TTTF;
    uint32_t FIBZEB_TTTFv[16384 /* FIBZEB_TTTF */];
    uint32_t FIBZEB_TRIGCM;
    uint32_t FIBZEB_TRIGCMI[4 /* FIBZEB_TRIGCM */];
    uint32_t FIBZEB_TRIGCME[4 /* FIBZEB_TRIGCM */];
    uint32_t FIBZEB_TRIGC;
    uint32_t FIBZEB_TRIGCv[4 /* FIBZEB_TRIGC */];
    uint32_t FIBZEB_TRIGFM;
    uint32_t FIBZEB_TRIGFMI[4 /* FIBZEB_TRIGFM */];
    uint32_t FIBZEB_TRIGFME[4 /* FIBZEB_TRIGFM */];
    uint32_t FIBZEB_TRIGF;
    uint32_t FIBZEB_TRIGFv[4 /* FIBZEB_TRIGF */];

} EXT_STR_h101_FIBZEB_onion;

/*******************************************************/

#define EXT_STR_h101_FIBZEB_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                    \
    do                                                                                                        \
    {                                                                                                         \
        ok = 1;                                                                                               \
        /* RAW */                                                                                             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TBLCM, UINT32, "FIBZEB_TBLCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBLCMI, UINT32, "FIBZEB_TBLCMI", "FIBZEB_TBLCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBLCME, UINT32, "FIBZEB_TBLCME", "FIBZEB_TBLCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TBLC, UINT32, "FIBZEB_TBLC", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBLCv, UINT32, "FIBZEB_TBLCv", "FIBZEB_TBLC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TBLFM, UINT32, "FIBZEB_TBLFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBLFMI, UINT32, "FIBZEB_TBLFMI", "FIBZEB_TBLFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBLFME, UINT32, "FIBZEB_TBLFME", "FIBZEB_TBLFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TBLF, UINT32, "FIBZEB_TBLF", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBLFv, UINT32, "FIBZEB_TBLFv", "FIBZEB_TBLF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TBTCM, UINT32, "FIBZEB_TBTCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBTCMI, UINT32, "FIBZEB_TBTCMI", "FIBZEB_TBTCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBTCME, UINT32, "FIBZEB_TBTCME", "FIBZEB_TBTCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TBTC, UINT32, "FIBZEB_TBTC", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBTCv, UINT32, "FIBZEB_TBTCv", "FIBZEB_TBTC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TBTFM, UINT32, "FIBZEB_TBTFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBTFMI, UINT32, "FIBZEB_TBTFMI", "FIBZEB_TBTFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBTFME, UINT32, "FIBZEB_TBTFME", "FIBZEB_TBTFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TBTF, UINT32, "FIBZEB_TBTF", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TBTFv, UINT32, "FIBZEB_TBTFv", "FIBZEB_TBTF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TTLCM, UINT32, "FIBZEB_TTLCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTLCMI, UINT32, "FIBZEB_TTLCMI", "FIBZEB_TTLCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTLCME, UINT32, "FIBZEB_TTLCME", "FIBZEB_TTLCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TTLC, UINT32, "FIBZEB_TTLC", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTLCv, UINT32, "FIBZEB_TTLCv", "FIBZEB_TTLC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TTLFM, UINT32, "FIBZEB_TTLFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTLFMI, UINT32, "FIBZEB_TTLFMI", "FIBZEB_TTLFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTLFME, UINT32, "FIBZEB_TTLFME", "FIBZEB_TTLFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TTLF, UINT32, "FIBZEB_TTLF", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTLFv, UINT32, "FIBZEB_TTLFv", "FIBZEB_TTLF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TTTCM, UINT32, "FIBZEB_TTTCM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTTCMI, UINT32, "FIBZEB_TTTCMI", "FIBZEB_TTTCM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTTCME, UINT32, "FIBZEB_TTTCME", "FIBZEB_TTTCM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TTTC, UINT32, "FIBZEB_TTTC", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTTCv, UINT32, "FIBZEB_TTTCv", "FIBZEB_TTTC");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TTTFM, UINT32, "FIBZEB_TTTFM", 256); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTTFMI, UINT32, "FIBZEB_TTTFMI", "FIBZEB_TTTFM");      \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTTFME, UINT32, "FIBZEB_TTTFME", "FIBZEB_TTTFM");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TTTF, UINT32, "FIBZEB_TTTF", 16384); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TTTFv, UINT32, "FIBZEB_TTTFv", "FIBZEB_TTTF");         \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TRIGCM, UINT32, "FIBZEB_TRIGCM", 4); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TRIGCMI, UINT32, "FIBZEB_TRIGCMI", "FIBZEB_TRIGCM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TRIGCME, UINT32, "FIBZEB_TRIGCME", "FIBZEB_TRIGCM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TRIGC, UINT32, "FIBZEB_TRIGC", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TRIGCv, UINT32, "FIBZEB_TRIGCv", "FIBZEB_TRIGC");      \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TRIGFM, UINT32, "FIBZEB_TRIGFM", 4); \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TRIGFMI, UINT32, "FIBZEB_TRIGFMI", "FIBZEB_TRIGFM");   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TRIGFME, UINT32, "FIBZEB_TRIGFME", "FIBZEB_TRIGFM");   \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, FIBZEB_TRIGF, UINT32, "FIBZEB_TRIGF", 4);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                \
            ok, si, offset, struct_t, printerr, FIBZEB_TRIGFv, UINT32, "FIBZEB_TRIGFv", "FIBZEB_TRIGF");      \
                                                                                                              \
    } while (0);

#endif /*__GUARD_H101_FIBZEB_EXT_H101_FIB23B_H__*/

/*******************************************************/
