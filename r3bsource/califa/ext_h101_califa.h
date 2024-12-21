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

#ifndef __GUARD_H101_CALIFA_EXT_H101_CALIFA_H__
#define __GUARD_H101_CALIFA_EXT_H101_CALIFA_H__

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

typedef struct EXT_STR_h101_CALIFA_t
{
    /* RAW */
    uint32_t CALIFA_TRGENE /* [0,2] */;
    uint32_t CALIFA_TRGENEI[2 EXT_STRUCT_CTRL(CALIFA_TRGENE)] /* [1,2] */;
    uint32_t CALIFA_TRGENEv[2 EXT_STRUCT_CTRL(CALIFA_TRGENE)] /* [0,65535] */;
    uint32_t CALIFA_ENE /* [0,5088] */;
    uint32_t CALIFA_ENEI[5088 EXT_STRUCT_CTRL(CALIFA_ENE)] /* [1,5088] */;
    uint32_t CALIFA_ENEv[5088 EXT_STRUCT_CTRL(CALIFA_ENE)] /* [0,65535] */;
    uint32_t CALIFA_TOT /* [0,5088] */;
    uint32_t CALIFA_TOTI[5088 EXT_STRUCT_CTRL(CALIFA_TOT)] /* [1,5088] */;
    uint32_t CALIFA_TOTv[5088 EXT_STRUCT_CTRL(CALIFA_TOT)] /* [0,65535] */;
    uint32_t CALIFA_TSLSB /* [0,5088] */;
    uint32_t CALIFA_TSLSBI[5088 EXT_STRUCT_CTRL(CALIFA_TSLSB)] /* [1,5088] */;
    uint32_t CALIFA_TSLSBv[5088 EXT_STRUCT_CTRL(CALIFA_TSLSB)] /* [-1,-1] */;
    uint32_t CALIFA_TSMSB /* [0,5088] */;
    uint32_t CALIFA_TSMSBI[5088 EXT_STRUCT_CTRL(CALIFA_TSMSB)] /* [1,5088] */;
    uint32_t CALIFA_TSMSBv[5088 EXT_STRUCT_CTRL(CALIFA_TSMSB)] /* [-1,-1] */;
    uint32_t CALIFA_NF /* [0,5088] */;
    uint32_t CALIFA_NFI[5088 EXT_STRUCT_CTRL(CALIFA_NF)] /* [1,5088] */;
    uint32_t CALIFA_NFv[5088 EXT_STRUCT_CTRL(CALIFA_NF)] /* [0,65535] */;
    uint32_t CALIFA_NS /* [0,5088] */;
    uint32_t CALIFA_NSI[5088 EXT_STRUCT_CTRL(CALIFA_NS)] /* [1,5088] */;
    uint32_t CALIFA_NSv[5088 EXT_STRUCT_CTRL(CALIFA_NS)] /* [0,65535] */;
    uint32_t CALIFA_TRG /* [0,5088] */;
    uint32_t CALIFA_TRGI[5088 EXT_STRUCT_CTRL(CALIFA_TRG)] /* [1,5088] */;
    uint32_t CALIFA_TRGv[5088 EXT_STRUCT_CTRL(CALIFA_TRG)] /* [0,65535] */;
    uint32_t CALIFA_OV /* [0,5088] */;
    uint32_t CALIFA_OVI[5088 EXT_STRUCT_CTRL(CALIFA_OV)] /* [1,5088] */;
    uint32_t CALIFA_OVv[5088 EXT_STRUCT_CTRL(CALIFA_OV)] /* [-1,-1] */;
    uint32_t CALIFA_WRTS_T1 /* [0,5088] */;
    uint32_t CALIFA_WRTS_T1I[5088 EXT_STRUCT_CTRL(CALIFA_WRTS_T1)] /* [1,5088] */;
    uint32_t CALIFA_WRTS_T1v[5088 EXT_STRUCT_CTRL(CALIFA_WRTS_T1)] /* [0,65535] */;
    uint32_t CALIFA_WRTS_T2 /* [0,5088] */;
    uint32_t CALIFA_WRTS_T2I[5088 EXT_STRUCT_CTRL(CALIFA_WRTS_T2)] /* [1,5088] */;
    uint32_t CALIFA_WRTS_T2v[5088 EXT_STRUCT_CTRL(CALIFA_WRTS_T2)] /* [0,65535] */;
    uint32_t CALIFA_WRTS_T3 /* [0,5088] */;
    uint32_t CALIFA_WRTS_T3I[5088 EXT_STRUCT_CTRL(CALIFA_WRTS_T3)] /* [1,5088] */;
    uint32_t CALIFA_WRTS_T3v[5088 EXT_STRUCT_CTRL(CALIFA_WRTS_T3)] /* [0,65535] */;
    uint32_t CALIFA_WRTS_T4 /* [0,5088] */;
    uint32_t CALIFA_WRTS_T4I[5088 EXT_STRUCT_CTRL(CALIFA_WRTS_T4)] /* [1,5088] */;
    uint32_t CALIFA_WRTS_T4v[5088 EXT_STRUCT_CTRL(CALIFA_WRTS_T4)] /* [0,65535] */;
    uint32_t CALIFA_PILEUP /* [0,5088] */;
    uint32_t CALIFA_PILEUPI[5088 EXT_STRUCT_CTRL(CALIFA_PILEUP)] /* [1,5088] */;
    uint32_t CALIFA_PILEUPv[5088 EXT_STRUCT_CTRL(CALIFA_PILEUP)] /* [0,65535] */;
    uint32_t CALIFA_DISCARD /* [0,5088] */;
    uint32_t CALIFA_DISCARDI[5088 EXT_STRUCT_CTRL(CALIFA_DISCARD)] /* [1,5088] */;
    uint32_t CALIFA_DISCARDv[5088 EXT_STRUCT_CTRL(CALIFA_DISCARD)] /* [0,65535] */;

} EXT_STR_h101_CALIFA;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_CALIFA_onion_t
{
    /* RAW */
    uint32_t CALIFA_TRGENE;
    uint32_t CALIFA_TRGENEI[2 /* CALIFA_TRGENE */];
    uint32_t CALIFA_TRGENEv[2 /* CALIFA_TRGENE */];
    uint32_t CALIFA_ENE;
    uint32_t CALIFA_ENEI[5088 /* CALIFA_ENE */];
    uint32_t CALIFA_ENEv[5088 /* CALIFA_ENE */];
    uint32_t CALIFA_TOT;
    uint32_t CALIFA_TOTI[5088 /* CALIFA_TOT */];
    uint32_t CALIFA_TOTv[5088 /* CALIFA_TOT */];
    uint32_t CALIFA_TSLSB;
    uint32_t CALIFA_TSLSBI[5088 /* CALIFA_TSLSB */];
    uint32_t CALIFA_TSLSBv[5088 /* CALIFA_TSLSB */];
    uint32_t CALIFA_TSMSB;
    uint32_t CALIFA_TSMSBI[5088 /* CALIFA_TSMSB */];
    uint32_t CALIFA_TSMSBv[5088 /* CALIFA_TSMSB */];
    uint32_t CALIFA_NF;
    uint32_t CALIFA_NFI[5088 /* CALIFA_NF */];
    uint32_t CALIFA_NFv[5088 /* CALIFA_NF */];
    uint32_t CALIFA_NS;
    uint32_t CALIFA_NSI[5088 /* CALIFA_NS */];
    uint32_t CALIFA_NSv[5088 /* CALIFA_NS */];
    uint32_t CALIFA_TRG;
    uint32_t CALIFA_TRGI[5088 /* CALIFA_TRG */];
    uint32_t CALIFA_TRGv[5088 /* CALIFA_TRG */];
    uint32_t CALIFA_OV;
    uint32_t CALIFA_OVI[5088 /* CALIFA_OV */];
    uint32_t CALIFA_OVv[5088 /* CALIFA_OV */];
    struct
    {
        uint32_t _;
        uint32_t I[5088 /*  */];
        uint32_t v[5088 /*  */];
    } CALIFA_WRTS_T[4];
    uint32_t CALIFA_PILEUP;
    uint32_t CALIFA_PILEUPI[5088 /* CALIFA_PILEUP */];
    uint32_t CALIFA_PILEUPv[5088 /* CALIFA_PILEUP */];
    uint32_t CALIFA_DISCARD;
    uint32_t CALIFA_DISCARDI[5088 /* CALIFA_DISCARD */];
    uint32_t CALIFA_DISCARDv[5088 /* CALIFA_DISCARD */];

} EXT_STR_h101_CALIFA_onion;

/*******************************************************/

#define EXT_STR_h101_CALIFA_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                           \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_TRGENE, UINT32, "CALIFA_TRGENE", 2);        \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_TRGENEI, UINT32, "CALIFA_TRGENEI", "CALIFA_TRGENE");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_TRGENEv, UINT32, "CALIFA_TRGENEv", "CALIFA_TRGENE");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_ENE, UINT32, "CALIFA_ENE", 5088);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_ENEI, UINT32, "CALIFA_ENEI", "CALIFA_ENE"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_ENEv, UINT32, "CALIFA_ENEv", "CALIFA_ENE"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_TOT, UINT32, "CALIFA_TOT", 5088);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_TOTI, UINT32, "CALIFA_TOTI", "CALIFA_TOT"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_TOTv, UINT32, "CALIFA_TOTv", "CALIFA_TOT"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_TSLSB, UINT32, "CALIFA_TSLSB", 5088);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_TSLSBI, UINT32, "CALIFA_TSLSBI", "CALIFA_TSLSB");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_TSLSBv, UINT32, "CALIFA_TSLSBv", "CALIFA_TSLSB");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_TSMSB, UINT32, "CALIFA_TSMSB", 5088);       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_TSMSBI, UINT32, "CALIFA_TSMSBI", "CALIFA_TSMSB");             \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_TSMSBv, UINT32, "CALIFA_TSMSBv", "CALIFA_TSMSB");             \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_NF, UINT32, "CALIFA_NF", 5088);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_NFI, UINT32, "CALIFA_NFI", "CALIFA_NF");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_NFv, UINT32, "CALIFA_NFv", "CALIFA_NF");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_NS, UINT32, "CALIFA_NS", 5088);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_NSI, UINT32, "CALIFA_NSI", "CALIFA_NS");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_NSv, UINT32, "CALIFA_NSv", "CALIFA_NS");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_TRG, UINT32, "CALIFA_TRG", 5088);           \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_TRGI, UINT32, "CALIFA_TRGI", "CALIFA_TRG"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_TRGv, UINT32, "CALIFA_TRGv", "CALIFA_TRG"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_OV, UINT32, "CALIFA_OV", 5088);             \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_OVI, UINT32, "CALIFA_OVI", "CALIFA_OV");    \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, CALIFA_OVv, UINT32, "CALIFA_OVv", "CALIFA_OV");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_WRTS_T1, UINT32, "CALIFA_WRTS_T1", 5088);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_WRTS_T1I, UINT32, "CALIFA_WRTS_T1I", "CALIFA_WRTS_T1");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_WRTS_T1v, UINT32, "CALIFA_WRTS_T1v", "CALIFA_WRTS_T1");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_WRTS_T2, UINT32, "CALIFA_WRTS_T2", 5088);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_WRTS_T2I, UINT32, "CALIFA_WRTS_T2I", "CALIFA_WRTS_T2");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_WRTS_T2v, UINT32, "CALIFA_WRTS_T2v", "CALIFA_WRTS_T2");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_WRTS_T3, UINT32, "CALIFA_WRTS_T3", 5088);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_WRTS_T3I, UINT32, "CALIFA_WRTS_T3I", "CALIFA_WRTS_T3");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_WRTS_T3v, UINT32, "CALIFA_WRTS_T3v", "CALIFA_WRTS_T3");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_WRTS_T4, UINT32, "CALIFA_WRTS_T4", 5088);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_WRTS_T4I, UINT32, "CALIFA_WRTS_T4I", "CALIFA_WRTS_T4");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_WRTS_T4v, UINT32, "CALIFA_WRTS_T4v", "CALIFA_WRTS_T4");       \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_PILEUP, UINT32, "CALIFA_PILEUP", 5088);     \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_PILEUPI, UINT32, "CALIFA_PILEUPI", "CALIFA_PILEUP");          \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_PILEUPv, UINT32, "CALIFA_PILEUPv", "CALIFA_PILEUP");          \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, CALIFA_DISCARD, UINT32, "CALIFA_DISCARD", 5088);   \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_DISCARDI, UINT32, "CALIFA_DISCARDI", "CALIFA_DISCARD");       \
        EXT_STR_ITEM_INFO_ZZP(                                                                                       \
            ok, si, offset, struct_t, printerr, CALIFA_DISCARDv, UINT32, "CALIFA_DISCARDv", "CALIFA_DISCARD");       \
                                                                                                                     \
    } while (0);

#endif /*__GUARD_H101_CALIFA_EXT_H101_CALIFA_H__*/

/*******************************************************/
