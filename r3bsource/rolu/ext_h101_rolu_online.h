/********************************************************
 *
 * Structure for ext_data_fetch_event() filling.
 *
 * Do not edit - automatically generated.
 */

#ifndef __GUARD_H101_ROLU_EXT_H101_ROLU_H__
#define __GUARD_H101_ROLU_EXT_H101_ROLU_H__

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

typedef struct EXT_STR_h101_ROLU_NOTRIG_t
{
    /* RAW */
    uint32_t ROLU1TTFLM /* [1,4] */;
    uint32_t ROLU1TTFLMI[4 EXT_STRUCT_CTRL(ROLU1TTFLM)] /* [1,4] */;
    uint32_t ROLU1TTFLME[4 EXT_STRUCT_CTRL(ROLU1TTFLM)] /* [1,128] */;
    uint32_t ROLU1TTFL /* [0,128] */;
    uint32_t ROLU1TTFLv[128 EXT_STRUCT_CTRL(ROLU1TTFL)] /* [0,65535] */;
    uint32_t ROLU1TTFTM /* [1,4] */;
    uint32_t ROLU1TTFTMI[4 EXT_STRUCT_CTRL(ROLU1TTFTM)] /* [1,4] */;
    uint32_t ROLU1TTFTME[4 EXT_STRUCT_CTRL(ROLU1TTFTM)] /* [1,128] */;
    uint32_t ROLU1TTFT /* [0,128] */;
    uint32_t ROLU1TTFTv[128 EXT_STRUCT_CTRL(ROLU1TTFT)] /* [0,65535] */;
    uint32_t ROLU1TTCLM /* [1,4] */;
    uint32_t ROLU1TTCLMI[4 EXT_STRUCT_CTRL(ROLU1TTCLM)] /* [1,4] */;
    uint32_t ROLU1TTCLME[4 EXT_STRUCT_CTRL(ROLU1TTCLM)] /* [1,128] */;
    uint32_t ROLU1TTCL /* [0,128] */;
    uint32_t ROLU1TTCLv[128 EXT_STRUCT_CTRL(ROLU1TTCL)] /* [0,65535] */;
    uint32_t ROLU1TTCTM /* [1,4] */;
    uint32_t ROLU1TTCTMI[4 EXT_STRUCT_CTRL(ROLU1TTCTM)] /* [1,4] */;
    uint32_t ROLU1TTCTME[4 EXT_STRUCT_CTRL(ROLU1TTCTM)] /* [1,128] */;
    uint32_t ROLU1TTCT /* [0,128] */;
    uint32_t ROLU1TTCTv[128 EXT_STRUCT_CTRL(ROLU1TTCT)] /* [0,65535] */;

} EXT_STR_h101_ROLU_NOTRIG;

/********************************************************
 *
 * Structure with multiple levels of arrays (partially)
 * recovered (recommended):
 */

typedef struct EXT_STR_h101_ROLU_NOTRIG_onion_t
{
    /* RAW */
    struct
    {
        uint32_t TTFLM;
        uint32_t TTFLMI[4 /* TTFLM */];
        uint32_t TTFLME[4 /* TTFLM */];
        uint32_t TTFL;
        uint32_t TTFLv[128 /* TTFL */];
        uint32_t TTFTM;
        uint32_t TTFTMI[4 /* TTFTM */];
        uint32_t TTFTME[4 /* TTFTM */];
        uint32_t TTFT;
        uint32_t TTFTv[128 /* TTFT */];
        uint32_t TTCLM;
        uint32_t TTCLMI[4 /* TTCLM */];
        uint32_t TTCLME[4 /* TTCLM */];
        uint32_t TTCL;
        uint32_t TTCLv[128 /* TTCL */];
        uint32_t TTCTM;
        uint32_t TTCTMI[4 /* TTCTM */];
        uint32_t TTCTME[4 /* TTCTM */];
        uint32_t TTCT;
        uint32_t TTCTv[128 /* TTCT */];
    } ROLU[1];

} EXT_STR_h101_ROLU_NOTRIG_onion;

/*******************************************************/

#define EXT_STR_h101_ROLU_NOTRIG_ITEMS_INFO(ok, si, offset, struct_t, printerr)                                      \
    do                                                                                                               \
    {                                                                                                                \
        ok = 1;                                                                                                      \
        /* RAW */                                                                                                    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFLM, UINT32, "ROLU1TTFLM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFLMI, UINT32, "ROLU1TTFLMI", "ROLU1TTFLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFLME, UINT32, "ROLU1TTFLME", "ROLU1TTFLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFL, UINT32, "ROLU1TTFL", 128);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFLv, UINT32, "ROLU1TTFLv", "ROLU1TTFL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFTM, UINT32, "ROLU1TTFTM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFTMI, UINT32, "ROLU1TTFTMI", "ROLU1TTFTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFTME, UINT32, "ROLU1TTFTME", "ROLU1TTFTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTFT, UINT32, "ROLU1TTFT", 128);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTFTv, UINT32, "ROLU1TTFTv", "ROLU1TTFT");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCLM, UINT32, "ROLU1TTCLM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCLMI, UINT32, "ROLU1TTCLMI", "ROLU1TTCLM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCLME, UINT32, "ROLU1TTCLME", "ROLU1TTCLM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCL, UINT32, "ROLU1TTCL", 128);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCLv, UINT32, "ROLU1TTCLv", "ROLU1TTCL");    \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCTM, UINT32, "ROLU1TTCTM", 4);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCTMI, UINT32, "ROLU1TTCTMI", "ROLU1TTCTM"); \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCTME, UINT32, "ROLU1TTCTME", "ROLU1TTCTM"); \
        EXT_STR_ITEM_INFO_LIM(ok, si, offset, struct_t, printerr, ROLU1TTCT, UINT32, "ROLU1TTCT", 128);              \
        EXT_STR_ITEM_INFO_ZZP(ok, si, offset, struct_t, printerr, ROLU1TTCTv, UINT32, "ROLU1TTCTv", "ROLU1TTCT");    \
                                                                                                                     \
    } while (0);
#endif /*__GUARD_H101_ROLU_EXT_H101_ROLU_H__*/

/*******************************************************/
