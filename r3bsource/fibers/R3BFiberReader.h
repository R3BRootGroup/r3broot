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

// This is generic class for the fiber detectors
// Boiler-plate code for bunched fiber detectors with multi-hit electronics on
// both ends.
// In the UCESB unpacker, signals should be mapped like:
//  SIGNAL(WONDERFULFIBER_TXYZ, ...);
// where:
//  X = M(Multi-anode PMT) side of the detector.
//  Y = L(eading) or T(railing) edge of signal.
//  Z = C(oarse) or F(ine) time from electronics.

#pragma once

#include "R3BLogger.h"
#include "R3BReader.h"

#include <Rtypes.h>
#include <TString.h>

class TClonesArray;

struct EXT_STR_h101_FIBZEA_t;
typedef struct EXT_STR_h101_FIBZEA_t EXT_STR_h101_FIBZEA;
typedef struct EXT_STR_h101_FIBZEA_onion_t EXT_STR_h101_FIBZEA_onion;

struct EXT_STR_h101_FIBZEB_t;
typedef struct EXT_STR_h101_FIBZEB_t EXT_STR_h101_FIBZEB;
typedef struct EXT_STR_h101_FIBZEB_onion_t EXT_STR_h101_FIBZEB_onion;

struct EXT_STR_h101_FIBTHREEA_t;
typedef struct EXT_STR_h101_FIBTHREEA_t EXT_STR_h101_FIBTHREEA;
typedef struct EXT_STR_h101_FIBTHREEA_onion_t EXT_STR_h101_FIBTHREEA_onion;

struct EXT_STR_h101_FIBTHREEB_t;
typedef struct EXT_STR_h101_FIBTHREEB_t EXT_STR_h101_FIBTHREEB;
typedef struct EXT_STR_h101_FIBTHREEB_onion_t EXT_STR_h101_FIBTHREEB_onion;

struct EXT_STR_h101_FIBSEVEN_t;
typedef struct EXT_STR_h101_FIBSEVEN_t EXT_STR_h101_FIBSEVEN;
typedef struct EXT_STR_h101_FIBSEVEN_onion_t EXT_STR_h101_FIBSEVEN_onion;

struct EXT_STR_h101_FIBEIGHT_t;
typedef struct EXT_STR_h101_FIBEIGHT_t EXT_STR_h101_FIBEIGHT;
typedef struct EXT_STR_h101_FIBEIGHT_onion_t EXT_STR_h101_FIBEIGHT_onion;

struct EXT_STR_h101_FIBTEN_t;
typedef struct EXT_STR_h101_FIBTEN_t EXT_STR_h101_FIBTEN;
typedef struct EXT_STR_h101_FIBTEN_onion_t EXT_STR_h101_FIBTEN_onion;

struct EXT_STR_h101_FIBELEVEN_t;
typedef struct EXT_STR_h101_FIBELEVEN_t EXT_STR_h101_FIBELEVEN;
typedef struct EXT_STR_h101_FIBELEVEN_onion_t EXT_STR_h101_FIBELEVEN_onion;

struct EXT_STR_h101_FIBTWELVE_t;
typedef struct EXT_STR_h101_FIBTWELVE_t EXT_STR_h101_FIBTWELVE;
typedef struct EXT_STR_h101_FIBTWELVE_onion_t EXT_STR_h101_FIBTWELVE_onion;

struct EXT_STR_h101_FIBTHIRTEEN_t;
typedef struct EXT_STR_h101_FIBTHIRTEEN_t EXT_STR_h101_FIBTHIRTEEN;
typedef struct EXT_STR_h101_FIBTHIRTEEN_onion_t EXT_STR_h101_FIBTHIRTEEN_onion;

struct EXT_STR_h101_FIBEO_t;
typedef struct EXT_STR_h101_FIBEO_t EXT_STR_h101_FIBEO;
typedef struct EXT_STR_h101_FIBEO_onion_t EXT_STR_h101_FIBEO_onion;

struct EXT_STR_h101_FIBEI_t;
typedef struct EXT_STR_h101_FIBEI_t EXT_STR_h101_FIBEI;
typedef struct EXT_STR_h101_FIBEI_onion_t EXT_STR_h101_FIBEI_onion;

struct EXT_STR_h101_FIBEZ_t;
typedef struct EXT_STR_h101_FIBEZ_t EXT_STR_h101_FIBEZ;
typedef struct EXT_STR_h101_FIBEZ_onion_t EXT_STR_h101_FIBEZ_onion;

struct EXT_STR_h101_FIBEE_t;
typedef struct EXT_STR_h101_FIBEE_t EXT_STR_h101_FIBEE;
typedef struct EXT_STR_h101_FIBEE_onion_t EXT_STR_h101_FIBEE_onion;

struct EXT_STR_h101_FIB_t;
typedef struct EXT_STR_h101_FIB_t EXT_STR_h101_FIB;
typedef struct EXT_STR_h101_FIB_onion_t EXT_STR_h101_FIB_onion;

class ext_data_struct_info;

#define R3BBUNCHED_FIBER_INIT_BEGIN_(NAME)                                                 \
    Int_t ok;                                                                              \
    EXT_STR_h101_##NAME##_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_##NAME, 0); \
    if (!ok)                                                                               \
    {                                                                                      \
        R3BLOG(error, "Failed to setup UCESB structure information.");                     \
        return kFALSE;                                                                     \
    }

#define R3BBUNCHED_FIBER_INIT(NAME, data)                                \
    do                                                                   \
    {                                                                    \
        R3BBUNCHED_FIBER_INIT_BEGIN_(NAME);                              \
        memset(data, 0, sizeof *data);                                   \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][0][0], data->NAME##_TMLC); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][0][1], data->NAME##_TMLF); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][1][0], data->NAME##_TMTC); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][1][1], data->NAME##_TMTF); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][0][0], data->NAME##_TSLC); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][0][1], data->NAME##_TSLF); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][1][0], data->NAME##_TSTC); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][1][1], data->NAME##_TSTF); \
    } while (0)

#define R3BBUNCHED_FIBER_INIT_MAPMT_ONLY(NAME, data)                     \
    do                                                                   \
    {                                                                    \
        R3BBUNCHED_FIBER_INIT_BEGIN_(NAME);                              \
        memset(data, 0, sizeof *data);                                   \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][0][0], data->NAME##_TMLC); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][0][1], data->NAME##_TMLF); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][1][0], data->NAME##_TMTC); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][1][1], data->NAME##_TMTF); \
    } while (0)

#define R3BBUNCHED_FIBER_INIT_MAPMT_TRIG(NAME, data)                        \
    do                                                                      \
    {                                                                       \
        R3BBUNCHED_FIBER_INIT_BEGIN_(NAME);                                 \
        memset(data, 0, sizeof *data);                                      \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][0][0], data->NAME##_TMLC);    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][0][1], data->NAME##_TMLF);    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][1][0], data->NAME##_TMTC);    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][1][1], data->NAME##_TMTF);    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][0][0], data->NAME##_TSLC);    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][0][1], data->NAME##_TSLF);    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][1][0], data->NAME##_TSTC);    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][1][1], data->NAME##_TSTF);    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[2][0][0], data->NAME##_TRIGMLC); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[2][0][1], data->NAME##_TRIGMLF); \
    } while (0)

#define R3B_MAPMT_INIT_TRIG(NAME, data)                                   \
    do                                                                    \
    {                                                                     \
        R3BBUNCHED_FIBER_INIT_BEGIN_(NAME);                               \
        memset(data, 0, sizeof *data);                                    \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][0][0], data->NAME##_TBLC);  \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][0][1], data->NAME##_TBLF);  \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][1][0], data->NAME##_TBTC);  \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[0][1][1], data->NAME##_TBTF);  \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][0][0], data->NAME##_TTLC);  \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][0][1], data->NAME##_TTLF);  \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][1][0], data->NAME##_TTTC);  \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[1][1][1], data->NAME##_TTTF);  \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[2][0][0], data->NAME##_TRIGC); \
        R3B_FIBER_UCESB_LINK_MULTIHIT(fMHL[2][0][1], data->NAME##_TRIGF); \
    } while (0)

#define R3B_FIBER_UCESB_LINK(dst, src)                      \
    do                                                      \
    {                                                       \
        dst._ = &src;                                       \
        dst._MI = src##I;                                   \
        dst._MI_len = LENGTH(src##I);                       \
        dst._v = src##v;                                    \
        dst._v_len = LENGTH(src##v);                        \
        if (dst._MI_len != dst._v_len)                      \
        {                                                   \
            R3BLOG(error, "UCESB single-hit array error."); \
            exit(EXIT_FAILURE);                             \
        }                                                   \
    } while (0)

#define LENGTH(x) (sizeof x / sizeof *x)
#define R3B_FIBER_UCESB_LINK_MULTIHIT(dst, src)            \
    do                                                     \
    {                                                      \
        dst._M = &src##M;                                  \
        dst._MI = src##MI;                                 \
        dst._MI_len = LENGTH(src##MI);                     \
        dst._ME = src##ME;                                 \
        dst._ME_len = LENGTH(src##ME);                     \
        if (dst._MI_len != dst._ME_len)                    \
        {                                                  \
            R3BLOG(error, "UCESB multi-hit array error."); \
            exit(EXIT_FAILURE);                            \
        }                                                  \
        dst._ = &src;                                      \
        dst._v = src##v;                                   \
        dst._v_len = LENGTH(src##v);                       \
    } while (0)

class R3BFiberReader : public R3BReader
{
  public:
    struct UCESBMultiHitLink
    {
        uint32_t const* _M;
        uint32_t const* _MI;
        size_t _MI_len;
        uint32_t const* _ME;
        size_t _ME_len;
        uint32_t const* _;
        uint32_t const* _v;
        size_t _v_len;
    };

    // Standard constructors
    R3BFiberReader(char const*, UInt_t, EXT_STR_h101_FIBZEA_onion*, size_t);
    R3BFiberReader(char const*, UInt_t, EXT_STR_h101_FIBZEB_onion*, size_t);
    R3BFiberReader(char const*, EXT_STR_h101_FIBTHREEA_onion*, size_t, UInt_t, UInt_t, UInt_t);
    R3BFiberReader(char const*, EXT_STR_h101_FIBTHREEB_onion*, size_t, UInt_t, UInt_t, UInt_t);
    R3BFiberReader(char const*, EXT_STR_h101_FIBSEVEN_onion*, size_t, UInt_t, UInt_t, UInt_t);
    R3BFiberReader(char const*, EXT_STR_h101_FIBEIGHT_onion*, size_t, UInt_t, UInt_t, UInt_t);
    R3BFiberReader(char const*, EXT_STR_h101_FIBTEN_onion*, size_t, UInt_t, UInt_t, UInt_t);
    R3BFiberReader(char const*, EXT_STR_h101_FIBELEVEN_onion*, size_t, UInt_t, UInt_t, UInt_t);
    R3BFiberReader(char const*, EXT_STR_h101_FIBTWELVE_onion*, size_t, UInt_t, UInt_t, UInt_t);
    R3BFiberReader(char const*, EXT_STR_h101_FIBTHIRTEEN_onion*, size_t, UInt_t, UInt_t, UInt_t);
    R3BFiberReader(char const*, UInt_t, EXT_STR_h101_FIBEO_onion*, size_t);
    R3BFiberReader(char const*, UInt_t, EXT_STR_h101_FIBEI_onion*, size_t);
    R3BFiberReader(char const*, UInt_t, EXT_STR_h101_FIBEZ_onion*, size_t);
    R3BFiberReader(char const*, UInt_t, EXT_STR_h101_FIBEE_onion*, size_t);

    // Function to set up the ucesb data structure of SPMT trigger info
    void SetUcesbSPMTTrigStructure(EXT_STR_h101_FIB_onion*, size_t);

    // Destructor
    virtual ~R3BFiberReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    inline void SetOnline(Bool_t option) { fOnline = option; }

  private:
    TString fShortName; // e.g. "Fi10", "Fi30", "Fi31", ...
    UInt_t fFiberNum;
    // [0=MAPMT,1=SPMT].
    UInt_t fChannelNum[2];
    // Output array
    TClonesArray* fMappedArray;
    // Reader specific data structure from ucesb
    EXT_STR_h101_FIBZEA_onion* fData23a;
    EXT_STR_h101_FIBZEB_onion* fData23b;
    EXT_STR_h101_FIBTHREEA_onion* fData3a;
    EXT_STR_h101_FIBTHREEB_onion* fData3b;
    EXT_STR_h101_FIBSEVEN_onion* fData7;
    EXT_STR_h101_FIBEIGHT_onion* fData8;
    EXT_STR_h101_FIBTEN_onion* fData10;
    EXT_STR_h101_FIBELEVEN_onion* fData11;
    EXT_STR_h101_FIBTWELVE_onion* fData12;
    EXT_STR_h101_FIBTHIRTEEN_onion* fData13;
    EXT_STR_h101_FIBEO_onion* fData30;
    EXT_STR_h101_FIBEI_onion* fData31;
    EXT_STR_h101_FIBEZ_onion* fData32;
    EXT_STR_h101_FIBEE_onion* fData33;
    // Reader specific data structure for SPMT trigger info
    // Only for fibers 10, 11, 12 and 13
    EXT_STR_h101_FIB_onion* fDataSPMTTrig;

    // Data offsets
    size_t fOffset;
    size_t fOffsetspmttrig;
    // Don't store data for online
    Bool_t fOnline;
    // [0=bottom,1=top,2=MAPMT-trig][0=leading,1=trailing][0=coarse,1=fine].
    UCESBMultiHitLink fMHL[3][2][2];

  public:
    ClassDefOverride(R3BFiberReader, 1);
};
