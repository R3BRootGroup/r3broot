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

// This is generic class for the fibers 30, 31, 32 33 used in s494 experiment
// modified from  R3BBunchedFiber code
// Boiler-plate code for bunched fiber detectors with multi-hit electronics on
// both ends.
// In the UCESB unpacker, signals should be mapped like:
//  SIGNAL(WONDERFULFIBER_TXYZ, ...);
// where:
//  X = M(Multi-anode PMT) or S(Single PMT) side of the detector.
//  Y = L(eading) or T(railing) edge of signal.
//  Z = C(oarse) or F(ine) time from electronics.
// For example:
//  SIGNAL(FIBONE_TMLC, ...);
// Note that sub-detectors (e.g. cloned parts side by side) are not mentioned
// in the mapping, since channels should anyhow not be intertwined, e.g. the
// first sub occupies the first 256/2 channels, the second the next 256/2
// channels and so on. The constructor however takes the number of subs to
// simplify the numbers in the user code a little.
// Layers are treated as separate detectors, e.g. Fi1a and Fi1b.
//
// Note that MappedData refers to a single multi-hit electronics channel
// connected to some bunch of fibers, NOT one fiber, so that Mapped2Cal
// calibrates only once per electronics channel.
//
// Cal2Hit however does handle electronics channels -> fibers.
//

#pragma once

#include "R3BLogger.h"
#include "R3BReader.h"

#include <TString.h>

class TClonesArray;

//
// Linking macros for multi-hit UCESB data for fiber detectors to generic
// pointers, makes it easier to generalize different detector layouts.
// Example:
//  R3B_FIBER_MAPMT_UCESB_LINK(fData->FIBONE);
//
// Something similar could spread to other detectors to generalize all/most
// unpackers to reduce copy-pasting readers, and to make sure error checking
// and error messages are consistent throughout the framework.
//

#define LENGTH(x) (sizeof x / sizeof *x)
#define R3B_FIBER_MAPMT_INIT_TRIG(NAME, data)                                                  \
    do                                                                                         \
    {                                                                                          \
        int ok;                                                                                \
        R3BLOG(info, "");                                                                      \
        EXT_STR_h101_##NAME##_ITEMS_INFO(ok, *a_struct_info, fOffset, EXT_STR_h101_##NAME, 0); \
        if (!ok)                                                                               \
        {                                                                                      \
            R3BLOG(error, "Failed to setup UCESB structure information.");                     \
            return kFALSE;                                                                     \
        }                                                                                      \
        memset(data, 0, sizeof *data);                                                         \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[0][0][0], data->NAME##_TBLC);                 \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[0][0][1], data->NAME##_TBLF);                 \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[0][1][0], data->NAME##_TBTC);                 \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[0][1][1], data->NAME##_TBTF);                 \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[1][0][0], data->NAME##_TTLC);                 \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[1][0][1], data->NAME##_TTLF);                 \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[1][1][0], data->NAME##_TTTC);                 \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[1][1][1], data->NAME##_TTTF);                 \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[2][0][0], data->NAME##_TRIGC);                \
        R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(fMHL[2][0][1], data->NAME##_TRIGF);                \
        return R3BFiberMAPMTReader::Init();                                                    \
    } while (0)

#define R3B_FIBER_MAPMT_UCESB_LINK_MULTIHIT(dst, src)     \
    do                                                    \
    {                                                     \
        dst._M = &src##M;                                 \
        dst._MI = src##MI;                                \
        dst._MI_len = LENGTH(src##MI);                    \
        dst._ME = src##ME;                                \
        dst._ME_len = LENGTH(src##ME);                    \
        if (dst._MI_len != dst._ME_len)                   \
        {                                                 \
            LOG(error) << "UCESB multi-hit array error."; \
            exit(EXIT_FAILURE);                           \
        }                                                 \
        dst._ = &src;                                     \
        dst._v = src##v;                                  \
        dst._v_len = LENGTH(src##v);                      \
    } while (0)

class R3BFiberMAPMTReader : public R3BReader
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

    R3BFiberMAPMTReader(char const*, UInt_t, size_t);
    Bool_t Init();
    Bool_t R3BRead();
    void Reset();

    // Accessor to select online mode
    inline void SetOnline(Bool_t option) { fOnline = option; }

  protected:
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    // [0=bottom,1=top,2=MAPMT-trig][0=leading,1=trailing][0=coarse,1=fine].
    UCESBMultiHitLink fMHL[3][2][2];

  private:
    TString fShortName; // e.g. "Fi1a", "Fi2", ...
    UInt_t fFiberNum;
    TClonesArray* fMappedArray;

  public:
    ClassDef(R3BFiberMAPMTReader, 1);
};
