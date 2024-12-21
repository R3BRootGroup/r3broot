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

#include "R3BFiberMAPMTReader.h"
#include "R3BFiberMAPMTMappedData.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <TClonesArray.h>

R3BFiberMAPMTReader::R3BFiberMAPMTReader(char const* a_name, UInt_t a_offset, size_t fiber_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fOffset(a_offset)
    , fShortName(a_name)
    , fFiberNum(fiber_num)
    , fMappedArray(new TClonesArray("R3BFiberMAPMTMappedData"))
    , fOnline(kFALSE)
{
}

Bool_t R3BFiberMAPMTReader::Init()
{
    for (size_t side_i = 0; side_i < 2; ++side_i)
    {
        for (size_t edge_i = 0; edge_i < 2; ++edge_i)
        {
            for (size_t prec_i = 0; prec_i < 2; ++prec_i)
            {
                auto const& ch = fMHL[side_i][edge_i][prec_i];
                if (fFiberNum != ch._MI_len)
                {
                    R3BLOG(fatal,
                           "Multi-hit array sizes mismatch (fFiberNum=" << fFiberNum << " != MI-len=" << ch._MI_len
                                                                        << ").");
                    return kFALSE;
                }
            }
        }
    }
    FairRootManager::Instance()->Register(fShortName + "Mapped", fShortName + " mapped data", fMappedArray, !fOnline);
    return kTRUE;
}

Bool_t R3BFiberMAPMTReader::R3BRead()
{
    R3BLOG(debug, "Read BEGIN for fib " << fShortName);
    for (size_t side_i = 0; side_i < 2; ++side_i)
    {
        for (size_t edge_i = 0; edge_i < 2; ++edge_i)
        {
            auto const& e = fMHL[side_i][edge_i];

            // Check that coarse and fine list are the same size.
            uint32_t c_M = *e[0]._M;
            uint32_t f_M = *e[1]._M;
            uint32_t c_ = *e[0]._;
            uint32_t f_ = *e[1]._;

            if (c_M != f_M || c_ != f_)
            {
                LOG(warn) << "Coarse and fine multi-hit list counts mismatch "
                             "(edge="
                          << edge_i << ";M{c=" << c_M << ",f=" << f_M << "};_{c=" << c_ << ",f=" << f_ << "}).";
                return kFALSE;
            }
            if (c_M > e[0]._MI_len || c_M > e[0]._ME_len || c_ > e[0]._v_len)
            {
                LOG(warn) << "Multi-hit indexing out of range "
                             "(edge="
                          << edge_i << ";(M=" << c_M << ")>=(MI=" << e[0]._MI_len << ",ME=" << e[0]._ME_len
                          << ");(_=" << c_ << ")>=(v=" << e[0]._v_len << ")).";
                return kFALSE;
            }

            // Simply dump edges, can only be sorted perfectly after calibration,
            // and hits should be somewhat sorted already, that helps many sorting
            // algos.
            uint32_t cur_entry = 0;
            for (uint32_t i = 0; i < c_M; i++)
            {
                uint32_t c_MI = e[0]._MI[i];
                uint32_t f_MI = e[1]._MI[i];
                uint32_t c_ME = e[0]._ME[i];
                uint32_t f_ME = e[1]._ME[i];

                //    LOG(warn)<<"Start: "<<side_i<<"; "<<edge_i<<", "<<c_MI<<", "<<f_MI<<", "<<c_ME<<", "<<f_ME;

                if (c_MI != f_MI || c_ME != f_ME)
                {
                    LOG(warn) << "Coarse and fine multi-hit data mismatch "
                                 "(edge="
                              << edge_i << ";MI{c=" << c_MI << ",f=" << f_MI << "};ME{c=" << c_ME << ",f=" << f_ME
                              << "}).";
                    return kFALSE;
                }
                for (; cur_entry < c_ME; cur_entry++)
                {
                    new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
                        R3BFiberMAPMTMappedData(side_i, c_MI, 0 == edge_i, e[0]._v[cur_entry], e[1]._v[cur_entry]);
                }
            }
        }
    }

    // MAPMT trigger times.
    {
        auto const& e = fMHL[2][0];

        // Check that coarse and fine list are the same size.
        uint32_t c_M = *e[0]._M;
        uint32_t f_M = *e[1]._M;
        uint32_t c_ = *e[0]._;
        uint32_t f_ = *e[1]._;

        if (c_ != f_)
        {
            LOG(warn) << "Coarse and fine single-hit list counts mismatch "
                         "(_{c="
                      << c_ << ",f=" << f_ << "}).";
            return kFALSE;
        }

        // Simply dump edges, can only be sorted perfectly after calibration,
        // and hits should be somewhat sorted already, that helps many sorting
        // algos.
        uint32_t cur_entry = 0;
        for (uint32_t i = 0; i < c_M; i++)
        {
            uint32_t c_MI = e[0]._MI[i];
            uint32_t f_MI = e[1]._MI[i];
            uint32_t c_ME = e[0]._ME[i];
            uint32_t f_ME = e[1]._ME[i];
            // Only take the first hit per trigger channel, seems like it can
            // have multi-hits.
            if (cur_entry < c_ME)
            {
                new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
                    R3BFiberMAPMTMappedData(2, c_MI, true, e[0]._v[cur_entry], e[1]._v[cur_entry]);
            }
            cur_entry = c_ME;
        }
    }
    // LOG(error) << "R3BFiberMAPMTReader::Read END";
    return kTRUE;
}

void R3BFiberMAPMTReader::Reset() { fMappedArray->Clear(); }

ClassImp(R3BFiberMAPMTReader)
