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

#include "R3BCalifaMappedData.h"
#include <fmt/core.h>

R3BCalifaMappedData::R3BCalifaMappedData(UShort_t crystalId,
                                         int16_t energy,
                                         int16_t nf,
                                         int16_t ns,
                                         uint64_t febextime,
                                         uint64_t wrts,
                                         uint32_t ov,
                                         uint16_t pu,
                                         uint16_t dc,
                                         uint16_t tot)
    : fCrystalId(crystalId)
    , fEnergy(energy)
    , fNf(nf)
    , fNs(ns)
    , fFebexTime(febextime)
    , fWrts(wrts)
    , fOverFlow(ov)
    , fPileup(pu)
    , fDiscard(dc)
    , fTot(tot)
{
}

std::string R3BCalifaMappedData::toString() const
{
    return fmt::format("CrystalId: {}, Energy(ch): {}, Nf: {}, Ns: {}, OverFlow: {}, Pileup: {}, Discard: {}",
                       GetCrystalId(),
                       GetEnergy(),
                       GetNf(),
                       GetNs(),
                       GetOverFlow(),
                       GetPileup(),
                       GetDiscard());
}

void R3BCalifaMappedData::Print(const Option_t*) const { std::cout << *this << std::endl; }

std::ostream& operator<<(std::ostream& os, const R3BCalifaMappedData& data)
{
    os << data.toString();
    return os;
}

ClassImp(R3BCalifaMappedData)
