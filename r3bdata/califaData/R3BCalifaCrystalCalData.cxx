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

#include "R3BCalifaCrystalCalData.h"
#include <fmt/core.h>

R3BCalifaCrystalCalData::R3BCalifaCrystalCalData(uint16_t cryid,
                                                 double energy,
                                                 double Nf,
                                                 double Ns,
                                                 ULong64_t time,
                                                 double tot_energy)
    : fCrystalId(cryid)
    , fEnergy(energy)
    , fNf(Nf)
    , fNs(Ns)
    , fTime(time)
    , fToT_Energy(tot_energy)
{
}

R3BCalifaCrystalCalData::R3BCalifaCrystalCalData(const R3BCalifaCrystalCalData& right)
    : R3BCalifaCrystalCalData(right.fCrystalId, right.fEnergy, right.fNf, right.fNs, right.fTime, right.fToT_Energy)
{
}

std::string R3BCalifaCrystalCalData::toString() const
{
    return fmt::format("CrystalId: {}, Energy(MeV): {:.3f}, Nf: {:.3f}, Ns: {:.3f}, Time: {}, ToT: {:.3f}",
                       GetCrystalId(),
                       GetEnergy(),
                       GetNf(),
                       GetNs(),
                       GetTime(),
                       GetToTEnergy());
}

void R3BCalifaCrystalCalData::Print(const Option_t*) const { std::cout << *this << std::endl; }

std::ostream& operator<<(std::ostream& os, const R3BCalifaCrystalCalData& data)
{
    os << data.toString();
    return os;
}

ClassImp(R3BCalifaCrystalCalData)
