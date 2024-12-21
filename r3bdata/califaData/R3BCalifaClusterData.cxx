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

#include "R3BCalifaClusterData.h"
#include "TMath.h"
#include <fmt/core.h>

R3BCalifaClusterData::R3BCalifaClusterData(std::vector<uint16_t> crystalList,
                                           double ene,
                                           double nf,
                                           double ns,
                                           double theta,
                                           double phi,
                                           ULong64_t time,
                                           uint8_t clusterType)
    : fCrystalList(crystalList)
    , fEnergy(ene)
    , fNf(nf)
    , fNs(ns)
    , fTheta(theta)
    , fPhi(phi)
    , fTime(time)
    , fClusterType(clusterType)
{
}

R3BCalifaClusterData::R3BCalifaClusterData(const R3BCalifaClusterData& right)
    : fCrystalList(right.fCrystalList)
    , fEnergy(right.fEnergy)
    , fNf(right.fNf)
    , fNs(right.fNs)
    , fTheta(right.fTheta)
    , fPhi(right.fPhi)
    , fTime(right.fTime)
    , fClusterType(right.fClusterType)
{
}

std::string R3BCalifaClusterData::toString() const
{
    return fmt::format(
        "Cluster type: {}, Energy(MeV): {:.3f}, Theta(deg): {:.3f}, Phi(deg): {:.3f}, Nf: {:.3f}, Ns: {:.3f}, Time: {}",
        GetClusterType(),
        GetEnergy(),
        GetTheta() * TMath::RadToDeg(),
        GetPhi() * TMath::RadToDeg(),
        GetNf(),
        GetNs(),
        GetTime());
}

void R3BCalifaClusterData::Print(const Option_t*) const { std::cout << *this << std::endl; }

std::ostream& operator<<(std::ostream& os, const R3BCalifaClusterData& data)
{
    os << data.toString();
    return os;
}

ClassImp(R3BCalifaClusterData)
