/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
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

R3BCalifaClusterData::R3BCalifaClusterData()
    : FairMultiLinkedData()
    , fEnergy(NAN)
    , fTheta(NAN)
    , fPhi(NAN)
    , fTime(0)
{
}

R3BCalifaClusterData::R3BCalifaClusterData(std::vector<Int_t> crystalList,
                                           std::vector<Double_t> energyList,
                                           std::vector<ULong64_t> timeList,
                                           Double_t ene,
                                           Double_t nf,
                                           Double_t ns,
                                           Double_t theta,
                                           Double_t phi,
                                           ULong64_t time,
                                           Int_t clusterType)
    : FairMultiLinkedData()
    , fCrystalList(crystalList)
    , fEnergyList(energyList)
    , fTimeList(timeList)
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
    : FairMultiLinkedData(right)
    , fCrystalList(right.fCrystalList)
    , fEnergyList(right.fEnergyList)
    , fTimeList(right.fTimeList)
    , fEnergy(right.fEnergy)
    , fNf(right.fNf)
    , fNs(right.fNs)
    , fTheta(right.fTheta)
    , fPhi(right.fPhi)
    , fTime(right.fTime)
    , fClusterType(right.fClusterType)
{
}

ClassImp(R3BCalifaClusterData);
