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

#include "R3BMCTrack.h"
#include "FairLogger.h"
#include <iostream>
#include <utility>

R3BMCTrack::R3BMCTrack()
    : fPdgCode(0)
    , fMotherId(0)
    , fStartVertex()
    , fMomentumMass()
    , fNPoints()
{
}

R3BMCTrack::R3BMCTrack(int pdgCode,
                       int motherId,
                       ROOT::Math::XYZTVector xyzt,
                       ROOT::Math::PxPyPzMVector pm,
                       std::array<int, kLAST + 1> nPoints)
    : fPdgCode(pdgCode)
    , fMotherId(motherId)
    , fStartVertex(std::move(xyzt))
    , fMomentumMass(std::move(pm))
    , fNPoints(nPoints)
{
}

R3BMCTrack::R3BMCTrack(TParticle* part, std::array<int, kLAST + 1> nPoints, int fMC)
    : fPdgCode(part->GetPdgCode())
    , fMotherId(part->GetMother(0))
    , fStartVertex(part->Vx(), part->Vy(), part->Vz(), fMC == 0 ? part->T() * 1e09 : part->T()) // G3 == 0 ?!
    , fMomentumMass(part->Px(), part->Py(), part->Pz(), part->GetMass())
    , fNPoints(nPoints)
{
}

void R3BMCTrack::Print(Option_t* option) const
{
    std::cout << "Track " << option << ", mother : " << fMotherId << ", Type " << fPdgCode << ", momentum ("
              << fMomentumMass.Px() << ", " << fMomentumMass.Py() << ", " << fMomentumMass.Pz() << ") GeV" << std::endl;
    std::cout << "Ref " << GetNPoints(kREF) << ", DCH " << GetNPoints(kDCH) << ", CAL " << GetNPoints(kCAL) << ", LAND "
              << GetNPoints(kLAND) << ", GFI " << GetNPoints(kGFI) << ", TOFd " << GetNPoints(kTOFD) << ", TOF "
              << GetNPoints(kTOF) << ", TRACKER " << GetNPoints(kTRA) << ", CALIFA " << GetNPoints(kCALIFA) << ", PSP "
              << GetNPoints(kPSP) << ", VETO " << GetNPoints(kVETO) << ", RPC " << GetNPoints(kRPC) << ", NeuLAND "
              << GetNPoints(kNEULAND) << std::endl;
#ifdef SOFIA
    std::cout << ", SCI " << GetNPoints(kSOFSCI) << ", AT " << GetNPoints(kSOFAT) << ", TRIM " << GetNPoints(kSOFTRIM)
              << ", MWPC1 " << GetNPoints(kSOFMWPC1) << ", TWIM " << GetNPoints(kSOFTWIM) << ", MWPC2 "
              << GetNPoints(kSOFMWPC2) << ", SOF ToF Wall " << GetNPoints(kSOFTofWall) << std::endl;
#endif
#ifdef GTPC
    std::cout << ", GTPC " << GetNPoints(kGTPC) << std::endl;
#endif
#ifdef ASYEOS
    std::cout << ", CHIMERA " << GetNPoints(kCHIMERA) << ", KRAB " << GetNPoints(kKRAB) << ", KRATTA "
              << GetNPoints(kKRATTA) << std::endl;
#endif
}

int R3BMCTrack::GetNPoints(DetectorId detId) const
{
    if (detId < 0 || detId >= fNPoints.size())
    {
        LOG(error) << "Unknown detector ID " << detId;
        return 0;
    }
    return fNPoints.at(detId);
}

ClassImp(R3BMCTrack)
