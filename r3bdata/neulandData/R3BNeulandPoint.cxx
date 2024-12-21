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

#include "R3BNeulandPoint.h"

using namespace std;

TVector3 R3BNeulandPoint::GetMomentum() const
{
    TVector3 v;
    // Momentum(v); // TODO: #fairroot: FairMCPoint not const
    v.SetXYZ(fPx, fPy, fPz);
    return v;
}

TVector3 R3BNeulandPoint::GetPosition() const
{
    TVector3 v;
    Position(v);
    return v;
}

ostream& operator<<(ostream& os, const R3BNeulandPoint& point)
{
    os << "R3BNeulandPoint: NeuLAND Point in Paddle " << point.GetPaddle() << " with light yield "
       << point.GetLightYield() << std::endl;
    // TODO: #fairroot: Move this to FairMCPoint
    os << "... with FairMCPoint for track " << point.GetTrackID() << std::endl;
    os << "    Position XYZ: " << point.GetPosition().X() << "    " << point.GetPosition().X() << "    "
       << point.GetPosition().Z() << std::endl;
    os << "    Momentum XYZ: " << point.GetMomentum().X() << "    " << point.GetMomentum().X() << "    "
       << point.GetMomentum().Z() << std::endl;
    os << "    Time: " << point.GetTime() << "    Length: " << point.GetLength()
       << "    Energy loss: " << point.GetEnergyLoss() << std::endl;
    return os;
}

void R3BNeulandPoint::Print(const Option_t*) const { std::cout << *this; }

ClassImp(R3BNeulandPoint)
