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

#include "R3BNeulandHit.h"

constexpr Double_t c = 29.97924580000000105;   // cm/ns
constexpr Double_t c2 = 898.75517873681758374; // cm²/ns²

R3BNeulandHit::R3BNeulandHit(Int_t paddle,
                             Double_t TdcL,
                             Double_t TdcR,
                             Double_t time,
                             Double_t QdcL,
                             Double_t QdcR,
                             Double_t energy,
                             const TVector3& pos,
                             const TVector3& pix)
    : fPaddle(paddle)
    , fTdcL(TdcL)
    , fTdcR(TdcR)
    , fT(time)
    , fQdcL(QdcL)
    , fQdcR(QdcR)
    , fE(energy)
    , fPosition(pos)
    , fPixel(pix)
{
}

Double_t R3BNeulandHit::GetBeta() const { return GetPosition().Mag() / (GetT() * c); }

Double_t R3BNeulandHit::GetEToF(const Double_t mass) const
{
    const Double_t v2 = GetPosition().Mag2() / std::pow(GetT(), 2); // cm²/ns²
    const Double_t gamma = 1. / std::sqrt(1. - (v2 / c2));
    const auto etof = (gamma - 1.) * mass;
    return 1.81522 + 0.984612 * etof; // TODO: EToF is ever so slightly off. Maybe some rounding/mass error?
}

std::ostream& operator<<(std::ostream& os, const R3BNeulandHit& hit)
{
    os << "R3BNeulandHit: NeuLAND Hit in Paddle " << hit.GetPaddle() << std::endl
       << "    TdcL: " << hit.GetTdcL() << "    TdcR: " << hit.GetTdcR() << "    Time: " << hit.GetT() << std::endl
       << "    QdcL: " << hit.GetQdcL() << "    QdcR: " << hit.GetQdcR() << "    Energy: " << hit.GetE() << std::endl
       << "    Position XYZ: " << hit.GetPosition().X() << "    " << hit.GetPosition().Y() << "    "
       << hit.GetPosition().Z() << std::endl;
    return os;
}

void R3BNeulandHit::Print(const Option_t*) const { std::cout << *this; }

ClassImp(R3BNeulandHit)
