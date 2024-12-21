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

// -------------------------------------------------------------------------
// -----            R3BAmsHitData source file                          -----
// -------------------------------------------------------------------------

#include "R3BAmsHitData.h"

// -----   Default constructor   -------------------------------------------
R3BAmsHitData::R3BAmsHitData()
    : FairMultiLinkedData()
    , fDetId(-1)
    , fNumHit(0)
    , fMulS(0)
    , fMulK(0)
    , fX(0)
    , fY(0)
    , fTheta(0)
    , fPhi(0)
    , fmaster(0, 0, 0)
    , fEnergyS(0)
    , fEnergyK(0)
{
}

// -----   Standard constructor   ------------------------------------------
R3BAmsHitData::R3BAmsHitData(Int_t detid,
                             Int_t numhit,
                             Double_t x,
                             Double_t y,
                             TVector3 master,
                             Double_t energy_s,
                             Double_t energy_k,
                             Int_t mulS,
                             Int_t mulK)
    : FairMultiLinkedData()
    , fDetId(detid)
    , fNumHit(numhit)
    , fMulS(mulS)
    , fMulK(mulK)
    , fX(x)
    , fY(y)
    , fTheta(master.Theta())
    , fPhi(master.Phi())
    , fmaster(master)
    , fEnergyS(energy_s)
    , fEnergyK(energy_k)
{
}
// -------------------------------------------------------------------------

R3BAmsHitData::R3BAmsHitData(const R3BAmsHitData& right)
    : FairMultiLinkedData(right)
    , fDetId(right.fDetId)
    , fNumHit(right.fNumHit)
    , fMulS(right.fMulS)
    , fMulK(right.fMulK)
    , fX(right.fX)
    , fY(right.fY)
    , fTheta(right.fmaster.Theta())
    , fPhi(right.fmaster.Phi())
    , fmaster(right.fmaster)
    , fEnergyS(right.fEnergyS)
    , fEnergyK(right.fEnergyK)
{
}

ClassImp(R3BAmsHitData);
