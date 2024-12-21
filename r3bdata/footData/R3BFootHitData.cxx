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

// ----------------------------------------------------------------------
// -----           R3BFootHitData source file                       -----
// ----------------------------------------------------------------------

#include "R3BFootHitData.h"

// -----   Standard constructor   ---------------------------------------
R3BFootHitData::R3BFootHitData(uint8_t detid,
                               uint16_t nbhit,
                               double pos,
                               TVector3 master,
                               double energy,
                               uint16_t mulstrip)
    : fDetId(detid)
    , fNbHit(nbhit)
    , fPos(pos)
    , fTheta(master.Theta())
    , fPhi(master.Phi())
    , fmaster(master)
    , fEnergy(energy)
    , fMulStrip(mulstrip)
{
}

// ----------------------------------------------------------------------
R3BFootHitData::R3BFootHitData(const R3BFootHitData& right)
    : fDetId(right.fDetId)
    , fNbHit(right.fNbHit)
    , fPos(right.fPos)
    , fTheta(right.fmaster.Theta())
    , fPhi(right.fmaster.Phi())
    , fmaster(right.fmaster)
    , fEnergy(right.fEnergy)
    , fMulStrip(right.fMulStrip)
{
}

ClassImp(R3BFootHitData)
