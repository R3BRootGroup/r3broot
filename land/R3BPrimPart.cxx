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

// -----------------------------------------------------------------------------
// -----                        R3BPrimPart                                -----
// -----                 Created 12-07-2011 by D.Kresan                    -----
// -----------------------------------------------------------------------------

#include <math.h>

#include "R3BPrimPart.h"

ClassImp(R3BPrimPart)

    // -----------------------------------------------------------------------------
    R3BPrimPart::R3BPrimPart()
{
    fPdg = 0;
    fPx = 0.;
    fPy = 0.;
    fPz = 0.;
    fX = 0.;
    fY = 0.;
    fZ = 0.;
    fT = 0.;
    fA = 0.;
    fM = 0.;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
R3BPrimPart::R3BPrimPart(Int_t pdg,
                         Double_t px,
                         Double_t py,
                         Double_t pz,
                         Double_t x,
                         Double_t y,
                         Double_t z,
                         Double_t t,
                         Double_t A,
                         Double_t M)
{
    fPdg = pdg;
    fPx = px;
    fPy = py;
    fPz = pz;
    fX = x;
    fY = y;
    fZ = z;
    fT = t;
    fA = A;
    fM = M;
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
R3BPrimPart::~R3BPrimPart() {}
// -----------------------------------------------------------------------------
