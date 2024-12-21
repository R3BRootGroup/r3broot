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

// ---------------------------------------------------------------
// -----                  R3BFrsData                         -----
// -----     Created 30/04/2021 by J.L. Rodriguez            -----
// ---------------------------------------------------------------

#include "R3BFrsData.h"

R3BFrsData::R3BFrsData()
    : fStaId(0)
    , fStoId(0)
    , fZ(0.)
    , fAq(0.)
    , fBeta(0.)
    , fTof(0.)
    , fBrho(0.)
    , fxs2(0.)
    , fxc(0.)
{
}

R3BFrsData::R3BFrsData(Int_t StaId,
                       Int_t StoId,
                       Double_t z,
                       Double_t aq,
                       Double_t beta,
                       Double_t brho,
                       Double_t xs2,
                       Double_t xc,
                       Double_t tof)
    : fStaId(StaId)
    , fStoId(StoId)
    , fZ(z)
    , fAq(aq)
    , fBeta(beta)
    , fTof(tof)
    , fBrho(brho)
    , fxs2(xs2)
    , fxc(xc)
{
}

ClassImp(R3BFrsData)
