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

#include "R3BPdcHitData.h"

// for the data analysis of the Pdc detector.
// Introduced by M.Heil, June 2020

R3BPdcHitData::R3BPdcHitData() {}

R3BPdcHitData::R3BPdcHitData(Double_t t, Double_t x, Double_t y, Double_t ELoss, Double_t ID, Int_t wire)
    : R3BHit(ID, x, y, ELoss, t)
    , fWireId(wire)

{
}
Int_t R3BPdcHitData::GetWireId() const { return fWireId; }

ClassImp(R3BPdcHitData)
