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

#include "R3BTofdHitData.h"

// for the data analysis of the Tofd detector.
// Introduced by M.Heil, May 2016

R3BTofdHitData::R3BTofdHitData()
    : fBarId(0)
    , fTof(0.)
{
}

R3BTofdHitData::R3BTofdHitData(Double_t t,
                               Double_t x,
                               Double_t y,
                               Double_t Z,
                               Double_t tdiff,
                               Double_t ELoss,
                               Double_t ID,
                               UInt_t iBar,
                               Double_t traw,
                               Double_t tof)
    : R3BHit(ID, x, y, ELoss, t)
    , fBarId(iBar)
    , fTimeRaw(traw)
    , fTof(tof)
{
}

UInt_t R3BTofdHitData::GetBarId() const { return fBarId; }

Double_t R3BTofdHitData::GetTimeRaw() const { return fTimeRaw; } // returm bar time without trigger, needed for Rolu Tof

Double_t R3BTofdHitData::GetTof() const { return fTof; }

ClassImp(R3BTofdHitData)
