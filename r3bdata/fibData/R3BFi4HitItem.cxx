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

#include "R3BFi4HitItem.h"

// for the data analysis of the s438b PSPX detectors.
// Introduced by Ralf, Jan 2016
// adapted by Joachim March 2016

R3BFi4HitItem::R3BFi4HitItem()
    : fDetector(0)
    , fFiber_n(0)
    , fEnergy_top(0)
    , fEnergy_bottom(0)
    , fTime_top(0)
    , fTime_bottom(0)
{
}

R3BFi4HitItem::R3BFi4HitItem(UInt_t detector,
                             UInt_t fiber,
                             Float_t xpos,
                             Float_t ypos,
                             Float_t zpos,
                             Float_t energy_top,
                             Float_t energy_bottom,
                             Float_t time_top,
                             Float_t time_bottom)
    : R3BHit(0, xpos, ypos, (energy_top + energy_bottom) / 2., 0.)
    , fDetector(detector)
    , fFiber_n(fiber)
    , fEnergy_top(energy_top)
    , fEnergy_bottom(energy_bottom)
    , fTime_top(time_top)
    , fTime_bottom(time_bottom)
{
}

ClassImp(R3BFi4HitItem)
