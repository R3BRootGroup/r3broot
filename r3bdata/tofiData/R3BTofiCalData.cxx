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

#include "R3BTofiCalData.h"
#include "FairLogger.h"
#include "TMath.h"
#include <cmath>
#include <iostream>

R3BTofiCalData::R3BTofiCalData()
    : fDetector(-1)
    , fBar(-1)
    , fSide(-1)
    , fLeading_ns(-1)
    , fTrailing_ns(-1)
{
}

R3BTofiCalData::R3BTofiCalData(UInt_t detector, UInt_t channel, UInt_t side, Double_t lead, Double_t trail)
    : fDetector(detector)
    , fBar(channel)
    , fSide(side)
    , fLeading_ns(lead)
    , fTrailing_ns(trail)
{
    //	if(fDetector<2) std::cout<<"Cal leading        : "<<fSide<<"; "<<fBar<<", "<<fLeading_ns<<",
    //"<<fTrailing_ns<<std::endl;
}

UInt_t R3BTofiCalData::GetDetectorId() const { return fDetector; }

UInt_t R3BTofiCalData::GetBarId() const { return fBar; }

UInt_t R3BTofiCalData::GetSideId() const { return fSide; }

Double_t R3BTofiCalData::GetTimeLeading_ns() const { return fLeading_ns; }

Double_t R3BTofiCalData::GetTimeTrailing_ns() const { return fTrailing_ns; }

ClassImp(R3BTofiCalData)
