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

#include "R3BTofiMappedData.h"
#include <iostream>

R3BTofiMappedData::R3BTofiMappedData()
    : fDetector(-1)
    , fSide(-1)
    , fBar(-1)
    , fTimeFine(-1)
    , fTimeCoarse(-1)
{
}

R3BTofiMappedData::R3BTofiMappedData(UInt_t detector,
                                     UInt_t side,
                                     UInt_t channel,
                                     UInt_t edge,
                                     UInt_t timeCoarse,
                                     UInt_t timeFine)
    : fDetector(detector)
    , fSide(side)
    , fBar(channel)
    , fEdge(edge)
    , fTimeCoarse(timeCoarse)
    , fTimeFine(timeFine)
{
    // if(fDetector < 2) std::cout<<"Mapped data        : "<<fSide<<"; "<<fBar<<"; "<<fEdge<<",
    // "<<fTimeCoarse<<std::endl;
}

UInt_t R3BTofiMappedData::GetDetectorId() const { return fDetector; }

UInt_t R3BTofiMappedData::GetSideId() const { return fSide; }

UInt_t R3BTofiMappedData::GetBarId() const { return fBar; }

UInt_t R3BTofiMappedData::GetEdgeId() const { return fEdge; }

UInt_t R3BTofiMappedData::GetTimeCoarse() const { return fTimeCoarse; }

UInt_t R3BTofiMappedData::GetTimeFine() const { return fTimeFine; }

ClassImp(R3BTofiMappedData)
