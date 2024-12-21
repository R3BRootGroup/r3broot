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

#include "R3BTofdMappedData.h"

R3BTofdMappedData::R3BTofdMappedData()
    : fDetector(-1)
    , fSide(-1)
    , fBar(-1)
    , fTimeFine(-1)
    , fTimeCoarse(-1)
{
}

R3BTofdMappedData::R3BTofdMappedData(UInt_t detector,
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
}

UInt_t R3BTofdMappedData::GetDetectorId() const { return fDetector; }

UInt_t R3BTofdMappedData::GetSideId() const { return fSide; }

UInt_t R3BTofdMappedData::GetBarId() const { return fBar; }

UInt_t R3BTofdMappedData::GetEdgeId() const { return fEdge; }

UInt_t R3BTofdMappedData::GetTimeCoarse() const { return fTimeCoarse; }

UInt_t R3BTofdMappedData::GetTimeFine() const { return fTimeFine; }

ClassImp(R3BTofdMappedData)
