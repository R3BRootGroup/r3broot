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

#include "R3BPdcMappedData.h"

R3BPdcMappedData::R3BPdcMappedData()
    : fPlane(-1)
    , fWire(-1)
    , fEdge(-1)
    , fTimeFine(-1)
    , fTimeCoarse(-1)
{
}

R3BPdcMappedData::R3BPdcMappedData(UInt_t plane, UInt_t wire, UInt_t edge, UInt_t timeCoarse, UInt_t timeFine)
    : fPlane(plane)
    , fWire(wire)
    , fEdge(edge)
    , fTimeCoarse(timeCoarse)
    , fTimeFine(timeFine)
{
}

UInt_t R3BPdcMappedData::GetPlaneId() const { return fPlane; }

UInt_t R3BPdcMappedData::GetWireId() const { return fWire; }

UInt_t R3BPdcMappedData::GetEdgeId() const { return fEdge; }

UInt_t R3BPdcMappedData::GetTimeCoarse() const { return fTimeCoarse; }

UInt_t R3BPdcMappedData::GetTimeFine() const { return fTimeFine; }

ClassImp(R3BPdcMappedData)
