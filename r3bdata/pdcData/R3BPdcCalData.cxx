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

#include "R3BPdcCalData.h"

R3BPdcCalData::R3BPdcCalData()
    : fPlane(-1)
    , fWire(-1)
    , fEdge(-1)
    , fTime_ns(-1)
{
}

R3BPdcCalData::R3BPdcCalData(UInt_t plane, UInt_t wire, UInt_t edge, Double_t time_ns)
    : fPlane(plane)
    , fWire(wire)
    , fEdge(edge)
    , fTime_ns(time_ns)
{
}

UInt_t R3BPdcCalData::GetPlaneId() const { return fPlane; }

UInt_t R3BPdcCalData::GetWireId() const { return fWire; }

UInt_t R3BPdcCalData::GetEdgeId() const { return fEdge; }

Double_t R3BPdcCalData::GetTime_ns() const { return fTime_ns; }

ClassImp(R3BPdcCalData)
