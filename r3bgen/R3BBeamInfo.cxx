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

#include "R3BBeamInfo.h"

R3BBeamInfo::R3BBeamInfo()
    : theta(0.)
    , phi(0.)
    , xPos(0.)
    , yPos(0.)
    , zPos(0.)
{
}

R3BBeamInfo::~R3BBeamInfo()
{
    //
    // Destructor. Makes nothing
    //
}

void R3BBeamInfo::print(void)
{
    //
    // Prints info
    //
}

void R3BBeamInfo::SetVertexPosition(Double_t x, Double_t y, Double_t z)
{
    //
    // Sets the position of the vertex (origin of the primary particle)
    //
    xPos = x;
    yPos = y;
    zPos = z;
}

void R3BBeamInfo::SetAngles(Double_t angleTheta, Double_t anglePhi)
{
    //
    // Sets the 4momentum
    //
    theta = angleTheta;
    phi = anglePhi;
}

ClassImp(R3BBeamInfo)
