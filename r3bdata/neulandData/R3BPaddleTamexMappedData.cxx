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

// ----------------------------------------------------------------
// -----              R3BPaddleTamexMappedData                -----
// -----             Created August 2nd 2016 by R.Plag        -----
// ----------------------------------------------------------------

#include "R3BPaddleTamexMappedData.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

R3BPaddleTamexMappedData::R3BPaddleTamexMappedData()
    : fPlane(0)
    , fBar(0)
    , fSide(0)
    , fEdge(0)
    , fCoarseTime1LE(-1)
    , fFineTime1LE(-1)
    , fCoarseTime1TE(-1)
    , fFineTime1TE(-1)
    , fCoarseTime2LE(-1)
    , fFineTime2LE(-1)
    , fCoarseTime2TE(-1)
    , fFineTime2TE(-1)
{
}

R3BPaddleTamexMappedData::R3BPaddleTamexMappedData(Int_t planeId, Int_t barId)
    : fPlane(planeId)
    , fBar(barId)
    , fSide(0)
    , fEdge(0)
    , fCoarseTime1LE(-1)
    , fFineTime1LE(-1)
    , fCoarseTime1TE(-1)
    , fFineTime1TE(-1)
    , fCoarseTime2LE(-1)
    , fFineTime2LE(-1)
    , fCoarseTime2TE(-1)
    , fFineTime2TE(-1)
{
}

ClassImp(R3BPaddleTamexMappedData)
