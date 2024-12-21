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
// -----              R3BFibMappedData                -----
// -----             Created Jan 2018 by M.Heil        -----
// ----------------------------------------------------------------

#include "R3BFibMappedData.h"

R3BFibMappedData::R3BFibMappedData()
    : fPlane(0)
    , fFiber(0)
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

R3BFibMappedData::R3BFibMappedData(Int_t planeId, Int_t fiberId)
    : fPlane(planeId)
    , fFiber(fiberId)
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

ClassImp(R3BFibMappedData)
