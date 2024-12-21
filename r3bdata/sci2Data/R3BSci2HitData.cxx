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

#include "R3BSci2HitData.h"

// for the data analysis of the Sci2 detector.

R3BSci2HitData::R3BSci2HitData()
    : fSciId(0)
    , fX(0.)
    , fTime(0.)
    , fTimeWithTref(0.)
{
}

R3BSci2HitData::R3BSci2HitData(Int_t sci, Double_t x, Double_t time, Double_t time_w_tref)
    : fSciId(sci)
    , fX(x)
    , fTime(time)
    , fTimeWithTref(time_w_tref)

{
}

ClassImp(R3BSci2HitData);
