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

#include "R3BRoluHitData.h"

R3BRoluHitData::R3BRoluHitData()
    : fDetector(0)
    , fChannel(0)
    , ftime(0.)
    , ftot(0.)
{
}

R3BRoluHitData::R3BRoluHitData(Int_t detector, Int_t channel, Double_t time, Double_t tot)

    : fDetector(detector)
    , fChannel(channel)
    , ftime(time)
    , ftot(tot)

{
}

ClassImp(R3BRoluHitData);
