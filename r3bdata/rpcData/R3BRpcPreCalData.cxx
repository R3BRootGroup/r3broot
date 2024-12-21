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

#include "R3BRpcPreCalData.h"

R3BRpcPreCalData::R3BRpcPreCalData()
    : FairMultiLinkedData()
    , fDetId(0)
    , fChannelId(0)
    , fTime(0)
    , fTot(0)
    , fSide(0)
{
}

R3BRpcPreCalData::R3BRpcPreCalData(UShort_t DetId, UShort_t channelId, double Time, double Tot, UShort_t Side)
    : FairMultiLinkedData()
    , fDetId(DetId)
    , fChannelId(channelId)
    , fTime(Time)
    , fTot(Tot)
    , fSide(Side)
{
}

ClassImp(R3BRpcPreCalData);
