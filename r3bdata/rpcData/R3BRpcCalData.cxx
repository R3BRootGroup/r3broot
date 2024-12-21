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

#include "R3BRpcCalData.h"

R3BRpcCalData::R3BRpcCalData()
    : FairMultiLinkedData()
    , fDetId(0)
    , fChannelId(0)
    , fTimeR_B(0)
    , fTimeL_T(0)
    , fTotR_B(0)
    , fTotL_T(0)
{
}

R3BRpcCalData::R3BRpcCalData(UShort_t detId,
                             UShort_t channelId,
                             double TimeR_B,
                             double TimeL_T,
                             double TotR_B,
                             double TotL_T)
    : FairMultiLinkedData()
    , fDetId(detId)
    , fChannelId(channelId)
    , fTimeR_B(TimeR_B)
    , fTimeL_T(TimeL_T)
    , fTotR_B(TotR_B)
    , fTotL_T(TotL_T)
{
}

ClassImp(R3BRpcCalData);
