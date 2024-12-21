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

#include "R3BNeulandCalData.h"

R3BNeulandCalData::R3BNeulandCalData(Int_t barId, Int_t side, Double_t time, Double_t triggertime, Int_t qdc)
    : fBarId(barId)
    , fSide(side)
    , fTime(time)
    , fTriggerTime(triggertime)
    , fQdc(qdc)
{
}

std::ostream& operator<<(std::ostream& os, const R3BNeulandCalData& calData)
{
    os << "R3BNeulandCalData: BarID " << calData.GetBarId() << "    Side " << calData.GetSide() << "    Time "
       << calData.GetTime() << "   Trigger time " << calData.GetTriggerTime() << "    QDC " << calData.GetQdc()
       << std::endl;
    return os;
}

void R3BNeulandCalData::Print(const Option_t*) const { std::cout << *this; }

ClassImp(R3BNeulandCalData)
