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

#include "R3BNeulandTacquilaMappedData.h"

R3BNeulandTacquilaMappedData::R3BNeulandTacquilaMappedData(UShort_t clock,
                                                           UShort_t tacData,
                                                           UShort_t stopT,
                                                           UShort_t qdcData,
                                                           Int_t plane,
                                                           Int_t paddle,
                                                           Int_t side)
    : fClock(clock)
    , fTacData(tacData)
    , fStopT(stopT)
    , fQdcData(qdcData)
    , fPlane(plane)
    , fPaddle(paddle)
    , fSide(side)
{
}

std::ostream& operator<<(std::ostream& os, const R3BNeulandTacquilaMappedData& mappedData)
{
    os << "R3BNeulandTacquilaMappedData:"
       << "    Clock " << mappedData.GetClock() << "    TacData " << mappedData.GetTacData() << "    StopT "
       << mappedData.GetStopT() << "    QdcData " << mappedData.GetQdcData() << "    Plane " << mappedData.GetPlane()
       << "    Paddle " << mappedData.GetPaddle() << "    Side " << mappedData.GetSide() << std::endl;
    return os;
}

void R3BNeulandTacquilaMappedData::Print(const Option_t*) const { std::cout << *this; }

ClassImp(R3BNeulandTacquilaMappedData)
