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

#include "R3BPspxPrecalData.h"

R3BPspxPrecalData::R3BPspxPrecalData()
    : fStrip(0)
    , fEnergy1(0)
    , fEnergy2(0)
{
}

R3BPspxPrecalData::R3BPspxPrecalData(Int_t strip, Float_t energy1, Float_t energy2)
    : fStrip(strip)
    , fEnergy1(energy1)
    , fEnergy2(energy2)
{
}

ClassImp(R3BPspxPrecalData)
