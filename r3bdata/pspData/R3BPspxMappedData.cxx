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

#include "R3BPspxMappedData.h"

R3BPspxMappedData::R3BPspxMappedData()
    : fStrip{ -1, -1 }
    , fEnergy{ -1, -1 }
{
}

R3BPspxMappedData::R3BPspxMappedData(Int_t strip1, Int_t strip2, Int_t energy1, Int_t energy2)
    : fStrip{ strip1, strip2 }
    , fEnergy{ energy1, energy2 }
{
}

void R3BPspxMappedData::SetValue(Int_t side, Int_t strip, Int_t energy)
{
    fStrip[side] = strip;
    fEnergy[side] = energy;
}

ClassImp(R3BPspxMappedData)
