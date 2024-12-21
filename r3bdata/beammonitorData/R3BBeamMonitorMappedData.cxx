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

#include "R3BBeamMonitorMappedData.h"
#include "FairLogger.h"
#include "TMath.h"
#define IS_NAN(x) TMath::IsNaN(x)
#include <iostream>

using namespace std;

R3BBeamMonitorMappedData::R3BBeamMonitorMappedData()
    : fIC(-1)      // IC
    , fSEETRAM(-1) // SEETRAM
    , fTOFDOR(-1)  // TOFD or
{
}

R3BBeamMonitorMappedData::R3BBeamMonitorMappedData(UInt_t IC, UInt_t SEETRAM, UInt_t TOFDOR)
    : fIC(IC)
    , fSEETRAM(SEETRAM)
    , fTOFDOR(TOFDOR)
{
    /*
     LOG(info)<<"LosBeamMonitorData: "<<IC<<",  "<<SEETRAM<<", "<<TOFDOR;

    */
}

UInt_t R3BBeamMonitorMappedData::GetIC() const { return fIC; }

UInt_t R3BBeamMonitorMappedData::GetSEETRAM() const { return fSEETRAM; }

UInt_t R3BBeamMonitorMappedData::GetTOFDOR() const { return fTOFDOR; }

ClassImp(R3BBeamMonitorMappedData)
