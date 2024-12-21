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

#include "R3BSci8MappedData.h"
#include "FairLogger.h"
#include "TMath.h"
#define IS_NAN(x) TMath::IsNaN(x)
#include <iostream>

using namespace std;

// for the data analysis of the Sci8 detector.

R3BSci8MappedData::R3BSci8MappedData()
    : fDetector(0)    // Detector number: 1...n
    , fChannel(0)     // Channel number:  1,2
    , fType(0)        // Type: 0 = VFTX time, 1 = TAMEX leading edge, 2 = TAMEX trailing edge
    , fTimeFine(-1)   // Fine time
    , fTimeCoarse(-1) // Coarse time
{
}

R3BSci8MappedData::R3BSci8MappedData(UInt_t detector, UInt_t channel, UInt_t type, UInt_t timeFine, UInt_t timeCoarse)
    : fDetector(detector)
    , fChannel(channel)
    , fType(type)
    , fTimeFine(timeFine)
    , fTimeCoarse(timeCoarse)
{
    // cout<<"R3BSci8MappedData: chann. "<< channel <<", type "<<type<<", fine/coarse Times: "<<timeFine<<", "
    //<<timeCoarse<<endl;
}

ClassImp(R3BSci8MappedData)
