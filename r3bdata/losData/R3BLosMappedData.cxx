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

#include "R3BLosMappedData.h"
#include "FairLogger.h"
#include "TMath.h"
#define IS_NAN(x) TMath::IsNaN(x)
#include <iostream>

using namespace std;

// for the data analysis of the Los detector.
// Introduced by Ralf, Jan 2016
// fTimeLCoarse, fTimeLFine, fTimeTCoarse, fTimeTFine added by Aleksandra, Sept 2016
// Structure changed to detector, channel, type, fineTime, coarseTime Nov 2016 (Ralf & Aleksandra)

R3BLosMappedData::R3BLosMappedData()
    : fDetector(-1)   // Detector number: 1...n
    , fChannel(-1)    // Channel number:  1...n (n = 4 (old LOS), 8 (new LOS))
    , fType(-1)       // Type: 0 = VFTX time, 1 = TAMEX leading edge, 2 = TAMEX trailing edge
    , fTimeFine(-1)   // Fine time
    , fTimeCoarse(-1) // Coarse time
{
}

R3BLosMappedData::R3BLosMappedData(UInt_t detector, UInt_t channel, UInt_t type, UInt_t timeFine, UInt_t timeCoarse)
    : fDetector(detector)
    , fChannel(channel)
    , fType(type)
    , fTimeFine(timeFine)
    , fTimeCoarse(timeCoarse)
{
    // cout<<"R3BLosMappedData: det. "<< detector<<", chann. "<< channel <<", type "<<type<<", fine/coarse Times:
    // "<<timeFine<<", "
    //<<timeCoarse<<endl;
    /* Everything fine here!

     if(timeFine <= 0 || IS_NAN(timeFine) ) LOG(info)<<"LosMappedData FineTime wrong value for channel: "<<channel<<"
     type: "<<type<<" value: "<<timeFine; if(timeCoarse < 0 || IS_NAN(timeCoarse)) LOG(info)<<"LosMappedData CoarseTime
     wrong value for channel: "<<channel<<" type: "<<type<<" value: "<<timeCoarse;

    */
}

UInt_t R3BLosMappedData::GetDetector() const { return fDetector; }

UInt_t R3BLosMappedData::GetChannel() const { return fChannel; }

UInt_t R3BLosMappedData::GetType() const { return fType; }
UInt_t R3BLosMappedData::GetTimeCoarse() const { return fTimeCoarse; }

UInt_t R3BLosMappedData::GetTimeFine() const { return fTimeFine; }

ClassImp(R3BLosMappedData)
