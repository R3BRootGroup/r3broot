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

#include "R3BRpcCal2Hit.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TRandom.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "TGeoManager.h"
#include "TGeoMatrix.h"

#include "R3BRpcCalData.h"
#include <list>
#include <vector>

R3BRpcCal2Hit::R3BRpcCal2Hit()
    : FairTask("R3B RPC Cal to Hit")
    , fParCont1(NULL)
    , fParCont2(NULL)
    , fParCont3(NULL)
    , fParCont4(NULL)
    , fRpcHitDataCA(NULL)
    , fRpcCalDataCA(NULL)
    , fR3BEventHeader(NULL)
    , fOnline(kFALSE)
{
}

R3BRpcCal2Hit::~R3BRpcCal2Hit()
{
    LOG(info) << "R3BRpcCal2Hit: Delete instance";
    if (fRpcHitDataCA)
        delete fRpcHitDataCA;
}

void R3BRpcCal2Hit::SetParContainers() {}

InitStatus R3BRpcCal2Hit::Init()
{

    // Parameter Container
    // Reading RPCHitPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BRpcCal2Hit:: FairRuntimeDb not opened";
    }

    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(fatal) << "R3BRpcCal2Hit::FairRootManager not found";
        return kFATAL;
    }

    fR3BEventHeader = dynamic_cast<R3BEventHeader*>(rootManager->GetObject("EventHeader."));
    if (!fR3BEventHeader)
    {
        LOG(error) << "R3BRpcCal2HitPar::Init() EventHeader. not found";
        return kFATAL;
    }

    fHitPar = dynamic_cast<R3BRpcHitPar*>(rtdb->getContainer("RpcHitPar"));
    if (!fHitPar)
    {
        LOG(error) << "R3BRpcCal2Hit::Init() Couldn't get handle on RPCHitPar container";
    }
    else
    {
        LOG(info) << "R3BRpcCal2Hit:: RPCHitPar container open";
    }

    fRpcCalDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("R3BRpcCalData"));
    if (!fRpcCalDataCA)
    {
        LOG(error) << "R3BRpcCal2HitPar::Init() R3BRpcCalData not found";
        return kFATAL;
    }

    // Register output array
    fRpcHitDataCA = new TClonesArray("R3BRpcHitData");
    rootManager->Register("R3BRpcHitData", "RPC Strip Hit", fRpcHitDataCA, !fOnline);

    // fill the TArray with Tot parameters!!!
    fParCont1 = fHitPar->GetCalParams1();
    fParCont2 = fHitPar->GetCalParams2();
    fParCont3 = fHitPar->GetCalParams3();
    fParCont4 = fHitPar->GetCalParams4();

    // Definition of a time stich object to correlate times coming from different systems
    fTimeStitch = new R3BCoarseTimeStitch();

    return kSUCCESS;
}

InitStatus R3BRpcCal2Hit::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BRpcCal2Hit::Exec(Option_t* opt)
{
    Reset();
    // loop over strip data
    Int_t nHits = fRpcCalDataCA->GetEntries();
    UInt_t iDetector = 0;
    double charge_left = -1000;
    double charge_right = -1000;
    double time_left = 0;
    double time_right = 0;
    double ichn_right = 0;
    double ichn_left = 0;
    UInt_t inum;

    for (Int_t i = 0; i < nHits; i++)
    {
        auto map1 = dynamic_cast<R3BRpcCalData*>(fRpcCalDataCA->At(i));
        iDetector = map1->GetDetId();
        inum = iDetector * 41 + map1->GetChannelId() - 1;

        if (iDetector == 0)
        {
            if (map1->GetTotR_B() >= charge_right)
            {
                charge_right = map1->GetTotR_B();
                time_right = map1->GetTimeR_B();
                ichn_right = map1->GetChannelId();
            }
            if (map1->GetTotL_T() >= charge_left)
            {

                charge_left = map1->GetTotL_T();
                time_left = map1->GetTimeL_T();
                ichn_left = map1->GetChannelId();
            }
            if (ichn_left == ichn_right)
            {

                double position = ((time_left - time_right) * CSTRIP / 2. - 10 * (fParCont1->GetAt(inum) - 200));

                double charge = (charge_left + charge_right) / 2.;

                double time = (time_left + time_right) / 2. - fParCont3->GetAt(inum);

                auto tof = fTimeStitch->GetTime(time - fR3BEventHeader->GetTStart());
                AddHitStrip(iDetector, ichn_right, time, position, charge, tof);
            }
        }
        // loop over Pmt data
        if (iDetector == 1)
        {
            double position = (map1->GetTimeR_B() - map1->GetTimeL_T()) * CSCINT / 2. - (fParCont1->GetAt(inum) - 2500);

            double charge = (map1->GetTotR_B() + map1->GetTotL_T()) / 2.;

            double time = (map1->GetTimeR_B() + map1->GetTimeL_T()) / 2.;

            auto tof = fTimeStitch->GetTime(time - fR3BEventHeader->GetTStart());

            AddHitStrip(iDetector, map1->GetChannelId(), time, position, charge, tof);
        }
    }
}

R3BRpcHitData* R3BRpcCal2Hit::AddHitStrip(UInt_t detId,
                                          UInt_t channel,
                                          double time,
                                          double pos,
                                          double charge,
                                          double tof)
{

    TClonesArray& clref = *fRpcHitDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BRpcHitData(detId, channel, time, pos, charge, tof);
}

void R3BRpcCal2Hit::Finish() {}

void R3BRpcCal2Hit::Reset()
{
    LOG(debug) << "Clearing RPCHItStructure Structure";
    if (fRpcHitDataCA)
    {
        fRpcHitDataCA->Clear();
    }
}

ClassImp(R3BRpcCal2Hit);
