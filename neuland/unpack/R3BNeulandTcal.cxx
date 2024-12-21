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

// ------------------------------------------------------------
// -----                  R3BNeuLandTcal                  -----
// -----            Created 27-01-2015 by M.Heil          -----
// ------------------------------------------------------------

#include "R3BNeulandTcal.h"

#include "R3BEventHeader.h"
#include "R3BNeulandPmt.h"
#include "R3BPaddleTamexMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TClonesArray.h"

R3BNeulandTcal::R3BNeulandTcal()
    : FairTask("LandTcal", 1)
    , fNEvents(0)
    , fMappedHit(NULL)
    , fPmt(new TClonesArray("R3BNeulandPmt"))
    , fNPmt(0)
    , fTcalPar(NULL)
    , fTrigger(-1)
    //    , fMap17Seen()
    //    , fMapStopTime()
    //    , fMapStopClock()
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{
}

R3BNeulandTcal::R3BNeulandTcal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNEvents(0)
    , fMappedHit(NULL)
    , fPmt(new TClonesArray("R3BNeulandPmt"))
    , fNPmt(0)
    , fTcalPar(NULL)
    , fTrigger(-1)
    //    , fMap17Seen()
    //    , fMapStopTime()
    //    , fMapStopClock()
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{
}

R3BNeulandTcal::~R3BNeulandTcal()
{
    if (fPmt)
    {
        delete fPmt;
        fPmt = NULL;
        fNPmt = 0;
    }
}

InitStatus R3BNeulandTcal::Init()
{
    LOG(info) << "R3BNeulandTcal::Init : read " << fTcalPar->GetNumModulePar() << " calibrated modules";
    // fTcalPar->printParams();

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(fatal) << "FairRootManager not found";
    }
    /*
        header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
        if (NULL == header)
        {
            FairLogger::GetLogger()->Fatal(MESSAGE_ORIGIN, "Branch R3BEventHeader not found");
        }
    */
    fMappedHit = dynamic_cast<TClonesArray*>(mgr->GetObject("NeulandMappedData"));
    if (NULL == fMappedHit)
    {
        LOG(fatal) << "Branch NeulandMappedData not found";
    }

    mgr->Register("NeulandPmt", "Land", fPmt, kTRUE);

    return kSUCCESS;
}

void R3BNeulandTcal::SetParContainers()
{
    FairRunAna* ana = FairRunAna::Instance();
    FairRuntimeDb* rtdb = ana->GetRuntimeDb();
    fTcalPar = dynamic_cast<R3BTCalPar*>(rtdb->getContainer("LandTCalPar"));
}

InitStatus R3BNeulandTcal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BNeulandTcal::Exec(Option_t*)
{
    /*
    if (fTrigger >= 0)
    {
        if (header->GetTrigger() != fTrigger)
        {
            return;
        }
    }
*/
    Int_t nHits = fMappedHit->GetEntriesFast();
    /*
        if (nHits > (fNofPMTs / 2))
        {
            return;
        }
    */
    R3BPaddleTamexMappedData* hit;
    Int_t iPlane;
    Int_t iBar;
    Int_t channel;
    Int_t tdc;
    R3BTCalModulePar* par;
    Double_t timeLE;
    Double_t timeTE;

    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        hit = dynamic_cast<R3BPaddleTamexMappedData*>(fMappedHit->At(ihit));
        if (NULL == hit)
        {
            continue;
        }
        iPlane = hit->GetPlaneId();
        iBar = hit->GetBarId();
        if (hit->Is17())
        {
            // 17-th channel
            continue;
        }

        for (Int_t iSide = 0; iSide < 2; iSide++)
        {

            // !!! There is the Edge missing GetModulePar( ... ) !!!
            // Convert TDC to [ns]
            if (!(par = fTcalPar->GetModuleParAt(iPlane, iBar, iSide)))
            {
                LOG(debug) << "R3BNeulandTcal::Exec : Tcal par not found, barId: " << iBar << ", side: " << iSide;
                continue;
            }

            tdc = hit->GetFineTime(iSide, 0); // PM, edge
            timeLE = par->GetTimeVFTX(tdc);

            if (!(par = fTcalPar->GetModuleParAt(iPlane, iBar, iSide + 2)))
            {
                LOG(debug) << "R3BNeulandTcal::Exec : Tcal par not found, barId: " << iBar << ", side: " << iSide;
                continue;
            }

            tdc = hit->GetFineTime(iSide, 1);
            timeTE = par->GetTimeVFTX(tdc);

            /*
                   if (timeLE < -1000.)
                   {
                       continue;
                   }
            */
            if (timeLE < 0. || timeLE > fClockFreq || timeTE < 0. || timeTE > fClockFreq)
            {
                LOG(error) << "R3BNeulandTcal::Exec : error in time calibration: ch= " << channel << ", tdc= " << tdc
                           << ", time leading edge = " << timeLE << ", time trailing edge = " << timeTE;
                continue;
            }

            timeLE = fClockFreq - timeLE + hit->GetCoarseTime(iSide, 0) * fClockFreq;
            timeTE = fClockFreq - timeTE + hit->GetCoarseTime(iSide, 1) * fClockFreq;

            new ((*fPmt)[fNPmt]) R3BNeulandPmt(iPlane, iBar, iSide, timeLE, timeTE - timeLE);
            fNPmt += 1;
        } // for Side
    }
}

void R3BNeulandTcal::FinishEvent()
{
    if (fVerbose && 0 == (fNEvents % 1))
    {
        LOG(info) << "R3BNeulandTcal::Exec : event=" << fNEvents << " nPMTs=" << fNPmt;
    }

    if (fPmt)
    {
        fPmt->Clear();
        fNPmt = 0;
    }

    fNEvents += 1;
}

void R3BNeulandTcal::FinishTask() {}

ClassImp(R3BNeulandTcal)
