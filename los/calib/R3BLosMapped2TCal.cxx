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
// -----                  R3BLosMapped2Cal                -----
// -----          Created Feb 4th 2016 by R.Plag          -----
// ------------------------------------------------------------

/* March 2016
 * Rewrote the Cal structure to provide individual leafs for the
 * left, top, right and bottom signals. This allows to plot
 * the time differences via cbmsim->Draw(...) interactively (aka without
 * looping over all channels) which is crucial for a quick check of the
 * detector status during the experiment.
 *
 *
 */

#include "R3BLosMapped2TCal.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BLosMappedData.h"
#include "R3BLosTCalData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TClonesArray.h"
#include "TMath.h"

#define IS_NAN(x) TMath::IsNaN(x)

R3BLosMapped2TCal::R3BLosMapped2TCal()
    : R3BLosMapped2TCal("R3BLosMapped2TCal", 1)
{
}

R3BLosMapped2TCal::R3BLosMapped2TCal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems(NULL)
    , fMappedTriggerItems(NULL)
    , fTCalItems(new TClonesArray("R3BLosTCalData"))
    , fTCalTriggerItems(new TClonesArray("R3BLosTCalData"))
    , fNofTCalItems(0)
    , fNofTcalPars(0)
    , fTcalPar(NULL)
    , fTrigger(-1) // trigger 1 - onspill, 2 - offspill, -1 - all events
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvent(0)
    , fOnline(kFALSE)
{
}

R3BLosMapped2TCal::~R3BLosMapped2TCal()
{
    R3BLOG(debug1, "Destructor");
    if (fTCalItems)
        delete fTCalItems;
    if (fTCalTriggerItems)
        delete fTCalTriggerItems;
}

InitStatus R3BLosMapped2TCal::Init()
{
    R3BLOG(info, "");
    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    // get access to Mapped data
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("LosMapped"));
    if (NULL == fMappedItems)
    {
        R3BLOG(fatal, "LosMapped not found");
        return kFATAL;
    }

    // get access to Trigger Mapped data
    fMappedTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("LosTriggerMapped"));
    R3BLOG_IF(warn, !fMappedTriggerItems, "LosTriggerMapped not found");

    // Request storage of TCal data in output tree
    mgr->Register("LosTCal", "LosTCal data", fTCalItems, !fOnline);
    fTCalItems->Clear();

    if (fMappedTriggerItems)
    {
        mgr->Register("LosTriggerTCal", "LosTriggerTCal data", fTCalTriggerItems, !fOnline);
        fTCalTriggerItems->Clear();
    }

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BLosMapped2TCal::SetParContainers()
{
    fTcalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("LosTCalPar"));
    if (!fTcalPar)
    {
        R3BLOG(fatal, "Could not get access to LosTCalPar-Container.");
        fNofTcalPars = 0;
        return;
    }
}

InitStatus R3BLosMapped2TCal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BLosMapped2TCal::Exec(Option_t* option)
{
    // check for requested trigger (Todo: should be done globablly / somewhere else)

    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    Int_t nHits = fMappedItems->GetEntriesFast();

    if (nHits == 0)
        return;

    for (Int_t ihit = 0; ihit < nHits; ihit++) // nHits = Nchannel_LOS * NTypes = 4 or 8 * 3
    {
        Double_t times_ns = 0. / 0.;
        Double_t times_raw_ns = 0. / 0.;

        R3BLosMappedData* hit = dynamic_cast<R3BLosMappedData*>(fMappedItems->At(ihit));
        if (!hit)
            continue;

        // channel numbers are stored 1-based (1..n)
        UInt_t iDet = hit->GetDetector(); // 1..
        UInt_t iCha = hit->GetChannel();  // 1..
        UInt_t iType = hit->GetType();    // 0,1,2,3

        if ((iDet < 1) || (iDet > fNofDetectors))
        {
            R3BLOG(warn, "Detector number out of range: " << iDet);
            continue;
        }

        if (hit->GetTimeCoarse() > 8192)
        {
            R3BLOG(warn, "Coarse counter > 8192: Det " << iDet << " , Ch: " << iCha << " , type: " << iType);
            continue;
        }

        // Fetch calib data for current channel
        // new:

        if (iType < 3)
        {
            R3BTCalModulePar* par = fTcalPar->GetModuleParAt(iDet, iCha, iType + 1);

            if (!par)
            {
                R3BLOG(warn, "Tcal par not found, Detector: " << iDet << ", Channel: " << iCha << ", Type: " << iType);
                continue;
            }

            // Convert TDC to [ns] ...

            times_raw_ns = par->GetTimeVFTX(hit->GetTimeFine());

            if (times_raw_ns < 0. || times_raw_ns > fClockFreq || IS_NAN(times_raw_ns))
            {

                R3BLOG(warn,
                       "Bad time in ns: det= " << iDet << ", ch= " << iCha << ", type= " << iType
                                               << ", time in channels = " << hit->GetTimeFine()
                                               << ", time in ns = " << times_raw_ns);
                continue;
            }

            // ... and add clock time
            times_ns = fClockFreq - times_raw_ns + hit->GetTimeCoarse() * fClockFreq;
        }
        else
        {
            // MTDC32 time need only a factor
            times_ns = hit->GetTimeFine() / 7.8 / 1000.; // range MTDC 3->7.8ps
        }

        AddTCalData(iDet, iCha, iType, times_ns);
        continue;
    }

    // Calibrate trigger channels -----------------------------------------------
    if (fMappedTriggerItems && fMappedTriggerItems->GetEntriesFast() > 0)
    {
        auto mapped_num = fMappedTriggerItems->GetEntriesFast();
        for (Int_t mapped_i = 0; mapped_i < mapped_num; mapped_i++)
        {
            auto mapped = dynamic_cast<R3BLosMappedData const*>(fMappedTriggerItems->At(mapped_i));
            UInt_t iDetector = mapped->GetDetector();
            UInt_t iChannel = mapped->GetChannel();
            UInt_t iType = mapped->GetType() + 1; // 1, 2, ... 4

            // Tcal parameters.
            auto par = fTcalPar->GetModuleParAt(2 + iDetector, iChannel, iType);
            if (!par)
            {
                R3BLOG(warn, "Trigger Tcal par not found.");
                continue;
            }

            // Convert TDC to [ns] ...
            Double_t time_ns = par->GetTimeVFTX(mapped->GetTimeFine());
            // ... and subtract it from the next clock cycle.
            if (time_ns > 0.)
            {
                time_ns = (mapped->GetTimeCoarse() + 1) * fClockFreq - time_ns;
                AddTriggerTCalData(iDetector, iChannel, iType - 1, time_ns);
            }
        }
    }

    ++fNEvent;
}

void R3BLosMapped2TCal::FinishEvent()
{
    if (fTCalItems)
    {
        fTCalItems->Clear();
    }
    if (fTCalTriggerItems)
    {
        fTCalTriggerItems->Clear();
    }
}

R3BLosTCalData* R3BLosMapped2TCal::AddTCalData(Int_t det, Int_t ch, Int_t typ, Double_t tns)
{
    // It fills the R3BLosTcalData
    TClonesArray& clref = *fTCalItems;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BLosTCalData(det, ch, typ, tns);
}

R3BLosTCalData* R3BLosMapped2TCal::AddTriggerTCalData(Int_t det, Int_t ch, Int_t typ, Double_t tns)
{
    // It fills the R3BLosTriggerTcalData
    TClonesArray& clref = *fTCalTriggerItems;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BLosTCalData(det, ch, typ, tns);
}
ClassImp(R3BLosMapped2TCal);
