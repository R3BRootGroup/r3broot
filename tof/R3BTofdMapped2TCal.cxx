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
// -----                  R3BTofdMapped2TCal              -----
// -----          Created Feb 4th 2016 by R.Plag          -----
// ------------------------------------------------------------

/* March 2016
 * Reconstructing data into PaddleItems. Each paddle holds two times,
 * (for PM1 and PM2). This allows to plot
 * the time differences e.g. via evt->Draw(...) interactively (aka without
 * looping over all channels) which is crucial for a quick check of the
 * detector status during the experiment.
 *
 */

#include "R3BTofdMapped2TCal.h"

#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "R3BTofdCalData.h"
#include "R3BTofdMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TClonesArray.h"
#include "TMath.h"

#define IS_NAN(x) TMath::IsNaN(x)

R3BTofdMapped2TCal::R3BTofdMapped2TCal()
    : R3BTofdMapped2TCal("R3BTofdMapped2TCal", 1)
{
}

R3BTofdMapped2TCal::R3BTofdMapped2TCal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems(nullptr)
    , fCalItems(new TClonesArray("R3BTofdCalData"))
    , fNofCalItems(0)
    , fTcalPar(0)
    , fNofTcalPars(0)
    , header(nullptr)
    , fTrigger(-1)
    , fNofPlanes(0)
    , fPaddlesPerPlane(0)
    , fNofEdges(0)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{
}

R3BTofdMapped2TCal::~R3BTofdMapped2TCal() { delete fCalItems; }

InitStatus R3BTofdMapped2TCal::Init()
{
    fNofTcalPars = fTcalPar->GetNumModulePar();
    if (fNofTcalPars == 0)
    {
        LOG(error) << "There are no TCal parameters in container TofdTCalPar";
        return kFATAL;
    }

    LOG(info) << "R3BTofdMapped2TCal::Init : read " << fNofTcalPars << " modules";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = (R3BEventHeader*)mgr->GetObject("EventHeader.");
    if (header == nullptr)
    {
        header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
        R3BLOG(warn, "R3BEventHeader was found instead of EventHeader.");
    }

    // get access to Mapped data
    fMappedItems = (TClonesArray*)mgr->GetObject("TofdMapped");
    R3BLOG_IF(fatal, NULL == fMappedItems, "TofdMapped not found");

    // request storage of Cal data in output tree
    mgr->Register("TofdCal", "Tofd cal data", fCalItems, kTRUE);

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BTofdMapped2TCal::SetParContainers()
{
    fTcalPar = (R3BTCalPar*)FairRuntimeDb::instance()->getContainer("TofdTCalPar");
    if (!fTcalPar)
    {
        LOG(error) << "Could not get access to TofdTCalPar-Container.";
        fNofTcalPars = 0;
        return;
    }
}

InitStatus R3BTofdMapped2TCal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BTofdMapped2TCal::Exec(Option_t* option)
{
    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    Int_t nHits = fMappedItems->GetEntriesFast();
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        R3BTofdMappedData* hit = (R3BTofdMappedData*)fMappedItems->At(ihit);

        Int_t iDetector = hit->GetDetector(); // 1..n
        Int_t iSide = hit->GetSide();         // 1/2
        Int_t iBar = hit->GetBar();           // 1..n
        Int_t iEdge = hit->GetEdge();         // 1/2

        if ((iDetector < 1) || (iDetector > fNofPlanes))
        {
            LOG(debug) << "R3BTofdMapped2TCal::Exec : Plane number out of range: " << iDetector;
            continue;
        }
        if ((iBar < 1) || (iBar > fPaddlesPerPlane))
        {
            LOG(debug) << "R3BTofdMapped2TCal::Exec : Bar number out of range: " << iBar << ", " << fPaddlesPerPlane;
            continue;
        }

        R3BTofdCalData* cal;

        // Look for a coincident CalDatum on this channel.
        Int_t nCals = fCalItems->GetEntriesFast();
        for (Int_t ical = 0; ical < nCals; ical++)
        {
            cal = (R3BTofdCalData*)fCalItems->At(ical);
            if (cal->GetBar() != iBar)
            {
                continue;
            }
            if (1 == iSide)
            {
            }
            if (IS_NAN(cal->GetBottom_ns()))
            {
            }
            else
            {
            }
        }

        R3BTofdCalData* calData = new ((*fCalItems)[fCalItems->GetEntriesFast()]) R3BTofdCalData(iDetector, iBar);

        R3BTCalModulePar const* par = fTcalPar->GetModuleParAt(iDetector, iBar, 2 * iSide + iEdge - 2);

        if (!par)
        {
            LOG(info) << "R3BTofdMapped2TCal::Exec : Tcal par not found, Plane: " << iDetector << ", Bar: " << iBar
                      << ", Edge: " << iEdge;
            continue;
        }

        // Convert TDC to [ns] ...
        Double_t time_ns = par->GetTimeVFTX(hit->GetTimeFine());
        // ... and subtract it from the next clock cycle.
        time_ns = (hit->GetTimeCoarse() + 1) * fClockFreq - time_ns;

        calData->Set(hit->GetSide() - 1, hit->GetEdge() - 1, time_ns);
    }
}

void R3BTofdMapped2TCal::FinishEvent()
{
    if (fCalItems)
    {
        fCalItems->Clear();
        fNofCalItems = 0;
    }
}

void R3BTofdMapped2TCal::FinishTask() {}

ClassImp(R3BTofdMapped2TCal);
