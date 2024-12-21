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
// -----                  R3BSci8Mapped2Cal                -----
// -----          Created August 7th 2018 by A. Kelic-Heil -----
// ------------------------------------------------------------

/*
 * Following R3BLosMapped2Cal
 */

#include "R3BSci8Mapped2Cal.h"

#include "R3BEventHeader.h"
#include "R3BSci8MappedData.h"
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

#define SCI8_COINC_WINDOW_V_NS 280
#define SCI8_COINC_WINDOW_T_NS 50 // Same as VFTX, as leading and trailing times are separately treated
#define IS_NAN(x) TMath::IsNaN(x)

R3BSci8Mapped2Cal::R3BSci8Mapped2Cal()
    : FairTask("Sci8Tcal", 1)
    , fMappedItems(NULL)
    , fCalItems(new TClonesArray("R3BSci8CalData"))
    , fNofCalItems(0)
    , fNofTcalPars(0)
    , fNofModules(0)
    , fTcalPar(NULL)
    , fTrigger(-1)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvent(0)
{
}

R3BSci8Mapped2Cal::R3BSci8Mapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems(NULL)
    , fCalItems(new TClonesArray("R3BSci8CalData"))
    , fNofCalItems(0)
    , fNofTcalPars(0)
    , fNofModules(0)
    , fTcalPar(NULL)
    , fTrigger(-1)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvent(0)
{
}

R3BSci8Mapped2Cal::~R3BSci8Mapped2Cal()
{
    if (fCalItems)
    {
        delete fCalItems;
        fCalItems = NULL;
        fNofCalItems = 0;
    }
}

InitStatus R3BSci8Mapped2Cal::Init()
{
    fNofTcalPars = fTcalPar->GetNumModulePar();
    if (fNofTcalPars == 0)
    {
        LOG(error) << "There are no TCal parameters in container Sci8TCalPar";
        return kFATAL;
    }

    LOG(info) << "R3BSci8Mapped2Cal::Init : read " << fNofModules << " modules";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    // get access to Mapped data
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("Sci8Mapped"));

    if (NULL == fMappedItems)
        LOG(fatal) << "Branch Sci8Mapped not found";

    // request storage of Cal data in output tree
    mgr->Register("Sci8Cal", "Land", fCalItems, kTRUE);
    fCalItems->Clear();

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BSci8Mapped2Cal::SetParContainers()
{
    fTcalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("Sci8TCalPar"));
    if (!fTcalPar)
    {
        LOG(error) << "Could not get access to Sci8TCalPar-Container.";
        fNofTcalPars = 0;
        return;
    }
}

InitStatus R3BSci8Mapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BSci8Mapped2Cal::Exec(Option_t* option)
{
    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    Int_t nHits = fMappedItems->GetEntriesFast();

    for (Int_t ihit = 0; ihit < nHits; ihit++) // nHits = Nchannel_Sci8 * NTypes = 2 * 3
    {

        R3BSci8MappedData* hit = dynamic_cast<R3BSci8MappedData*>(fMappedItems->At(ihit));
        if (!hit)
            continue;

        // channel numbers are stored 1-based (1..n)
        UInt_t iDet = hit->GetDetector(); // 1..
        UInt_t iCha = hit->GetChannel();  // 1..
        UInt_t iType = hit->GetType();    // 0,1,2

        if ((iDet < 1) || (iDet > fNofDetectors))
        {
            LOG(info) << "R3BSci8Mapped2Cal::Exec : Detector number out of range: " << iDet;
            continue;
        }

        // Fetch calib data for current channel
        // new:

        R3BTCalModulePar* par = fTcalPar->GetModuleParAt(iDet, iCha, iType + 1);

        if (!par)
        {
            LOG(info) << "R3BSci8Mapped2Cal::Exec : Tcal par not found, Detector: " << iDet << ", Channel: " << iCha
                      << ", Type: " << iType;
            continue;
        }

        // Convert TDC to [ns] ...

        Double_t times_raw_ns = par->GetTimeVFTX(hit->GetTimeFine());

        if (times_raw_ns < 0. || times_raw_ns > fClockFreq || IS_NAN(times_raw_ns))
        {

            LOG(info) << "R3BSci8Mapped2Cal::Exec : Bad time in ns: det= " << iDet << ", ch= " << iCha
                      << ", type= " << iType << ", time in channels = " << hit->GetTimeFine()
                      << ", time in ns = " << times_raw_ns;
            continue;
        }

        // ... and add clock time
        Double_t times_ns = fClockFreq - times_raw_ns + hit->GetTimeCoarse() * fClockFreq;

        /* Note: we have multi-hit data...
         *
         * So the map needs to have one item per detector and (multi-)hit
         * Then we need to establish a time window
         * Here, we have the hits unsorted in time and channel. If we
         * reconstruct a detector hit using a time window, So:
         *
         * For each single hit, search the list of detector hits. If a
         * matching hit is found (dt < window and item not yet set), add
         * item. Else create new detector hit.
         *
         * This way, we theoretically *might* end up with two calItems
         * which are actually just one. Hm... this should be very rare.
         * Care about that later if it becomes necessary.
         *
         */

        // see if there is already a detector hit around that time
        R3BSci8CalData* calItem = NULL;

        int iCal;
        for (iCal = 0; iCal < fNofCalItems; iCal++)
        {
            R3BSci8CalData* aCalItem = dynamic_cast<R3BSci8CalData*>(fCalItems->At(iCal));

            if (aCalItem->GetDetector() != iDet)
            {
                // Do not consider an item for another detector.
                continue;
            }

            Double_t SCI8_COINC_WINDOW_NS;
            Double_t Tdev;
            Bool_t SCI8_COINC = false;

            if (iType == 0)
            {
                SCI8_COINC_WINDOW_NS = SCI8_COINC_WINDOW_V_NS;
                Tdev = fabs(aCalItem->GetMeanTimeVFTX() - times_ns);
                if (Tdev < SCI8_COINC_WINDOW_NS)
                    SCI8_COINC = true;
            }
            if (iType == 1)
            {
                SCI8_COINC_WINDOW_NS = SCI8_COINC_WINDOW_T_NS;
                Tdev = fabs(aCalItem->GetMeanTimeTAMEXL() - times_ns);
                if (Tdev < SCI8_COINC_WINDOW_NS && aCalItem->GetTAMEXLNcha() > 0)
                    SCI8_COINC = true;
                if (IS_NAN(Tdev) && aCalItem->GetTAMEXLNcha() == 0)
                    SCI8_COINC = true; // First Tamex leading time
            }
            if (iType == 2)
            {
                SCI8_COINC_WINDOW_NS = SCI8_COINC_WINDOW_T_NS;
                Tdev = fabs(aCalItem->GetMeanTimeTAMEXT() - times_ns);
                if (Tdev < SCI8_COINC_WINDOW_NS && aCalItem->GetTAMEXTNcha() > 0)
                    SCI8_COINC = true;
                if (IS_NAN(Tdev) && aCalItem->GetTAMEXTNcha() == 0)
                    SCI8_COINC = true; // First Tamex trailing time
            }

            if (SCI8_COINC)
            {
                // check if item is already set. If so, we need to skip this event!
                switch (iCha)
                {

                    case 0:
                    {
                        if (iType == 0 && !IS_NAN(aCalItem->fTimeV_r_ns))
                            goto skip_event_pileup;
                        if (iType == 1 && !IS_NAN(aCalItem->fTimeL_r_ns))
                            goto skip_event_pileup;
                        if (iType == 2 && !IS_NAN(aCalItem->fTimeT_r_ns))
                            goto skip_event_pileup;
                    }
                    break;
                    case 1:
                    {
                        if (iType == 0 && !IS_NAN(aCalItem->fTimeV_l_ns))
                            goto skip_event_pileup;
                        if (iType == 1 && !IS_NAN(aCalItem->fTimeL_l_ns))
                            goto skip_event_pileup;
                        if (iType == 2 && !IS_NAN(aCalItem->fTimeT_l_ns))
                            goto skip_event_pileup;
                    }
                    break;
                }
                if (!calItem)
                    calItem = aCalItem;
            }
        }
        if (!calItem)
        {
            // there is no detector hit with matching time. Hence, create a new one.
            calItem = new ((*fCalItems)[fNofCalItems]) R3BSci8CalData(iDet);
            fNofCalItems += 1;
        }
        // set the time to the correct cal item
        if (iCha == 1)
        {

            if (iType == 0)
            {
                calItem->fTimeV_r_ns = times_ns;
                if (calItem->fTimeV_r_ns < 0. || IS_NAN(calItem->fTimeV_r_ns))
                    LOG(info) << "Problem with  fTimeV_r_ns: " << calItem->fTimeV_r_ns << " " << times_ns << " "
                              << endl;
            }

            if (iType == 1)
            {
                calItem->fTimeL_r_ns = times_ns;
                if (calItem->fTimeL_r_ns < 0. || IS_NAN(calItem->fTimeL_r_ns))
                    LOG(info) << "Problem with  fTimeL_r_ns: " << calItem->fTimeL_r_ns << " " << times_ns << " "
                              << endl;
            }

            if (iType == 2)
            {
                calItem->fTimeT_r_ns = times_ns;
                if (calItem->fTimeT_r_ns < 0. || IS_NAN(calItem->fTimeT_r_ns))
                    LOG(info) << "Problem with  fTimeT_r_ns: " << calItem->fTimeT_r_ns << " " << times_ns << " "
                              << endl;
            }
        }

        if (iCha == 2)
        {

            if (iType == 0)
            {
                calItem->fTimeV_l_ns = times_ns;
                if (calItem->fTimeV_l_ns < 0. || IS_NAN(calItem->fTimeV_l_ns))
                    LOG(info) << "Problem with  fTimeV_l_ns: " << calItem->fTimeV_l_ns << " " << times_ns << " "
                              << endl;
            }

            if (iType == 1)
            {
                calItem->fTimeL_l_ns = times_ns;
                if (calItem->fTimeL_l_ns < 0. || IS_NAN(calItem->fTimeL_l_ns))
                    LOG(info) << "Problem with  fTimeL_l_ns: " << calItem->fTimeL_l_ns << " " << times_ns << " "
                              << endl;
            }

            if (iType == 2)
            {
                calItem->fTimeT_l_ns = times_ns;
                if (calItem->fTimeT_l_ns < 0. || IS_NAN(calItem->fTimeT_l_ns))
                    LOG(info) << "Problem with  fTimeT_l_ns: " << calItem->fTimeT_l_ns << " " << times_ns << " "
                              << endl;
            }
        }

        continue;
    skip_event_pileup:
        LOG(warn) << "R3BSci8Mapped2Cal::Exec : " << fNEvent << " iCha: " << iCha << " iType: " << iType
                  << " iCal: " << iCal << " Skip event because of pileup.";
    }

    ++fNEvent;
}

void R3BSci8Mapped2Cal::FinishEvent()
{
    if (fCalItems)
    {
        fCalItems->Clear();
        fNofCalItems = 0;
    }
}

void R3BSci8Mapped2Cal::FinishTask() {}

ClassImp(R3BSci8Mapped2Cal)
