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

#include "R3BLosMapped2Cal.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"
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

#define LOS_COINC_WINDOW_V_NS 200;
#define LOS_COINC_WINDOW_TL_NS 200; // leading
#define LOS_COINC_WINDOW_TT_NS 400; // trailing, longer because of pileup
#define LOS_COINC_WINDOW_M_NS 400;  // 200  // ???
#define IS_NAN(x) TMath::IsNaN(x)

R3BLosMapped2Cal::R3BLosMapped2Cal()
    : R3BLosMapped2Cal("R3BLosMapped2Cal", 1)
{
}

R3BLosMapped2Cal::R3BLosMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems(NULL)
    , fMappedTriggerItems(NULL)
    , fCalItems(new TClonesArray("R3BLosCalData"))
    , fCalTriggerItems(new TClonesArray("R3BLosCalData"))
    , fNofCalItems(0)
    , fNofTcalPars(0)
    , fTcalPar(NULL)
    , fTrigger(-1) // trigger 1 - onspill, 2 - offspill, -1 - all events
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvent(0)
    , fOnline(kFALSE)
{
}

R3BLosMapped2Cal::~R3BLosMapped2Cal()
{
    R3BLOG(debug1, "Destructor");
    if (fCalItems)
        delete fCalItems;
    if (fCalTriggerItems)
        delete fCalTriggerItems;
}

InitStatus R3BLosMapped2Cal::Init()
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

    // Request storage of Cal data in output tree
    mgr->Register("LosCal", "LosCal data", fCalItems, !fOnline);
    fCalItems->Clear();

    if (fMappedTriggerItems)
    {
        mgr->Register("LosTriggerCal", "LosTriggerCal data", fCalTriggerItems, !fOnline);
        fCalTriggerItems->Clear();
    }

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BLosMapped2Cal::SetParContainers()
{
    fTcalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("LosTCalPar"));
    if (!fTcalPar)
    {
        R3BLOG(fatal, "Could not get access to LosTCalPar-Container.");
        fNofTcalPars = 0;
        return;
    }
}

InitStatus R3BLosMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BLosMapped2Cal::Exec(Option_t* option)
{
    // check for requested trigger (Todo: should be done globablly / somewhere else)

    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    Int_t nHits = fMappedItems->GetEntriesFast();

    if (nHits == 0)
        return;

    // if(nHits >0) cout<<"Mapped hits: "<<nHits<<", No det.: "<<fNofDetectors<<endl;

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

        //   cout<<"Mapped info: "<<ihit<<"; "<<iDet<<", "<<iCha<<"; "<<iType<<", timeFine "<<hit->GetTimeFine()<<endl;

        // if(fNEvent == 273 || fNEvent == 362 || fNEvent == 554)
        //   cout<<"R3BLosMapped2Cal: Channel "<<iCha<<", type "<<iType<<", nHits "<<nHits<<", ihit "<<ihit<<", timeFine
        //   "<<hit->GetTimeFine()<<endl;

        //   if(nHits%8 != 0) return;

        //   if(nHits != 24) return;

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

        // cout<<"Mapped2Cal :"<<ihit<<"; "<<iDet<<", "<<iCha<<", "<<iType<<", "<<times_ns<<",
        // "<<hit->GetTimeFine()<<endl;

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
        R3BLosCalData* calItem = NULL;

        int iCal;
        for (iCal = 0; iCal < fNofCalItems; iCal++)
        {
            R3BLosCalData* aCalItem = dynamic_cast<R3BLosCalData*>(fCalItems->At(iCal));

            //       cout<<"aCalItem->GetDetector() "<<aCalItem->GetDetector()<<"; "<<iDet<<endl;

            if (aCalItem->GetDetector() != iDet)
            {
                // Do not consider an item for another detector.
                continue;
            }

            Double_t LOS_COINC_WINDOW_NS;
            Double_t Tdev;
            Bool_t LOS_COINC = false;

            if (iType == 0)
            {
                LOS_COINC_WINDOW_NS = LOS_COINC_WINDOW_V_NS;
                Tdev = fabs(aCalItem->GetMeanTimeVFTX() - times_ns);
                if (Tdev < LOS_COINC_WINDOW_NS)
                    LOS_COINC = true;
                //	if(nHits%48 != 0 ) cout<<"Tdev VFTX: "<<Tdev<<", "<<aCalItem->GetMeanTimeVFTX()<<endl;
            }
            if (iType == 1)
            {
                LOS_COINC_WINDOW_NS = LOS_COINC_WINDOW_TL_NS;
                Tdev = fabs(aCalItem->GetMeanTimeTAMEXL() - times_ns);
                if (Tdev < LOS_COINC_WINDOW_NS && aCalItem->GetTAMEXLNcha() > 0)
                    LOS_COINC = true;
                if (IS_NAN(Tdev) && aCalItem->GetTAMEXLNcha() == 0 && aCalItem->GetVFTXNcha() == 8)
                    LOS_COINC = true; // First Tamex leading time
            }
            if (iType == 2)
            {
                LOS_COINC_WINDOW_NS = LOS_COINC_WINDOW_TT_NS;
                Tdev = fabs(aCalItem->GetMeanTimeTAMEXT() - times_ns);
                if (Tdev < LOS_COINC_WINDOW_NS && aCalItem->GetTAMEXTNcha() > 0)
                    LOS_COINC = true;
                if (IS_NAN(Tdev) && aCalItem->GetTAMEXTNcha() == 0 && aCalItem->GetVFTXNcha() == 8)
                    LOS_COINC = true; // First Tamex trailing time
            }
            if (iType == 3)
            {
                LOS_COINC_WINDOW_NS = LOS_COINC_WINDOW_M_NS;
                Tdev = fabs(aCalItem->GetMeanTimeMTDC32() - times_ns);
                if (Tdev < LOS_COINC_WINDOW_NS && aCalItem->GetMTDC32Ncha() > 0)
                    LOS_COINC = true;
                if (IS_NAN(Tdev) && aCalItem->GetMTDC32Ncha() == 0)
                    LOS_COINC = true; // First MTDC time
            }

            if (LOS_COINC)
            {
                // check if item is already set. If so, we need to skip this event!
                switch (iCha)
                {

                    case 1:
                    { // change to 1
                        /*if (iType == 0 && !IS_NAN(aCalItem->GetTimeV_ns(0)))
                            goto skip_event_pileup;
                        if (iType == 1 && !IS_NAN(aCalItem->GetTimeL_ns(0)))
                            goto skip_event_pileup;
                        if (iType == 2 && !IS_NAN(aCalItem->GetTimeT_ns(0)))
                            goto skip_event_pileup;*/
                    }
                    break;
                    case 3:
                    { // change to 3
                        /* if (iType == 0 && !IS_NAN(aCalItem->GetTimeV_ns(2)))
                             goto skip_event_pileup;
                         if (iType == 1 && !IS_NAN(aCalItem->GetTimeL_ns(2)))
                             goto skip_event_pileup;
                         if (iType == 2 && !IS_NAN(aCalItem->GetTimeT_ns(2)))
                             goto skip_event_pileup;*/
                    }
                    break;
                    case 5:
                    { // change to 5
                        /* if (iType == 0 && !IS_NAN(aCalItem->GetTimeV_ns(4)))
                             goto skip_event_pileup;
                         if (iType == 1 && !IS_NAN(aCalItem->GetTimeL_ns(4)))
                             goto skip_event_pileup;
                         if (iType == 2 && !IS_NAN(aCalItem->GetTimeT_ns(4)))
                             goto skip_event_pileup;*/
                    }
                    break;
                    case 7:
                    { // change to 7
                        /* if (iType == 0 && !IS_NAN(aCalItem->GetTimeV_ns(6)))
                             goto skip_event_pileup;
                         if (iType == 1 && !IS_NAN(aCalItem->GetTimeL_ns(6)))
                             goto skip_event_pileup;
                         if (iType == 2 && !IS_NAN(aCalItem->GetTimeT_ns(6)))
                             goto skip_event_pileup;*/
                    }
                    break;
                    case 2:
                    { // change to 2
                        /* if (iType == 0 && !IS_NAN(aCalItem->GetTimeV_ns(1)))
                             goto skip_event_pileup;
                         if (iType == 1 && !IS_NAN(aCalItem->GetTimeL_ns(1)))
                             goto skip_event_pileup;
                         if (iType == 2 && !IS_NAN(aCalItem->GetTimeT_ns(1)))
                             goto skip_event_pileup;*/
                    }
                    break;
                    case 4:
                    { // change to 4
                        /*if (iType == 0 && !IS_NAN(aCalItem->GetTimeV_ns(3)))
                            goto skip_event_pileup;
                        if (iType == 1 && !IS_NAN(aCalItem->GetTimeL_ns(3)))
                            goto skip_event_pileup;
                        if (iType == 2 && !IS_NAN(aCalItem->GetTimeT_ns(3)))
                            goto skip_event_pileup;*/
                    }
                    break;
                    case 6:
                    { // change to 6
                        /*if (iType == 0 && !IS_NAN(aCalItem->GetTimeV_ns(5)))
                            goto skip_event_pileup;
                        if (iType == 1 && !IS_NAN(aCalItem->GetTimeL_ns(5)))
                            goto skip_event_pileup;
                        if (iType == 2 && !IS_NAN(aCalItem->GetTimeT_ns(5)))
                            goto skip_event_pileup;*/
                    }
                    break;
                    case 8:
                    { // change to 8
                      // if (iType == 0 && !IS_NAN(aCalItem->GetTimeV_ns(7)))
                      //      goto skip_event_pileup;
                      //  if (iType == 1 && !IS_NAN(aCalItem->GetTimeL_ns(7)))
                      //      goto skip_event_pileup;
                      //  if (iType == 2 && !IS_NAN(aCalItem->GetTimeT_ns(7)))
                      //    goto skip_event_pileup;
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
            calItem = new ((*fCalItems)[fCalItems->GetEntriesFast()]) R3BLosCalData(iDet);
            fNofCalItems += 1;
        }
        // set the time to the correct cal item

        if (iType == 0)
        {
            calItem->fTimeV_ns[iCha - 1] = times_ns;
            if (calItem->fTimeV_ns[iCha - 1] < 0. || IS_NAN(calItem->fTimeV_ns[iCha - 1]))
                LOG(info) << "Problem with  fTimeV_ns: " << calItem->fTimeV_ns[iCha - 1] << " " << times_ns << " "
                          << endl;
        }

        if (iType == 1)
        {
            calItem->fTimeL_ns[iCha - 1] = times_ns;
            if (calItem->fTimeL_ns[iCha - 1] < 0. || IS_NAN(calItem->fTimeL_ns[iCha - 1]))
                LOG(info) << "Problem with  fTimeL_ns: " << calItem->fTimeL_ns[iCha - 1] << " " << times_ns << " "
                          << endl;
        }

        if (iType == 2)
        {
            calItem->fTimeT_ns[iCha - 1] = times_ns;
            if (calItem->fTimeT_ns[iCha - 1] < 0. || IS_NAN(calItem->fTimeT_ns[iCha - 1]))
                LOG(info) << "Problem with  fTimeT_ns: " << calItem->fTimeT_ns[iCha - 1] << " " << times_ns << " "
                          << endl;
        }

        if (iType == 3)
        {
            calItem->fTimeM_ns[iCha - 1] = times_ns;
            if (calItem->fTimeM_ns[iCha - 1] < 0. || IS_NAN(calItem->fTimeM_ns[iCha - 1]))
                LOG(info) << "Problem with  fTimeM_ns: " << calItem->fTimeM_ns[iCha - 1] << " " << times_ns << " "
                          << endl;
        }

        // if(fNEvent == 25383 || fNEvent == 322367 || fNEvent == 399481)

        // if(fNEvent == 24733) cout<<"Mapped2Cal "<<fNEvent<<"; "<<fNofCalItems<<", "<<nHits<<", "<<iCha<<",
        // "<<iType<<", "<<
        //                            times_ns<<", "<<hit->GetTimeFine()<<", "<<hit->GetTimeCoarse()<<endl;

        continue;
        // skip_event_pileup:
        //   LOG(warn) << "R3BLosMapped2Cal::Exec : " << fNEvent << " iCha: " << iCha << " iType: " << iType
        //              << " iCal: " << iCal << " Skip event because of pileup.";
    }

    // Calibrate trigger channels -----------------------------------------------
    if (fMappedTriggerItems && fMappedTriggerItems->GetEntriesFast() > 0)
    {
        auto mapped_num = fMappedTriggerItems->GetEntriesFast();
        UInt_t prevdet = 0;
        R3BLosCalData* caltrigger = NULL;
        for (Int_t mapped_i = 0; mapped_i < mapped_num; mapped_i++)
        {
            auto mapped = dynamic_cast<R3BLosMappedData const*>(fMappedTriggerItems->At(mapped_i));
            UInt_t iDetector = mapped->GetDetector();
            UInt_t iChannel = mapped->GetChannel();
            UInt_t iType = mapped->GetType() + 1; // 1, 2, ... 4

            if (!caltrigger || prevdet != mapped->GetDetector())
            {

                caltrigger = new ((*fCalTriggerItems)[fCalTriggerItems->GetEntriesFast()]) R3BLosCalData(iDetector);
                prevdet = mapped->GetDetector();
            }

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
                if (iType == 1)
                    caltrigger->fTimeV_ns[0] = time_ns;
                else if (iType == 2)
                    caltrigger->fTimeL_ns[0] = time_ns;
                else if (iType == 3)
                    caltrigger->fTimeT_ns[0] = time_ns;
                else if (iType == 4)
                    caltrigger->fTimeM_ns[0] = time_ns;
            }
        }
    }

    ++fNEvent;
}

void R3BLosMapped2Cal::FinishEvent()
{
    if (fCalItems)
    {
        fCalItems->Clear();
        fNofCalItems = 0;
    }
    if (fCalTriggerItems)
    {
        fCalTriggerItems->Clear();
    }
}

ClassImp(R3BLosMapped2Cal);
