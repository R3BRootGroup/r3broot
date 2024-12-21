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
// -----               R3BFiberMAPMTMapped2TCal           -----
// -----          Created May 13th 2021 by V.Panin        -----
// ------------------------------------------------------------

#include "R3BFiberMAPMTMapped2Cal.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "R3BFiberMAPMTCalData.h"
#include "R3BFiberMappedData.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"
#include "TClonesArray.h"
#include <FairRootManager.h>
#include <cassert>

R3BFiberMAPMTMapped2Cal::R3BFiberMAPMTMapped2Cal(const char* a_name, Int_t a_verbose)
    : FairTask(TString("R3B") + a_name + "Mapped2Cal", a_verbose)
    , fName(a_name)
    , fMAPMTTCalPar(nullptr)
    , fMAPMTTrigTCalPar(nullptr)
    , fMappedItems(nullptr)
    , fCalItems(new TClonesArray("R3BFiberMAPMTCalData"))
    , fCalTriggerItems(new TClonesArray("R3BFiberMAPMTCalData"))
    , fClockFreq(1000. / 150)
    , fnEvents(0)
    , fOnline(kFALSE)
{
}

R3BFiberMAPMTMapped2Cal::~R3BFiberMAPMTMapped2Cal()
{
    if (fCalItems)
        delete fCalItems;
    if (fCalTriggerItems)
        delete fCalTriggerItems;
}

InitStatus R3BFiberMAPMTMapped2Cal::Init()
{
    R3BLOG(info, "For firber " << fName);
    auto mgr = FairRootManager::Instance();
    if (!mgr)
    {
        R3BLOG(fatal, "FairRootManager not found.");
        return kFATAL;
    }
    auto name = fName + "Mapped";
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name));
    if (!fMappedItems)
    {
        R3BLOG(fatal, "Branch " << name << " not found.");
        return kFATAL;
    }
    mgr->Register(fName + "Cal", "Fiber Cal Data", fCalItems, !fOnline);
    mgr->Register(fName + "TriggerCal", "Fiber TriggerCal Data", fCalTriggerItems, !fOnline);

    return kSUCCESS;
}

void R3BFiberMAPMTMapped2Cal::SetParContainers()
{
#define GET_TCALPAR(NAME)                                                                            \
    do                                                                                               \
    {                                                                                                \
        auto name = fName + #NAME "TCalPar";                                                         \
        f##NAME##TCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer(name)); \
        if (!f##NAME##TCalPar)                                                                       \
        {                                                                                            \
            R3BLOG(error, "Could not get access to " << name << " container.");                      \
        }                                                                                            \
    } while (0)
    GET_TCALPAR(MAPMT);
    GET_TCALPAR(MAPMTTrig);
    return;
}

InitStatus R3BFiberMAPMTMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BFiberMAPMTMapped2Cal::Exec(Option_t* option)
{
    auto mapped_num = fMappedItems->GetEntriesFast();
    R3BLOG(debug, "fMappedItems=" << fMappedItems->GetName() << '.');

    if (mapped_num == 0)
        return;

    for (auto i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BFiberMappedData*>(fMappedItems->At(i));
        assert(mapped);

        auto channel = mapped->GetChannel();
        R3BLOG(debug,
               "Channel=" << channel << ":Side=" << mapped->GetSide()
                          << ":Edge=" << (mapped->IsLeading() ? "Leading" : "Trailing") << '.');

        // Fetch tcal parameters.
        R3BTCalModulePar* par;
        if (mapped->IsTrigger())
        {
            par = fMAPMTTrigTCalPar->GetModuleParAt(1, channel, 1);
        }
        else
        {
            auto tcal_channel_i = channel * 2 - (mapped->IsLeading() ? 1 : 0);
            par = fMAPMTTCalPar->GetModuleParAt(1, tcal_channel_i, mapped->GetSide());
        }
        if (!par)
        {
            R3BLOG(warn, "(" << fName << "): Channel=" << channel << ": TCal par not found.");
            continue;
        }

        // Calibrate fine time.
        auto fine_raw = mapped->GetFine();
        if (-1 == fine_raw)
        {
            // TODO: Is this really ok?
            continue;
        }
        auto fine_ns = par->GetTimeClockTDC(fine_raw);
        R3BLOG(debug, "Fine raw=" << fine_raw << " -> ns=" << fine_ns << '.');

        // we have to differ between single PMT which is on Tamex and MAPMT which is on clock TDC
        Double_t time_ns = -1;
        if (fine_ns < 0. || fine_ns > fClockFreq)
        {
            R3BLOG_IF(debug,
                      fOnline == false,
                      "(" << fName << "): Channel=" << channel << ": Bad CTDC fine time (raw=" << fine_raw
                          << ",ns=" << fine_ns << ").");
            continue;
        }

        // Calculate final time with clock cycles.
        //		time_ns = mapped->GetCoarse() * fClockFreq +
        //		(mapped->IsLeading() ? -fine_ns : fine_ns);
        // new clock TDC firmware need here a minus
        time_ns = mapped->GetCoarse() * fClockFreq - fine_ns;

        R3BLOG(debug, "(" << fName << "): Channel=" << channel << ": Time=" << time_ns << "ns.");

        if (fName == "Fi30" || fName == "Fi31" || fName == "Fi32" || fName == "Fi33")
        {
            if (3 == mapped->GetSide())
            {
                new ((*fCalTriggerItems)[fCalTriggerItems->GetEntriesFast()])
                    R3BFiberMAPMTCalData(mapped->GetSide(), channel, mapped->IsLeading(), time_ns);
            }
            else
            {
                new ((*fCalItems)[fCalItems->GetEntriesFast()])
                    R3BFiberMAPMTCalData(mapped->GetSide(), channel, mapped->IsLeading(), time_ns);
            }
        }
        if (fName == "Fi23a" || fName == "Fi23b")
        {
            // for Fib23a and Fib23b some anode channels have two fibers. Here the channel number is not the fiber
            // number!
            if (2 == mapped->GetSide())
            {
                new ((*fCalTriggerItems)[fCalTriggerItems->GetEntriesFast()])
                    R3BFiberMAPMTCalData(mapped->GetSide(), channel, mapped->IsLeading(), time_ns);
            }
            else
            {
                Int_t iFib = 0;
                if (channel < 65)
                {
                    iFib = channel * 2;
                    new ((*fCalItems)[fCalItems->GetEntriesFast()])
                        R3BFiberMAPMTCalData(mapped->GetSide(), iFib - 1, mapped->IsLeading(), time_ns);
                    new ((*fCalItems)[fCalItems->GetEntriesFast()])
                        R3BFiberMAPMTCalData(mapped->GetSide(), iFib, mapped->IsLeading(), time_ns);
                }
                else if (channel > 64 && channel < 193)
                {
                    iFib = 64 + channel;
                    new ((*fCalItems)[fCalItems->GetEntriesFast()])
                        R3BFiberMAPMTCalData(mapped->GetSide(), iFib, mapped->IsLeading(), time_ns);
                }
                else if (channel > 192 && channel < 257)
                {
                    iFib = 258 + (channel - 193) * 2;
                    new ((*fCalItems)[fCalItems->GetEntriesFast()])
                        R3BFiberMAPMTCalData(mapped->GetSide(), iFib - 1, mapped->IsLeading(), time_ns);
                    new ((*fCalItems)[fCalItems->GetEntriesFast()])
                        R3BFiberMAPMTCalData(mapped->GetSide(), iFib, mapped->IsLeading(), time_ns);
                }
            }
        }
    }
    fnEvents++;
}

void R3BFiberMAPMTMapped2Cal::FinishEvent()
{
    if (fCalItems)
    {
        fCalItems->Clear();
    }
    if (fCalTriggerItems)
    {
        fCalTriggerItems->Clear();
    }
}

ClassImp(R3BFiberMAPMTMapped2Cal);
