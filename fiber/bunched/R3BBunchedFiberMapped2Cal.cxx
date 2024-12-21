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
// -----             R3BBunchedFiberMapped2TCal           -----
// -----          Created Feb 13th 2018 by M.Heil         -----
// ------------------------------------------------------------

#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include <FairRootManager.h>

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberMapped2Cal.h"
#include "R3BFiberMappedData.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

#include "TClonesArray.h"
#include <cassert>

R3BBunchedFiberMapped2Cal::R3BBunchedFiberMapped2Cal(const char* a_name,
                                                     Int_t a_verbose,
                                                     enum R3BTCalEngine::CTDCVariant a_variant,
                                                     Bool_t a_skip_spmt)
    : FairTask(TString("R3B") + a_name + "Mapped2Cal", a_verbose)
    , fName(a_name)
    , fSPMTElectronics(Electronics::CTDC)
    , fCTDCVariant(a_variant)
    , fSkipSPMT(a_skip_spmt)
    , fMAPMTTCalPar(nullptr)
    , fMAPMTTrigTCalPar(nullptr)
    , fSPMTTCalPar(nullptr)
    , fMappedItems(nullptr)
    , fCalItems(new TClonesArray("R3BBunchedFiberCalData"))
    , fCalTriggerItems(new TClonesArray("R3BBunchedFiberCalData"))
    , fClockFreq(1000. / (R3BTCalEngine::CTDC_16_BWD_150 == a_variant ? 150 : 250))
    , fTamexFreq(1000. / VFTX_CLOCK_MHZ)
    , fnEvents(0)
    , fOnline(kFALSE)
{
}

R3BBunchedFiberMapped2Cal::R3BBunchedFiberMapped2Cal(const char* a_name,
                                                     Int_t a_verbose,
                                                     enum Electronics a_spmt_electronics,
                                                     enum R3BTCalEngine::CTDCVariant a_variant,
                                                     Bool_t a_skip_spmt)
    : FairTask(TString("R3B") + a_name + "Mapped2Cal", a_verbose)
    , fName(a_name)
    , fSPMTElectronics(a_spmt_electronics)
    , fCTDCVariant(a_variant)
    , fSkipSPMT(a_skip_spmt)
    , fMAPMTTCalPar(nullptr)
    , fMAPMTTrigTCalPar(nullptr)
    , fSPMTTCalPar(nullptr)
    , fMappedItems(nullptr)
    , fCalItems(new TClonesArray("R3BBunchedFiberCalData"))
    , fCalTriggerItems(new TClonesArray("R3BBunchedFiberCalData"))
    , fnEvents(0)
    , fOnline(kFALSE)
{
}

R3BBunchedFiberMapped2Cal::~R3BBunchedFiberMapped2Cal()
{
    if (fCalItems)
        delete fCalItems;
    if (fCalTriggerItems)
        delete fCalTriggerItems;
}

InitStatus R3BBunchedFiberMapped2Cal::Init()
{
    R3BLOG(info, "");
    if (fSPMTElectronics)
    {
        switch (fSPMTElectronics)
        {
            case CTDC:
                fTamexFreq = 1000. / (R3BTCalEngine::CTDC_16_BWD_150 == fCTDCVariant ? 150 : 250);
                break;
            case TAMEX:
                fTamexFreq = 1000. / VFTX_CLOCK_MHZ;
                break;
            default:
                assert(0 && "This should not happen!");
        }
        fClockFreq = 1000. / (R3BTCalEngine::CTDC_16_BWD_150 == fCTDCVariant ? 150 : 250);
    }

    auto mgr = FairRootManager::Instance();
    if (!mgr)
    {
        R3BLOG(fatal, "FairRootManager not found.");
        return kERROR;
    }
    auto name = fName + "Mapped";
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name));
    if (!fMappedItems)
    {
        R3BLOG(error, "Branch " << name << " not found.");
        return kERROR;
    }
    mgr->Register(fName + "Cal", fName + " cal data", fCalItems, !fOnline);
    mgr->Register(fName + "TriggerCal", fName + " trigger data", fCalTriggerItems, !fOnline);

    return kSUCCESS;
}

void R3BBunchedFiberMapped2Cal::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");
#define GET_TCALPAR(NAME)                                                       \
    do                                                                          \
    {                                                                           \
        auto name = fName + #NAME "TCalPar";                                    \
        f##NAME##TCalPar = dynamic_cast<R3BTCalPar*>(rtdb->getContainer(name)); \
        if (!f##NAME##TCalPar)                                                  \
        {                                                                       \
            R3BLOG(error, "Could not get access to " << name << " container");  \
        }                                                                       \
    } while (0)

    GET_TCALPAR(MAPMT);
    GET_TCALPAR(MAPMTTrig);
    GET_TCALPAR(SPMT);

    if (0 == fMAPMTTCalPar->GetNumModulePar())
    {
        R3BLOG(error, "No TCal parameters in containers " << fMAPMTTCalPar->GetName());
    }

    if ((!fSkipSPMT || 0 == fSPMTTCalPar->GetNumModulePar()))
    {
        R3BLOG(error, "No TCal parameters in containers " << fSPMTTCalPar->GetName());
    }
}

InitStatus R3BBunchedFiberMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BBunchedFiberMapped2Cal::Exec(Option_t* option)
{
    auto mapped_num = fMappedItems->GetEntriesFast();
    R3BLOG(debug, "fMappedItems=" << fMappedItems->GetName() << '.');
    if (mapped_num == 0)
    {
        return;
    }

    for (auto i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BFiberMappedData*>(fMappedItems->At(i));
        auto channel = mapped->GetChannel();
        R3BLOG(debug, "Channel=" << channel << ":Edge=" << (mapped->IsLeading() ? "Leading" : "Trailing") << '.');

        // Fetch tcal parameters.
        R3BTCalModulePar* par;
        if (mapped->IsMAPMTTrigger())
        {
            par = fMAPMTTrigTCalPar->GetModuleParAt(1, channel, 1);
        }
        else
        {
            auto tcal_channel_i = channel * 2 - (mapped->IsLeading() ? 1 : 0);
            par = mapped->IsMAPMT() ? fMAPMTTCalPar->GetModuleParAt(1, tcal_channel_i, 1)
                                    : fSPMTTCalPar->GetModuleParAt(1, tcal_channel_i, 1);
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

        // we have to differ between single PMT which is on Tamex and MAPMT which is on clock TDC
        Double_t time_ns = -1;
        if (mapped->IsMAPMT() || mapped->IsMAPMTTrigger())
        {
            auto fine_ns = par->GetTimeClockTDC(fine_raw);
            R3BLOG(debug, "Fine raw=" << fine_raw << " -> ns=" << fine_ns << '.');
            if (fine_ns < 0. || fine_ns > fClockFreq)
            {
                R3BLOG(error,
                       "(" << fName << "): Channel=" << channel << ": Bad CTDC fine time (raw=" << fine_raw
                           << ",ns=" << fine_ns << ").");
                continue;
            }

            // Calculate final time with clock cycles.
            //		time_ns = mapped->GetCoarse() * fClockFreq +
            //		(mapped->IsLeading() ? -fine_ns : fine_ns);
            // new clock TDC firmware need here a minus
            time_ns = (mapped->GetCoarse() + 1) * fClockFreq - fine_ns;

            R3BLOG(debug, "(" << fName << "): Channel=" << channel << ": Time=" << time_ns << "ns.");
        }
        else
        {
            // Calculate final time with clock cycles.
            //		time_ns = mapped->GetCoarse() * fTamexFreq +
            //		(mapped->IsLeading() ? -fine_ns : fine_ns);
            // new clock TDC firmware need here a minus
            // time_ns = mapped->GetCoarse() * fTamexFreq - fine_ns;
            auto fine_ns = par->GetTimeVFTX(fine_raw);
            if (fine_ns < 0. || fine_ns > fTamexFreq)
            {
                R3BLOG(error,
                       "(" << fName << "): Channel=" << channel << ": Bad Tamex fine time (raw=" << fine_raw
                           << ",ns=" << fine_ns << ").");
                continue;
            }
            time_ns = (mapped->GetCoarse() + 1) * fTamexFreq - fine_ns;

            R3BLOG(debug, "Channel=" << channel << ": Time=" << time_ns << "ns.");
        }
        if (3 == mapped->GetSide())
        {
            new ((*fCalTriggerItems)[fCalTriggerItems->GetEntriesFast()])
                R3BBunchedFiberCalData(mapped->GetSide(), channel, mapped->IsLeading(), time_ns);
        }
        else
        {
            new ((*fCalItems)[fCalItems->GetEntriesFast()])
                R3BBunchedFiberCalData(mapped->GetSide(), channel, mapped->IsLeading(), time_ns);
        }
    }
    fnEvents++;
}

void R3BBunchedFiberMapped2Cal::FinishEvent()
{
    if (fCalItems)
        fCalItems->Clear();
    if (fCalTriggerItems)
        fCalTriggerItems->Clear();
}

ClassImp(R3BBunchedFiberMapped2Cal);
