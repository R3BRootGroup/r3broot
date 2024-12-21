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
// -----                  R3BBunchedFiberSPMTTrigMapped2TCal                -----
// -----          Created Feb 13th 2018 by M.Heil          -----
// ------------------------------------------------------------

#include "R3BBunchedFiberSPMTTrigMapped2Cal.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "R3BBunchedFiberCalData.h"
#include "R3BFiberMappedData.h"
#include "R3BTCalEngine.h"
#include "TClonesArray.h"
#include <FairRootManager.h>
#include <cassert>

R3BBunchedFiberSPMTTrigMapped2Cal::R3BBunchedFiberSPMTTrigMapped2Cal(Int_t a_verbose)
    : FairTask("R3BBunchedFiberSPMTTrigMapped2Cal", a_verbose)
    , fTCalPar(nullptr)
    , fMappedItems(nullptr)
    , fCalItems(new TClonesArray("R3BBunchedFiberCalData"))
{
}

R3BBunchedFiberSPMTTrigMapped2Cal::~R3BBunchedFiberSPMTTrigMapped2Cal()
{
    if (fCalItems)
        delete fCalItems;
}

InitStatus R3BBunchedFiberSPMTTrigMapped2Cal::Init()
{
    if (!fTCalPar)
    {
        LOG(error) << "TCal parameter containers missing, did you forget SetParContainers?";
        return kERROR;
    }
    if (0 == fTCalPar->GetNumModulePar())
    {
        LOG(error) << "No TCal parameters in containers " << fTCalPar->GetName() << ".";
        return kERROR;
    }
    auto mgr = FairRootManager::Instance();
    if (!mgr)
    {
        LOG(error) << "FairRootManager not found.";
        return kERROR;
    }
    auto name = "BunchedFiberSPMTTrigMapped";
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name));
    if (!fMappedItems)
    {
        LOG(error) << "Branch " << name << " not found.";
        return kERROR;
    }
    mgr->Register("BunchedFiberSPMTTrigCal", "Land", fCalItems, kTRUE);
    return kSUCCESS;
}

void R3BBunchedFiberSPMTTrigMapped2Cal::SetParContainers()
{
    auto name = "BunchedFiberSPMTTrigTCalPar";
    fTCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer(name));
    if (!fTCalPar)
    {
        LOG(error) << "Could not get access to " << name << " container.";
    }
}

InitStatus R3BBunchedFiberSPMTTrigMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BBunchedFiberSPMTTrigMapped2Cal::Exec(Option_t* option)
{
    auto mapped_num = fMappedItems->GetEntriesFast();
    LOG(debug) << "R3BBunchedFiberSPMTTrigMapped2Cal::Exec:fMappedItems=" << fMappedItems->GetName() << '.';
    for (auto i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BFiberMappedData*>(fMappedItems->At(i));
        auto channel = mapped->GetChannel();
        auto par = (R3BTCalModulePar*)fTCalPar->GetModuleParAt(1, channel, 1);
        if (!par)
        {
            LOG(warn) << "R3BBunchedFiberSPMTTrigMapped2Cal::Exec (" << fName << "): Channel=" << channel
                      << ": TCal par not found.";
            continue;
        }

        // Calibrate fine time.
        auto fine_raw = mapped->GetFine();
        if (-1 == fine_raw)
        {
            // TODO: Is this really ok?
            continue;
        }
        auto fine_ns = par->GetTimeVFTX(fine_raw);
        LOG(debug) << " R3BBunchedFiberSPMTTrigMapped2Cal::Exec: Fine raw=" << fine_raw << " -> ns=" << fine_ns << '.';

        if (fine_ns < 0. || fine_ns > 5.)
        {
            LOG(error) << "R3BBunchedFiberSPMTTrigMapped2Cal::Exec (" << fName << "): Channel=" << channel
                       << ": Bad Tamex fine time (raw=" << fine_raw << ",ns=" << fine_ns << ").";
            continue;
        }

        // Calculate final time with clock cycles.
        //		time_ns = mapped->GetCoarse() * fTamexFreq +
        //		(mapped->IsLeading() ? -fine_ns : fine_ns);
        // new clock TDC firmware need here a minus
        auto time_ns = mapped->GetCoarse() * 5. - fine_ns;

        LOG(debug) << " R3BBunchedFiberSPMTTrigMapped2Cal::Exec:Channel=" << channel << ": Time=" << time_ns << "ns.";
        new ((*fCalItems)[fCalItems->GetEntriesFast()])
            R3BBunchedFiberCalData(mapped->GetSide(), channel, mapped->IsLeading(), time_ns);
    }
}

void R3BBunchedFiberSPMTTrigMapped2Cal::FinishEvent()
{
    if (fCalItems)
        fCalItems->Clear();
}

ClassImp(R3BBunchedFiberSPMTTrigMapped2Cal);
