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
// -----                  R3BSfibMapped2TCal                -----
// -----          Created Feb 13th 2018 by M.Heil          -----
// ------------------------------------------------------------

#include "R3BSfibMapped2Cal.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "R3BSfibCalData.h"
#include "R3BSfibMappedData.h"
#include "R3BTCalEngine.h"
#include "TClonesArray.h"
#include <FairRootManager.h>
#include <cassert>

R3BSfibMapped2Cal::R3BSfibMapped2Cal(Int_t a_verbose, enum R3BTCalEngine::CTDCVariant a_variant)
    : FairTask("R3BSfibMapped2Cal", a_verbose)
    , fTCalPar(nullptr)
    , fMappedItems(nullptr)
    , fCalItems(new TClonesArray("R3BSfibCalData"))
    , fClockFreq(1000. / (R3BTCalEngine::CTDC_16_BWD_150 == a_variant ? 150 : 250))
{
}

R3BSfibMapped2Cal::~R3BSfibMapped2Cal() { delete fCalItems; }

InitStatus R3BSfibMapped2Cal::Init()
{
    if (!fTCalPar)
    {
        LOG(error) << "TCal parameter containers missing, "
                      "did you forget SetParContainers?";
        return kERROR;
    }
    auto mgr = FairRootManager::Instance();
    if (!mgr)
    {
        LOG(error) << "FairRootManager not found.";
        return kERROR;
    }
    auto name = "SfibMapped";
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name));
    if (!fMappedItems)
    {
        LOG(error) << "Branch " << name << " not found.";
        return kERROR;
    }
    mgr->Register("SfibCal", "Land", fCalItems, kTRUE);
    return kSUCCESS;
}

void R3BSfibMapped2Cal::SetParContainers()
{
    auto name = "SfibTCalPar";
    fTCalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer(name));
    if (!fTCalPar)
    {
        LOG(error) << "Could not get access to " << name << " container.";
    }
}

InitStatus R3BSfibMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BSfibMapped2Cal::Exec(Option_t* option)
{
    auto mapped_num = fMappedItems->GetEntriesFast();
    LOG(debug) << "R3BSfibMapped2Cal::Exec:fMappedItems=" << fMappedItems->GetName() << '.';
    for (auto i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BSfibMappedData*>(fMappedItems->At(i));
        assert(mapped);

        auto channel = mapped->GetChannel();
        LOG(debug) << " R3BSfibMapped2Cal::Exec:Channel=" << channel
                   << ":Edge=" << (mapped->IsLeading() ? "Leading" : "Trailing") << '.';

        // Fetch tcal parameters.
        R3BTCalModulePar* par;
        auto tcal_channel_i = channel * 2 - (mapped->IsLeading() ? 1 : 0);
        par = fTCalPar->GetModuleParAt(1 + mapped->IsTop(), tcal_channel_i, 1);
        // std::cout << 1 + mapped->IsTop() << ' ' << channel << ' ' << tcal_channel_i << ' ' << par << std::endl;
        if (!par)
        {
            LOG(warn) << "R3BSfibMapped2Cal::Exec: Channel=" << channel << ": TCal par not found.";
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
        LOG(debug) << " R3BSfibMapped2Cal::Exec: Fine raw=" << fine_raw << " -> ns=" << fine_ns << '.';

        Double_t time_ns = -1;
        if (fine_ns < 0. || fine_ns >= fClockFreq)
        {
            LOG(error) << "R3BSfibMapped2Cal::Exec (" << fName << "): Channel=" << channel
                       << ": Bad CTDC fine time (raw=" << fine_raw << ",ns=" << fine_ns << ").";
            continue;
        }

        // Calculate final time with clock cycles.
        //		time_ns = mapped->GetCoarse() * fClockFreq +
        //		(mapped->IsLeading() ? -fine_ns : fine_ns);
        // new clock TDC firmware need here a minus
        time_ns = mapped->GetCoarse() * fClockFreq - fine_ns;

        LOG(debug) << " R3BSfibMapped2Cal::Exec: Channel=" << channel << ": Time=" << time_ns << "ns.";
        new ((*fCalItems)[fCalItems->GetEntriesFast()])
            R3BSfibCalData(1 + mapped->IsTop(), channel, mapped->IsLeading(), time_ns);
    }
}

void R3BSfibMapped2Cal::FinishEvent() { fCalItems->Clear(); }

void R3BSfibMapped2Cal::FinishTask() {}

ClassImp(R3BSfibMapped2Cal)
