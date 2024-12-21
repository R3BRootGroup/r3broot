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
// -----                  R3BTofiMapped2Cal                -----
// -----          Created March 2021 by M.Heil          -----
// ------------------------------------------------------------

/*
 * March 2016
 * Reconstructing data into PaddleItems. Each paddle holds two times,
 * (for PM1 and PM2). This allows to plot
 * the time differences e.g. via evt->Draw(...) interactively (aka without
 * looping over all channels) which is crucial for a quick check of the
 * detector status during the experiment.
 */

#include <FairRootManager.h>
#include <assert.h>

#include "R3BTofiMapped2Cal.h"

#include "TClonesArray.h"
#include "TMath.h"

#include "FairLogger.h"
#include "FairRuntimeDb.h"

#include "R3BTCalEngine.h"
#include "R3BTofiCalData.h"
#include "R3BTofiMappedData.h"

#include "R3BEventHeader.h"

#define IS_NAN(x) TMath::IsNaN(x)

namespace
{
    double const c_range_ns = 2048 * 5;
};

R3BTofiMapped2Cal::R3BTofiMapped2Cal()
    : FairTask("R3BTofiMapped2Cal", 1)
    , fMappedItems(nullptr)
    , fMappedTriggerItems(nullptr)
    , fCalItems(new TClonesArray("R3BTofiCalData"))
    , fCalTriggerItems(new TClonesArray("R3BTofiCalData"))
    , fTcalPar(0)
    , fNofTcalPars(0)
    , fNofPlanes(0)
    , fPaddlesPerPlane(0)
    , fTrigger(-1)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fCalLookup()
{
}

R3BTofiMapped2Cal::R3BTofiMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems(nullptr)
    , fMappedTriggerItems(nullptr)
    , fCalItems(new TClonesArray("R3BTofiCalData"))
    , fCalTriggerItems(new TClonesArray("R3BTofiCalData"))
    , fTcalPar(0)
    , fNofTcalPars(0)
    , fNofPlanes(0)
    , fPaddlesPerPlane(0)
    , fTrigger(-1)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fCalLookup()
{
}

R3BTofiMapped2Cal::~R3BTofiMapped2Cal()
{
    delete fCalItems;
    delete fCalTriggerItems;
}

size_t R3BTofiMapped2Cal::GetCalLookupIndex(R3BTofiMappedData const& a_mapped) const
{
    size_t i = ((a_mapped.GetDetectorId() - 1) * fPaddlesPerPlane + (a_mapped.GetBarId() - 1)) * 2 +
               (a_mapped.GetSideId() - 1);

    //  cout<<i<<";"<<fCalLookup.size()<<"; "<<a_mapped.GetDetectorId()<<"; "<<a_mapped.GetBarId()<<endl   ;

    assert(i < fCalLookup.size());
    return i;
}

InitStatus R3BTofiMapped2Cal::Init()
{
    fNofTcalPars = fTcalPar->GetNumModulePar();
    if (fNofTcalPars == 0)
    {
        LOG(error) << "There are no TCal parameters in container TofiTCalPar";
        return kFATAL;
    }

    LOG(info) << "R3BTofiMapped2Cal::Init : read " << fNofTcalPars << " modules";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    // get access to Mapped data
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofiMapped"));
    if (NULL == fMappedItems)
        LOG(fatal) << "Branch TofiMapped not found";
    fMappedTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofiTriggerMapped"));
    if (NULL == fMappedTriggerItems)
        LOG(fatal) << "Branch TofiTriggerMapped not found";

    // request storage of Cal data in output tree
    mgr->Register("TofiCal", "Land", fCalItems, kTRUE);
    mgr->Register("TofiTriggerCal", "Land", fCalTriggerItems, kTRUE);

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BTofiMapped2Cal::SetParContainers()
{
    fTcalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("TofiTCalPar"));
    if (!fTcalPar)
    {
        LOG(error) << "Could not get access to TofiTCalPar-Container.";
        fNofTcalPars = 0;
    }
}

InitStatus R3BTofiMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BTofiMapped2Cal::Exec(Option_t* option)
{
    // check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    for (auto it = fCalLookup.begin(); fCalLookup.end() != it; ++it)
    {
        auto& vec = *it;
        vec.clear();
    }

    Int_t mapped_num = fMappedItems->GetEntriesFast();

    // Calibrate time to nanoseconds.
    struct Cal
    {
        Cal(R3BTofiMappedData const* a_mapped, double a_time_ns)
            : mapped(a_mapped)
            , time_ns(a_time_ns)
        {
        }
        R3BTofiMappedData const* mapped;
        double time_ns;
    };
    std::vector<std::vector<Cal>> cal_vec(fNofPlanes * fPaddlesPerPlane * 2);
    for (Int_t mapped_i = 0; mapped_i < mapped_num; mapped_i++)
    {
        auto mapped = dynamic_cast<R3BTofiMappedData const*>(fMappedItems->At(mapped_i));

        if ((mapped->GetDetectorId() < 1) || (mapped->GetDetectorId() > fNofPlanes))
        {
            LOG(debug) << "R3BTofiMapped2Cal::Exec : Plane number out of range: " << mapped->GetDetectorId();
            continue;
        }
        if ((mapped->GetBarId() < 1) || (mapped->GetBarId() > fPaddlesPerPlane))
        {
            LOG(debug) << "R3BTofiMapped2Cal::Exec : Bar number out of range: " << mapped->GetBarId() << ", "
                       << fPaddlesPerPlane;
            continue;
        }

        // Tcal parameters.
        auto* par = fTcalPar->GetModuleParAt(
            mapped->GetDetectorId(), mapped->GetBarId(), 2 * mapped->GetSideId() + mapped->GetEdgeId() - 2);
        if (!par)
        {
            LOG(error) << "R3BTofiMapped2Cal::Exec : Tcal par not found, Plane: " << mapped->GetDetectorId()
                       << ", Bar: " << mapped->GetBarId() << ", Side: " << mapped->GetSideId()
                       << ", Edge: " << mapped->GetEdgeId();
            continue;
        }

        // Convert TDC to [ns] ...
        Double_t time_ns = par->GetTimeVFTX(mapped->GetTimeFine());
        // ... and subtract it from the next clock cycle.
        time_ns = (mapped->GetTimeCoarse() + 1) * fClockFreq - time_ns;

        //       cout<<"time_ns: "<<time_ns<<". "<<mapped->GetBarId()<<"; "<<mapped->GetSideId()<<",
        //       "<<mapped->GetEdgeId()<<endl;

        auto& entry = cal_vec.at(((mapped->GetDetectorId() - 1) * fPaddlesPerPlane + mapped->GetBarId() - 1) * 2 +
                                 mapped->GetSideId() - 1);
        entry.push_back(Cal(mapped, time_ns));
    }

    // Iterate through calibrated times and match leading/trailing pairs.
    // Note that the reader saves all leading edges first, then appends the trailing.
    //
    // TODO: This algo would not properly handle e.g. double leading before a trailing,
    //       should be fixed to be future proof!
    //
    // puts("---");
    // puts("Event");
    for (auto it = cal_vec.begin(); cal_vec.end() != it; ++it)
    {
        auto& ch = *it;
        if (ch.empty())
            continue;
        // for (auto it2 = ch.begin(); ch.end() != it2; ++it2) {
        // std::cout << it2->mapped->GetEdgeId() << ": " << it2->time_ns << '\n';
        //}
        size_t lead_i = 0;
        size_t trail_i;
        for (trail_i = 0; trail_i < ch.size() && 2 != ch.at(trail_i).mapped->GetEdgeId(); ++trail_i)
            ;
        // std::cout << "Trail=" << trail_i << '\n';
        for (;;)
        {
            for (; trail_i < ch.size(); ++trail_i)
            {
                /*
                                cout << "lead: " << lead_i << " trail: " << trail_i << endl;
                                cout << "data in Det:" << ch.at(trail_i).mapped->GetDetectorId() << " Bar: " <<
                   ch.at(trail_i).mapped->GetBarId()
                                 << " Side: " << ch.at(trail_i).mapped->GetSideId() <<" lead: " << ch.at(lead_i).time_ns
                                 << " trail: " << ch.at(trail_i).time_ns
                                 << " ToT: " << ch.at(trail_i).time_ns - ch.at(lead_i).time_ns << endl;
                */
                auto dt = ch.at(trail_i).time_ns - ch.at(lead_i).time_ns;
                if (dt < -c_range_ns / 2)
                {
                    // Wrap-around.
                    // cout << "wrap-around"<<endl;
                    break;
                }
                if (dt > c_range_ns / 2)
                {
                    // std::cout << lead_i << ' ' << trail_i << ": " << ch.at(lead_i).time_ns << ' ' <<
                    // ch.at(trail_i).time_ns << '\n';
                    // Missing leading edge with wrap-around at the same time.
                    // cout << "Missing leading edge with wrap-around at the same time." << endl;
                    continue;
                }
                if (dt > 0)
                {
                    break;
                }
            }
            if (trail_i == ch.size() || 1 != ch.at(lead_i).mapped->GetEdgeId())
            {
                break;
            }

            auto lead = ch.at(lead_i).mapped;
            auto cal = new ((*fCalItems)[fCalItems->GetEntriesFast()]) R3BTofiCalData(lead->GetDetectorId(),
                                                                                      lead->GetBarId(),
                                                                                      lead->GetSideId(),
                                                                                      ch.at(lead_i).time_ns,
                                                                                      ch.at(trail_i).time_ns);

            //    cout << "stored Det:" << lead->GetDetectorId() << " Bar: " << lead->GetBarId()
            //         << " Side: " << lead->GetSideId() <<" lead: " << ch.at(lead_i).time_ns
            //         << " trail: " << ch.at(trail_i).time_ns
            //         << " ToT: " << ch.at(trail_i).time_ns - ch.at(lead_i).time_ns << endl;

            ++lead_i;
            ++trail_i;
        }
    }

    // Calibrate trigger channels.
    mapped_num = fMappedTriggerItems->GetEntriesFast();
    for (Int_t mapped_i = 0; mapped_i < mapped_num; mapped_i++)
    {
        auto mapped = dynamic_cast<R3BTofiMappedData const*>(fMappedTriggerItems->At(mapped_i));

        if (mapped->GetDetectorId() != fNofPlanes + 1)
        {
            LOG(debug) << "R3BTofiMapped2Cal::Exec : Trigger plane number out of range: " << mapped->GetDetectorId();
            continue;
        }

        // Tcal parameters.
        auto* par = fTcalPar->GetModuleParAt(mapped->GetDetectorId(), mapped->GetBarId(), 1);
        if (!par)
        {
            LOG(info) << "R3BTofiMapped2Cal::Exec : Trigger Tcal par not found, Plane: " << mapped->GetDetectorId()
                      << ", Bar: " << mapped->GetBarId();
            continue;
        }

        // Convert TDC to [ns] ...
        Double_t time_ns = par->GetTimeVFTX(mapped->GetTimeFine());
        // ... and subtract it from the next clock cycle.
        time_ns = (mapped->GetTimeCoarse() + 1) * fClockFreq - time_ns;

        new ((*fCalTriggerItems)[fCalTriggerItems->GetEntriesFast()])
            R3BTofiCalData(mapped->GetDetectorId(), mapped->GetBarId(), 1, time_ns, 0);
    }
}

void R3BTofiMapped2Cal::FinishEvent()
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

void R3BTofiMapped2Cal::FinishTask() {}

void R3BTofiMapped2Cal::SetNofModules(Int_t planes, Int_t ppp)
{
    fNofPlanes = planes;
    fPaddlesPerPlane = ppp;
    // #planes * #bars * 2, each entry pair is a full bar, i.e. merged sides +
    // edges.
    fCalLookup.resize(planes * ppp * 2);
}

ClassImp(R3BTofiMapped2Cal)
