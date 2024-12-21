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

#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include <FairRootManager.h>

#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "R3BTofDMapped2Cal.h"
#include "R3BTofDMappingPar.h"
#include "R3BTofdCalData.h"
#include "R3BTofdMappedData.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <assert.h>

#define IS_NAN(x) TMath::IsNaN(x)

namespace
{
    double const c_range_ns = 2048 * 5;
};

R3BTofDMapped2Cal::R3BTofDMapped2Cal()
    : R3BTofDMapped2Cal("R3BTofDMapped2Cal", 1)
{
}

R3BTofDMapped2Cal::R3BTofDMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems(NULL)
    , fMappedTriggerItems(NULL)
    , fCalItems(new TClonesArray("R3BTofdCalData"))
    , fCalTriggerItems(new TClonesArray("R3BTofdCalData"))
    , fNofTcalPars(0)
    , fNofPlanes(0)
    , fPaddlesPerPlane(0)
    , fMapPar(NULL)
    , fTcalPar(NULL)
    , fTrigger(-1)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fCalLookup()
    , fOnline(kFALSE)
{
}

R3BTofDMapped2Cal::~R3BTofDMapped2Cal()
{
    if (fCalItems)
        delete fCalItems;
    if (fCalTriggerItems)
        delete fCalTriggerItems;
}

size_t R3BTofDMapped2Cal::GetCalLookupIndex(R3BTofdMappedData const& a_mapped) const
{
    size_t i = ((a_mapped.GetDetectorId() - 1) * fPaddlesPerPlane + (a_mapped.GetBarId() - 1)) * 2 +
               (a_mapped.GetSideId() - 1);

    assert(i < fCalLookup.size());
    return i;
}

void R3BTofDMapped2Cal::SetParContainers()
{
    // Parameter Container
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(error, !rtdb, "FairRuntimeDb not found");

    fMapPar = dynamic_cast<R3BTofDMappingPar*>(FairRuntimeDb::instance()->getContainer("tofdMappingPar"));
    R3BLOG_IF(warn, !fMapPar, "Could not get access to tofdMappingPar container");

    fTcalPar = dynamic_cast<R3BTCalPar*>(rtdb->getContainer("TofdTCalPar"));
    if (!fTcalPar)
    {
        R3BLOG(error, "Could not get access to TofdTCalPar-Container.");
        fNofTcalPars = 0;
    }
    return;
}

void R3BTofDMapped2Cal::SetParameter()
{
    //--- Parameter Container ---
    R3BLOG_IF(info, fMapPar, "Nb of planes " << fMapPar->GetNbPlanes() << " and paddles " << fMapPar->GetNbPaddles());
    if (fMapPar)
    {
        fNofPlanes = fMapPar->GetNbPlanes();
        fPaddlesPerPlane = fMapPar->GetNbPaddles();
    }

    fNofTcalPars = fTcalPar->GetNumModulePar();
    R3BLOG_IF(fatal, fNofTcalPars == 0, "There are no TCal parameters in container TofdTCalPar");
    return;
}

InitStatus R3BTofDMapped2Cal::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();
    if (!mgr)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    // get access to Mapped data
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdMapped"));
    if (!fMappedItems)
    {
        R3BLOG(fatal, "TofdMapped not found");
        return kFATAL;
    }

    fMappedTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdTriggerMapped"));
    R3BLOG_IF(warn, !fMappedTriggerItems, "TofdTriggerMapped not found");

    // request storage of TCal data in output tree
    mgr->Register("TofdCal", "TofdCal data", fCalItems, !fOnline);
    if (fMappedTriggerItems)
    {
        mgr->Register("TofdTriggerCal", "TofdTriggerCal data", fCalTriggerItems, !fOnline);
    }

    SetParameter();
    R3BLOG(info, "Read " << fNofTcalPars << " modules");
    return kSUCCESS;
}

InitStatus R3BTofDMapped2Cal::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BTofDMapped2Cal::Exec(Option_t* option)
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
        Cal(R3BTofdMappedData const* a_mapped, double a_time_ns)
            : mapped(a_mapped)
            , time_ns(a_time_ns)
        {
        }
        R3BTofdMappedData const* mapped;
        double time_ns;
    };

    std::vector<std::vector<Cal>> cal_vec(fNofPlanes * fPaddlesPerPlane * 2);
    for (Int_t mapped_i = 0; mapped_i < mapped_num; mapped_i++)
    {
        auto mapped = dynamic_cast<R3BTofdMappedData const*>(fMappedItems->At(mapped_i));

        if ((mapped->GetDetectorId() < 1) || (mapped->GetDetectorId() > fNofPlanes))
        {
            R3BLOG(debug, "Plane number out of range: " << mapped->GetDetectorId());
            continue;
        }
        if ((mapped->GetBarId() < 1) || (mapped->GetBarId() > fPaddlesPerPlane))
        {
            R3BLOG(debug, "Bar number out of range: " << mapped->GetBarId() << ", " << fPaddlesPerPlane);
            continue;
        }

        // Tcal parameters.
        auto* par = fTcalPar->GetModuleParAt(
            mapped->GetDetectorId(), mapped->GetBarId(), 2 * mapped->GetSideId() + mapped->GetEdgeId() - 2);
        if (!par)
        {
            R3BLOG(error,
                   "Tcal par not found, Plane: " << mapped->GetDetectorId() << ", Bar: " << mapped->GetBarId()
                                                 << ", Side: " << mapped->GetSideId()
                                                 << ", Edge: " << mapped->GetEdgeId());
            continue;
        }

        // Convert TDC to [ns] ...
        Double_t time_ns = par->GetTimeVFTX(mapped->GetTimeFine());
        // ... and subtract it from the next clock cycle.
        time_ns = (mapped->GetTimeCoarse() + 1) * fClockFreq - time_ns;

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

        for (;;)
        {
            for (; trail_i < ch.size(); ++trail_i)
            {
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
            AddTCalData(lead->GetDetectorId(),
                        lead->GetBarId(),
                        lead->GetSideId(),
                        ch.at(lead_i).time_ns,
                        ch.at(trail_i).time_ns);
            ++lead_i;
            ++trail_i;
        }
    }

    // Calibrate trigger channels.
    if (fMappedTriggerItems)
    {
        Int_t trigger_hits = fMappedTriggerItems->GetEntriesFast();
        for (Int_t mapped_i = 0; mapped_i < trigger_hits; mapped_i++)
        {
            auto mapped = dynamic_cast<R3BTofdMappedData const*>(fMappedTriggerItems->At(mapped_i));

            if (mapped->GetDetectorId() != fNofPlanes + 1)
            {
                R3BLOG(debug, "Trigger plane number out of range: " << mapped->GetDetectorId());
                continue;
            }

            // Tcal parameters.
            auto* par = fTcalPar->GetModuleParAt(mapped->GetDetectorId(), mapped->GetBarId(), 1);
            if (!par)
            {
                R3BLOG(warn,
                       "Trigger Tcal par not found, Plane: " << mapped->GetDetectorId()
                                                             << ", Bar: " << mapped->GetBarId());
                continue;
            }

            // Convert TDC to [ns] ...
            Double_t time_ns = par->GetTimeVFTX(mapped->GetTimeFine());
            // ... and subtract it from the next clock cycle.
            time_ns = (mapped->GetTimeCoarse() + 1) * fClockFreq - time_ns;

            AddTriggerTCalData(mapped->GetDetectorId(), mapped->GetBarId(), time_ns);
        }
    }
}

void R3BTofDMapped2Cal::FinishEvent()
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

void R3BTofDMapped2Cal::SetNofModules(Int_t planes, Int_t ppp)
{
    fNofPlanes = planes;
    fPaddlesPerPlane = ppp;
    // #planes * #bars * 2, each entry pair is a full bar,
    // i.e. merged sides + edges.
    fCalLookup.resize(planes * ppp * 2);
}

// -----   Private method AddTCalData  --------------------------------------------
R3BTofdCalData* R3BTofDMapped2Cal::AddTCalData(UInt_t detid,
                                               UInt_t barid,
                                               UInt_t sideid,
                                               Double_t lead_time,
                                               Double_t trail_time)
{
    // It fills the R3BTofdCalData
    TClonesArray& clref = *fCalItems;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTofdCalData(detid, barid, sideid, lead_time, trail_time);
}

// -----   Private method AddTriggerTCalData  --------------------------------------------
R3BTofdCalData* R3BTofDMapped2Cal::AddTriggerTCalData(UInt_t detid, UInt_t barid, Double_t lead_time)
{
    // It fills the R3BTofdCalData
    TClonesArray& clref = *fCalTriggerItems;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTofdCalData(detid, barid, 1, lead_time, 0.);
}

ClassImp(R3BTofDMapped2Cal)
