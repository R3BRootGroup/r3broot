/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BLosProvideTStart.h"
#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"

#include <FairRootManager.h>

R3BLosProvideTStart::R3BLosProvideTStart()
    : FairTask("R3BLosProvideTStart", 0)
    , fLosCalData("LosCal")
    , fLosTriggerCalData("LosTriggerCal")
    , fLosHitData("LosHit")
    , fLosTriggerData("LosTriggerTCal")
    , fEventHeader(nullptr)
    , edgeL(0.)
    , edgeR(0.)
    , fUseTrigHit(kFALSE)
{
}

InitStatus R3BLosProvideTStart::Init()
{
    R3BLOG(info, "");
    fLosCalData.Init();
    fLosTriggerCalData.Init();

    if (fUseTrigHit)
    {
        fLosHitData.Init();
        fLosTriggerData.Init();
    }

    auto ioman = FairRootManager::Instance();
    if (ioman == nullptr)
    {
        throw std::runtime_error("R3BLosProvideTStart: No FairRootManager");
    }

    fEventHeader = dynamic_cast<R3BEventHeader*>(ioman->GetObject("EventHeader."));
    if (fEventHeader == nullptr)
    {
        R3BLOG(fatal, "EventHeader. not found");
    }
    // Definition of a time stich object to correlate times coming from different systems
    fTimeStitch = new R3BCoarseTimeStitch();

    return kSUCCESS;
}

void R3BLosProvideTStart::Exec(Option_t*)
{
    if (fUseTrigHit)
    {
        fEventHeader->SetTStart(GetTStartTrigHit());
        return;
    }

    fEventHeader->SetTStart(GetTStart());
    fEventHeader->SetTStartSimple(GetTStart_without_trigger());
}

Double_t R3BLosProvideTStart::GetTStart_without_trigger() const
{
    const auto losCalData = fLosCalData.Retrieve();
    if (losCalData.empty())
    {
        return std::numeric_limits<Double_t>::quiet_NaN();
    }
    else
    {
        return losCalData.back()->GetMeanTimeTAMEXL();
    }
}

Double_t R3BLosProvideTStart::GetTStart() const
{
    const auto losCalData = fLosCalData.Retrieve();
    const auto losTriggerCalData = fLosTriggerCalData.Retrieve();
    Float_t diff;
    if (losCalData.empty())
    {
        return std::numeric_limits<Double_t>::quiet_NaN();
    }
    else if (losTriggerCalData.empty())
    {
        R3BLOG(debug1, "CalData info for LOS");
        return losCalData.back()->GetMeanTimeVFTX();
    }
    else
    {
        if (losTriggerCalData.back()->GetTimeV_ns(0) > 0.)
        {
            R3BLOG(debug1, "CalData with VFTX trigger info for LOS");

            return fTimeStitch->GetTime(
                losCalData.back()->GetMeanTimeVFTX() - losTriggerCalData.back()->GetTimeV_ns(0), "vftx", "vftx");
        }
        else
        {
            R3BLOG(debug1, "CalData with Tamex trigger info for LOS");
            return fTimeStitch->GetTime(
                losCalData.back()->GetMeanTimeVFTX() - losTriggerCalData.back()->GetTimeL_ns(0), "vftx", "tamex");
        }
    }
}

Double_t R3BLosProvideTStart::GetTStartTrigHit() const
{
    const auto losHitData = fLosHitData.Retrieve();
    const auto losTriggerData = fLosTriggerData.Retrieve();
    if (losHitData.empty())
    {
        return std::numeric_limits<Double_t>::quiet_NaN();
    }
    else if (losTriggerData.empty())
    {
        return std::numeric_limits<Double_t>::quiet_NaN();
    }
    else
    {
        for (auto it = losHitData.rbegin(); it != losHitData.rend(); ++it)
        {
            Double_t tref_t =
                fTimeStitch->GetTime((*it)->GetTime() - losTriggerData.front()->GetRawTimeNs(), "vftx", "vftx");

            if (tref_t > edgeL && tref_t < edgeR)
                return tref_t;
        }
    }
    return std::numeric_limits<Double_t>::quiet_NaN();
}

bool R3BLosProvideTStart::IsBeam() const { return !std::isnan(GetTStart()); }

ClassImp(R3BLosProvideTStart)
