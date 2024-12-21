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

#include "R3BStrawtubesMapped2Cal.h"

#include "R3BStrawtubesCalData.h"
#include "R3BStrawtubesMappedData.h"
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

#define COINC_NS 20
#define IS_NAN(x) TMath::IsNaN(x)
#define VFTX2_COARSE_NS (1000.0 / VFTX_CLOCK_MHZ)

R3BStrawtubesMapped2Cal::R3BStrawtubesMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems()
    , fCalItems(new TClonesArray("R3BStrawtubesCalData"))
    , fNofCalItems(0)
    , fTcalPar()
{
}

R3BStrawtubesMapped2Cal::~R3BStrawtubesMapped2Cal()
{
    if (fCalItems)
    {
        delete fCalItems;
    }
}

void R3BStrawtubesMapped2Cal::Exec(Option_t* option)
{
    Int_t mapped_num = fMappedItems->GetEntriesFast();
    for (Int_t mapped_i = 0; mapped_i < mapped_num; mapped_i++)
    {
        R3BStrawtubesMappedData* mapped = dynamic_cast<R3BStrawtubesMappedData*>(fMappedItems->At(mapped_i));
        if (!mapped)
        {
            LOG(error) << "R3BStrawtubesMapped2Cal::Exec: NULL mapped item?";
            continue;
        }

        R3BTCalModulePar* tcal = fTcalPar->GetModuleParAt(mapped->GetPlane(), mapped->GetStraw(), mapped->GetSide());
        if (!tcal)
        {
            LOG(error) << "R3BStrawtubesMapped2Cal::Exec: Tcal par not found, "
                          "mapped=(plane="
                       << (int)mapped->GetPlane() << ",side=" << (int)mapped->GetSide()
                       << ",straw=" << (int)mapped->GetStraw() << ")";
            continue;
        }

        Double_t time_ns = tcal->GetTimeVFTX(mapped->GetTimeFine());
        if (time_ns < 0. || time_ns > VFTX2_COARSE_NS)
        {
            LOG(error) << "R3BStrawtubesMapped2Cal::Exec: Bad fine time, "
                          "mapped=(plane="
                       << mapped->GetPlane() << ",side=" << mapped->GetSide() << ",straw=" << mapped->GetStraw()
                       << ", raw=" << mapped->GetTimeFine() << ", ns=" << time_ns << ".";
            continue;
        }

        time_ns = VFTX2_COARSE_NS - time_ns + mapped->GetTimeCoarse() * VFTX2_COARSE_NS;

        R3BStrawtubesCalData* cal = NULL;
        for (int cal_i = 0; cal_i < fCalItems->GetEntriesFast(); ++cal_i)
        {
            R3BStrawtubesCalData* cal_test = dynamic_cast<R3BStrawtubesCalData*>(fCalItems->At(cal_i));
            // We want the same straw.
            if (cal_test->GetPlane() != mapped->GetPlane() || cal_test->GetStraw() != mapped->GetStraw())
            {
                continue;
            }
            // Pile-up galore?
            if (COINC_NS > fabs(cal_test->GetTime(mapped->GetSide()) - time_ns))
            {
                LOG(error) << "R3BStrawtubesMapped2Cal::Exec: Crazy pileup.";
                continue;
            }
            // Out of range.
            if (COINC_NS < fabs(cal_test->GetTime(1 - mapped->GetSide()) - time_ns))
            {
                continue;
            }
            cal = cal_test;
            break;
        }
        if (!cal)
        {
            cal = new ((*fCalItems)[fNofCalItems]) R3BStrawtubesCalData(mapped->GetPlane(), mapped->GetStraw());
            ++fNofCalItems;
        }
        cal->SetTime(mapped->GetSide(), time_ns);
    }
}

void R3BStrawtubesMapped2Cal::FinishEvent()
{
    if (fCalItems)
    {
        fCalItems->Clear();
    }
    fNofCalItems = 0;
}

void R3BStrawtubesMapped2Cal::FinishTask() {}

InitStatus R3BStrawtubesMapped2Cal::Init()
{
    if (0 == fTcalPar->GetNumModulePar())
    {
        LOG(error) << "There are no TCal parameters in container StrawtubesTCalPar";
        return kFATAL;
    }

    FairRootManager* mgr = FairRootManager::Instance();
    if (!mgr)
    {
        LOG(fatal) << "FairRootManager not found";
    }

    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("StrawtubesMapped"));
    if (!fMappedItems)
    {
        LOG(fatal) << "Branch StrawtubesMapped not found";
    }

    mgr->Register("StrawtubesCal", "Land", fCalItems, kTRUE);

    return kSUCCESS;
}

InitStatus R3BStrawtubesMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BStrawtubesMapped2Cal::SetParContainers()
{
    fTcalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("StrawtubesTCalPar"));
    if (!fTcalPar)
    {
        LOG(error) << "No StrawtubesTCalPar container.";
        return;
    }
}

ClassImp(R3BStrawtubesMapped2Cal)
