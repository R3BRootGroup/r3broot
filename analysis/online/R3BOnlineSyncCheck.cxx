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

#include "R3BOnlineSyncCheck.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BShared.h"
#include "R3BSyncCheckData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TH2F.h"
#include "THttpServer.h"

R3BOnlineSyncCheck::R3BOnlineSyncCheck()
    : R3BOnlineSyncCheck("OnlineSyncCheck", 1)
{
}

R3BOnlineSyncCheck::R3BOnlineSyncCheck(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCA(nullptr)
    , header(nullptr)
    , fTpat(-1)
    , fTrig(-1)
    , canvas(nullptr)
{
}

TString R3BOnlineSyncCheck::EnumName(int value)
{
    switch (value)
    {
        case MASTER_SC:
            return { TString("MASTER") };
        case MASTERREF_SC:
            return { TString("MASTER_REF") };
        case MUSIC_SC:
            return { TString("MUSIC") };
        case RPC_SC:
            return { TString("RPC") };
        case S2_SC:
            return { TString("S2") };
        case FOOT1_SC:
            return { TString("FOOT1") };
        case FOOT2_SC:
            return { TString("FOOT2") };
        case FOOT13_SC:
            return { TString("FOOT13") };
        case FOOT4_SC:
            return { TString("FOOT4") };
        case FOOT11_SC:
            return { TString("FOOT11") };
        case FOOT6_SC:
            return { TString("FOOT6") };
        case FOOT7_SC:
            return { TString("FOOT7") };
        case FOOT12_SC:
            return { TString("FOOT12") };
        default:
            return { TString("Invalid enum value") };
    }
}

InitStatus R3BOnlineSyncCheck::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, nullptr == mgr, "FairRootManager not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(fatal, !header, "Branch EventHeader. not found");

    fCA = dynamic_cast<TClonesArray*>(mgr->GetObject("SyncCheckData"));
    R3BLOG_IF(fatal, !fCA, "Branch SyncCheckData not found");

    canvas = new TCanvas("Sync check", "Sync check");
    canvas->Divide(NXPADS, NYPADS);

    for (int i = 0; i < NUM_SC; ++i)
    {
        canvas->cd(i + 1);
        TString hname = EnumName(i);
        auto* fh2 = R3B::root_owned<TH2F>(hname.Data(), hname.Data(), XBINS, XMIN, XMAX, YBINS, YMIN, YMAX);
        fh2->GetXaxis()->SetTitle(EnumName(MASTER_SC).Data());
        fh2->GetYaxis()->SetTitle(EnumName(i).Data());
        fh2->GetXaxis()->CenterTitle(true);
        fh2->GetYaxis()->CenterTitle(true);
        fh2->Draw("colz");
        fh2_array.push_back(fh2);
    }
    auto* mainfol = new TFolder("SYNC_CHECK", "SYNC_CHECK");
    mainfol->Add(canvas);
    run->AddObject(mainfol);
    run->GetHttpServer()->RegisterCommand("Reset_Sync_Check", Form("/Objects/%s/->Reset_Histo()", GetName()));
    return kSUCCESS;
}

void R3BOnlineSyncCheck::Reset_Histo()
{
    R3BLOG(info, "");
    for (const auto& hist : fh2_array)
    {
        hist->Reset();
    }
}

void R3BOnlineSyncCheck::Exec(Option_t* option)
{
    if ((fTpat >= 0) && (header != nullptr) && ((header->GetTpat() & fTpat) != fTpat))
    {
        return;
    }
    if ((fTrig >= 0) && (header != nullptr) && (header->GetTrigger() != fTrig))
    {
        return;
    }
    if (fCA->GetEntriesFast() == 0)
    {
        return;
    }

    auto* sdata = dynamic_cast<R3BSyncCheckData*>(fCA->At(0));
    R3BLOG_IF(fatal, !sdata, "Sync check data is empty!");

    uint32_t val = 0;
    for (int i = 0; i < NUM_SC; ++i)
    {
        switch (i)
        {
            case MUSIC_SC:
                val = sdata->GetMusic();
                break;
            case RPC_SC:
                val = sdata->GetRpc();
                break;
            case S2_SC:
                val = sdata->GetS2();
                break;
            case FOOT1_SC:
                val = sdata->GetFoot(1);
                break;
            case FOOT2_SC:
                val = sdata->GetFoot(2);
                break;
            case FOOT13_SC:
                val = sdata->GetFoot(13);
                break;
            case FOOT4_SC:
                val = sdata->GetFoot(4);
                break;
            case FOOT11_SC:
                val = sdata->GetFoot(11);
                break;
            case FOOT6_SC:
                val = sdata->GetFoot(6);
                break;
            case FOOT7_SC:
                val = sdata->GetFoot(7);
                break;
            case FOOT12_SC:
                val = sdata->GetFoot(12);
                break;
            default:
                val = 0;
                break;
        }
        fh2_array.at(i)->Fill(sdata->GetMaster(), val);
    }
}

void R3BOnlineSyncCheck::FinishEvent() { fCA->Clear(); }

void R3BOnlineSyncCheck::FinishTask() { canvas->Write(); }

ClassImp(R3BOnlineSyncCheck);
