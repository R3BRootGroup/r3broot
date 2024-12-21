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
// -----             R3BFootOnlineSpectra                 -----
// -----    Created 16/07/21 by J.L. Rodriguez-Sanchez    -----
// -----          Fill FOOT online histograms             -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with FOOT online data
 */

#include "R3BFootOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BFootCalData.h"
#include "R3BFootHitData.h"
#include "R3BFootMappedData.h"
#include "R3BLogger.h"
#include "R3BShared.h"

#include <THttpServer.h>

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>
#include <TCanvas.h>
#include <TFolder.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLine.h>

#include <TClonesArray.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

R3BFootOnlineSpectra::R3BFootOnlineSpectra()
    : R3BFootOnlineSpectra("FootOnlineSpectra", 1)
{
}

R3BFootOnlineSpectra::R3BFootOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

InitStatus R3BFootOnlineSpectra::Init()
{
    R3BLOG(info, "");
    // Looking for FairRootManager
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, mgr == nullptr, "FairRootManager not found");

    // Look for the R3BEventHeader
    fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (fEventHeader == nullptr)
    {
        R3BLOG(error, "EventHeader. not found");
    }
    else
    {
        R3BLOG(info, "EventHeader. found");
    }

    // Get access to Mapped data
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("FootMappedData"));
    R3BLOG_IF(fatal, fMappedItems == nullptr, "FootMappedData not found");

    // Get access to Cal data
    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("FootCalData"));
    R3BLOG_IF(warn, fCalItems == nullptr, "FootCalData not found");

    // Get access to Hit data
    fHitItems = dynamic_cast<TClonesArray*>(mgr->GetObject("FootHitData"));
    R3BLOG_IF(warn, fHitItems == nullptr, "FootHitData not found");

    // Create histograms for all the detectors
    // Energy range for strips
    Double_t binsE = 3000;
    Double_t minE = 0;
    Double_t maxE = 3000;

    char Name1[255];
    char Name2[255];

    // MAIN FOLDER-FOOT
    TFolder* mainfol = new TFolder("Foot", "Foot-Si info");
    TFolder* mapfol = new TFolder("Map", "Map FOOT info");
    TFolder* calfol = new TFolder("Cal", "Cal FOOT info");
    TFolder* hitfol = new TFolder("Hit", "Hit FOOT info");
    mainfol->Add(mapfol);

    auto cMap = new TCanvas("FOOT_mapped", "mapped info", 10, 10, 500, 500);
    cMap->Divide(4, 3);
    mapfol->Add(cMap);

    //================  Mapped data =====================
    fh2_EnergyVsStrip.resize(fNbDet);
    int i_pad = 1; // pad iterator
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh2_energy_vs_strip_det_%d", i + 1);
        sprintf(Name2, "Mapped energy vs strip number for FOOT Det: %d", i + 1);
        fh2_EnergyVsStrip[i] = R3B::root_owned<TH2F>(Name1, Name2, 640, 1, 641, binsE, minE, maxE);
        fh2_EnergyVsStrip[i]->GetXaxis()->SetTitle("Strip number");
        fh2_EnergyVsStrip[i]->GetYaxis()->SetTitle("Energy [channels]");
        fh2_EnergyVsStrip[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_EnergyVsStrip[i]->GetXaxis()->CenterTitle(true);
        int foot_num = i + 1;
        if (foot_num == 1 || foot_num == 2 || foot_num == 13 || foot_num == 4 || foot_num == 11 || foot_num == 6 ||
            foot_num == 7 || foot_num == 12 || foot_num == 9 || foot_num == 10 || foot_num == 15 || foot_num == 16)
        {
            cMap->cd(i_pad);
            fh2_EnergyVsStrip[i]->Draw("col");
            for (int i_asic = 1; i_asic < 10; i_asic++)
            {
                TLine* l = new TLine(64.5 * i_asic, minE, 64.5 * i_asic, maxE);
                l->Draw("same");
                l->SetLineStyle(7);
                l->SetLineWidth(1);
                l->SetLineColor(13);
            }
            i_pad++;
        }
    }

    //============ CAL data ==================

    auto cCal = new TCanvas("FOOT_cal", "cal info", 10, 10, 500, 500);
    cCal->Divide(4, 3);
    calfol->Add(cCal);

    if (fCalItems)
    {
        i_pad = 1; // pad iterator
        fh2_EnergyVsStrip_cal.resize(fNbDet);
        for (Int_t i = 0; i < fNbDet; i++)
        {
            sprintf(Name1, "fh2_energy_vs_strip_cal_det_%d", i + 1);
            sprintf(Name2, "Cal-energy vs strip number for FOOT Det: %d", i + 1);
            fh2_EnergyVsStrip_cal[i] = R3B::root_owned<TH2F>(Name1, Name2, 640, 1, 641, binsE, minE, maxE);
            fh2_EnergyVsStrip_cal[i]->GetXaxis()->SetTitle("Strip number");
            fh2_EnergyVsStrip_cal[i]->GetYaxis()->SetTitle("Energy [channels]");
            fh2_EnergyVsStrip_cal[i]->GetYaxis()->SetTitleOffset(1.4);
            fh2_EnergyVsStrip_cal[i]->GetXaxis()->CenterTitle(true);
            fh2_EnergyVsStrip_cal[i]->GetYaxis()->CenterTitle(true);
            int foot_num = i + 1;
            if (foot_num == 1 || foot_num == 2 || foot_num == 13 || foot_num == 4 || foot_num == 11 || foot_num == 6 ||
                foot_num == 7 || foot_num == 12 || foot_num == 9 || foot_num == 10 || foot_num == 15 || foot_num == 16)
            {
                cCal->cd(i_pad);
                fh2_EnergyVsStrip_cal[i]->Draw("col");
                for (int i_asic = 1; i_asic < 10; i_asic++)
                {
                    TLine* l = new TLine(64.5 * i_asic, minE, 64.5 * i_asic, maxE);
                    l->Draw("same");
                    l->SetLineStyle(7);
                    l->SetLineWidth(1);
                    l->SetLineColor(13);
                }
                i_pad++;
            }
        }
    }

    //================ HIT data ==========================

    auto cHit = new TCanvas("FOOT_hit", "hit info", 10, 10, 500, 500);
    cHit->Divide(3, 2);
    hitfol->Add(cHit);
    if (fHitItems)
    {
        i_pad = 1; // pad iterator
        fh1_pos.resize(fNbDet);
        fh1_ene.resize(fNbDet);
        for (Int_t i = 0; i < fNbDet; i++)
        { // one histo per detector
            sprintf(Name1, "fh1_pos_det_%d", i + 1);
            sprintf(Name2, "Cluster position for FOOT Det: %d", i + 1);
            fh1_pos[i] = R3B::root_owned<TH1F>(Name1, Name2, 600, -50., 50.);
            fh1_pos[i]->GetXaxis()->SetTitle("Position [mm]");
            fh1_pos[i]->GetYaxis()->SetTitle("Counts");
            fh1_pos[i]->GetYaxis()->SetTitleOffset(1.4);
            fh1_pos[i]->GetXaxis()->CenterTitle(true);
            fh1_pos[i]->GetYaxis()->CenterTitle(true);
            sprintf(Name1, "fh1_ene_det_%d", i + 1);
            sprintf(Name2, "Cluster energy for FOOT Det: %d", i + 1);
            fh1_ene[i] = R3B::root_owned<TH1F>(Name1, Name2, binsE, minE, maxE);
            fh1_ene[i]->GetXaxis()->SetTitle("Energy");
            fh1_ene[i]->GetYaxis()->SetTitle("Counts");
            fh1_ene[i]->GetYaxis()->SetTitleOffset(1.4);
            fh1_ene[i]->GetXaxis()->CenterTitle(true);
            fh1_ene[i]->GetYaxis()->CenterTitle(true);

            int foot_num = i + 1;
            if (foot_num == 2 || foot_num == 4 || foot_num == 11)
            {
                cHit->cd(i_pad);
                fh1_pos[i]->Draw();
                hitfol->Add(fh1_pos[i]);
                i_pad++;
                cHit->cd(i_pad);
                fh1_ene[i]->Draw();
                hitfol->Add(fh1_ene[i]);
                i_pad++;
                hitfol->Add(fh1_pos[i]);
            }
        }
    }

    auto cInBeam = new TCanvas("FOOT_inBeam", "inBeam info", 10, 10, 500, 500);
    cInBeam->Divide(3, 2);

    if (fHitItems)
    {
        fh1_mult.resize(fNbDet);
        fh2_BeamSpot = R3B::root_owned<TH2F>("BeamSpot", "BeamSpot", 600, -50., 50., 600, -50., 50.);
        fh2_BeamSpot->GetYaxis()->SetTitle("Position [mm]");
        fh2_BeamSpot->GetYaxis()->SetTitleOffset(1.4);
        fh2_BeamSpot->GetXaxis()->CenterTitle(true);
        fh2_BeamSpot->GetYaxis()->CenterTitle(true);
        cInBeam->cd(1);
        fh2_BeamSpot->Draw("colz");

        fh2_BeamSpotE = R3B::root_owned<TH2F>("BeamSpotEnergy", "BeamSpotEnergy", binsE, minE, maxE, binsE, minE, maxE);
        fh2_BeamSpotE->GetXaxis()->SetTitle("Energy [ch]");
        fh2_BeamSpotE->GetYaxis()->SetTitle("Energy [ch]");
        fh2_BeamSpotE->GetYaxis()->SetTitleOffset(1.4);
        fh2_BeamSpotE->GetXaxis()->CenterTitle(true);
        fh2_BeamSpotE->GetYaxis()->CenterTitle(true);
        cInBeam->cd(4);
        fh2_BeamSpotE->Draw("colz");
        for (Int_t i = 0; i < fNbDet; i++)
        {
            sprintf(Name1, "fh1_mult_strip_%d", i + 1);
            sprintf(Name2, "Cluster strip Mult for FOOT Det: %d", i + 1);
            fh1_mult[i] = R3B::root_owned<TH1F>(Name1, Name2, 10, 0, 10);
            fh1_mult[i]->GetXaxis()->SetTitle("Multplicity");
            fh1_mult[i]->GetYaxis()->SetTitle("Counts");
            fh1_mult[i]->GetYaxis()->SetTitleOffset(1.4);
            fh1_mult[i]->GetXaxis()->CenterTitle(true);
            fh1_mult[i]->GetYaxis()->CenterTitle(true);
        }
        cInBeam->cd(2);
        fh1_mult[1]->Draw();
        cInBeam->cd(3);
        fh1_mult[3]->Draw();
    }
    if (fCalItems)
    {
        mainfol->Add(calfol);
    }
    if (fHitItems)
    {
        mainfol->Add(hitfol);
    }

    // Looking for FairRunOnline
    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);
    run->AddObject(mainfol);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_FOOT", Form("/Objects/%s/->Reset_FOOT_Histo()", GetName()));

    return kSUCCESS;
}

void R3BFootOnlineSpectra::Reset_FOOT_Histo()
{
    LOG(info) << "R3BFootOnlineSpectra::Reset_FOOT_Histo";

    // Mapped data
    for (Int_t i = 0; i < fNbDet; i++)
    {
        fh2_EnergyVsStrip[i]->Reset();
    }

    // Cal data
    if (fCalItems)
        for (Int_t i = 0; i < fNbDet; i++)
        {
            fh2_EnergyVsStrip_cal[i]->Reset();
        }

    // Hit data
    if (fHitItems)
    {
        for (Int_t i = 0; i < fNbDet; i++)
        {
            fh1_pos[i]->Reset();
            fh1_ene[i]->Reset();
            fh1_mult[i]->Reset();
        }
        fh2_BeamSpot->Reset();
        fh2_BeamSpotE->Reset();
    }
    return;
}

void R3BFootOnlineSpectra::Exec(Option_t* option)
{

    if (fEventHeader->GetTrigger() != fTrigger && fTrigger > -1)
        return;

    if (fTpat > 0 && (fEventHeader->GetTpat() & fTpat) != fTpat)
        return;

    // Fill mapped data
    if (fMappedItems && fMappedItems->GetEntriesFast() > 0)
    {
        auto nHits = fMappedItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BFootMappedData* hit = dynamic_cast<R3BFootMappedData*>(fMappedItems->At(ihit));
            if (!hit)
                continue;
            fh2_EnergyVsStrip[hit->GetDetId() - 1]->Fill(hit->GetStripId(), hit->GetEnergy());
        }
    }

    // Fill cal data
    if (fCalItems && fCalItems->GetEntriesFast() > 0)
    {
        auto nHits = fCalItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BFootCalData* hit = dynamic_cast<R3BFootCalData*>(fCalItems->At(ihit));
            if (!hit)
                continue;
            fh2_EnergyVsStrip_cal[hit->GetDetId() - 1]->Fill(hit->GetStripId(), hit->GetEnergy());
        }
    }

    // Fill hit data
    if (fHitItems && fHitItems->GetEntriesFast() > 0)
    {
        auto nHits = fHitItems->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BFootHitData* hit = dynamic_cast<R3BFootHitData*>(fHitItems->At(ihit));
            if (!hit)
                continue;
            fh1_pos[hit->GetDetId() - 1]->Fill(hit->GetPos());
            fh1_ene[hit->GetDetId() - 1]->Fill(hit->GetEnergy());
            fh1_mult[hit->GetDetId() - 1]->Fill(hit->GetNbHit());
        }
        for (Int_t ihit = 0; ihit < nHits - 1; ihit++)
        {
            for (Int_t jhit = ihit + 1; jhit < nHits; jhit++)
            {
                R3BFootHitData* hitI = dynamic_cast<R3BFootHitData*>(fHitItems->At(ihit));
                R3BFootHitData* hitJ = dynamic_cast<R3BFootHitData*>(fHitItems->At(jhit));
                if (!hitI)
                    continue;
                if (!hitJ)
                    continue;
                if ((hitI->GetDetId() == 4 && hitJ->GetDetId() == 2) ||
                    (hitJ->GetDetId() == 4 && hitI->GetDetId() == 2))
                {
                    if ((hitI->GetDetId() == 4 && hitJ->GetDetId() == 2))
                    {
                        fh2_BeamSpot->Fill(hitI->GetPos(), hitJ->GetPos());
                        fh2_BeamSpotE->Fill(hitI->GetEnergy(), hitJ->GetEnergy());
                    }
                    else
                    {
                        fh2_BeamSpot->Fill(hitJ->GetPos(), hitI->GetPos());
                        fh2_BeamSpotE->Fill(hitJ->GetEnergy(), hitI->GetEnergy());
                        if ((hitI->GetDetId() == 4 && hitJ->GetDetId() == 2) ||
                            (hitJ->GetDetId() == 4 && hitI->GetDetId() == 2))
                        {
                            if ((hitI->GetDetId() == 4 && hitJ->GetDetId() == 2))
                            {
                                fh2_BeamSpot->Fill(hitI->GetPos(), hitJ->GetPos());
                                fh2_BeamSpotE->Fill(hitI->GetEnergy(), hitJ->GetEnergy());
                            }
                            else
                            {
                                fh2_BeamSpot->Fill(hitJ->GetPos(), hitI->GetPos());
                                fh2_BeamSpotE->Fill(hitJ->GetEnergy(), hitI->GetEnergy());
                            }
                        }
                    }
                }
            }
        }
    }
    fNEvents += 1;
    return;
}

void R3BFootOnlineSpectra::FinishEvent()
{
    if (fMappedItems)
    {
        fMappedItems->Clear();
    }
    if (fCalItems)
    {
        fCalItems->Clear();
    }
    if (fHitItems)
    {
        fHitItems->Clear();
    }
    return;
}

void R3BFootOnlineSpectra::FinishTask()
{
    if (fMappedItems)
    {
        for (Int_t i = 0; i < fNbDet; i++)
            fh2_EnergyVsStrip[i]->Write();
    }
    if (fCalItems)
        for (Int_t i = 0; i < fNbDet; i++)
        {
            fh2_EnergyVsStrip_cal[i]->Write();
        }

    if (fHitItems)
    {
        for (Int_t i = 0; i < fNbDet; i++)
        {
            fh1_pos[i]->Write();
            fh1_ene[i]->Write();
            fh1_mult[i]->Write();
        }
        fh2_BeamSpot->Write();
        fh2_BeamSpotE->Write();
    }
    return;
}

ClassImp(R3BFootOnlineSpectra)
