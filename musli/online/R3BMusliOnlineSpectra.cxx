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

/*
 * This task should fill histograms with MUSIC online data
 */

#include "R3BMusliOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BMusliCalData.h"
#include "R3BMusliHitData.h"
#include "R3BMusliMappedData.h"
#include "THttpServer.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TVector3.h"

#include "TClonesArray.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TRandom.h"
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

R3BMusliOnlineSpectra::R3BMusliOnlineSpectra()
    : R3BMusliOnlineSpectra("MusliOnlineSpectra", 1)
{
}

R3BMusliOnlineSpectra::R3BMusliOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMapItemsMusli(NULL)
    , fCalItemsMusli(NULL)
    , fHitItemsMusli(NULL)
    , fNEvents(0)
{
}

R3BMusliOnlineSpectra::~R3BMusliOnlineSpectra()
{
    LOG(info) << "R3BMusliOnlineSpectra::Delete instance";
    if (fMapItemsMusli)
        delete fMapItemsMusli;
}

InitStatus R3BMusliOnlineSpectra::Init()
{
    LOG(info) << "R3BMusliOnlineSpectra::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(fatal) << "R3BMusliOnlineSpectra::Init FairRootManager not found";
        return kFATAL;
    }
    // header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // get access to the MAPPED data of the MUSIC
    fMapItemsMusli = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliMappedData"));
    if (!fMapItemsMusli)
    {
        return kFATAL;
    }

    // get access to the CAL data of the MUSIC
    fCalItemsMusli = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliCalData"));
    if (!fCalItemsMusli)
    {
        LOG(info) << "R3BMusliOnlineSpectra::Init() No MusliCalData found";
    }

    // get access to the HIT data of the MUSIC
    fHitItemsMusli = dynamic_cast<TClonesArray*>(mgr->GetObject("MusliHitData"));
    if (!fHitItemsMusli)
    {
        LOG(info) << "R3BMusliOnlineSpectra::Init() No MusliHitData found";
    }

    // Create histograms for detectors
    char Name1[255];
    char Name2[255];

    // === =========== === //
    // === MAPPED DATA === //
    // === =========== === //

    // --- Multiplicities --- //
    cMusliMap_Mult = new TCanvas("MusliMap_Mult", "", 10, 10, 800, 700);
    cMusliMap_Mult->Divide(2, 1);

    fh1_Muslimap_mult = new TH1I("MultMapTot", "Total Multiplicity at MAPPED level", 20, -0.5, 19.5);
    fh1_Muslimap_mult->GetXaxis()->SetTitle("fSignal");
    fh1_Muslimap_mult->GetYaxis()->SetTitle("Total Counts");
    fh1_Muslimap_mult->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_mult->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_mult->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_mult->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_mult->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_mult->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_mult->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_mult->SetFillColor(31);
    cMusliMap_Mult->cd(1);
    fh1_Muslimap_mult->Draw("");

    fh2_Muslimap_mult = new TH2I(
        "MultMapDist", "Multiplicity Distribution per fSignal at MAPPED level", 20, -0.5, 19.5, 50, -0.5, 49.5);
    fh2_Muslimap_mult->GetXaxis()->SetTitle("fSignal");
    fh2_Muslimap_mult->GetYaxis()->SetTitle("Counts per Event");
    fh2_Muslimap_mult->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_mult->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_mult->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_mult->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_mult->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_mult->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_mult->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_mult->SetFillColor(31);
    cMusliMap_Mult->cd(2);
    fh2_Muslimap_mult->Draw("COL");

    // --- Energies --- //
    cMusliMap_E2 = new TCanvas("Musli_Emap2", "", 10, 10, 800, 700);
    cMusliMap_E2->Divide(4, 2);
    for (Int_t j = 0; j < 8; j++)
    {
        sprintf(Name1, "fh1_Musli_Emap_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, " Eraw A%02d and A%02d if mult==1", 2 * j + 1, 2 * j + 2);
        fh1_Muslimap_E[j] = new TH1I(Name1, Name2, 64000, 0, 64000);
        fh1_Muslimap_E[j]->GetXaxis()->SetTitle("Raw Energy [channels]");
        fh1_Muslimap_E[j]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_E[j]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_E[j]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_E[j]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_E[j]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_E[j]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_E[j]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_E[j]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_E[j]->SetFillColor(31);
        cMusliMap_E2->cd(j + 1);
        fh1_Muslimap_E[j]->Draw("");
    }

    cMusliMap_E4 = new TCanvas("Musli_Emap4", "", 10, 10, 800, 700);
    cMusliMap_E4->Divide(2, 2);
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh1_Musli_Emap_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "Eraw A%02d to A%02d if mult==1", j * 4 + 1, j * 4 + 5);
        fh1_Muslimap_E[j + 8] = new TH1I(Name1, Name2, 64000, 0, 64000);
        fh1_Muslimap_E[j + 8]->GetXaxis()->SetTitle("Raw Energy [channels]");
        fh1_Muslimap_E[j + 8]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_E[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_E[j + 8]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_E[j + 8]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_E[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_E[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_E[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_E[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_E[j + 8]->SetFillColor(31);
        cMusliMap_E4->cd(j + 1);
        fh1_Muslimap_E[j + 8]->Draw("");
    }

    cMusliMap_E8 = new TCanvas("Musli_Emap8", "", 10, 10, 800, 700);
    cMusliMap_E8->Divide(2, 1);
    for (Int_t j = 0; j < 2; j++)
    {
        sprintf(Name1, "fh1_Musli_Emap_mean_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "Eraw A%02d to A%02d if mult==1", j * 8 + 1, j * 8 + 8);
        fh1_Muslimap_E[j + 12] = new TH1I(Name1, Name2, 64000, 0, 64000);
        fh1_Muslimap_E[j + 12]->GetXaxis()->SetTitle("Raw Energy [channels]");
        fh1_Muslimap_E[j + 12]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_E[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_E[j + 12]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_E[j + 12]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_E[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_E[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_E[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_E[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_E[j + 12]->SetFillColor(31);
        cMusliMap_E8->cd(j + 1);
        fh1_Muslimap_E[j + 12]->Draw("");
    }

    cMusliMap_E16 = new TCanvas("Musli_Emap16", "", 10, 10, 800, 700);
    sprintf(Name1, "fh1_Musli_Emap_mean_a1_to_a16");
    sprintf(Name2, "Eraw A01 to A16 if mult==1");
    fh1_Muslimap_E[14] = new TH1I(Name1, Name2, 64000, 0, 64000);
    fh1_Muslimap_E[14]->GetXaxis()->SetTitle("Raw Energy [channels]");
    fh1_Muslimap_E[14]->GetYaxis()->SetTitle("Counts");
    fh1_Muslimap_E[14]->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_E[14]->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_E[14]->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_E[14]->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_E[14]->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_E[14]->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_E[14]->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_E[14]->SetFillColor(31);
    cMusliMap_E16->cd();
    fh1_Muslimap_E[14]->Draw("");

    // --- Energy vs Energy --- //
    cMusliMap_EvsE2 = new TCanvas("Musli_EvsEmap2", "", 10, 10, 800, 700);
    cMusliMap_EvsE2->Divide(4, 2);
    for (Int_t j = 0; j < 7; j++)
    {
        sprintf(Name1, "fh2_Musli_Ea%da%d_vs_Ea%da%d", 2 * j + 1, 2 * j + 2, 2 * j + 3, 2 * j + 4);
        sprintf(Name2,
                "Eraw A%02d and A%02d vs Eraw A%02d and A%02d if mult==1",
                2 * j + 1,
                2 * j + 2,
                2 * j + 3,
                2 * j + 4);
        fh2_Muslimap_EvsE[j] = new TH2I(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
        fh2_Muslimap_EvsE[j]->GetXaxis()->SetTitle("Raw Energy [channel]");
        fh2_Muslimap_EvsE[j]->GetYaxis()->SetTitle("Raw Energy [channel]");
        fh2_Muslimap_EvsE[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslimap_EvsE[j]->GetXaxis()->CenterTitle(true);
        fh2_Muslimap_EvsE[j]->GetYaxis()->CenterTitle(true);
        fh2_Muslimap_EvsE[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsE[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsE[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsE[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsE[j]->SetFillColor(31);
        cMusliMap_EvsE2->cd(j + 1);
        fh2_Muslimap_EvsE[j]->Draw("COL");
    }

    cMusliMap_EvsE4 = new TCanvas("Musli_EvsEmap4", "", 10, 10, 800, 700);
    cMusliMap_EvsE4->Divide(3, 1);
    for (Int_t j = 0; j < 3; j++)
    {
        sprintf(Name1, "fh2_Musli_Ea%da%d_vs_Ea%da%d", 4 * j + 1, 4 * j + 4, 4 * j + 5, 4 * j + 8);
        sprintf(
            Name2, "Eraw A%02d to A%02d vs Eraw A%02d to A%02d if mult==1", 4 * j + 1, 4 * j + 4, 4 * j + 5, 4 * j + 8);
        fh2_Muslimap_EvsE[j + 8] = new TH2I(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
        fh2_Muslimap_EvsE[j + 8]->GetXaxis()->SetTitle("Raw Energy [channel]");
        fh2_Muslimap_EvsE[j + 8]->GetYaxis()->SetTitle("Raw Energy [channel]");
        fh2_Muslimap_EvsE[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslimap_EvsE[j + 8]->GetXaxis()->CenterTitle(true);
        fh2_Muslimap_EvsE[j + 8]->GetYaxis()->CenterTitle(true);
        fh2_Muslimap_EvsE[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsE[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsE[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsE[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsE[j + 8]->SetFillColor(31);
        cMusliMap_EvsE4->cd(j + 1);
        fh2_Muslimap_EvsE[j + 8]->Draw("COL");
    }

    cMusliMap_EvsE8 = new TCanvas("Musli_EvsEmap8", "", 10, 10, 800, 700);
    sprintf(Name1, "fh2_Musli_Ea1a8_vs_Ea9a16");
    sprintf(Name2, "Eraw A01 to A08 vs Eraw A09 to A16 if mult==1");
    fh2_Muslimap_EvsE[12] = new TH2I(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
    fh2_Muslimap_EvsE[12]->GetXaxis()->SetTitle("Raw Energy [channel]");
    fh2_Muslimap_EvsE[12]->GetYaxis()->SetTitle("Raw Energy [channel]");
    fh2_Muslimap_EvsE[12]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_EvsE[12]->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_EvsE[12]->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_EvsE[12]->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_EvsE[12]->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_EvsE[12]->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_EvsE[12]->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_EvsE[12]->SetFillColor(31);
    cMusliMap_EvsE8->cd();
    fh2_Muslimap_EvsE[12]->Draw("COL");

    // --- Drift times --- //
    cMusliMap_DT2 = new TCanvas("Musli_DTmap2", "", 10, 10, 800, 700);
    cMusliMap_DT2->Divide(4, 2);
    for (Int_t j = 0; j < 8; j++)
    {
        sprintf(Name1, "fh1_Musli_DTmap_mean_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, "DTraw A%02d and A%02d if mult==1", 2 * j + 1, 2 * j + 2);
        fh1_Muslimap_DT[j] = new TH1I(Name1, Name2, 2500, -15000, 35000);
        fh1_Muslimap_DT[j]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
        fh1_Muslimap_DT[j]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_DT[j]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_DT[j]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_DT[j]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_DT[j]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_DT[j]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_DT[j]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_DT[j]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_DT[j]->SetFillColor(31);
        cMusliMap_DT2->cd(j + 1);
        fh1_Muslimap_DT[j]->Draw("");
    }

    cMusliMap_DT4 = new TCanvas("Musli_DTmap4", "", 10, 10, 800, 700);
    cMusliMap_DT4->Divide(2, 2);
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh1_Musli_DTmap_mean_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "DTraw A%02d to A%02d if mult==1", j * 4 + 1, j * 4 + 5);
        fh1_Muslimap_DT[j + 8] = new TH1I(Name1, Name2, 2500, -15000, 35000);
        fh1_Muslimap_DT[j + 8]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
        fh1_Muslimap_DT[j + 8]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_DT[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_DT[j + 8]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_DT[j + 8]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_DT[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_DT[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_DT[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_DT[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_DT[j + 8]->SetFillColor(31);
        cMusliMap_DT4->cd(j + 1);
        fh1_Muslimap_DT[j + 8]->Draw("");
    }

    cMusliMap_DT8 = new TCanvas("Musli_DTmap8", "", 10, 10, 800, 700);
    cMusliMap_DT8->Divide(2, 1);
    for (Int_t j = 0; j < 2; j++)
    {
        sprintf(Name1, "fh1_Musli_DTmap_mean_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "DTraw A%02d to A%02d if mult==1", j * 8 + 1, j * 8 + 8);
        fh1_Muslimap_DT[j + 12] = new TH1I(Name1, Name2, 2500, -15000, 35000);
        fh1_Muslimap_DT[j + 12]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
        fh1_Muslimap_DT[j + 12]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_DT[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_DT[j + 12]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_DT[j + 12]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_DT[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_DT[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_DT[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_DT[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_DT[j + 12]->SetFillColor(31);
        cMusliMap_DT8->cd(j + 1);
        fh1_Muslimap_DT[j + 12]->Draw("");
    }

    cMusliMap_DT16 = new TCanvas("Musli_DTmap16", "", 10, 10, 800, 700);
    sprintf(Name1, "fh1_Musli_DTmap_mean_a1_to_a16");
    sprintf(Name2, "DTraw A01 to A16 if mult==1");
    fh1_Muslimap_DT[14] = new TH1I(Name1, Name2, 2500, -15000, 35000);
    fh1_Muslimap_DT[14]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
    fh1_Muslimap_DT[14]->GetYaxis()->SetTitle("Counts");
    fh1_Muslimap_DT[14]->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_DT[14]->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_DT[14]->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_DT[14]->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_DT[14]->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_DT[14]->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_DT[14]->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_DT[14]->SetFillColor(31);
    cMusliMap_DT16->cd();
    fh1_Muslimap_DT[14]->Draw("");

    // --- Drift Time vs Drift Time --- //
    cMusliMap_DTvsDT2 = new TCanvas("Musli_DTvsDTmap2", "", 10, 10, 800, 700);
    cMusliMap_DTvsDT2->Divide(4, 2);
    for (Int_t j = 0; j < 7; j++)
    {
        sprintf(Name1, "fh2_Musli_DTa%da%d_vs_DTa%da%d", 2 * j + 1, 2 * j + 2, 2 * j + 3, 2 * j + 4);
        sprintf(Name2,
                "DTraw A%02d and A%02d vs DTraw A%02d and A%02d if mult==1",
                2 * j + 1,
                2 * j + 2,
                2 * j + 3,
                2 * j + 4);
        fh2_Muslimap_DTvsDT[j] = new TH2I(Name1, Name2, 1000, -5000, 30000, 1000, -5000, 30000);
        fh2_Muslimap_DTvsDT[j]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_Muslimap_DTvsDT[j]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_Muslimap_DTvsDT[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslimap_DTvsDT[j]->GetXaxis()->CenterTitle(true);
        fh2_Muslimap_DTvsDT[j]->GetYaxis()->CenterTitle(true);
        fh2_Muslimap_DTvsDT[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslimap_DTvsDT[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslimap_DTvsDT[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslimap_DTvsDT[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslimap_DTvsDT[j]->SetFillColor(31);
        cMusliMap_DTvsDT2->cd(j + 1);
        fh2_Muslimap_DTvsDT[j]->Draw("COL");
    }

    cMusliMap_DTvsDT4 = new TCanvas("Musli_DTvsDTmap4", "", 10, 10, 800, 700);
    cMusliMap_DTvsDT4->Divide(3, 1);
    for (Int_t j = 0; j < 3; j++)
    {
        sprintf(Name1, "fh2_Musli_DTa%da%d_vs_DTa%da%d", 4 * j + 1, 4 * j + 4, 4 * j + 5, 4 * j + 8);
        sprintf(Name2,
                "DTraw A%02d to A%02d vs DTraw A%02d to A%02d if mult==1",
                4 * j + 1,
                4 * j + 4,
                4 * j + 5,
                4 * j + 8);
        fh2_Muslimap_DTvsDT[j + 8] = new TH2I(Name1, Name2, 1000, -5000, 30000, 1000, -5000, 30000);
        fh2_Muslimap_DTvsDT[j + 8]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_Muslimap_DTvsDT[j + 8]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_Muslimap_DTvsDT[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslimap_DTvsDT[j + 8]->GetXaxis()->CenterTitle(true);
        fh2_Muslimap_DTvsDT[j + 8]->GetYaxis()->CenterTitle(true);
        fh2_Muslimap_DTvsDT[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslimap_DTvsDT[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslimap_DTvsDT[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslimap_DTvsDT[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslimap_DTvsDT[j + 8]->SetFillColor(31);
        cMusliMap_DTvsDT4->cd(j + 1);
        fh2_Muslimap_DTvsDT[j + 8]->Draw("COL");
    }

    cMusliMap_DTvsDT8 = new TCanvas("Musli_DTvsDTmap8", "", 10, 10, 800, 700);
    sprintf(Name1, "fh2_Musli_DTa1a8_vs_DTa9a16");
    sprintf(Name2, "DTraw A01 to A08 vs DTraw A09 to A16 if mult==1");
    fh2_Muslimap_DTvsDT[12] = new TH2I(Name1, Name2, 1000, -5000, 30000, 1000, -5000, 30000);
    fh2_Muslimap_DTvsDT[12]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
    fh2_Muslimap_DTvsDT[12]->GetYaxis()->SetTitle("Raw Drift Time [100ps]");
    fh2_Muslimap_DTvsDT[12]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_DTvsDT[12]->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_DTvsDT[12]->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_DTvsDT[12]->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_DTvsDT[12]->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_DTvsDT[12]->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_DTvsDT[12]->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_DTvsDT[12]->SetFillColor(31);
    cMusliMap_DTvsDT8->cd();
    fh2_Muslimap_DTvsDT[12]->Draw("COL");

    // --- Energy vs Drift time --- //
    cMusliMap_EvsDT2 = new TCanvas("Musli_EvsDTmap2", "", 10, 10, 800, 700);
    cMusliMap_EvsDT2->Divide(4, 2);
    for (Int_t j = 0; j < 8; j++)
    {
        sprintf(Name1, "fh2_Musli_EvsDTmap_mean_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, "Eraw vs DTraw A%02d and A%02d if mult==1", 2 * j + 1, 2 * j + 2);
        fh2_Muslimap_EvsDT[j] = new TH2I(Name1, Name2, 1000, -5000, 35000, 1280, 0, 64000);
        fh2_Muslimap_EvsDT[j]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_Muslimap_EvsDT[j]->GetYaxis()->SetTitle("Raw Energy [channel]");
        fh2_Muslimap_EvsDT[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslimap_EvsDT[j]->GetXaxis()->CenterTitle(true);
        fh2_Muslimap_EvsDT[j]->GetYaxis()->CenterTitle(true);
        fh2_Muslimap_EvsDT[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsDT[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsDT[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsDT[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsDT[j]->SetFillColor(31);
        cMusliMap_EvsDT2->cd(j + 1);
        fh2_Muslimap_EvsDT[j]->Draw("COL");
    }

    cMusliMap_EvsDT4 = new TCanvas("Musli_EvsDTmap4", "", 10, 10, 800, 700);
    cMusliMap_EvsDT4->Divide(2, 2);
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh2_Musli_EvsDTmap_mean_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "Eraw vs DTraw A%02d to A%02d if mult==1", j * 4 + 1, j * 4 + 5);
        fh2_Muslimap_EvsDT[j + 8] = new TH2I(Name1, Name2, 1000, -5000, 35000, 1280, 0, 64000);
        fh2_Muslimap_EvsDT[j + 8]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_Muslimap_EvsDT[j + 8]->GetYaxis()->SetTitle("Raw Energy [channel]");
        fh2_Muslimap_EvsDT[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslimap_EvsDT[j + 8]->GetXaxis()->CenterTitle(true);
        fh2_Muslimap_EvsDT[j + 8]->GetYaxis()->CenterTitle(true);
        fh2_Muslimap_EvsDT[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsDT[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsDT[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsDT[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsDT[j + 8]->SetFillColor(31);
        cMusliMap_EvsDT4->cd(j + 1);
        fh2_Muslimap_EvsDT[j + 8]->Draw("COL");
    }

    cMusliMap_EvsDT8 = new TCanvas("Musli_EvsDTmap8", "", 10, 10, 800, 700);
    cMusliMap_EvsDT8->Divide(2, 1);
    for (Int_t j = 0; j < 2; j++)
    {
        sprintf(Name1, "fh2_Musli_EvsDTmap_mean_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "Eraw vs DTraw A%02d to A%02d if mult==1", j * 8 + 1, j * 8 + 8);
        fh2_Muslimap_EvsDT[j + 12] = new TH2I(Name1, Name2, 1000, -5000, 35000, 1280, 0, 64000);
        fh2_Muslimap_EvsDT[j + 12]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
        fh2_Muslimap_EvsDT[j + 12]->GetYaxis()->SetTitle("Raw Energy [channel]");
        fh2_Muslimap_EvsDT[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslimap_EvsDT[j + 12]->GetXaxis()->CenterTitle(true);
        fh2_Muslimap_EvsDT[j + 12]->GetYaxis()->CenterTitle(true);
        fh2_Muslimap_EvsDT[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsDT[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsDT[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslimap_EvsDT[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslimap_EvsDT[j + 12]->SetFillColor(31);
        cMusliMap_EvsDT8->cd(j + 1);
        fh2_Muslimap_EvsDT[j + 12]->Draw("COL");
    }

    cMusliMap_EvsDT16 = new TCanvas("Musli_EvsDTmap16", "", 10, 10, 800, 700);
    sprintf(Name1, "fh2_Musli_EvsDTmap_mean_a1_to_a16");
    sprintf(Name2, "Eraw vs DTraw A01 to A16 if mult==1");
    fh2_Muslimap_EvsDT[14] = new TH2I(Name1, Name2, 1000, -5000, 35000, 1280, 0, 64000);
    fh2_Muslimap_EvsDT[14]->GetXaxis()->SetTitle("Raw Drift Time [100ps]");
    fh2_Muslimap_EvsDT[14]->GetYaxis()->SetTitle("Raw Energy [channel]");
    fh2_Muslimap_EvsDT[14]->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_EvsDT[14]->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_EvsDT[14]->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_EvsDT[14]->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_EvsDT[14]->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_EvsDT[14]->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_EvsDT[14]->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_EvsDT[14]->SetFillColor(31);
    cMusliMap_EvsDT16->cd();
    fh2_Muslimap_EvsDT[14]->Draw("COL");

    // --- Mean Energy --- //

    cMusliMap_Emean = new TCanvas("MusliMap_Emean", "", 10, 10, 800, 700);
    cMusliMap_Emean->Divide(3, 3);

    sprintf(Name1, "fh1_MusliMap_Emean_from_Emap2");
    sprintf(Name2, "Mean Energy from Emap2");
    fh1_Muslimap_Emean2 = new TH1F(Name1, Name2, 64000, 0, 64000);
    fh1_Muslimap_Emean2->GetXaxis()->SetTitle("Energy [channels]");
    fh1_Muslimap_Emean2->GetYaxis()->SetTitle("Counts");
    fh1_Muslimap_Emean2->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_Emean2->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_Emean2->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_Emean2->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_Emean2->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_Emean2->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_Emean2->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_Emean2->SetFillColor(31);
    cMusliMap_Emean->cd(1);
    fh1_Muslimap_Emean2->Draw("");

    sprintf(Name1, "fh1_MusliMap_Emean_from_Emap4");
    sprintf(Name2, "Mean Energy from Emap4");
    fh1_Muslimap_Emean4 = new TH1F(Name1, Name2, 64000, 0, 64000);
    fh1_Muslimap_Emean4->GetXaxis()->SetTitle("Energy [channels]");
    fh1_Muslimap_Emean4->GetYaxis()->SetTitle("Counts");
    fh1_Muslimap_Emean4->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_Emean4->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_Emean4->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_Emean4->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_Emean4->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_Emean4->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_Emean4->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_Emean4->SetFillColor(31);
    cMusliMap_Emean->cd(2);
    fh1_Muslimap_Emean4->Draw("");

    sprintf(Name1, "fh1_MusliMap_Emean_from_Emap8");
    sprintf(Name2, "Mean Energy from Emap8");
    fh1_Muslimap_Emean8 = new TH1F(Name1, Name2, 64000, 0, 64000);
    fh1_Muslimap_Emean8->GetXaxis()->SetTitle("Energy [channels]");
    fh1_Muslimap_Emean8->GetYaxis()->SetTitle("Counts");
    fh1_Muslimap_Emean8->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_Emean8->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_Emean8->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_Emean8->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_Emean8->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_Emean8->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_Emean8->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_Emean8->SetFillColor(31);
    cMusliMap_Emean->cd(3);
    fh1_Muslimap_Emean8->Draw("");

    sprintf(Name1, "fh2_MusliMap_Emean2_vs_Emean4");
    sprintf(Name2, "Emean2 vs Emean4");
    fh2_Muslimap_Emean2vsEmean4 = new TH2F(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
    fh2_Muslimap_Emean2vsEmean4->GetXaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean2vsEmean4->GetYaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean2vsEmean4->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_Emean2vsEmean4->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_Emean2vsEmean4->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_Emean2vsEmean4->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean2vsEmean4->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean2vsEmean4->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean2vsEmean4->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean2vsEmean4->SetFillColor(31);
    cMusliMap_Emean->cd(4);
    fh2_Muslimap_Emean2vsEmean4->Draw("COL");

    sprintf(Name1, "fh2_MusliMap_Emean2_vs_Emean8");
    sprintf(Name2, "Emean2 vs Emean8");
    fh2_Muslimap_Emean2vsEmean8 = new TH2F(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
    fh2_Muslimap_Emean2vsEmean8->GetXaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean2vsEmean8->GetYaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean2vsEmean8->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_Emean2vsEmean8->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_Emean2vsEmean8->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_Emean2vsEmean8->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean2vsEmean8->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean2vsEmean8->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean2vsEmean8->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean2vsEmean8->SetFillColor(31);
    cMusliMap_Emean->cd(5);
    fh2_Muslimap_Emean2vsEmean8->Draw("COL");

    sprintf(Name1, "fh2_MusliMap_Emean2_vs_E16");
    sprintf(Name2, "Emean2 vs E16");
    fh2_Muslimap_Emean2vsE16 = new TH2F(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
    fh2_Muslimap_Emean2vsE16->GetXaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean2vsE16->GetYaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean2vsE16->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_Emean2vsE16->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_Emean2vsE16->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_Emean2vsE16->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean2vsE16->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean2vsE16->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean2vsE16->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean2vsE16->SetFillColor(31);
    cMusliMap_Emean->cd(6);
    fh2_Muslimap_Emean2vsE16->Draw("COL");

    sprintf(Name1, "fh2_MusliMap_Emean4_vs_Emean8");
    sprintf(Name2, "Emean4 vs Emean8");
    fh2_Muslimap_Emean4vsEmean8 = new TH2F(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
    fh2_Muslimap_Emean4vsEmean8->GetXaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean4vsEmean8->GetYaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean4vsEmean8->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_Emean4vsEmean8->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_Emean4vsEmean8->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_Emean4vsEmean8->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean4vsEmean8->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean4vsEmean8->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean4vsEmean8->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean4vsEmean8->SetFillColor(31);
    cMusliMap_Emean->cd(7);
    fh2_Muslimap_Emean4vsEmean8->Draw("COL");

    sprintf(Name1, "fh2_MusliMap_Emean4_vs_E16");
    sprintf(Name2, "Emean4 vs E16");
    fh2_Muslimap_Emean4vsE16 = new TH2F(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
    fh2_Muslimap_Emean4vsE16->GetXaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean4vsE16->GetYaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean4vsE16->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_Emean4vsE16->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_Emean4vsE16->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_Emean4vsE16->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean4vsE16->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean4vsE16->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean4vsE16->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean4vsE16->SetFillColor(31);
    cMusliMap_Emean->cd(8);
    fh2_Muslimap_Emean4vsE16->Draw("COL");

    sprintf(Name1, "fh2_MusliMap_Emean8_vs_E16");
    sprintf(Name2, "Emean8 vs E16");
    fh2_Muslimap_Emean8vsE16 = new TH2F(Name1, Name2, 1280, 0, 64000, 1280, 0, 64000);
    fh2_Muslimap_Emean8vsE16->GetXaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean8vsE16->GetYaxis()->SetTitle("Energy [channels]");
    fh2_Muslimap_Emean8vsE16->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslimap_Emean8vsE16->GetXaxis()->CenterTitle(true);
    fh2_Muslimap_Emean8vsE16->GetYaxis()->CenterTitle(true);
    fh2_Muslimap_Emean8vsE16->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean8vsE16->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean8vsE16->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslimap_Emean8vsE16->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslimap_Emean8vsE16->SetFillColor(31);
    cMusliMap_Emean->cd(9);
    fh2_Muslimap_Emean8vsE16->Draw("COL");

    // --- DeltaDT --- //

    cMusliMap_DeltaDT = new TCanvas("MusliMap_DeltaDT", "", 10, 10, 800, 700);

    sprintf(Name1, "fh1_MusliMap_DeltaDT");
    sprintf(Name2, "Delta DTmap2 between first and last pair ");
    fh1_Muslimap_DeltaDT = new TH1I(Name1, Name2, 3000, -1500, 1500);
    fh1_Muslimap_DeltaDT->GetXaxis()->SetTitle("Delta Drift Time [100ps]");
    fh1_Muslimap_DeltaDT->GetYaxis()->SetTitle("Counts");
    fh1_Muslimap_DeltaDT->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_DeltaDT->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_DeltaDT->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_DeltaDT->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_DeltaDT->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_DeltaDT->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_DeltaDT->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_DeltaDT->SetFillColor(31);
    cMusliMap_DeltaDT->cd();
    fh1_Muslimap_DeltaDT->Draw("");

    // --- T[ch] - T[17]=Ttrig --- //
    cMusliMap_TchTtrig = new TCanvas("MusliMap_TchTtrig", "", 10, 10, 800, 700);
    cMusliMap_TchTtrig->Divide(4, 4);
    for (Int_t j = 0; j < 8; j++)
    {
        sprintf(Name1, "fh1_MusliMap_Traw_a%d_a%d_Ttrig", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, "DeltaT with respect to Ttrig, A%02d and A%02d if mult==1", 2 * j + 1, 2 * j + 2);
        fh1_Muslimap_TchTtrig[j] = new TH1D(Name1, Name2, 20000, -1000000, 1000000);
        fh1_Muslimap_TchTtrig[j]->GetXaxis()->SetTitle("Delta T [100ps]");
        fh1_Muslimap_TchTtrig[j]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_TchTtrig[j]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_TchTtrig[j]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_TchTtrig[j]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_TchTtrig[j]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_TchTtrig[j]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_TchTtrig[j]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_TchTtrig[j]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_TchTtrig[j]->SetFillColor(31);
        cMusliMap_TchTtrig->cd(j + 1);
        fh1_Muslimap_TchTtrig[j]->Draw("");
    }
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh1_MusliMap_Traw_a%d_to_a%d_Ttrig", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "DeltaT with respect to Ttrig, A%02d to A%02d if mult==1", j * 4 + 1, j * 4 + 5);
        fh1_Muslimap_TchTtrig[j + 8] = new TH1D(Name1, Name2, 20000, -1000000, 1000000);
        fh1_Muslimap_TchTtrig[j + 8]->GetXaxis()->SetTitle("Delta T [100ps]");
        fh1_Muslimap_TchTtrig[j + 8]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_TchTtrig[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_TchTtrig[j + 8]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_TchTtrig[j + 8]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_TchTtrig[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_TchTtrig[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_TchTtrig[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_TchTtrig[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_TchTtrig[j + 8]->SetFillColor(31);
        cMusliMap_TchTtrig->cd(j + 9);
        fh1_Muslimap_TchTtrig[j + 8]->Draw("");
    }
    for (Int_t j = 0; j < 2; j++)
    {
        sprintf(Name1, "fh1_MusliMap_Traw_a%d_to_a%d_Ttrig", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "DeltaT with respect to Ttrig, A%02d to A%02d if mult==1", j * 8 + 1, j * 8 + 8);
        fh1_Muslimap_TchTtrig[j + 12] = new TH1D(Name1, Name2, 20000, -1000000, 1000000);
        fh1_Muslimap_TchTtrig[j + 12]->GetXaxis()->SetTitle("Delta T [100ps]");
        fh1_Muslimap_TchTtrig[j + 12]->GetYaxis()->SetTitle("Counts");
        fh1_Muslimap_TchTtrig[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslimap_TchTtrig[j + 12]->GetXaxis()->CenterTitle(true);
        fh1_Muslimap_TchTtrig[j + 12]->GetYaxis()->CenterTitle(true);
        fh1_Muslimap_TchTtrig[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslimap_TchTtrig[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslimap_TchTtrig[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslimap_TchTtrig[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslimap_TchTtrig[j + 12]->SetFillColor(31);
        cMusliMap_TchTtrig->cd(j + 13);
        fh1_Muslimap_TchTtrig[j + 12]->Draw("");
    }
    sprintf(Name1, "fh1_MusliMap_Traw_a1_to_a16_Ttrig");
    sprintf(Name2, "DeltaT with respect to Ttrig, A01 to A16 if mult==1");
    fh1_Muslimap_TchTtrig[14] = new TH1D(Name1, Name2, 20000, -1000000, 1000000);
    fh1_Muslimap_TchTtrig[14]->GetXaxis()->SetTitle("Delta T [100ps]");
    fh1_Muslimap_TchTtrig[14]->GetYaxis()->SetTitle("Counts");
    fh1_Muslimap_TchTtrig[14]->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_TchTtrig[14]->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_TchTtrig[14]->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_TchTtrig[14]->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_TchTtrig[14]->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_TchTtrig[14]->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_TchTtrig[14]->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_TchTtrig[14]->SetFillColor(31);
    cMusliMap_TchTtrig->cd(15);
    fh1_Muslimap_TchTtrig[14]->Draw("");

    sprintf(Name1, "fh1_MusliMap_Tref_Ttrig");
    sprintf(Name2, "DeltaT = Tref - Ttrig if mult==1");
    fh1_Muslimap_TchTtrig[15] = new TH1D(Name1, Name2, 20000, -1000000, 1000000);
    fh1_Muslimap_TchTtrig[15]->GetXaxis()->SetTitle("Delta T [100ps]");
    fh1_Muslimap_TchTtrig[15]->GetYaxis()->SetTitle("Counts");
    fh1_Muslimap_TchTtrig[15]->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslimap_TchTtrig[15]->GetXaxis()->CenterTitle(true);
    fh1_Muslimap_TchTtrig[15]->GetYaxis()->CenterTitle(true);
    fh1_Muslimap_TchTtrig[15]->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslimap_TchTtrig[15]->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslimap_TchTtrig[15]->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslimap_TchTtrig[15]->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslimap_TchTtrig[15]->SetFillColor(31);
    cMusliMap_TchTtrig->cd(16);
    fh1_Muslimap_TchTtrig[15]->Draw("");

    // === ======== === //
    // === CAL DATA === //
    // === ======== === //

    // --- Multiplicities --- //
    cMusliCal_Mult = new TCanvas("MusliCal_Mult", "", 10, 10, 800, 700);
    cMusliCal_Mult->Divide(2, 1);

    fh1_Muslical_mult = new TH1I("MultCalTot", "Total Multiplicity at CAL level", 17, -0.5, 16.5);
    fh1_Muslical_mult->GetXaxis()->SetTitle("fSignal");
    fh1_Muslical_mult->GetYaxis()->SetTitle("Total Counts");
    fh1_Muslical_mult->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslical_mult->GetXaxis()->CenterTitle(true);
    fh1_Muslical_mult->GetYaxis()->CenterTitle(true);
    fh1_Muslical_mult->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslical_mult->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslical_mult->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslical_mult->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslical_mult->SetFillColor(31);
    cMusliCal_Mult->cd(1);
    fh1_Muslical_mult->Draw("");

    fh2_Muslical_mult =
        new TH2I("MultCalDist", "Multiplicity Distribution per fSignal at CAL level", 17, -0.5, 16.5, 50, -0.5, 49.5);
    fh2_Muslical_mult->GetXaxis()->SetTitle("fSignal");
    fh2_Muslical_mult->GetYaxis()->SetTitle("Counts per Event");
    fh2_Muslical_mult->GetYaxis()->SetTitleOffset(1.1);
    fh2_Muslical_mult->GetXaxis()->CenterTitle(true);
    fh2_Muslical_mult->GetYaxis()->CenterTitle(true);
    fh2_Muslical_mult->GetXaxis()->SetLabelSize(0.045);
    fh2_Muslical_mult->GetXaxis()->SetTitleSize(0.045);
    fh2_Muslical_mult->GetYaxis()->SetLabelSize(0.045);
    fh2_Muslical_mult->GetYaxis()->SetTitleSize(0.045);
    fh2_Muslical_mult->SetFillColor(31);
    cMusliCal_Mult->cd(2);
    fh2_Muslical_mult->Draw("COL");

    // --- Energies --- //
    cMusliCal_E2 = new TCanvas("Musli_Ecal2", "", 10, 10, 800, 700);
    cMusliCal_E2->Divide(4, 2);
    for (Int_t j = 0; j < 8; j++)
    {
        sprintf(Name1, "fh1_Musli_Ecal_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, " Ecal A%02d and A%02d if mult==1", 2 * j + 1, 2 * j + 2);
        fh1_Muslical_E[j] = new TH1D(Name1, Name2, 64000, 0, 64000);
        fh1_Muslical_E[j]->GetXaxis()->SetTitle("Energy [channels]");
        fh1_Muslical_E[j]->GetYaxis()->SetTitle("Counts");
        fh1_Muslical_E[j]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslical_E[j]->GetXaxis()->CenterTitle(true);
        fh1_Muslical_E[j]->GetYaxis()->CenterTitle(true);
        fh1_Muslical_E[j]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslical_E[j]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslical_E[j]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslical_E[j]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslical_E[j]->SetFillColor(31);
        cMusliCal_E2->cd(j + 1);
        fh1_Muslical_E[j]->Draw("");
    }

    cMusliCal_E4 = new TCanvas("Musli_Ecal4", "", 10, 10, 800, 700);
    cMusliCal_E4->Divide(2, 2);
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh1_Musli_Ecal_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "Ecal A%02d to A%02d if mult==1", j * 4 + 1, j * 4 + 5);
        fh1_Muslical_E[j + 8] = new TH1D(Name1, Name2, 64000, 0, 64000);
        fh1_Muslical_E[j + 8]->GetXaxis()->SetTitle("Energy [channels]");
        fh1_Muslical_E[j + 8]->GetYaxis()->SetTitle("Counts");
        fh1_Muslical_E[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslical_E[j + 8]->GetXaxis()->CenterTitle(true);
        fh1_Muslical_E[j + 8]->GetYaxis()->CenterTitle(true);
        fh1_Muslical_E[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslical_E[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslical_E[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslical_E[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslical_E[j + 8]->SetFillColor(31);
        cMusliCal_E4->cd(j + 1);
        fh1_Muslical_E[j + 8]->Draw("");
    }

    cMusliCal_E8 = new TCanvas("Musli_Ecal8", "", 10, 10, 800, 700);
    cMusliCal_E8->Divide(2, 1);
    for (Int_t j = 0; j < 2; j++)
    {
        sprintf(Name1, "fh1_Musli_Ecal_mean_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "Ecal A%02d to A%02d if mult==1", j * 8 + 1, j * 8 + 8);
        fh1_Muslical_E[j + 12] = new TH1D(Name1, Name2, 64000, 0, 64000);
        fh1_Muslical_E[j + 12]->GetXaxis()->SetTitle("Energy [channels]");
        fh1_Muslical_E[j + 12]->GetYaxis()->SetTitle("Counts");
        fh1_Muslical_E[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslical_E[j + 12]->GetXaxis()->CenterTitle(true);
        fh1_Muslical_E[j + 12]->GetYaxis()->CenterTitle(true);
        fh1_Muslical_E[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslical_E[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslical_E[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslical_E[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslical_E[j + 12]->SetFillColor(31);
        cMusliCal_E8->cd(j + 1);
        fh1_Muslical_E[j + 12]->Draw("");
    }

    cMusliCal_E16 = new TCanvas("Musli_Ecal16", "", 10, 10, 800, 700);
    sprintf(Name1, "fh1_Musli_Ecal_mean_a1_to_a16");
    sprintf(Name2, "Ecal A01 to A16 if mult==1");
    fh1_Muslical_E[14] = new TH1D(Name1, Name2, 64000, 0, 64000);
    fh1_Muslical_E[14]->GetXaxis()->SetTitle("Energy [channels]");
    fh1_Muslical_E[14]->GetYaxis()->SetTitle("Counts");
    fh1_Muslical_E[14]->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslical_E[14]->GetXaxis()->CenterTitle(true);
    fh1_Muslical_E[14]->GetYaxis()->CenterTitle(true);
    fh1_Muslical_E[14]->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslical_E[14]->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslical_E[14]->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslical_E[14]->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslical_E[14]->SetFillColor(31);
    cMusliCal_E16->cd();
    fh1_Muslical_E[14]->Draw("");

    // --- Drift times --- //

    cMusliCal_DT2 = new TCanvas("Musli_DTcal2", "", 10, 10, 800, 700);
    cMusliCal_DT2->Divide(4, 2);
    for (Int_t j = 0; j < 8; j++)
    {
        sprintf(Name1, "fh1_Musli_DTcal_mean_a%d_a%d", 2 * j + 1, 2 * j + 2);
        sprintf(Name2, "DTcal A%02d and A%02d if mult==1", 2 * j + 1, 2 * j + 2);
        fh1_Muslical_DT[j] = new TH1D(Name1, Name2, 2000, -100, 100);
        fh1_Muslical_DT[j]->GetXaxis()->SetTitle("Calibrated Drift Time in [mm]");
        fh1_Muslical_DT[j]->GetYaxis()->SetTitle("Counts");
        fh1_Muslical_DT[j]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslical_DT[j]->GetXaxis()->CenterTitle(true);
        fh1_Muslical_DT[j]->GetYaxis()->CenterTitle(true);
        fh1_Muslical_DT[j]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslical_DT[j]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslical_DT[j]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslical_DT[j]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslical_DT[j]->SetFillColor(31);
        cMusliCal_DT2->cd(j + 1);
        fh1_Muslical_DT[j]->Draw("");
    }

    cMusliCal_DT4 = new TCanvas("Musli_DTcal4", "", 10, 10, 800, 700);
    cMusliCal_DT4->Divide(2, 2);
    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh1_Musli_DTcal_mean_a%d_to_a%d", j * 4 + 1, j * 4 + 4);
        sprintf(Name2, "DTcal A%02d to A%02d if mult==1", j * 4 + 1, j * 4 + 5);
        fh1_Muslical_DT[j + 8] = new TH1D(Name1, Name2, 2000, -100, 100);
        fh1_Muslical_DT[j + 8]->GetXaxis()->SetTitle("Calibrated Drift Time in [mm]");
        fh1_Muslical_DT[j + 8]->GetYaxis()->SetTitle("Counts");
        fh1_Muslical_DT[j + 8]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslical_DT[j + 8]->GetXaxis()->CenterTitle(true);
        fh1_Muslical_DT[j + 8]->GetYaxis()->CenterTitle(true);
        fh1_Muslical_DT[j + 8]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslical_DT[j + 8]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslical_DT[j + 8]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslical_DT[j + 8]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslical_DT[j + 8]->SetFillColor(31);
        cMusliCal_DT4->cd(j + 1);
        fh1_Muslical_DT[j + 8]->Draw("");
    }

    cMusliCal_DT8 = new TCanvas("Musli_DTcal8", "", 10, 10, 800, 700);
    cMusliCal_DT8->Divide(2, 1);
    for (Int_t j = 0; j < 2; j++)
    {
        sprintf(Name1, "fh1_Musli_DTcal_mean_a%d_to_a%d", j * 8 + 1, j * 8 + 8);
        sprintf(Name2, "DTcal A%02d to A%02d if mult==1", j * 8 + 1, j * 8 + 8);
        fh1_Muslical_DT[j + 12] = new TH1D(Name1, Name2, 2000, -100, 100);
        fh1_Muslical_DT[j + 12]->GetXaxis()->SetTitle("Calibrated Drift Time in [mm]");
        fh1_Muslical_DT[j + 12]->GetYaxis()->SetTitle("Counts");
        fh1_Muslical_DT[j + 12]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslical_DT[j + 12]->GetXaxis()->CenterTitle(true);
        fh1_Muslical_DT[j + 12]->GetYaxis()->CenterTitle(true);
        fh1_Muslical_DT[j + 12]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslical_DT[j + 12]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslical_DT[j + 12]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslical_DT[j + 12]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslical_DT[j + 12]->SetFillColor(31);
        cMusliCal_DT8->cd(j + 1);
        fh1_Muslical_DT[j + 12]->Draw("");
    }

    cMusliCal_DT16 = new TCanvas("Musli_DTcal16", "", 10, 10, 800, 700);
    sprintf(Name1, "fh1_Musli_DTcal_mean_a1_to_a16");
    sprintf(Name2, "DTcal A01 to A16 if mult==1");
    fh1_Muslical_DT[14] = new TH1D(Name1, Name2, 5000, 10000, 35000);
    fh1_Muslical_DT[14]->GetXaxis()->SetTitle("Calibrated Drift Time in [mm]");
    fh1_Muslical_DT[14]->GetYaxis()->SetTitle("Counts");
    fh1_Muslical_DT[14]->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslical_DT[14]->GetXaxis()->CenterTitle(true);
    fh1_Muslical_DT[14]->GetYaxis()->CenterTitle(true);
    fh1_Muslical_DT[14]->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslical_DT[14]->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslical_DT[14]->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslical_DT[14]->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslical_DT[14]->SetFillColor(31);
    cMusliCal_DT16->cd();
    fh1_Muslical_DT[14]->Draw("");

    // === ======== === //
    // === HIT DATA === //
    // === ======== === //
    Int_t nb_anodes[4] = { 2, 4, 8, 16 };

    cMusliHit_E = new TCanvas("Musli_Ehit", "", 10, 10, 800, 700);
    cMusliHit_E->Divide(2, 2);

    cMusliHit_Z = new TCanvas("Musli_Zcharge", "", 10, 10, 800, 700);
    cMusliHit_Z->Divide(2, 2);

    cMusliHit_EvsX = new TCanvas("Musli_EhitVsX", "", 10, 10, 800, 700);
    cMusliHit_EvsX->Divide(2, 2);

    cMusliHit_ZvsX = new TCanvas("Musli_ZchargeVsX", "", 10, 10, 800, 700);
    cMusliHit_ZvsX->Divide(2, 2);

    cMusliHit_EvsTheta = new TCanvas("Musli_EhitVsTheta", "", 10, 10, 800, 700);
    cMusliHit_EvsTheta->Divide(2, 2);

    cMusliHit_ZvsTheta = new TCanvas("Musli_ZchargeVsTheta", "", 10, 10, 800, 700);
    cMusliHit_ZvsTheta->Divide(2, 2);

    for (Int_t j = 0; j < 4; j++)
    {
        sprintf(Name1, "fh1_Musli_Ehit_%i", nb_anodes[j]);
        sprintf(Name2, " Eaverage for data with nb anodes %02d  if mult==1", nb_anodes[j]);
        fh1_Muslihit_E[j] = new TH1D(Name1, Name2, 64000, 0, 64000);
        fh1_Muslihit_E[j]->GetXaxis()->SetTitle("Energy [channels]");
        fh1_Muslihit_E[j]->GetYaxis()->SetTitle("Counts");
        fh1_Muslihit_E[j]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslihit_E[j]->GetXaxis()->CenterTitle(true);
        fh1_Muslihit_E[j]->GetYaxis()->CenterTitle(true);
        fh1_Muslihit_E[j]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslihit_E[j]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslihit_E[j]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslihit_E[j]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslihit_E[j]->SetFillColor(31);
        cMusliHit_E->cd(j + 1);
        fh1_Muslihit_E[j]->Draw("");

        sprintf(Name1, "fh1_Musli_Zcharge_%i", nb_anodes[j]);
        sprintf(Name2, " Z for data with nb anodes %02d  if mult==1", nb_anodes[j]);
        fh1_Muslihit_Z[j] = new TH1D(Name1, Name2, 3000, 0, 15);
        fh1_Muslihit_Z[j]->GetXaxis()->SetTitle("Charge  [atomic number]");
        fh1_Muslihit_Z[j]->GetYaxis()->SetTitle("Counts");
        fh1_Muslihit_Z[j]->GetYaxis()->SetTitleOffset(1.1);
        fh1_Muslihit_Z[j]->GetXaxis()->CenterTitle(true);
        fh1_Muslihit_Z[j]->GetYaxis()->CenterTitle(true);
        fh1_Muslihit_Z[j]->GetXaxis()->SetLabelSize(0.045);
        fh1_Muslihit_Z[j]->GetXaxis()->SetTitleSize(0.045);
        fh1_Muslihit_Z[j]->GetYaxis()->SetLabelSize(0.045);
        fh1_Muslihit_Z[j]->GetYaxis()->SetTitleSize(0.045);
        fh1_Muslihit_Z[j]->SetFillColor(31);
        cMusliHit_Z->cd(j + 1);
        fh1_Muslihit_Z[j]->Draw("");

        sprintf(Name1, "fh2_Musli_EhitVsX%i", nb_anodes[j]);
        sprintf(Name2, " Eaverage for data with nb anodes %02d  if mult==1", nb_anodes[j]);
        fh2_Muslihit_EvsX[j] = new TH2D(Name1, Name2, 1000, -50, 50, 1600, 0, 64000);
        fh2_Muslihit_EvsX[j]->GetXaxis()->SetTitle("X [mm]");
        fh2_Muslihit_EvsX[j]->GetYaxis()->SetTitle("Energy [channels]");
        fh2_Muslihit_EvsX[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslihit_EvsX[j]->GetXaxis()->CenterTitle(true);
        fh2_Muslihit_EvsX[j]->GetYaxis()->CenterTitle(true);
        fh2_Muslihit_EvsX[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslihit_EvsX[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslihit_EvsX[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslihit_EvsX[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslihit_EvsX[j]->SetFillColor(31);
        cMusliHit_EvsX->cd(j + 1);
        fh2_Muslihit_EvsX[j]->Draw("COL");

        sprintf(Name1, "fh2_Musli_Zcharge_%i_vs_X", nb_anodes[j]);
        sprintf(Name2, " Z for data with nb anodes %02d  if mult==1", nb_anodes[j]);
        fh2_Muslihit_ZvsX[j] = new TH2D(Name1, Name2, 1000, -50, 50, 750, 0, 15);
        fh2_Muslihit_ZvsX[j]->GetXaxis()->SetTitle("X [mm]");
        fh2_Muslihit_ZvsX[j]->GetYaxis()->SetTitle("Charge  [atomic number]");
        fh2_Muslihit_ZvsX[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslihit_ZvsX[j]->GetXaxis()->CenterTitle(true);
        fh2_Muslihit_ZvsX[j]->GetYaxis()->CenterTitle(true);
        fh2_Muslihit_ZvsX[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslihit_ZvsX[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslihit_ZvsX[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslihit_ZvsX[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslihit_ZvsX[j]->SetFillColor(31);
        cMusliHit_ZvsX->cd(j + 1);
        fh2_Muslihit_ZvsX[j]->Draw("COL");

        sprintf(Name1, "fh2_Musli_EhitVsTheta%i", nb_anodes[j]);
        sprintf(Name2, " Eaverage for data with nb anodes %02d  if mult==1", nb_anodes[j]);
        fh2_Muslihit_EvsTheta[j] = new TH2D(Name1, Name2, 1000, -0.05, 0.05, 1600, 0, 64000);
        fh2_Muslihit_EvsTheta[j]->GetXaxis()->SetTitle("Theta [mm]");
        fh2_Muslihit_EvsTheta[j]->GetYaxis()->SetTitle("Energy [channels]");
        fh2_Muslihit_EvsTheta[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslihit_EvsTheta[j]->GetXaxis()->CenterTitle(true);
        fh2_Muslihit_EvsTheta[j]->GetYaxis()->CenterTitle(true);
        fh2_Muslihit_EvsTheta[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslihit_EvsTheta[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslihit_EvsTheta[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslihit_EvsTheta[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslihit_EvsTheta[j]->SetFillColor(31);
        cMusliHit_EvsTheta->cd(j + 1);
        fh2_Muslihit_EvsTheta[j]->Draw("COL");

        sprintf(Name1, "fh2_Musli_Zcharge_%i_vs_Theta", nb_anodes[j]);
        sprintf(Name2, " Z for data with nb anodes %02d  if mult==1", nb_anodes[j]);
        fh2_Muslihit_ZvsTheta[j] = new TH2D(Name1, Name2, 1000, -0.05, 0.05, 750, 0, 15);
        fh2_Muslihit_ZvsTheta[j]->GetXaxis()->SetTitle("Theta [mrad]");
        fh2_Muslihit_ZvsTheta[j]->GetYaxis()->SetTitle("Charge  [atomic number]");
        fh2_Muslihit_ZvsTheta[j]->GetYaxis()->SetTitleOffset(1.1);
        fh2_Muslihit_ZvsTheta[j]->GetXaxis()->CenterTitle(true);
        fh2_Muslihit_ZvsTheta[j]->GetYaxis()->CenterTitle(true);
        fh2_Muslihit_ZvsTheta[j]->GetXaxis()->SetLabelSize(0.045);
        fh2_Muslihit_ZvsTheta[j]->GetXaxis()->SetTitleSize(0.045);
        fh2_Muslihit_ZvsTheta[j]->GetYaxis()->SetLabelSize(0.045);
        fh2_Muslihit_ZvsTheta[j]->GetYaxis()->SetTitleSize(0.045);
        fh2_Muslihit_ZvsTheta[j]->SetFillColor(31);
        cMusliHit_ZvsTheta->cd(j + 1);
        fh2_Muslihit_ZvsTheta[j]->Draw("COL");
    }

    cMusliHit_Theta = new TCanvas("Musli_Theta", "", 10, 10, 800, 700);
    fh1_Muslihit_Theta = new TH1D("fh1_Musli_Theta", "Theta for data with nb_anodes = 2 if mult==1", 1000, -0.05, 0.05);
    fh1_Muslihit_Theta->GetXaxis()->SetTitle("Theta [mrad]");
    fh1_Muslihit_Theta->GetYaxis()->SetTitle("Counts");
    fh1_Muslihit_Theta->GetYaxis()->SetTitleOffset(1.1);
    fh1_Muslihit_Theta->GetXaxis()->CenterTitle(true);
    fh1_Muslihit_Theta->GetYaxis()->CenterTitle(true);
    fh1_Muslihit_Theta->GetXaxis()->SetLabelSize(0.045);
    fh1_Muslihit_Theta->GetXaxis()->SetTitleSize(0.045);
    fh1_Muslihit_Theta->GetYaxis()->SetLabelSize(0.045);
    fh1_Muslihit_Theta->GetYaxis()->SetTitleSize(0.045);
    fh1_Muslihit_Theta->SetFillColor(31);
    cMusliHit_Theta->cd();
    fh1_Muslihit_Theta->Draw("COL");

    // === ================= === //
    // === MAIN FOLDER-Musli === //
    // === ================= === //

    TFolder* mainfolMusli = new TFolder("MusicLightIon", "Music Light Ion info");

    mainfolMusli->Add(cMusliMap_Mult);

    mainfolMusli->Add(cMusliMap_E2);
    mainfolMusli->Add(cMusliMap_E4);
    mainfolMusli->Add(cMusliMap_E8);
    mainfolMusli->Add(cMusliMap_E16);

    mainfolMusli->Add(cMusliMap_EvsE2);
    mainfolMusli->Add(cMusliMap_EvsE4);
    mainfolMusli->Add(cMusliMap_EvsE8);

    mainfolMusli->Add(cMusliMap_DT2);
    mainfolMusli->Add(cMusliMap_DT4);
    mainfolMusli->Add(cMusliMap_DT8);
    mainfolMusli->Add(cMusliMap_DT16);

    mainfolMusli->Add(cMusliMap_DTvsDT2);
    mainfolMusli->Add(cMusliMap_DTvsDT4);
    mainfolMusli->Add(cMusliMap_DTvsDT8);

    mainfolMusli->Add(cMusliMap_EvsDT2);
    mainfolMusli->Add(cMusliMap_EvsDT4);
    mainfolMusli->Add(cMusliMap_EvsDT8);
    mainfolMusli->Add(cMusliMap_EvsDT16);

    mainfolMusli->Add(cMusliMap_Emean);

    mainfolMusli->Add(cMusliMap_DeltaDT);

    mainfolMusli->Add(cMusliMap_TchTtrig);

    mainfolMusli->Add(cMusliCal_Mult);

    mainfolMusli->Add(cMusliCal_E2);
    mainfolMusli->Add(cMusliCal_E4);
    mainfolMusli->Add(cMusliCal_E8);
    mainfolMusli->Add(cMusliCal_E16);

    mainfolMusli->Add(cMusliCal_DT2);
    mainfolMusli->Add(cMusliCal_DT4);
    mainfolMusli->Add(cMusliCal_DT8);
    mainfolMusli->Add(cMusliCal_DT16);

    mainfolMusli->Add(cMusliHit_E);
    mainfolMusli->Add(cMusliHit_Z);
    mainfolMusli->Add(cMusliHit_Theta);
    mainfolMusli->Add(cMusliHit_EvsX);
    mainfolMusli->Add(cMusliHit_ZvsX);
    mainfolMusli->Add(cMusliHit_EvsTheta);
    mainfolMusli->Add(cMusliHit_ZvsTheta);

    run->AddObject(mainfolMusli);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_MusLI_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BMusliOnlineSpectra::Reset_Histo()
{
    LOG(info) << "R3BMusliOnlineSpectra::Reset_Histo";

    // --- Reset Histo of mapped level --- //
    fh1_Muslimap_mult->Reset();
    fh2_Muslimap_mult->Reset();
    for (Int_t j = 0; j < 15; j++)
    {
        fh1_Muslimap_E[j]->Reset();
        fh1_Muslimap_DT[j]->Reset();
        fh2_Muslimap_EvsDT[j]->Reset();
        fh1_Muslimap_TchTtrig[j]->Reset();
    }
    fh1_Muslimap_TchTtrig[15]->Reset();
    for (Int_t j = 0; j < 13; j++)
    {
        if (j == 7 || j == 11)
            continue;
        fh2_Muslimap_EvsE[j]->Reset();
        fh2_Muslimap_DTvsDT[j]->Reset();
    }
    fh1_Muslimap_Emean2->Reset();
    fh1_Muslimap_Emean4->Reset();
    fh1_Muslimap_Emean8->Reset();
    fh2_Muslimap_Emean2vsEmean4->Reset();
    fh2_Muslimap_Emean2vsEmean8->Reset();
    fh2_Muslimap_Emean2vsE16->Reset();
    fh2_Muslimap_Emean4vsEmean8->Reset();
    fh2_Muslimap_Emean4vsE16->Reset();
    fh2_Muslimap_Emean8vsE16->Reset();
    fh1_Muslimap_DeltaDT->Reset();

    // --- Reset Histo of cal level --- //
    fh1_Muslical_mult->Reset();
    fh2_Muslical_mult->Reset();
    for (Int_t j = 0; j < 15; j++)
    {
        fh1_Muslical_E[j]->Reset();
        fh1_Muslical_DT[j]->Reset();
    }

    // --- Reset Histo of hit level --- //
    for (Int_t j = 0; j < 4; j++)
    {
        fh1_Muslihit_E[j]->Reset();
        fh1_Muslihit_Z[j]->Reset();
        fh2_Muslihit_EvsX[j]->Reset();
        fh2_Muslihit_ZvsX[j]->Reset();
        fh2_Muslihit_EvsTheta[j]->Reset();
        fh2_Muslihit_ZvsTheta[j]->Reset();
    }
    fh1_Muslihit_Theta->Reset();
}

void R3BMusliOnlineSpectra::Exec(Option_t* option)
{

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BMusliOnlineSpectra::Exec FairRootManager not found";

    Int_t nHits;
    UInt_t rank;

    // --- ----------- --- //
    // --- MAPPED DATA --- //
    // --- ----------- --- //

    if (fMapItemsMusli && fMapItemsMusli->GetEntriesFast() > 0)
    {
        // Local variables for mapped data
        UInt_t mult[18], e[15], t[18];
        UInt_t tref = 0;
        UInt_t ttrig = 0;
        Int_t dt[15];
        Float_t emean2 = 0;
        Float_t emean4 = 0;
        Float_t emean8 = 0;
        UShort_t cptemean2 = 0;
        UShort_t cptemean4 = 0;
        UShort_t cptemean8 = 0;
        for (UShort_t i = 0; i < 15; i++)
        {
            mult[i] = 0;
            e[i] = 0;
            t[i] = 0;
            dt[i] = 0;
        }
        for (UShort_t i = 15; i < 18; i++)
        {
            mult[i] = 0;
            t[i] = 0;
        }

        // read mapped data
        nHits = fMapItemsMusli->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusliMappedData* hit = dynamic_cast<R3BMusliMappedData*>(fMapItemsMusli->At(ihit));
            if (!hit)
                continue;
            rank = hit->GetSignal() - 1;
            fh1_Muslimap_mult->Fill(rank + 1);
            if (mult[rank] == 0) // we take the first hit
            {
                e[rank] = hit->GetEnergy();
                t[rank] = hit->GetTime();
            }
            mult[rank]++;
        }

        // hot fix for s091 and s118 with a bad cable for Tref
        if (0 < mult[16] && mult[16] <= 2)
            tref = t[16];

        if (mult[17] == 1)
            ttrig = t[17];

        // calculate drift time (and DeltaT to check searching window) and fill mapped data histo
        for (UShort_t i = 0; i < 15; i++)
        {
            if (mult[i] == 1)
            {
                fh1_Muslimap_E[i]->Fill(e[i]);
                if (ttrig > 0)
                {
                    fh1_Muslimap_TchTtrig[i]->Fill((double)(t[i] - ttrig));
                }
                if (tref > 0)
                {
                    dt[i] = t[i] - tref;
                    fh1_Muslimap_DT[i]->Fill(dt[i]);
                    fh2_Muslimap_EvsDT[i]->Fill(dt[i], e[i]);
                    if (ttrig > 0)
                        fh1_Muslimap_TchTtrig[15]->Fill(tref - ttrig);
                }
            }
            fh2_Muslimap_mult->Fill(i + 1, mult[i]);
        }
        for (UShort_t i = 15; i < 18; i++)
        {
            if (mult[i] == 1)
            {
                fh1_Muslimap_mult->Fill(i + 1);
            }
            fh2_Muslimap_mult->Fill(i + 1, mult[i]);
        }
        for (UShort_t i = 0; i < 13; i++)
        {
            if (i == 7 || i == 11)
                continue;
            if (mult[i] == 1 && mult[i + 1] == 1)
            {
                fh2_Muslimap_EvsE[i]->Fill(e[i + 1], e[i]);
                if (tref > 0)
                {
                    fh2_Muslimap_DTvsDT[i]->Fill(dt[i + 1], dt[i]);
                }
            }
        }
        if (mult[0] == 1 && mult[7] == 1 && dt[0] > 0 && dt[7] > 0)
        {
            fh1_Muslimap_DeltaDT->Fill(dt[7] - dt[0]);
        }

        // calculate the mean energy from the different type of data
        // fType = 2
        for (UShort_t i = 0; i < 8; i++)
        {
            if (mult[i] == 1)
            {
                emean2 += e[i];
                cptemean2++;
            }
        }
        if (cptemean2 == 8)
        {
            emean2 = emean2 / 8.;
            fh1_Muslimap_Emean2->Fill(emean2);
            if (mult[14] == 1)
                fh2_Muslimap_Emean2vsE16->Fill(e[14], emean2);
        }
        else
            emean2 = 0.;
        // fType = 4
        for (UShort_t i = 8; i < 12; i++)
        {
            if (mult[i] == 1)
            {
                emean4 += e[i];
                cptemean4++;
            }
        }
        if (cptemean4 == 4)
        {
            emean4 = emean4 / 4.;
            fh1_Muslimap_Emean4->Fill(emean4);
            if (emean2 > 0)
                fh2_Muslimap_Emean2vsEmean4->Fill(emean4, emean2);
            if (mult[14] == 1)
                fh2_Muslimap_Emean4vsE16->Fill(e[14], emean4);
        }
        else
            emean4 = 0.;
        // fType = 8
        for (UShort_t i = 12; i < 14; i++)
        {
            if (mult[i] == 1)
            {
                emean8 += e[i];
                cptemean8++;
            }
        }
        if (cptemean8 == 2)
        {
            emean8 = emean8 / 2.;
            fh1_Muslimap_Emean8->Fill(emean8);
            if (emean2 > 0)
                fh2_Muslimap_Emean2vsEmean8->Fill(emean8, emean2);
            if (emean4 > 0)
                fh2_Muslimap_Emean4vsEmean8->Fill(emean8, emean4);
            if (mult[14] == 1)
                fh2_Muslimap_Emean8vsE16->Fill(e[14], emean8);
        }
        else
            emean8 = 0.;

    } // end of if(fMapItemsMusli && fMapItemsMusli->GetEntriesFast()>0)

    // --- -------- --- //
    // --- CAL DATA --- //
    // --- -------- --- //

    if (fCalItemsMusli && fCalItemsMusli->GetEntriesFast() > 0)
    {
        // Local variables for cal data
        Double_t multcal[15], ecal[15], dtcal[18];
        for (UShort_t i = 0; i < 15; i++)
        {
            multcal[i] = 0;
            ecal[i] = 0;
            dtcal[i] = 0;
        }

        // read cal data
        nHits = fCalItemsMusli->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusliCalData* hit = dynamic_cast<R3BMusliCalData*>(fCalItemsMusli->At(ihit));
            if (!hit)
                continue;
            rank = hit->GetSignal() - 1;
            multcal[rank]++;
            fh1_Muslical_mult->Fill(rank + 1);
            ecal[rank] = hit->GetEnergy();
            dtcal[rank] = hit->GetDT();
        }
        // implement the spectra only if mult[signal] == 1
        for (Int_t i = 0; i < 15; i++)
        {
            if (multcal[i] == 1)
            {
                fh1_Muslical_E[i]->Fill(ecal[i]);
                fh1_Muslical_DT[i]->Fill(dtcal[i]);
            }
            fh2_Muslical_mult->Fill(i + 1, multcal[i]);
        }

    } // end of if (fCalItemsMusli)

    // --- -------- --- //
    // --- HIT DATA --- //
    // --- -------- --- //

    if (fHitItemsMusli && fHitItemsMusli->GetEntriesFast() > 0)
    {
        // Local variables for hit data
        Double_t multhit[4], ehit[4], zhit[4], xhit[4], theta[4];
        for (UShort_t i = 0; i < 4; i++)
        {
            multhit[i] = 0;
            ehit[i] = 0;
            zhit[i] = 0;
            xhit[i] = -1000.;
            theta[i] = -1000.;
        }

        // read  data
        nHits = fHitItemsMusli->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusliHitData* hit = dynamic_cast<R3BMusliHitData*>(fHitItemsMusli->At(ihit));
            if (!hit)
                continue;
            rank = hit->GetType() - 1;
            multhit[rank]++;
            ehit[rank] = hit->GetEave();
            zhit[rank] = hit->GetZcharge();
            xhit[rank] = hit->GetX();
            theta[rank] = hit->GetTheta();
        }
        // implement the spectra only if mult[signal] == 1
        for (Int_t i = 0; i < 4; i++)
        {
            if (multhit[i] == 1)
            {
                fh1_Muslihit_E[i]->Fill(ehit[i]);
                fh1_Muslihit_Z[i]->Fill(zhit[i]);
                fh2_Muslihit_EvsX[i]->Fill(xhit[i], ehit[i]);
                fh2_Muslihit_ZvsX[i]->Fill(xhit[i], zhit[i]);
                fh2_Muslihit_EvsTheta[i]->Fill(theta[0], ehit[i]);
                fh2_Muslihit_ZvsTheta[i]->Fill(theta[0], zhit[i]);
                if (i == 0)
                    fh1_Muslihit_Theta->Fill(theta[i]);
            }
        }
    }
    fNEvents += 1;
}

void R3BMusliOnlineSpectra::FinishEvent()
{
    if (fMapItemsMusli)
    {
        fMapItemsMusli->Clear();
    }
    if (fCalItemsMusli)
    {
        fCalItemsMusli->Clear();
    }
    if (fHitItemsMusli)
    {
        fHitItemsMusli->Clear();
    }
}

void R3BMusliOnlineSpectra::FinishTask()
{

    if (fMapItemsMusli)
    {
        cMusliMap_Mult->Write();

        cMusliMap_E2->Write();
        cMusliMap_E4->Write();
        cMusliMap_E8->Write();
        cMusliMap_E16->Write();

        cMusliMap_EvsE2->Write();
        cMusliMap_EvsE4->Write();
        cMusliMap_EvsE8->Write();

        cMusliMap_DT2->Write();
        cMusliMap_DT4->Write();
        cMusliMap_DT8->Write();
        cMusliMap_DT16->Write();

        cMusliMap_DTvsDT2->Write();
        cMusliMap_DTvsDT4->Write();
        cMusliMap_DTvsDT8->Write();

        cMusliMap_EvsDT2->Write();
        cMusliMap_EvsDT4->Write();
        cMusliMap_EvsDT8->Write();
        cMusliMap_EvsDT16->Write();

        cMusliMap_Emean->Write();

        cMusliMap_DeltaDT->Write();

        cMusliMap_TchTtrig->Write();
    }

    if (fCalItemsMusli)
    {
        cMusliCal_Mult->Write();

        cMusliCal_E2->Write();
        cMusliCal_E4->Write();
        cMusliCal_E8->Write();
        cMusliCal_E16->Write();

        cMusliCal_DT2->Write();
        cMusliCal_DT4->Write();
        cMusliCal_DT8->Write();
        cMusliCal_DT16->Write();
    }

    if (fHitItemsMusli)
    {
        cMusliHit_E->Write();
        cMusliHit_Z->Write();
        cMusliHit_Theta->Write();
    }
}

ClassImp(R3BMusliOnlineSpectra);
