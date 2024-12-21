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

#include "R3BCalifaDemoOnlineSpectra.h"
#include "R3BCalifaClusterData.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaMappedData.h"
#include "R3BEventHeader.h"
#include "R3BWRData.h"
#include "THttpServer.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

R3BCalifaDemoOnlineSpectra::R3BCalifaDemoOnlineSpectra()
    : FairTask("CALIFAOnlineSpectra", 1)
    , fMappedItemsCalifa(NULL)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
    , fWRItemsCalifa(NULL)
    , fWRItemsMaster(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fCalifaNumPetals(1)
    , fNumCrystalPetal(64)
    , fMapHistos_bins(3000)
    , fMapHistos_max(3000)
    , fRaw2Cal(kFALSE)
    , fLogScale(kTRUE)
    , fFebex2Preamp(kTRUE)
    , fCalON(kTRUE)
{
}

R3BCalifaDemoOnlineSpectra::R3BCalifaDemoOnlineSpectra(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsCalifa(NULL)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
    , fWRItemsCalifa(NULL)
    , fWRItemsMaster(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fCalifaNumPetals(1)
    , fNumCrystalPetal(64)
    , fMapHistos_bins(3000)
    , fMapHistos_max(3000)
    , fRaw2Cal(kFALSE)
    , fLogScale(kTRUE)
    , fFebex2Preamp(kTRUE)
    , fCalON(kTRUE)
{
}

R3BCalifaDemoOnlineSpectra::~R3BCalifaDemoOnlineSpectra()
{
    LOG(info) << "R3BCalifaDemoOnlineSpectra: Delete instance";
    delete fMappedItemsCalifa;
    delete fCalItemsCalifa;
    delete fHitItemsCalifa;
    delete fWRItemsCalifa;
    delete fWRItemsMaster;
}

InitStatus R3BCalifaDemoOnlineSpectra::Init()
{
    LOG(info) << "R3BCalifaDemoOnlineSpectra::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BCalifaDemoOnlineSpectra::Init FairRootManager not found";

    // Define Preamp sequence for histograms
    Int_t PreampOrder[16] = { 6, 5, 4, 3, 2, 1, 0, 7, 8, 15, 14, 13, 12, 11, 10, 9 };
    for (Int_t i = 0; i < 16; i++)
    {
        fOrderFebexPreamp[i] = PreampOrder[i];
    }
    Int_t BinsChannelFebex = 65535;

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    FairRunOnline* run = FairRunOnline::Instance();

    run->GetHttpServer()->Register("", this);

    // get access to Mapped data
    fMappedItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaMappedData"));
    if (!fMappedItemsCalifa)
    {
        return kFATAL;
    }

    // get access to Cal data
    fCalItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaCrystalCalData"));
    if (!fCalItemsCalifa)
    {
        LOG(info) << "R3BCalifaDemoOnlineSpectra::Init CalifaCrystalCalData not found";
        fCalON = kFALSE;
    }

    // get access to Hit data
    fHitItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaClusterData"));
    if (!fHitItemsCalifa)
    {
        LOG(info) << "R3BCalifaDemoOnlineSpectra::Init CalifaClusterData not found";
    }

    // get access to WR-Califa data
    fWRItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("WRCalifaData"));
    if (!fWRItemsCalifa)
    {
        LOG(info) << "R3BCalifaDemoOnlineSpectra::Init WRCalifaData not found";
    }

    // get access to WR-Master data
    fWRItemsMaster = dynamic_cast<TClonesArray*>(mgr->GetObject("WRMasterData"));
    if (!fWRItemsMaster)
    {
        LOG(info) << "R3BCalifaDemoOnlineSpectra::Init WRMasterData not found";
    }

    // reading the file
    Bool_t noFile = kFALSE;
    ifstream* FileHistos = new ifstream(fCalifaFile);
    if (!FileHistos->is_open())
    {
        LOG(warn) << "R3BCalifaDemoOnlineSpectra:  No Histogram definition file";
        noFile = kTRUE;
    }

    Double_t arry_bins[fCalifaNumPetals * fNumCrystalPetal];
    Double_t arry_maxE[fCalifaNumPetals * fNumCrystalPetal];
    Double_t arry_minE[fCalifaNumPetals * fNumCrystalPetal];

    Double_t bins = fMapHistos_bins;
    Double_t maxE = fMapHistos_max;
    Double_t minE = 0.;

    if (!noFile)
    {
        for (Int_t i = 0; i < fCalifaNumPetals * fNumCrystalPetal; i++)
        {
            *FileHistos >> arry_bins[i] >> arry_minE[i] >> arry_maxE[i];
        }
    }
    else
    {
        for (Int_t i = 0; i < fCalifaNumPetals * fNumCrystalPetal; i++)
        {
            arry_bins[i] = fMapHistos_bins;
            arry_minE[i] = 0;
            arry_maxE[i] = fMapHistos_max;
        }
    }

    // CANVAS 1
    cCalifa1 = new TCanvas("Califa_petal_vs_cryNb", "Califa_petal_vs_cryNb", 10, 10, 500, 500);

    fh_Califa_cryId_petal = new TH2F("fh_Califa_cryNb_petal",
                                     "Califa petal number vs crystal number",
                                     fNumCrystalPetal,
                                     0.5,
                                     fNumCrystalPetal + 0.5,
                                     fCalifaNumPetals,
                                     0.5,
                                     fCalifaNumPetals + 0.5);
    fh_Califa_cryId_petal->GetXaxis()->SetTitle("Crystal number [1-64]");
    fh_Califa_cryId_petal->GetYaxis()->SetTitle("Petal number [1-8]");
    fh_Califa_cryId_petal->GetYaxis()->SetTitleOffset(1.2);
    fh_Califa_cryId_petal->GetXaxis()->CenterTitle(true);
    fh_Califa_cryId_petal->GetYaxis()->CenterTitle(true);
    fh_Califa_cryId_petal->Draw("COLZ");

    TString Name1;
    TString Name2;
    TString Yaxis1;

    TString Name3;
    TString Name4;
    TString Yaxis2;

    // raw data
    Name1 = "Califa_Map_cryNb_energy";
    Name2 = "Califa_Map energy vs cryNb";
    Yaxis1 = "Energy (channels)";

    // cal data
    Name3 = "Califa_Cal_cryNb_energy";
    Name4 = "Califa_Cal energy vs cryNb";
    Yaxis2 = "Energy (keV)";

    // CANVAS 2
    cCalifa2 = new TCanvas(Name1, Name2, 10, 10, 500, 500);
    fh_Califa_cryId_energy = new TH2F(Name1,
                                      Name2,
                                      fNumCrystalPetal * fCalifaNumPetals,
                                      1.,
                                      fNumCrystalPetal * fCalifaNumPetals + 1.,
                                      BinsChannelFebex,
                                      0,
                                      BinsChannelFebex);
    fh_Califa_cryId_energy->GetXaxis()->SetTitle("Crystal number");
    fh_Califa_cryId_energy->GetYaxis()->SetTitle(Yaxis1);
    fh_Califa_cryId_energy->GetYaxis()->SetTitleOffset(1.4);
    fh_Califa_cryId_energy->GetXaxis()->CenterTitle(true);
    fh_Califa_cryId_energy->GetYaxis()->CenterTitle(true);
    fh_Califa_cryId_energy->Draw("COLZ");

    fh_Califa_cryId_energy_cal = new TH2F(Name3,
                                          Name4,
                                          fNumCrystalPetal * fCalifaNumPetals,
                                          1.,
                                          fNumCrystalPetal * fCalifaNumPetals + 1.,
                                          bins,
                                          minE,
                                          maxE);
    fh_Califa_cryId_energy_cal->GetXaxis()->SetTitle("Crystal number");
    fh_Califa_cryId_energy_cal->GetYaxis()->SetTitle(Yaxis2);
    fh_Califa_cryId_energy_cal->GetYaxis()->SetTitleOffset(1.4);
    fh_Califa_cryId_energy_cal->GetXaxis()->CenterTitle(true);
    fh_Califa_cryId_energy_cal->GetYaxis()->CenterTitle(true);

    // CANVAS 3
    cCalifa3 = new TCanvas("Califa Energy_per_petal", "Califa Energy per petal", 10, 10, 500, 500);
    cCalifa3->Divide(4, 2);
    char Name5[255];
    char Name6[255];
    TString Xaxis1;
    char Name7[255];
    char Name8[255];
    TString Xaxis2;
    for (Int_t i = 0; i < fCalifaNumPetals; i++)
    {
        if (i < 7)
        { // FIXME:for experiment s444
            sprintf(Name5, "h_Califa_Map_energy_for_petal_%d", i + 1);
            sprintf(Name6, "Califa Map Energy for petal %d", i + 1);
            Xaxis1 = "Energy (channels)";
            fh_Califa_energy_per_petal[i] = new TH1F(Name5, Name6, BinsChannelFebex, 0, BinsChannelFebex);
            fh_Califa_energy_per_petal[i]->GetXaxis()->SetTitle(Xaxis1);

            sprintf(Name7, "h_Califa_Cal_energy_for_petal_%d", i + 1);
            sprintf(Name8, "Califa Cal Energy for petal %d", i + 1);
            Xaxis2 = "Energy (keV)";
            fh_Califa_energy_per_petal_cal[i] = new TH1F(Name7, Name8, bins, minE, maxE);
            fh_Califa_energy_per_petal_cal[i]->GetXaxis()->SetTitle(Xaxis2);

            cCalifa3->cd(i + 1);
            gPad->SetLogy();
            fh_Califa_energy_per_petal[i]->Draw();
        }
    }

    // CANVAS 3bis for proton range: experiment s444
    cCalifa3b =
        new TCanvas("Califa Energy_per_petal_proton", "Califa Energy per petal for proton range", 10, 10, 500, 500);
    cCalifa3b->Divide(2, 1);

    sprintf(Name5, "h_Califa_Map_energy_for_petal__%d_proton", 1);
    sprintf(Name6, "Califa Map Energy for petal %d proton range", 1);
    Xaxis1 = "Energy (channels)";
    fh_Califa_energy_per_petal[7] = new TH1F(Name5, Name6, BinsChannelFebex, 0, BinsChannelFebex);
    fh_Califa_energy_per_petal[7]->GetXaxis()->SetTitle(Xaxis1);

    sprintf(Name7, "h_Califa_Cal_energy_for_petal_%d_proton", 1);
    sprintf(Name8, "Califa Cal Energy for petal %d proton range", 1);
    Xaxis2 = "Energy (keV)";
    fh_Califa_energy_per_petal_cal[7] = new TH1F(Name7, Name8, bins, minE, maxE);
    fh_Califa_energy_per_petal_cal[7]->GetXaxis()->SetTitle(Xaxis2);

    cCalifa3b->cd(1);
    gPad->SetLogy();
    fh_Califa_energy_per_petal[7]->Draw();

    sprintf(Name5, "h_Califa_Map_energy_for_petal__%d_proton", 4);
    sprintf(Name6, "Califa Map Energy for petal %d proton range", 4);
    Xaxis1 = "Energy (channels)";
    fh_Califa_energy_per_petal[8] = new TH1F(Name5, Name6, BinsChannelFebex, 0, BinsChannelFebex);
    fh_Califa_energy_per_petal[8]->GetXaxis()->SetTitle(Xaxis1);

    sprintf(Name7, "h_Califa_Cal_energy_for_petal_%d_proton", 4);
    sprintf(Name8, "Califa Cal Energy for petal %d proton range", 4);
    Xaxis2 = "Energy (keV)";
    fh_Califa_energy_per_petal_cal[8] = new TH1F(Name7, Name8, bins, minE, maxE);
    fh_Califa_energy_per_petal_cal[8]->GetXaxis()->SetTitle(Xaxis2);

    cCalifa3b->cd(2);
    gPad->SetLogy();
    fh_Califa_energy_per_petal[8]->Draw();

    // CANVAS 4
    char Name9[255];
    char Name10[255];
    char Name11[255];
    TString Xaxis3;
    char Name12[255];
    char Name13[255];
    TString Xaxis4;
    for (Int_t i = 0; i < fCalifaNumPetals; i++)
    {
        for (Int_t k = 0; k < 4; k++)
        {

            if (i < 7)
            {
                sprintf(Name9, "Califa_Petal_%d_section_%d", i + 1, k + 1);
                cCalifa4[i][k] = new TCanvas(Name9, Name9, 10, 10, 500, 500);
                cCalifa4[i][k]->Divide(4, 4);
            }
            else
            {
                if (k < 2)
                {
                    sprintf(Name9, "Califa_Petal_%d_section_%d_proton", 1, k + 1);
                    cCalifa4[i][k] = new TCanvas(Name9, Name9, 10, 10, 500, 500);
                    cCalifa4[i][k]->Divide(4, 4);
                }
                else
                {
                    sprintf(Name9, "Califa_Petal_%d_section_%d_proton", 4, k - 1);
                    cCalifa4[i][k] = new TCanvas(Name9, Name9, 10, 10, 500, 500);
                    cCalifa4[i][k]->Divide(4, 4);
                }
            }

            for (Int_t j = 0; j < 16; j++)
            {

                if (i < 7)
                {
                    sprintf(Name10, "h_Califa_Map_Petal_%d_%d_Crystal_%d_energy", i + 1, k + 1, j + 1);
                    sprintf(Name11, "Map level, Petal %d, Preamp %d, Febex ch. %d", i + 1, k + 1, j);
                    sprintf(Name12, "h_Califa_Cal_Petal_%d_%d_Crystal_%d_energy", i + 1, k + 1, j + 1);
                    sprintf(Name13, "Cal level, Petal %d, Preamp %d, Febex ch. %d", i + 1, k + 1, j);
                }
                else
                {
                    if (k < 2)
                    {
                        sprintf(Name10, "h_Califa_Map_Petal_%d_%d_Crystal_%d_energy_p", 1, k + 1, j + 1);
                        sprintf(Name11, "Map level, Petal %d, Preamp %d, Febex ch. %d", 1, k + 1, j);
                        sprintf(Name12, "h_Califa_Cal_Petal_%d_%d_Crystal_%d_energy_p", 1, k + 1, j + 1);
                        sprintf(Name13, "Cal level, Petal %d, Preamp %d, Febex ch. %d", 1, k + 1, j);
                    }
                    else
                    {
                        sprintf(Name10, "h_Califa_Map_Petal_%d_%d_Crystal_%d_energy_p", 4, k - 1, j + 1);
                        sprintf(Name11, "Map level, Petal %d, Preamp %d, Febex ch. %d", 4, k - 1, j);
                        sprintf(Name12, "h_Califa_Cal_Petal_%d_%d_Crystal_%d_energy_p", 4, k - 1, j + 1);
                        sprintf(Name13, "Cal level, Petal %d, Preamp %d, Febex ch. %d", 4, k - 1, j);
                    }
                }
                Xaxis3 = "Energy (channels)";
                Xaxis4 = "Energy (keV)";

                fh_Califa_crystals[i][j + 16 * k] = new TH1F(Name10, Name11, BinsChannelFebex, 0, BinsChannelFebex);
                fh_Califa_crystals[i][j + 16 * k]->SetTitleSize(1.5, "t");
                fh_Califa_crystals[i][j + 16 * k]->GetXaxis()->SetTitle(Xaxis3);
                fh_Califa_crystals[i][j + 16 * k]->GetXaxis()->SetLabelSize(0.06);
                fh_Califa_crystals[i][j + 16 * k]->GetYaxis()->SetLabelSize(0.07);
                fh_Califa_crystals[i][j + 16 * k]->GetXaxis()->SetTitleSize(0.05);
                fh_Califa_crystals[i][j + 16 * k]->GetXaxis()->CenterTitle(true);
                fh_Califa_crystals[i][j + 16 * k]->GetYaxis()->CenterTitle(true);
                fh_Califa_crystals[i][j + 16 * k]->GetXaxis()->SetTitleOffset(0.95);
                fh_Califa_crystals[i][j + 16 * k]->SetFillColor(kGreen);

                fh_Califa_crystals_cal[i][j + 16 * k] =
                    new TH1F(Name12, Name13, arry_bins[j + 16 * k], arry_minE[j + 16 * k], arry_maxE[j + 16 * k]);
                fh_Califa_crystals_cal[i][j + 16 * k]->SetTitleSize(1.5, "t");
                fh_Califa_crystals_cal[i][j + 16 * k]->GetXaxis()->SetTitle(Xaxis4);
                fh_Califa_crystals_cal[i][j + 16 * k]->GetXaxis()->SetLabelSize(0.06);
                fh_Califa_crystals_cal[i][j + 16 * k]->GetYaxis()->SetLabelSize(0.07);
                fh_Califa_crystals_cal[i][j + 16 * k]->GetXaxis()->SetTitleSize(0.05);
                fh_Califa_crystals_cal[i][j + 16 * k]->GetXaxis()->CenterTitle(true);
                fh_Califa_crystals_cal[i][j + 16 * k]->GetYaxis()->CenterTitle(true);
                fh_Califa_crystals_cal[i][j + 16 * k]->GetXaxis()->SetTitleOffset(0.95);
                fh_Califa_crystals_cal[i][j + 16 * k]->SetFillColor(kGreen);

                cCalifa4[i][k]->cd(j + 1);
                gPad->SetLogy();
                fh_Califa_crystals[i][j + 16 * k]->Draw();
            }
        }
    }

    // CANVAS 5
    cCalifa5 = new TCanvas("Califa_petal_1_vs_petal_3and4", "Coinc. 1 vs 3and4", 10, 10, 500, 500);

    fh_Califa_coinc_petal1 = new TH2F(
        "fh_Califa_petal_1_vs_petal_3and4", "Califa coinc. petals 1 vs 3 and 4", bins, minE, maxE, bins, minE, maxE);
    fh_Califa_coinc_petal1->GetXaxis()->SetTitle("Energy petal 1  (keV)");
    fh_Califa_coinc_petal1->GetYaxis()->SetTitle("Energy petal 3 and 4 (keV)");
    fh_Califa_coinc_petal1->GetYaxis()->SetTitleOffset(1.2);
    fh_Califa_coinc_petal1->GetXaxis()->CenterTitle(true);
    fh_Califa_coinc_petal1->GetYaxis()->CenterTitle(true);
    fh_Califa_coinc_petal1->Draw("COLZ");

    // CANVAS 6
    cCalifa6 = new TCanvas("Califa_petal_2_vs_petal_5", "Coinc. 2 vs 5", 10, 10, 500, 500);

    fh_Califa_coinc_petal2 =
        new TH2F("fh_Califa_petal_2_vs_petal_5", "Califa coinc. petals 2 vs 5", bins, minE, maxE, bins, minE, maxE);
    fh_Califa_coinc_petal2->GetXaxis()->SetTitle("Energy petal 2 (keV)");
    fh_Califa_coinc_petal2->GetYaxis()->SetTitle("Energy petal 5 (keV)");
    fh_Califa_coinc_petal2->GetYaxis()->SetTitleOffset(1.2);
    fh_Califa_coinc_petal2->GetXaxis()->CenterTitle(true);
    fh_Califa_coinc_petal2->GetYaxis()->CenterTitle(true);
    fh_Califa_coinc_petal2->Draw("COLZ");

    // CANVAS 7
    cCalifa7 = new TCanvas("Califa_petal_1_vs_petal_4_proton", "Coinc. 1 vs 4 for proton range", 10, 10, 500, 500);

    fh_Califa_coinc_petal3 = new TH2F("fh_Califa_petal_1_vs_petal_4_proton_range",
                                      "Califa coinc. petals 1 vs 4 for proton range",
                                      bins,
                                      minE,
                                      maxE,
                                      bins,
                                      minE,
                                      maxE);
    fh_Califa_coinc_petal3->GetXaxis()->SetTitle("Energy petal 1  (keV)");
    fh_Califa_coinc_petal3->GetYaxis()->SetTitle("Energy petal 4 (keV)");
    fh_Califa_coinc_petal3->GetYaxis()->SetTitleOffset(1.2);
    fh_Califa_coinc_petal3->GetXaxis()->CenterTitle(true);
    fh_Califa_coinc_petal3->GetYaxis()->CenterTitle(true);
    fh_Califa_coinc_petal3->Draw("COLZ");

    // CANVAS 8
    cCalifa8 = new TCanvas("Califa_Theta_vs_Phi", "Theta vs Phi", 10, 10, 500, 500);
    fh_Califa_theta_phi = new TH2F("fh_Califa_theta_vs_phi", "Califa theta vs phi", 90, 0, 90, 180, -180, 180);
    fh_Califa_theta_phi->GetXaxis()->SetTitle("Theta (degrees)");
    fh_Califa_theta_phi->GetYaxis()->SetTitle("Phi (degrees)");
    fh_Califa_theta_phi->GetYaxis()->SetTitleOffset(1.2);
    fh_Califa_theta_phi->GetXaxis()->CenterTitle(true);
    fh_Califa_theta_phi->GetYaxis()->CenterTitle(true);
    fh_Califa_theta_phi->Draw("COLZ");

    // CANVAS 9
    char Name14[255], Name15[255], Name16[255];
    for (Int_t k = 0; k < 7; k++)
    { // 7 petals: s444
        sprintf(Name14, "Califa_Petal_%d_theta_vs_energy", k + 1);
        sprintf(Name15, "fh_Califa_Petal_%d_theta_vs_energy", k + 1);
        sprintf(Name16, "Califa theta vs energy for petal %d", k + 1);
        cCalifa_hitpetal[k] = new TCanvas(Name14, Name14, 10, 10, 500, 500);
        fh_Califa_theta_energy[k] = new TH2F(Name15, Name16, 45, 0, 90, bins, minE, maxE);
        fh_Califa_theta_energy[k]->GetXaxis()->SetTitle("Theta (degrees)");
        fh_Califa_theta_energy[k]->GetYaxis()->SetTitle("Energy (keV)");
        fh_Califa_theta_energy[k]->GetYaxis()->SetTitleOffset(1.4);
        fh_Califa_theta_energy[k]->GetXaxis()->CenterTitle(true);
        fh_Califa_theta_energy[k]->GetYaxis()->CenterTitle(true);
        fh_Califa_theta_energy[k]->Draw("COLZ");
    }
    // Theta vs energy
    sprintf(Name14, "Califa_calorimeter_theta_vs_energy");
    sprintf(Name15, "fh_Califa_Petal_total_theta_vs_energy");
    sprintf(Name16, "Califa theta vs energy for full calorimeter");
    cCalifa_hitpetal[7] = new TCanvas(Name14, Name14, 10, 10, 500, 500);
    fh_Califa_theta_energy[7] = new TH2F(Name15, Name16, 45, 0, 90, bins, minE, maxE);
    fh_Califa_theta_energy[7]->GetXaxis()->SetTitle("Theta (degrees)");
    fh_Califa_theta_energy[7]->GetYaxis()->SetTitle("Energy (keV)");
    fh_Califa_theta_energy[7]->GetYaxis()->SetTitleOffset(1.4);
    fh_Califa_theta_energy[7]->GetXaxis()->CenterTitle(true);
    fh_Califa_theta_energy[7]->GetYaxis()->CenterTitle(true);
    fh_Califa_theta_energy[7]->Draw("COLZ");

    // Total
    sprintf(Name14, "Califa_calorimeter_total_energy_per_hit");
    sprintf(Name15, "fh_Califa_Petal_total_energy");
    sprintf(Name16, "Califa total energy per hit for full calorimeter");
    cCalifa_hitenergy = new TCanvas(Name14, Name14, 10, 10, 500, 500);
    fh_Califa_total_energy = new TH1F(Name15, Name16, bins, minE, maxE);
    fh_Califa_total_energy->GetXaxis()->SetTitle("Energy (keV)");
    fh_Califa_total_energy->GetYaxis()->SetTitle("Counts");
    fh_Califa_total_energy->GetYaxis()->SetTitleOffset(1.4);
    fh_Califa_total_energy->GetXaxis()->CenterTitle(true);
    fh_Califa_total_energy->GetYaxis()->CenterTitle(true);
    fh_Califa_total_energy->Draw("");
    gPad->SetLogy();

    // Difference between WRs
    sprintf(Name14, "WR_Califa_Master");
    sprintf(Name15, "fh_WR_Califa_Master");
    sprintf(Name16, "WR-Califa - WR-Master");
    cCalifa_wr = new TCanvas(Name14, Name14, 10, 10, 500, 500);
    fh_Califa_wr = new TH1F(Name15, Name16, 600, -1, 4100);
    fh_Califa_wr->GetXaxis()->SetTitle("WR-Califa - WR-Master");
    fh_Califa_wr->GetYaxis()->SetTitle("Counts");
    fh_Califa_wr->GetYaxis()->SetTitleOffset(1.4);
    fh_Califa_wr->GetXaxis()->CenterTitle(true);
    fh_Califa_wr->GetYaxis()->CenterTitle(true);
    fh_Califa_wr->Draw("");

    // CANVAS 10
    char Name17[255];
    char Name18[255];
    char Name19[255];
    char Name20[255];
    cCalifa10 = new TCanvas("Califa_Multiplicities", "Califa_Multiplicities", 10, 10, 500, 500);
    cCalifa10->Divide(3, 3);
    fh_Califa_Mult[0] = new TH1F("fh_Califa_Mult", "Califa multiplicities", 20, 0, 20);
    fh_Califa_MultHit[0] = new TH1F("fh_Califa_MultHit", "Califa Hit multiplicities", 20, 0, 20);
    fh_Califa_Mult[0]->GetXaxis()->SetTitle("Mult (crystal:blue, cluster:red)");
    fh_Califa_Mult[0]->GetXaxis()->CenterTitle(true);
    fh_Califa_Mult[0]->GetYaxis()->CenterTitle(true);
    fh_Califa_Mult[0]->GetYaxis()->SetTitleOffset(1.2);
    fh_Califa_Mult[0]->GetXaxis()->SetTitleOffset(1.1);
    cCalifa10->cd(1);
    gPad->SetLogy();
    fh_Califa_Mult[0]->Draw();
    fh_Califa_MultHit[0]->SetLineColor(kRed);
    fh_Califa_MultHit[0]->Draw("SAME");

    for (Int_t i = 0; i < fCalifaNumPetals; i++)
    {
        sprintf(Name17, "h_Califa_Mult_Petal_%d", i + 1);
        sprintf(Name18, "Califa multiplicities for petal %d", i + 1);
        sprintf(Name19, "h_Califa_HitMult_Petal_%d", i + 1);
        sprintf(Name20, "Califa Hit multiplicities for petal %d", i + 1);
        fh_Califa_Mult[i + 1] = new TH1F(Name17, Name18, 20, 0, 20);
        if (i < fCalifaNumPetals - 1)
            fh_Califa_MultHit[i + 1] = new TH1F(Name19, Name20, 20, 0, 20);
        fh_Califa_Mult[i + 1]->GetXaxis()->SetTitle("Mult (crystal:blue, cluster:red)");
        fh_Califa_Mult[i + 1]->GetXaxis()->CenterTitle(true);
        fh_Califa_Mult[i + 1]->GetYaxis()->CenterTitle(true);
        fh_Califa_Mult[i + 1]->GetYaxis()->SetTitleOffset(1.2);
        fh_Califa_Mult[i + 1]->GetXaxis()->SetTitleOffset(1.1);
        cCalifa10->cd(i + 2);
        gPad->SetLogy();
        fh_Califa_Mult[i + 1]->Draw();
        if (i < fCalifaNumPetals - 1)
        {
            fh_Califa_MultHit[i + 1]->SetLineColor(kRed);
            fh_Califa_MultHit[i + 1]->Draw("SAME");
        }
    }

    // MAIN FOLDER-Califa
    TFolder* mainfolCalifa = new TFolder("CALIFA", "CALIFA info");
    mainfolCalifa->Add(cCalifa_wr);
    mainfolCalifa->Add(cCalifa1);
    mainfolCalifa->Add(cCalifa2);
    mainfolCalifa->Add(cCalifa3);
    mainfolCalifa->Add(cCalifa3b); // proton range, s444
    for (Int_t i = 0; i < fCalifaNumPetals; i++)
    {
        for (Int_t k = 0; k < 4; k++)
        {
            mainfolCalifa->Add(cCalifa4[i][k]);
        }
    }
    if (fCalON)
    {
        mainfolCalifa->Add(cCalifa5);
        mainfolCalifa->Add(cCalifa6);
        mainfolCalifa->Add(cCalifa7);
        for (Int_t i = 0; i < fCalifaNumPetals; i++)
            mainfolCalifa->Add(cCalifa_hitpetal[i]);
        mainfolCalifa->Add(cCalifa8);
        mainfolCalifa->Add(cCalifa_hitenergy);
        mainfolCalifa->Add(cCalifa10);
    }
    run->AddObject(mainfolCalifa);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_Califa", Form("/Objects/%s/->Reset_CALIFA_DEMO_Histo()", GetName()));
    // Register command to change the histogram scales (Log/Lineal)
    run->GetHttpServer()->RegisterCommand("Log_Califa", Form("/Objects/%s/->Log_CALIFA_Histo()", GetName()));
    // Register command for moving between Febex and Preamp channels
    run->GetHttpServer()->RegisterCommand("Febex2Preamp_Califa",
                                          Form("/Objects/%s/->Febex2Preamp_CALIFA_Histo()", GetName()));
    // Register command for moving between histograms of mapped and cal levels
    if (fCalON)
    {
        run->GetHttpServer()->RegisterCommand("Map2Cal_Califa",
                                              Form("/Objects/%s/->Map2Cal_CALIFA_Histo()", GetName()));
    }

    return kSUCCESS;
}

void R3BCalifaDemoOnlineSpectra::Reset_CALIFA_DEMO_Histo()
{
    LOG(info) << "R3BCalifaDemoOnlineSpectra::Reset_CALIFA_Histo";

    fh_Califa_cryId_petal->Reset();
    fh_Califa_wr->Reset();
    fh_Califa_cryId_energy->Reset();
    fh_Califa_coinc_petal1->Reset();
    fh_Califa_coinc_petal2->Reset();
    fh_Califa_coinc_petal3->Reset();
    fh_Califa_energy_per_petal[8]->Reset(); // s444
    fh_Califa_Mult[0]->Reset();
    fh_Califa_theta_phi->Reset();

    for (Int_t i = 0; i < fCalifaNumPetals; i++)
    {
        fh_Califa_theta_energy[i]->Reset();
        fh_Califa_energy_per_petal[i]->Reset();
        fh_Califa_Mult[i + 1]->Reset();
        fh_Califa_MultHit[i]->Reset();
        for (Int_t k = 0; k < 4; k++)
        {
            for (Int_t j = 0; j < 16; j++)
            {
                fh_Califa_crystals[i][j + 16 * k]->Reset();
            }
        }
    }

    if (fCalON)
    {
        fh_Califa_cryId_energy_cal->Reset();
        fh_Califa_energy_per_petal_cal[8]->Reset(); // s444
        fh_Califa_total_energy->Reset();
        for (Int_t i = 0; i < fCalifaNumPetals; i++)
        {
            fh_Califa_energy_per_petal_cal[i]->Reset();
            for (Int_t k = 0; k < 4; k++)
            {
                for (Int_t j = 0; j < 16; j++)
                {
                    fh_Califa_crystals_cal[i][j + 16 * k]->Reset();
                }
            }
        }
    }
}

void R3BCalifaDemoOnlineSpectra::Log_CALIFA_Histo()
{
    LOG(info) << "R3BCalifaDemoOnlineSpectra::Log_CALIFA_Histo";

    for (Int_t i = 0; i < fCalifaNumPetals; i++)
    {
        for (Int_t k = 0; k < 4; k++)
        {
            for (Int_t j = 0; j < 16; j++)
            {
                cCalifa4[i][k]->cd(j + 1);
                if (fLogScale)
                {
                    gPad->SetLogy(0);
                }
                else
                {
                    gPad->SetLogy(1);
                }
            }
        }
    }
    for (Int_t i = 0; i < fCalifaNumPetals - 1; i++)
    {
        cCalifa3->cd(i + 1);
        if (fLogScale)
        {
            gPad->SetLogy(0);
        }
        else
        {
            gPad->SetLogy(1);
        }
    }
    for (Int_t i = 0; i < 2; i++)
    {
        cCalifa3b->cd(i + 1);
        if (fLogScale)
        {
            gPad->SetLogy(0);
        }
        else
        {
            gPad->SetLogy(1);
        }
    }

    cCalifa_hitenergy->cd();
    if (fLogScale)
    {
        gPad->SetLogy(0);
    }
    else
    {
        gPad->SetLogy(1);
    }

    if (fLogScale)
        fLogScale = kFALSE;
    else
        fLogScale = kTRUE;
}

void R3BCalifaDemoOnlineSpectra::Map2Cal_CALIFA_Histo()
{
    LOG(info) << "R3BCalifaDemoOnlineSpectra::Map2Cal_CALIFA_Histo";

    char Name[255];

    if (fFebex2Preamp)
    { // Febex sequence

        if (fRaw2Cal)
        {
            cCalifa2->cd();
            fh_Califa_cryId_energy->Draw("COLZ");

            for (Int_t i = 0; i < fCalifaNumPetals - 1; i++)
            {
                cCalifa3->cd(i + 1);
                fh_Califa_energy_per_petal[i]->Draw();
            }

            for (Int_t i = 0; i < 2; i++)
            {
                cCalifa3b->cd(i + 1);
                fh_Califa_energy_per_petal[7 + i]->Draw();
            }

            for (Int_t i = 0; i < fCalifaNumPetals; i++)
            {
                for (Int_t k = 0; k < 4; k++)
                {
                    for (Int_t j = 0; j < 16; j++)
                    {
                        cCalifa4[i][k]->cd(j + 1);
                        gPad->SetLogy();
                        fh_Califa_crystals[i][j + 16 * k]->SetFillColor(kGreen);
                        if (i == 7 && k < 2)
                            sprintf(Name, "Map level, Petal %d, Preamp %d, Febex ch. %d", 1, k + 1, j);
                        else if (i == 7 && k > 1)
                            sprintf(Name, "Map level, Petal %d, Preamp %d, Febex ch. %d", 4, k - 1, j);
                        else
                            sprintf(Name, "Map level, Petal %d, Preamp %d, Febex ch. %d", i + 1, k + 1, j);
                        fh_Califa_crystals[i][j + 16 * k]->SetTitle(Name);
                        fh_Califa_crystals[i][j + 16 * k]->Draw();
                    }
                }
            }

            fRaw2Cal = kFALSE;
        }
        else
        {

            cCalifa2->cd();
            fh_Califa_cryId_energy_cal->Draw("COLZ");

            for (Int_t i = 0; i < fCalifaNumPetals - 1; i++)
            {
                cCalifa3->cd(i + 1);
                fh_Califa_energy_per_petal_cal[i]->Draw();
            }

            for (Int_t i = 0; i < 2; i++)
            {
                cCalifa3b->cd(i + 1);
                fh_Califa_energy_per_petal_cal[7 + i]->Draw();
            }

            for (Int_t i = 0; i < fCalifaNumPetals; i++)
            {
                for (Int_t k = 0; k < 4; k++)
                {
                    for (Int_t j = 0; j < 16; j++)
                    {
                        cCalifa4[i][k]->cd(j + 1);
                        gPad->SetLogy();
                        fh_Califa_crystals_cal[i][j + 16 * k]->SetFillColor(kGreen);
                        if (i == 7 && k < 2)
                            sprintf(Name, "Cal level, Petal %d, Preamp %d, Febex ch. %d", 1, k + 1, j);
                        else if (i == 7 && k > 1)
                            sprintf(Name, "Cal level, Petal %d, Preamp %d, Febex ch. %d", 4, k - 1, j);
                        else
                            sprintf(Name, "Cal level, Petal %d, Preamp %d, Febex ch. %d", i + 1, k + 1, j);
                        fh_Califa_crystals_cal[i][j + 16 * k]->SetTitle(Name);
                        fh_Califa_crystals_cal[i][j + 16 * k]->Draw();
                    }
                }
            }

            fRaw2Cal = kTRUE;
        }
    }
    else
    { // Preamp sequence

        if (fRaw2Cal)
        {
            cCalifa2->cd();
            fh_Califa_cryId_energy->Draw("COLZ");

            for (Int_t i = 0; i < fCalifaNumPetals - 1; i++)
            {
                cCalifa3->cd(i + 1);
                fh_Califa_energy_per_petal[i]->Draw();
            }

            for (Int_t i = 0; i < 2; i++)
            {
                cCalifa3b->cd(i + 1);
                fh_Califa_energy_per_petal[7 + i]->Draw();
            }

            for (Int_t i = 0; i < fCalifaNumPetals; i++)
            {
                for (Int_t k = 0; k < 4; k++)
                {
                    for (Int_t j = 0; j < 16; j++)
                    {
                        cCalifa4[i][k]->cd(fOrderFebexPreamp[j] + 1);
                        gPad->SetLogy();
                        fh_Califa_crystals[i][j + 16 * k]->SetFillColor(45);
                        if (i == 7 && k < 2)
                            sprintf(Name,
                                    "Map level, Petal %d, Preamp %d, Preamp ch. %d",
                                    1,
                                    k + 1,
                                    fOrderFebexPreamp[j] + 1);
                        else if (i == 7 && k > 1)
                            sprintf(Name,
                                    "Map level, Petal %d, Preamp %d, Preamp ch. %d",
                                    4,
                                    k - 1,
                                    fOrderFebexPreamp[j] + 1);
                        else
                            sprintf(Name,
                                    "Map level, Petal %d, Preamp %d, Preamp ch. %d",
                                    i + 1,
                                    k + 1,
                                    fOrderFebexPreamp[j] + 1);
                        fh_Califa_crystals[i][j + 16 * k]->SetTitle(Name);
                        fh_Califa_crystals[i][j + 16 * k]->Draw();
                    }
                }
            }

            fRaw2Cal = kFALSE;
        }
        else
        {

            cCalifa2->cd();
            fh_Califa_cryId_energy_cal->Draw("COLZ");

            for (Int_t i = 0; i < fCalifaNumPetals - 1; i++)
            {
                cCalifa3->cd(i + 1);
                fh_Califa_energy_per_petal_cal[i]->Draw();
            }

            for (Int_t i = 0; i < 2; i++)
            {
                cCalifa3b->cd(i + 1);
                fh_Califa_energy_per_petal_cal[7 + i]->Draw();
            }

            for (Int_t i = 0; i < fCalifaNumPetals; i++)
            {
                for (Int_t k = 0; k < 4; k++)
                {
                    for (Int_t j = 0; j < 16; j++)
                    {
                        cCalifa4[i][k]->cd(fOrderFebexPreamp[j] + 1);
                        gPad->SetLogy();
                        fh_Califa_crystals_cal[i][j + 16 * k]->SetFillColor(45);
                        if (i == 7 && k < 2)
                            sprintf(Name,
                                    "Cal level, Petal %d, Preamp %d, Preamp ch. %d",
                                    1,
                                    k + 1,
                                    fOrderFebexPreamp[j] + 1);
                        else if (i == 7 && k > 1)
                            sprintf(Name,
                                    "Cal level, Petal %d, Preamp %d, Preamp ch. %d",
                                    4,
                                    k - 1,
                                    fOrderFebexPreamp[j] + 1);
                        else
                            sprintf(Name,
                                    "Cal level, Petal %d, Preamp %d, Preamp ch. %d",
                                    i + 1,
                                    k + 1,
                                    fOrderFebexPreamp[j] + 1);
                        fh_Califa_crystals_cal[i][j + 16 * k]->SetTitle(Name);
                        fh_Califa_crystals_cal[i][j + 16 * k]->Draw();
                    }
                }
            }

            fRaw2Cal = kTRUE;
        }
    }
}

void R3BCalifaDemoOnlineSpectra::Febex2Preamp_CALIFA_Histo()
{
    LOG(info) << "R3BCalifaDemoOnlineSpectra::Febex2Preamp_CALIFA_Histo";

    char Name[255];

    // Raw data
    if (!fRaw2Cal)
    {
        if (fFebex2Preamp)
        { // Febex to Preamp sequence

            for (Int_t i = 0; i < fCalifaNumPetals; i++)
            {
                for (Int_t k = 0; k < 4; k++)
                {
                    for (Int_t j = 0; j < 16; j++)
                    {
                        cCalifa4[i][k]->cd(fOrderFebexPreamp[j] + 1);
                        gPad->SetLogy();
                        fh_Califa_crystals[i][j + 16 * k]->SetFillColor(45);
                        if (i == 7 && k < 2)
                            sprintf(Name,
                                    "Map level, Petal %d, Preamp %d, Preamp ch. %d",
                                    1,
                                    k + 1,
                                    fOrderFebexPreamp[j] + 1);
                        else if (i == 7 && k > 1)
                            sprintf(Name,
                                    "Map level, Petal %d, Preamp %d, Preamp ch. %d",
                                    4,
                                    k - 1,
                                    fOrderFebexPreamp[j] + 1);
                        else
                            sprintf(Name,
                                    "Map level, Petal %d, Preamp %d, Preamp ch. %d",
                                    i + 1,
                                    k + 1,
                                    fOrderFebexPreamp[j] + 1);
                        fh_Califa_crystals[i][j + 16 * k]->SetTitle(Name);
                        fh_Califa_crystals[i][j + 16 * k]->Draw();
                    }
                }
            }
            fFebex2Preamp = kFALSE;
        }
        else
        {

            for (Int_t i = 0; i < fCalifaNumPetals; i++)
            {
                for (Int_t k = 0; k < 4; k++)
                {
                    for (Int_t j = 0; j < 16; j++)
                    {
                        cCalifa4[i][k]->cd(j + 1);
                        gPad->SetLogy();
                        fh_Califa_crystals[i][j + 16 * k]->SetFillColor(kGreen);
                        if (i == 7 && k < 2)
                            sprintf(Name, "Map level, Petal %d, Preamp %d, Febex ch. %d", 1, k + 1, j);
                        else if (i == 7 && k > 1)
                            sprintf(Name, "Map level, Petal %d, Preamp %d, Febex ch. %d", 4, k - 1, j);
                        else
                            sprintf(Name, "Map level, Petal %d, Preamp %d, Febex ch. %d", i + 1, k + 1, j);
                        fh_Califa_crystals[i][j + 16 * k]->SetTitle(Name);
                        fh_Califa_crystals[i][j + 16 * k]->Draw();
                    }
                }
            }
            fFebex2Preamp = kTRUE;
        }
    }
    else
    { // Cal data

        if (fFebex2Preamp)
        { // Febex to Preamp sequence

            for (Int_t i = 0; i < fCalifaNumPetals; i++)
            {
                for (Int_t k = 0; k < 4; k++)
                {
                    for (Int_t j = 0; j < 16; j++)
                    {
                        cCalifa4[i][k]->cd(fOrderFebexPreamp[j] + 1);
                        gPad->SetLogy();
                        fh_Califa_crystals_cal[i][j + 16 * k]->SetFillColor(45);
                        if (i == 7 && k < 2)
                            sprintf(Name,
                                    "Cal level, Petal %d, Preamp %d, Preamp ch. %d",
                                    1,
                                    k + 1,
                                    fOrderFebexPreamp[j] + 1);
                        else if (i == 7 && k > 1)
                            sprintf(Name,
                                    "Cal level, Petal %d, Preamp %d, Preamp ch. %d",
                                    4,
                                    k - 1,
                                    fOrderFebexPreamp[j] + 1);
                        else
                            sprintf(Name,
                                    "Cal level, Petal %d, Preamp %d, Preamp ch. %d",
                                    i + 1,
                                    k + 1,
                                    fOrderFebexPreamp[j] + 1);
                        fh_Califa_crystals_cal[i][j + 16 * k]->SetTitle(Name);
                        fh_Califa_crystals_cal[i][j + 16 * k]->Draw();
                    }
                }
            }
            fFebex2Preamp = kFALSE;
        }
        else
        {

            for (Int_t i = 0; i < fCalifaNumPetals; i++)
            {
                for (Int_t k = 0; k < 4; k++)
                {
                    for (Int_t j = 0; j < 16; j++)
                    {
                        cCalifa4[i][k]->cd(j + 1);
                        gPad->SetLogy();
                        fh_Califa_crystals_cal[i][j + 16 * k]->SetFillColor(kGreen);
                        if (i == 7 && k < 2)
                            sprintf(Name, "Cal level, Petal %d, Preamp %d, Febex ch. %d", 1, k + 1, j);
                        else if (i == 7 && k > 1)
                            sprintf(Name, "Cal level, Petal %d, Preamp %d, Febex ch. %d", 4, k - 1, j);
                        else
                            sprintf(Name, "Cal level, Petal %d, Preamp %d, Febex ch. %d", i + 1, k + 1, j);
                        fh_Califa_crystals_cal[i][j + 16 * k]->SetTitle(Name);
                        fh_Califa_crystals_cal[i][j + 16 * k]->Draw();
                    }
                }
            }
            fFebex2Preamp = kTRUE;
        }
    }
}

void R3BCalifaDemoOnlineSpectra::Exec(Option_t* option)
{

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BCalifaDemoOnlineSpectra::Exec FairRootManager not found";

    uint64_t wrc = 0;
    if (fWRItemsCalifa && fWRItemsCalifa->GetEntriesFast())
    {
        Int_t nHits = fWRItemsCalifa->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsCalifa->At(ihit));
            if (!hit)
                continue;
            wrc = hit->GetTimeStamp();
        }
    }
    uint64_t wrm = 0;
    if (fWRItemsMaster && fWRItemsMaster->GetEntriesFast())
    {
        Int_t nHits = fWRItemsMaster->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsMaster->At(ihit));
            if (!hit)
                continue;
            wrm = hit->GetTimeStamp();
        }
    }
    if (wrc != 0 || wrm != 0)
    {
        fh_Califa_wr->Fill(wrc - wrm);
    }

    Int_t counter[8];
    if (fMappedItemsCalifa && fMappedItemsCalifa->GetEntriesFast())
    {
        Int_t nHits = fMappedItemsCalifa->GetEntriesFast();

        fh_Califa_Mult[0]->Fill(nHits);

        for (Int_t ccoo = 0; ccoo < 8; ccoo++)
            counter[ccoo] = 0;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaMappedData* hit = dynamic_cast<R3BCalifaMappedData*>(fMappedItemsCalifa->At(ihit));
            if (!hit)
                continue;

            Int_t cryId = Map_For_s444(hit->GetCrystalId());
            Int_t petal = (Int_t)(cryId) / fNumCrystalPetal;        // from 0 to 6 (s444)
            Int_t cryId_petal = cryId - fNumCrystalPetal * (petal); // from 0 to 63

            counter[petal]++;

            fh_Califa_cryId_energy->Fill(cryId + 1, hit->GetEnergy());
            // individual energy histo for each crystalId
            fh_Califa_crystals[petal][cryId_petal]->Fill(hit->GetEnergy());
            // energy(channels) sum for each petal
            if (petal < 7)
                fh_Califa_energy_per_petal[petal]->Fill(hit->GetEnergy()); // s444
            if (petal == 7 && cryId > 447 && cryId < 480)
            {
                fh_Califa_energy_per_petal[petal]->Fill(hit->GetEnergy());
            }
            if (petal == 7 && cryId > 479)
            {
                fh_Califa_energy_per_petal[8]->Fill(hit->GetEnergy());
            }
            // crystalId vs petal number
            fh_Califa_cryId_petal->Fill(cryId_petal + 1, petal + 1);
        }
        for (Int_t coo = 0; coo < 8; coo++)
            fh_Califa_Mult[coo + 1]->Fill(counter[coo]);
    }

    if (fCalON == kTRUE && fCalItemsCalifa && fCalItemsCalifa->GetEntriesFast())
    {
        Int_t nHits = fCalItemsCalifa->GetEntriesFast();

        Double_t E0 = 0., E1 = 0., E2 = 0., E3 = 0., E4 = 0., E5 = 0., E6 = 0.;
        Double_t E0p = 0., E1p = 0.; // proton range for petal 1 vs 3and4, dual preamps

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaCrystalCalData* hit = dynamic_cast<R3BCalifaCrystalCalData*>(fCalItemsCalifa->At(ihit));
            if (!hit)
                continue;

            Int_t cryId = Map_For_s444(hit->GetCrystalId());
            Int_t petal = (Int_t)(cryId) / fNumCrystalPetal;        // from 0 to 7
            Int_t cryId_petal = cryId - fNumCrystalPetal * (petal); // from 0 to 63

            fh_Califa_cryId_energy_cal->Fill(cryId + 1, hit->GetEnergy());
            // individual energy histo for each crystalId
            fh_Califa_crystals_cal[petal][cryId_petal]->Fill(hit->GetEnergy());
            // energy(channels) sum for each petal
            if (petal < 7)
                fh_Califa_energy_per_petal_cal[petal]->Fill(hit->GetEnergy()); // s444
            // s444
            if (petal == 0)
                E0 = E0 + hit->GetEnergy();
            if (petal == 1)
                E1 = E1 + hit->GetEnergy();
            if (petal == 2)
                E2 = E2 + hit->GetEnergy();
            if (petal == 3)
                E3 = E3 + hit->GetEnergy();
            if (petal == 4)
                E4 = E4 + hit->GetEnergy();
            if (petal == 5)
                E5 = E5 + hit->GetEnergy();
            if (petal == 6)
                E6 = E6 + hit->GetEnergy();
            if (petal == 7 && cryId > 447 && cryId < 480)
            {
                E0p = E0p + hit->GetEnergy();
                fh_Califa_energy_per_petal_cal[petal]->Fill(hit->GetEnergy());
            }
            if (petal == 7 && cryId > 479)
            {
                E1p = E1p + hit->GetEnergy();
                fh_Califa_energy_per_petal_cal[8]->Fill(hit->GetEnergy());
            }
        }
        if (E0 > 0 && E2 + E3 > 0)
            fh_Califa_coinc_petal1->Fill(E0, E2 + E3);
        if (E1 > 0 && E4 > 0)
            fh_Califa_coinc_petal2->Fill(E1, E4);
        if (E0p > 0 && E1p > 0)
            fh_Califa_coinc_petal3->Fill(E0p, E1p);
    }

    Int_t counterHits[8];
    if (fHitItemsCalifa && fHitItemsCalifa->GetEntriesFast())
    {
        Int_t nHits = fHitItemsCalifa->GetEntriesFast();
        fh_Califa_MultHit[0]->Fill(nHits);
        for (Int_t ccoo = 0; ccoo < 8; ccoo++)
            counterHits[ccoo] = 0;

        Double_t theta = 0., phi = 0.;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaClusterData* hit = dynamic_cast<R3BCalifaClusterData*>(fHitItemsCalifa->At(ihit));
            if (!hit)
                continue;
            theta = hit->GetTheta() / TMath::Pi() * 180.;
            phi = hit->GetPhi() / TMath::Pi() * 180.;
            if (theta < 90.)
            {
                fh_Califa_theta_phi->Fill(theta, phi);
            }
            else
            {
                fh_Califa_theta_phi->Fill(theta - 90., phi);
            }
            if (phi > 120)
            { // petal 1
                if (theta < 90)
                    fh_Califa_theta_energy[0]->Fill(theta, hit->GetEnergy());
                else
                    fh_Califa_theta_energy[0]->Fill(theta - 90., hit->GetEnergy());
                counterHits[0]++;
            }
            else if (phi < -120)
            { // petal 2
                fh_Califa_theta_energy[1]->Fill(theta, hit->GetEnergy());
                counterHits[1]++;
            }
            else if (phi > -100 && phi < -52. && theta < 90)
            { // petal 3
                fh_Califa_theta_energy[2]->Fill(theta, hit->GetEnergy());
                counterHits[2]++;
            }
            else if (phi >= -52 && phi < -10 && theta < 90)
            { // petal 4
                fh_Califa_theta_energy[3]->Fill(theta, hit->GetEnergy());
                counterHits[3]++;
            }
            else if (phi >= -70 && phi < -10 && theta > 90)
            { // petal 4 dual preamp
                fh_Califa_theta_energy[3]->Fill(theta - 90., hit->GetEnergy());
            }
            else if (phi < 45 && phi > 10)
            { // petal 5
                fh_Califa_theta_energy[4]->Fill(theta, hit->GetEnergy());
                counterHits[4]++;
            }
            else if (phi < 73 && phi > 45)
            { // petal 6
                fh_Califa_theta_energy[5]->Fill(theta, hit->GetEnergy());
                counterHits[5]++;
            }
            else if (phi > 73 && phi < 120)
            { // petal 7
                fh_Califa_theta_energy[6]->Fill(theta, hit->GetEnergy());
                counterHits[6]++;
            }
            // total
            fh_Califa_theta_energy[7]->Fill(theta, hit->GetEnergy());
            fh_Califa_total_energy->Fill(hit->GetEnergy());
        }
        for (Int_t coo = 0; coo < 7; coo++)
            fh_Califa_MultHit[coo + 1]->Fill(counterHits[coo]);
    }

    fNEvents += 1;
}

void R3BCalifaDemoOnlineSpectra::FinishEvent()
{
    if (fMappedItemsCalifa)
    {
        fMappedItemsCalifa->Clear();
    }
    if (fCalItemsCalifa)
    {
        fCalItemsCalifa->Clear();
    }
    if (fHitItemsCalifa)
    {
        fHitItemsCalifa->Clear();
    }
    if (fWRItemsCalifa)
    {
        fWRItemsCalifa->Clear();
    }
    if (fWRItemsMaster)
    {
        fWRItemsMaster->Clear();
    }
}

void R3BCalifaDemoOnlineSpectra::FinishTask()
{

    if (fWRItemsCalifa && fWRItemsMaster)
    {
        cCalifa_wr->Write();
    }

    if (fMappedItemsCalifa)
    {
        fh_Califa_cryId_energy->Write();
        cCalifa1->Write();
        cCalifa2->Write();
        cCalifa3->Write();
        cCalifa3b->Write();
        for (Int_t i = 0; i < fCalifaNumPetals; i++)
        {
            for (Int_t k = 0; k < 4; k++)
            {
                cCalifa4[i][k]->Write();
            }
        }
        cCalifa5->Write();
        cCalifa6->Write();
        cCalifa7->Write();
        cCalifa10->Write();
    }

    // mapped data
    if (fMappedItemsCalifa)
    {
        for (Int_t i = 0; i < fCalifaNumPetals; i++)
        {
            for (Int_t k = 0; k < 4; k++)
            {
                for (Int_t j = 0; j < 16; j++)
                {
                    fh_Califa_crystals[i][j + 16 * k]->Write();
                }
            }
        }
    }
    // cal data
    if (fCalItemsCalifa)
    {
        for (Int_t i = 0; i < fCalifaNumPetals; i++)
        {
            for (Int_t k = 0; k < 4; k++)
            {
                for (Int_t j = 0; j < 16; j++)
                {
                    fh_Califa_crystals_cal[i][j + 16 * k]->Write();
                }
            }
        }
    }

    // hit data
    if (fHitItemsCalifa)
    {
        cCalifa8->Write();
        cCalifa_hitenergy->Write();
        for (Int_t i = 0; i < fCalifaNumPetals - 1; i++)
            cCalifa_hitpetal[i]->Write();
    }
}

Int_t R3BCalifaDemoOnlineSpectra::Map_For_s444(Int_t val)
{
    // With this mapping all the code can remain unchanged
    // with the new software mapping for s444 (after 12 Feb 2019).
    // This is only valid for the s444 experiment
    // (but this is the case for all the code in this class!!)
    Int_t id_1[512] = {
        1142, 1144, 1138, 1140, 1141, 1143, 1137, 1139, 1014, 1016, 1010, 1012, 1013, 1015, 1009, 1011, 1398, 1400,
        1394, 1396, 1397, 1399, 1393, 1395, 1270, 1272, 1266, 1268, 1269, 1271, 1265, 1267, 1654, 1655, 1650, 1651,
        1653, 1656, 1649, 1652, 1526, 1528, 1522, 1524, 1525, 1527, 1521, 1523, 6654, 6655, 6650, 6651, 6653, 6656,
        6649, 6652, 6526, 6528, 6522, 6524, 6525, 6527, 6521, 6523, 1910, 1911, 1906, 1907, 1909, 1912, 1905, 1908,
        1782, 1783, 1778, 1779, 1781, 1784, 1777, 1780, 6910, 6911, 6906, 6907, 6909, 6912, 6905, 6908, 6782, 6783,
        6778, 6779, 6781, 6784, 6777, 6780, 1126, 1128, 1122, 1124, 1125, 1127, 1121, 1123, 998,  1000, 994,  996,
        997,  999,  993,  995,  1382, 1384, 1378, 1380, 1381, 1383, 1377, 1379, 1254, 1256, 1250, 1252, 1253, 1255,
        1249, 1251, 1638, 1639, 1634, 1635, 1637, 1640, 1633, 1636, 1510, 1512, 1506, 1508, 1509, 1511, 1505, 1507,
        1894, 1895, 1890, 1891, 1893, 1896, 1889, 1892, 1766, 1767, 1762, 1763, 1765, 1768, 1761, 1764, 1102, 1104,
        1098, 1100, 1101, 1103, 1097, 1099, 974,  976,  970,  972,  973,  975,  969,  971,  1358, 1360, 1354, 1356,
        1357, 1359, 1353, 1355, 1230, 1232, 1226, 1228, 1229, 1231, 1225, 1227, 1614, 1615, 1610, 1611, 1613, 1616,
        1609, 1612, 1486, 1488, 1482, 1484, 1485, 1487, 1481, 1483, 1870, 1871, 1866, 1867, 1869, 1872, 1865, 1868,
        1742, 1743, 1738, 1739, 1741, 1744, 1737, 1740, 1094, 1096, 1090, 1092, 1093, 1095, 1089, 1091, 966,  968,
        962,  964,  965,  967,  961,  963,  1350, 1352, 1346, 1348, 1349, 1351, 1345, 1347, 1222, 1224, 1218, 1220,
        1221, 1223, 1217, 1219, 1606, 1607, 1602, 1603, 1605, 1608, 1601, 1604, 1478, 1480, 1474, 1476, 1477, 1479,
        1473, 1475, 6606, 6607, 6602, 6603, 6605, 6608, 6601, 6604, 6478, 6480, 6474, 6476, 6477, 6479, 6473, 6475,
        1862, 1863, 1858, 1859, 1861, 1864, 1857, 1860, 1734, 1735, 1730, 1731, 1733, 1736, 1729, 1732, 6862, 6863,
        6858, 6859, 6861, 6864, 6857, 6860, 6734, 6735, 6730, 6731, 6733, 6736, 6729, 6732, 1078, 1080, 1074, 1076,
        1077, 1079, 1073, 1075, 950,  952,  946,  948,  949,  951,  945,  947,  1334, 1336, 1330, 1332, 1333, 1335,
        1329, 1331, 1206, 1208, 1202, 1204, 1205, 1207, 1201, 1203, 1590, 1591, 1586, 1587, 1589, 1592, 1585, 1588,
        1462, 1464, 1458, 1460, 1461, 1463, 1457, 1459, 1846, 1847, 1842, 1843, 1845, 1848, 1841, 1844, 1718, 1719,
        1714, 1715, 1717, 1720, 1713, 1716, 1070, 1072, 1066, 1068, 1069, 1071, 1065, 1067, 942,  944,  938,  940,
        941,  943,  937,  939,  1326, 1328, 1322, 1324, 1325, 1327, 1321, 1323, 1198, 1200, 1194, 1196, 1197, 1199,
        1193, 1195, 1582, 1583, 1578, 1579, 1581, 1584, 1577, 1580, 1454, 1456, 1450, 1452, 1453, 1455, 1449, 1451,
        1838, 1839, 1834, 1835, 1837, 1840, 1833, 1836, 1710, 1711, 1706, 1707, 1709, 1712, 1705, 1708, 1062, 1064,
        1058, 1060, 1061, 1063, 1057, 1059, 934,  936,  930,  932,  933,  935,  929,  931,  1318, 1320, 1314, 1316,
        1317, 1319, 1313, 1315, 1190, 1192, 1186, 1188, 1189, 1191, 1185, 1187, 1574, 1575, 1570, 1571, 1573, 1576,
        1569, 1572, 1446, 1448, 1442, 1444, 1445, 1447, 1441, 1443, 1830, 1831, 1826, 1827, 1829, 1832, 1825, 1828,
        1702, 1703, 1698, 1699, 1701, 1704, 1697, 1700
    };
    Int_t id_2[512] = {
        246, 245, 244, 243, 242, 241, 240, 247, 248, 255, 254, 253, 252, 251, 250, 249, 230, 229, 228, 227, 226, 225,
        224, 231, 232, 239, 238, 237, 236, 235, 234, 233, 214, 213, 212, 211, 210, 209, 208, 215, 216, 223, 222, 221,
        220, 219, 218, 217, 198, 197, 196, 195, 194, 193, 192, 199, 200, 207, 206, 205, 204, 203, 202, 201, 182, 181,
        180, 179, 178, 177, 176, 183, 184, 191, 190, 189, 188, 187, 186, 185, 166, 165, 164, 163, 162, 161, 160, 167,
        168, 175, 174, 173, 172, 171, 170, 169, 502, 501, 500, 499, 498, 497, 496, 503, 504, 511, 510, 509, 508, 507,
        506, 505, 486, 485, 484, 483, 482, 481, 480, 487, 488, 495, 494, 493, 492, 491, 490, 489, 470, 469, 468, 467,
        466, 465, 464, 471, 472, 479, 478, 477, 476, 475, 474, 473, 454, 453, 452, 451, 450, 449, 448, 455, 456, 463,
        462, 461, 460, 459, 458, 457, 54,  53,  52,  51,  50,  49,  48,  55,  56,  63,  62,  61,  60,  59,  58,  57,
        38,  37,  36,  35,  34,  33,  32,  39,  40,  47,  46,  45,  44,  43,  42,  41,  22,  21,  20,  19,  18,  17,
        16,  23,  24,  31,  30,  29,  28,  27,  26,  25,  6,   5,   4,   3,   2,   1,   0,   7,   8,   15,  14,  13,
        12,  11,  10,  9,   150, 149, 148, 147, 146, 145, 144, 151, 152, 159, 158, 157, 156, 155, 154, 153, 134, 133,
        132, 131, 130, 129, 128, 135, 136, 143, 142, 141, 140, 139, 138, 137, 118, 117, 116, 115, 114, 113, 112, 119,
        120, 127, 126, 125, 124, 123, 122, 121, 102, 101, 100, 99,  98,  97,  96,  103, 104, 111, 110, 109, 108, 107,
        106, 105, 86,  85,  84,  83,  82,  81,  80,  87,  88,  95,  94,  93,  92,  91,  90,  89,  70,  69,  68,  67,
        66,  65,  64,  71,  72,  79,  78,  77,  76,  75,  74,  73,  310, 309, 308, 307, 306, 305, 304, 311, 312, 319,
        318, 317, 316, 315, 314, 313, 294, 293, 292, 291, 290, 289, 288, 295, 296, 303, 302, 301, 300, 299, 298, 297,
        278, 277, 276, 275, 274, 273, 272, 279, 280, 287, 286, 285, 284, 283, 282, 281, 262, 261, 260, 259, 258, 257,
        256, 263, 264, 271, 270, 269, 268, 267, 266, 265, 374, 373, 372, 371, 370, 369, 368, 375, 376, 383, 382, 381,
        380, 379, 378, 377, 358, 357, 356, 355, 354, 353, 352, 359, 360, 367, 366, 365, 364, 363, 362, 361, 342, 341,
        340, 339, 338, 337, 336, 343, 344, 351, 350, 349, 348, 347, 346, 345, 326, 325, 324, 323, 322, 321, 320, 327,
        328, 335, 334, 333, 332, 331, 330, 329, 438, 437, 436, 435, 434, 433, 432, 439, 440, 447, 446, 445, 444, 443,
        442, 441, 422, 421, 420, 419, 418, 417, 416, 423, 424, 431, 430, 429, 428, 427, 426, 425, 406, 405, 404, 403,
        402, 401, 400, 407, 408, 415, 414, 413, 412, 411, 410, 409, 390, 389, 388, 387, 386, 385, 384, 391, 392, 399,
        398, 397, 396, 395, 394, 393
    };
    Int_t nbcry = 0;
    for (Int_t i = 0; i < 512; i++)
    {
        if (id_1[i] == val)
            nbcry = id_2[i];
    }

    // New arrangement due to dual preamps
    if (nbcry >= 0 && nbcry < 64)
    {
        nbcry = nbcry + 128; // USC1: petal 3
    }
    else if (nbcry > 63 && nbcry < 80)
    {
        nbcry = nbcry - 64 + 448 + 32; // USC2: proton range, preamp 1, section 1
    }
    else if (nbcry > 79 && nbcry < 96)
    {
        nbcry = nbcry - 80 + 192; // USC2: gamma range, preamp 1, section 1
    }
    else if (nbcry > 95 && nbcry < 112)
    {
        nbcry = nbcry - 96 + 448 + 48; // USC2: proton range, preamp 1, section 2
    }
    else if (nbcry > 111 && nbcry < 128)
    {
        nbcry = nbcry - 112 + 192 + 16; // USC2: gamma range, preamp 1, section 2
    }
    else if (nbcry > 127 && nbcry < 144)
    {
        nbcry = nbcry - 128 + 192 + 32; // USC2: preamp 2, section 1
    }
    else if (nbcry > 143 && nbcry < 160)
    {
        nbcry = nbcry - 144 + 192 + 48; // USC2: preamp 2, section 2
    }
    else if (nbcry > 159 && nbcry < 176)
    {
        nbcry = nbcry - 160 + 448; // TUDA1: proton range, preamp 1, section 1
    }
    else if (nbcry > 175 && nbcry < 192)
    {
        nbcry = nbcry - 176; // TUDA1: gamma range, preamp 1, section 1
    }
    else if (nbcry > 191 && nbcry < 208)
    {
        nbcry = nbcry - 192 + 448 + 16; // TUDA1: proton range, preamp 1, section 2
    }
    else if (nbcry > 207 && nbcry < 224)
    {
        nbcry = nbcry - 208 + 16; // TUDA1: gamma range, preamp 1, section 2
    }
    else if (nbcry > 223 && nbcry < 240)
    {
        nbcry = nbcry - 224 + 32; // TUDA1: preamp 2, section 1
    }
    else if (nbcry > 239 && nbcry < 256)
    {
        nbcry = nbcry - 240 + 48; // TUDA1: preamp 2, section 2
    }
    else if (nbcry > 447)
    {
        nbcry = nbcry - 448 + 64; // TUDA2: petal 2
    }

    return nbcry;
}

ClassImp(R3BCalifaDemoOnlineSpectra)
