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
// -----                  R3BAmsOnlineSpectra             -----
// -----    Created 22/07/18  by J.L. Rodriguez-Sanchez   -----
// -----           Fill AMS online histograms             -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with AMS online data
 */

#include "R3BAmsOnlineSpectra.h"
#include "R3BAmsHitData.h"
#include "R3BAmsMappedData.h"
#include "R3BAmsMappingPar.h"
#include "R3BAmsStripCalData.h"
#include "R3BEventHeader.h"
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
#include "TRandom.h"

#include "TClonesArray.h"
#include "TMath.h"
#include "TVector3.h"
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

R3BAmsOnlineSpectra::R3BAmsOnlineSpectra()
    : FairTask("AmsOnlineSpectra", 1)
    , fMappedItemsAms(NULL)
    , fCalItemsAms(NULL)
    , fHitItemsAms(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(6)
{
}

R3BAmsOnlineSpectra::R3BAmsOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItemsAms(NULL)
    , fCalItemsAms(NULL)
    , fHitItemsAms(NULL)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(6)
{
}

R3BAmsOnlineSpectra::~R3BAmsOnlineSpectra()
{
    LOG(info) << "R3BAmsOnlineSpectra: Delete instance";
    if (fMappedItemsAms)
        delete fMappedItemsAms;
    if (fCalItemsAms)
        delete fCalItemsAms;
    if (fHitItemsAms)
        delete fHitItemsAms;
}

void R3BAmsOnlineSpectra::SetParContainers()
{
    // Parameter Container
    // Reading amsStripCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "FairRuntimeDb not opened!";
    }

    fMap_Par = dynamic_cast<R3BAmsMappingPar*>(rtdb->getContainer("amsMappingPar"));
    if (!fMap_Par)
    {
        LOG(error) << "R3BAmsOnlineSpectra::Couldn't get handle on amsMappingPar container";
    }
    else
    {
        LOG(info) << "R3BAmsOnlineSpectra::amsMappingPar container open";
    }
}

void R3BAmsOnlineSpectra::SetParameter()
{
    if (!fMap_Par)
    {
        LOG(error) << "R3BAmsOnlineSpectra::amsMappingPar container not found";
    }
    //--- Parameter Container ---
    fNbDet = fMap_Par->GetNumDets(); // Number of detectors
    LOG(info) << "R3BAmsOnlineSpectra::NumDets " << fNbDet;
}

InitStatus R3BAmsOnlineSpectra::Init()
{
    LOG(info) << "R3BAmsOnlineSpectra::Init()";

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BAmsOnlineSpectra::Init FairRootManager not found";

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // create histograms of all detectors

    // get access to Mapped data
    fMappedItemsAms = dynamic_cast<TClonesArray*>(mgr->GetObject("AmsMappedData"));
    if (!fMappedItemsAms)
    {
        LOG(fatal) << "R3BAmsOnlineSpectra::AmsMappedData not found";
        return kFATAL;
    }

    // get access to Cal data
    fCalItemsAms = dynamic_cast<TClonesArray*>(mgr->GetObject("AmsStripCalData"));
    if (!fCalItemsAms)
        LOG(warn) << "R3BAmsOnlineSpectra::AmsStripCalData not found";

    // get access to Hit data
    fHitItemsAms = dynamic_cast<TClonesArray*>(mgr->GetObject("AmsHitData"));
    if (!fHitItemsAms)
        LOG(warn) << "R3BAmsOnlineSpectra::AmsHitData not found";

    SetParameter();

    // Energy range for strips
    Double_t bins = 200;
    Double_t maxE = 3500;
    Double_t minE = 0;

    //  CANVAS 1  -------------------------------
    cMap = new TCanvas("AMS_mapped", "mapped info", 10, 10, 500, 500);
    cMap->Divide(fNbDet / 2, 2);

    //  CANVAS 2  -------------------------------
    cCalR = new TCanvas("AMS_cal_Right", "Cal info right-side", 10, 10, 500, 500);
    cCalR->Divide(2, fNbDet / 2);
    cCalL = new TCanvas("AMS_cal_Left", "Cal info left-side", 10, 10, 500, 500);
    cCalL->Divide(2, fNbDet / 2);

    if (fNbDet == 6)
    {
        //  CANVAS 3  -------------------------------
        cHit[0] = new TCanvas("AMS_hit_right_layer1", "Hit info", 10, 10, 500, 500);
        cHit[0]->Divide(2, 2);
        cHit[1] = new TCanvas("AMS_hit_right_layer2_top", "Hit info", 10, 10, 500, 500);
        cHit[1]->Divide(2, 2);
        cHit[2] = new TCanvas("AMS_hit_right_layer2_bottom", "Hit info", 10, 10, 500, 500);
        cHit[2]->Divide(2, 2);
        cHit[3] = new TCanvas("AMS_hit_left_layer1", "Hit info", 10, 10, 500, 500);
        cHit[3]->Divide(2, 2);
        cHit[4] = new TCanvas("AMS_hit_left_layer2_top", "Hit info", 10, 10, 500, 500);
        cHit[4]->Divide(2, 2);
        cHit[5] = new TCanvas("AMS_hit_left_layer2_bottom", "Hit info", 10, 10, 500, 500);
        cHit[5]->Divide(2, 2);

        //  CANVAS 4  -------------------------------
        cHitAngles = new TCanvas("AMS_Theta_vs_Phi", "Theta vs Phi per layer", 10, 10, 500, 500);
        cHitAngles->Divide(2, 1);

        //  CANVAS 5  -------------------------------
        cHitEnergyCor =
            new TCanvas("AMS_E-inner_vs_E-outer", "Energy-inner vs Energy-outer per side", 10, 10, 500, 500);
        cHitEnergyCor->Divide(2, 1);

        //  CANVAS 6  -------------------------------
        cHitThetaCor =
            new TCanvas("AMS_Theta-inner_vs_Theta-outer", "Theta-inner vs Theta-outer per side", 10, 10, 500, 500);
        cHitThetaCor->Divide(2, 2);

        //  CANVAS 7  -------------------------------
        cHitPhiCor = new TCanvas("AMS_Phi-inner_vs_Phi-outer", "Phi-inner vs Phi-outer per side", 10, 10, 500, 500);
        cHitPhiCor->Divide(2, 1);
    }
    else
    {
        //  CANVAS 3  -------------------------------
        cHit[0] = new TCanvas("AMS_hit_left", "Hit info", 10, 10, 500, 500);
        cHit[0]->Divide(2, 2);
        cHit[1] = new TCanvas("AMS_hit_top", "Hit info", 10, 10, 500, 500);
        cHit[1]->Divide(2, 2);
        cHit[2] = new TCanvas("AMS_hit_bottom", "Hit info", 10, 10, 500, 500);
        cHit[2]->Divide(2, 2);
        cHit[3] = new TCanvas("AMS_hit_right", "Hit info", 10, 10, 500, 500);
        cHit[3]->Divide(2, 2);
    }

    char Name1[255];
    char Name2[255];

    // Mapped data
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh_Ams_energy_allStrips_%d", i + 1);
        sprintf(Name2, "Energy vs strip number for AMS Det: %d", i + 1);
        fh_Ams_energy_allStrips[i] = new TH2F(Name1, Name2, 1024, 0, 1024, bins, minE, maxE);
        fh_Ams_energy_allStrips[i]->GetXaxis()->SetTitle("Strip number");
        fh_Ams_energy_allStrips[i]->GetYaxis()->SetTitle("Energy [channels]");
        fh_Ams_energy_allStrips[i]->GetYaxis()->SetTitleOffset(1.4);
        fh_Ams_energy_allStrips[i]->GetXaxis()->CenterTitle(true);
        fh_Ams_energy_allStrips[i]->GetYaxis()->CenterTitle(true);
        cMap->cd(i + 1);
        fh_Ams_energy_allStrips[i]->Draw("col");
    }

    // Cal data right
    char Side[50];
    for (Int_t i = 0; i < fNbDet / 2; i++)
    { // two histo per detector
        for (Int_t j = 0; j < 2; j++)
        {
            if (j == 0)
                sprintf(Side, "S");
            else
                sprintf(Side, "K");
            sprintf(Name1, "fh_Ams_energy_allCalStrips_%d", i * 2 + j);
            sprintf(Name2, "Cal_Energy vs strip number for AMS Det %d and Side-%s", i + 1, Side);
            if (j == 0)
                fh_Ams_energy_allCalStrips[i * 2 + j] = new TH2F(Name1, Name2, 640, 0, 640, bins, minE, maxE);
            else
                fh_Ams_energy_allCalStrips[i * 2 + j] = new TH2F(Name1, Name2, 384, 0, 384, bins, minE, maxE);
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetXaxis()->SetTitle("Strip number");
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetYaxis()->SetTitle("Energy [channels]");
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetYaxis()->SetTitleSize(0.06);
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetYaxis()->SetTitleOffset(0.85);
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetXaxis()->CenterTitle(true);
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetYaxis()->CenterTitle(true);
            cCalR->cd(i * 2 + 1 + j);
            fh_Ams_energy_allCalStrips[i * 2 + j]->Draw("col");
        }
    }

    // Cal data left
    for (Int_t i = fNbDet / 2; i < fNbDet; i++)
    { // two histo per detector
        for (Int_t j = 0; j < 2; j++)
        {
            if (j == 0)
                sprintf(Side, "S");
            else
                sprintf(Side, "K");
            sprintf(Name1, "fh_Ams_energy_allCalStrips_%d", i * 2 + j);
            sprintf(Name2, "Cal_Energy vs strip number for AMS Det %d and Side-%s", i + 1, Side);
            if (j == 0)
                fh_Ams_energy_allCalStrips[i * 2 + j] = new TH2F(Name1, Name2, 640, 0, 640, bins, minE, maxE);
            else
                fh_Ams_energy_allCalStrips[i * 2 + j] = new TH2F(Name1, Name2, 384, 0, 384, bins, minE, maxE);
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetXaxis()->SetTitle("Strip number");
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetYaxis()->SetTitle("Energy [channels]");
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetYaxis()->SetTitleSize(0.06);
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetYaxis()->SetTitleOffset(0.85);
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetXaxis()->CenterTitle(true);
            fh_Ams_energy_allCalStrips[i * 2 + j]->GetYaxis()->CenterTitle(true);
            cCalL->cd((i - fNbDet / 2) * 2 + 1 + j);
            fh_Ams_energy_allCalStrips[i * 2 + j]->Draw("col");
        }
    }

    // Hit data
    for (Int_t i = 0; i < fNbDet; i++)
    { // one histo per detector
        sprintf(Name1, "fh_Ams_hit_Mul_%d", i + 1);
        sprintf(Name2, "Multiplicity Det %d (blue: S-side, red: K-side)", i + 1);
        fh_Ams_hit_Mul[i][0] = new TH1F(Name1, Name2, 50, -0.5, 49.5);
        fh_Ams_hit_Mul[i][0]->GetXaxis()->SetTitle("Multiplicity");
        fh_Ams_hit_Mul[i][0]->GetYaxis()->SetTitle("Counts");
        fh_Ams_hit_Mul[i][0]->GetXaxis()->CenterTitle(true);
        fh_Ams_hit_Mul[i][0]->GetYaxis()->CenterTitle(true);
        fh_Ams_hit_Mul[i][0]->GetYaxis()->SetTitleOffset(1.1);
        fh_Ams_hit_Mul[i][0]->GetXaxis()->SetTitleOffset(1.1);
        cHit[i]->cd(1);
        fh_Ams_hit_Mul[i][0]->Draw("");
        sprintf(Name1, "fh_Ams_hit_MulK_%d", i + 1);
        fh_Ams_hit_Mul[i][1] = new TH1F(Name1, "", 50, -0.5, 49.5);
        fh_Ams_hit_Mul[i][1]->SetLineColor(2);
        fh_Ams_hit_Mul[i][1]->Draw("same");

        sprintf(Name1, "fh_Ams_hit_Pos_%d", i + 1);
        sprintf(Name2, "S-side vs K-side for AMS Det %d", i + 1);
        fh_Ams_hit_Pos[i] = new TH2F(Name1, Name2, 350, 0, 70., 200, 0, 40.);
        fh_Ams_hit_Pos[i]->GetXaxis()->SetTitle("S-side [mm]");
        fh_Ams_hit_Pos[i]->GetYaxis()->SetTitle("K-side [mm]");
        fh_Ams_hit_Pos[i]->GetXaxis()->CenterTitle(true);
        fh_Ams_hit_Pos[i]->GetYaxis()->CenterTitle(true);
        fh_Ams_hit_Pos[i]->GetYaxis()->SetTitleOffset(1.1);
        fh_Ams_hit_Pos[i]->GetXaxis()->SetTitleOffset(1.1);
        cHit[i]->cd(2);
        fh_Ams_hit_Pos[i]->Draw("col");

        sprintf(Name1, "fh_Ams_hit_Energies_%d", i + 1);
        sprintf(Name2, "Energy_S vs Energy_K for AMS Det %d", i + 1);
        fh_Ams_hit_E[i] = new TH2F(Name1, Name2, 2000, 0, 1000., 2000, 0, 1000.);
        fh_Ams_hit_E[i]->GetXaxis()->SetTitle("Energy_S [ADC units]");
        fh_Ams_hit_E[i]->GetYaxis()->SetTitle("Energy_K [ADC units]");
        fh_Ams_hit_E[i]->GetXaxis()->CenterTitle(true);
        fh_Ams_hit_E[i]->GetYaxis()->CenterTitle(true);
        fh_Ams_hit_E[i]->GetYaxis()->SetTitleOffset(1.2);
        fh_Ams_hit_E[i]->GetXaxis()->SetTitleOffset(1.1);
        cHit[i]->cd(3);
        fh_Ams_hit_E[i]->Draw("col");

        sprintf(Name1, "fh_Ams_hit_Energy_vs_theta_%d", i + 1);
        sprintf(Name2, "Energy_S vs Theta for AMS Det %d", i + 1);
        fh_Ams_hit_E_theta[i] = new TH2F(Name1, Name2, 300, 0., 90., 500, 0., 1000.);
        fh_Ams_hit_E_theta[i]->GetXaxis()->SetTitle("Theta [degrees]");
        fh_Ams_hit_E_theta[i]->GetYaxis()->SetTitle("Energy_S [ADC units]");
        fh_Ams_hit_E_theta[i]->GetXaxis()->CenterTitle(true);
        fh_Ams_hit_E_theta[i]->GetYaxis()->CenterTitle(true);
        fh_Ams_hit_E_theta[i]->GetYaxis()->SetTitleOffset(1.1);
        fh_Ams_hit_E_theta[i]->GetXaxis()->SetTitleOffset(1.1);
        cHit[i]->cd(4);
        fh_Ams_hit_E_theta[i]->Draw("col");
    }

    // Experiment S455
    if (fNbDet == 6)
    {
        for (Int_t i = 0; i < 2; i++) // two layers: inner (1) and outer (2)
        {
            if (i == 0)
            {
                sprintf(Name1, "fh2_Ams_theta_vs_phi_layer_%d", i + 1);
                sprintf(Name2, "AMS: theta vs phi for layer %d (inner)", i + 1);
            }
            else
            {
                sprintf(Name1, "fh2_Ams_theta_vs_phi_layer_%d", i + 1);
                sprintf(Name2, "AMS: theta vs phi for layer %d (outer)", i + 1);
            }
            fh2_ams_theta_phi[i] = new TH2F(Name1, Name2, 90, 0., 90., 180, -180., 180.);
            fh2_ams_theta_phi[i]->GetXaxis()->SetTitle("Theta [degrees]");
            fh2_ams_theta_phi[i]->GetYaxis()->SetTitle("Phi [degrees]");
            fh2_ams_theta_phi[i]->GetXaxis()->SetTitleOffset(1.1);
            fh2_ams_theta_phi[i]->GetYaxis()->SetTitleOffset(1.1);
            fh2_ams_theta_phi[i]->GetXaxis()->CenterTitle(true);
            fh2_ams_theta_phi[i]->GetYaxis()->CenterTitle(true);
            cHitAngles->cd(i + 1);
            fh2_ams_theta_phi[i]->Draw("COLZ");
        }

        for (Int_t i = 0; i < 2; i++) // two sides: right (1) and left (2)
        {
            // Energy correlations per side
            if (i == 0)
            {
                sprintf(Name1, "fh2_Ams_Einner_vs_Eouter_right");
                sprintf(Name2, "AMS: E-inner vs E-outer for right side");
            }
            else
            {
                sprintf(Name1, "fh2_Ams_Einner_vs_Eouter_left");
                sprintf(Name2, "AMS: E-inner vs E-outer for left side");
            }
            fh2_ams_e1_e2[i] = new TH2F(Name1, Name2, 500, 0., 1000., 500, 0., 1000.);
            fh2_ams_e1_e2[i]->GetXaxis()->SetTitle("E-inner [ADC units]");
            fh2_ams_e1_e2[i]->GetYaxis()->SetTitle("E-outer [ADC units]");
            fh2_ams_e1_e2[i]->GetXaxis()->SetTitleOffset(1.1);
            fh2_ams_e1_e2[i]->GetYaxis()->SetTitleOffset(1.4);
            fh2_ams_e1_e2[i]->GetXaxis()->CenterTitle(true);
            fh2_ams_e1_e2[i]->GetYaxis()->CenterTitle(true);
            cHitEnergyCor->cd(i + 1);
            fh2_ams_e1_e2[i]->Draw("COLZ");

            // Theta correlations per side
            if (i == 0)
            {
                sprintf(Name1, "fh2_Ams_#thetainner_vs_#thetaouterup_right");
                sprintf(Name2, "AMS: #theta-inner vs #theta-outer-up for right side");
            }
            else
            {
                sprintf(Name1, "fh2_Ams_#thetainner_vs_#thetaouterup_left");
                sprintf(Name2, "AMS: #theta-inner vs #theta-outer-up for left side");
            }
            fh2_ams_theta1_theta2[i] = new TH2F(Name1, Name2, 90, 0., 90., 90, 0., 90.);
            fh2_ams_theta1_theta2[i]->GetXaxis()->SetTitle("#theta-inner [degrees]");
            fh2_ams_theta1_theta2[i]->GetYaxis()->SetTitle("#theta-outer [degrees]");
            fh2_ams_theta1_theta2[i]->GetXaxis()->SetTitleOffset(1.);
            fh2_ams_theta1_theta2[i]->GetYaxis()->SetTitleOffset(1.5);
            fh2_ams_theta1_theta2[i]->GetXaxis()->CenterTitle(true);
            fh2_ams_theta1_theta2[i]->GetYaxis()->CenterTitle(true);
            cHitThetaCor->cd(i * 2 + 1);
            fh2_ams_theta1_theta2[i]->Draw("COLZ");

            // Theta correlations per side
            if (i == 0)
            {
                sprintf(Name1, "fh2_Ams_#thetainner_vs_#thetaouterdown_right");
                sprintf(Name2, "AMS: #theta-inner vs #theta-outer-down for right side");
            }
            else
            {
                sprintf(Name1, "fh2_Ams_#thetainner_vs_#thetaouterdown_left");
                sprintf(Name2, "AMS: #theta-inner vs #theta-outer-down for left side");
            }
            fh2_ams_theta1_theta2[i + 2] = new TH2F(Name1, Name2, 90, 0., 90., 90, 0., 90.);
            fh2_ams_theta1_theta2[i + 2]->GetXaxis()->SetTitle("#theta-inner [degrees]");
            fh2_ams_theta1_theta2[i + 2]->GetYaxis()->SetTitle("#theta-outer [degrees]");
            fh2_ams_theta1_theta2[i + 2]->GetXaxis()->SetTitleOffset(1.);
            fh2_ams_theta1_theta2[i + 2]->GetYaxis()->SetTitleOffset(1.5);
            fh2_ams_theta1_theta2[i + 2]->GetXaxis()->CenterTitle(true);
            fh2_ams_theta1_theta2[i + 2]->GetYaxis()->CenterTitle(true);
            cHitThetaCor->cd(i * 2 + 2);
            fh2_ams_theta1_theta2[i + 2]->Draw("COLZ");

            // Phi correlations per side
            if (i == 0)
            {
                sprintf(Name1, "fh2_Ams_#phiinner_vs_#phiouter_right");
                sprintf(Name2, "AMS: #phi-inner vs #phi-outer for right side");
                fh2_ams_phi1_phi2[i] = new TH2F(Name1, Name2, 180, 90., 270., 180, 90., 270.);
            }
            else
            {
                sprintf(Name1, "fh2_Ams_#phiinner_vs_#phiouter_left");
                sprintf(Name2, "AMS: #phi-inner vs #phi-outer for left side");
                fh2_ams_phi1_phi2[i] = new TH2F(Name1, Name2, 180, -180., 180., 180, -180., 180.);
            }
            fh2_ams_phi1_phi2[i]->GetXaxis()->SetTitle("#phi-inner [degrees]");
            fh2_ams_phi1_phi2[i]->GetYaxis()->SetTitle("#phi-outer [degrees]");
            fh2_ams_phi1_phi2[i]->GetXaxis()->SetTitleOffset(1.);
            fh2_ams_phi1_phi2[i]->GetYaxis()->SetTitleOffset(1.5);
            fh2_ams_phi1_phi2[i]->GetXaxis()->CenterTitle(true);
            fh2_ams_phi1_phi2[i]->GetYaxis()->CenterTitle(true);
            cHitPhiCor->cd(i + 1);
            fh2_ams_phi1_phi2[i]->Draw("COLZ");
        }
    }

    // MAIN FOLDER-AMS
    TFolder* mainfolAms = new TFolder("AMS", "AMS info");

    // Folder for mapped data
    TFolder* mapfolAms = new TFolder("Map", "Map AMS info");
    mapfolAms->Add(cMap);
    mainfolAms->Add(mapfolAms);

    // Folder for cal data
    if (fCalItemsAms)
    {
        TFolder* calfolAms = new TFolder("Cal", "Cal AMS info");
        calfolAms->Add(cCalL);
        calfolAms->Add(cCalR);
        mainfolAms->Add(calfolAms);
    }

    // Folder for hit data
    if (fHitItemsAms)
    {
        TFolder* hitfolAms = new TFolder("Hit", "Hit AMS info");
        for (Int_t i = 0; i < fNbDet; i++)
            hitfolAms->Add(cHit[i]);
        mainfolAms->Add(hitfolAms);

        // Experiment s455
        if (fNbDet == 6)
        {
            TFolder* CorSidefolAms = new TFolder("Correlations_per_side", "Hit AMS info for correlations per side");
            CorSidefolAms->Add(cHitAngles);
            CorSidefolAms->Add(cHitThetaCor);
            CorSidefolAms->Add(cHitPhiCor);
            CorSidefolAms->Add(cHitEnergyCor);
            mainfolAms->Add(CorSidefolAms);
        }
    }
    run->AddObject(mainfolAms);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_AMS", Form("/Objects/%s/->Reset_AMS_Histo()", GetName()));

    return kSUCCESS;
}

void R3BAmsOnlineSpectra::Reset_AMS_Histo()
{
    LOG(info) << "R3BAmsOnlineSpectra::Reset_AMS_Histo";

    // Mapped data
    for (Int_t i = 0; i < fNbDet; i++)
    {
        fh_Ams_energy_allStrips[i]->Reset();
    }

    // Cal data
    if (fCalItemsAms)
        for (Int_t i = 0; i < fNbDet; i++)
        {
            for (Int_t j = 0; j < 2; j++)
            {
                fh_Ams_energy_allCalStrips[i * 2 + j]->Reset();
            }
        }

    // Hit data
    if (fHitItemsAms)
    {
        for (Int_t i = 0; i < fNbDet; i++)
        {
            fh_Ams_hit_Mul[i][0]->Reset();
            fh_Ams_hit_Mul[i][1]->Reset();
            fh_Ams_hit_Pos[i]->Reset();
            fh_Ams_hit_E[i]->Reset();
            fh_Ams_hit_E_theta[i]->Reset();
        }

        if (fNbDet == 6)
            for (Int_t i = 0; i < 2; i++)
            {
                fh2_ams_theta_phi[i]->Reset();
                fh2_ams_theta1_theta2[i * 2]->Reset();
                fh2_ams_theta1_theta2[i * 2 + 1]->Reset();
                fh2_ams_phi1_phi2[i]->Reset();
                fh2_ams_e1_e2[i]->Reset();
            }
    }
}

void R3BAmsOnlineSpectra::Exec(Option_t* option)
{
    // Fill mapped data
    if (fMappedItemsAms && fMappedItemsAms->GetEntriesFast() > 0)
    {
        Int_t nHits = fMappedItemsAms->GetEntriesFast();
        // std::cout << "hit:"<<nHits << std::endl;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BAmsMappedData* hit = dynamic_cast<R3BAmsMappedData*>(fMappedItemsAms->At(ihit));
            if (!hit)
                continue;
            fh_Ams_energy_allStrips[hit->GetDetectorId()]->Fill(hit->GetStripId(), hit->GetEnergy());
        }
    }

    // Fill cal data
    if (fCalItemsAms && fCalItemsAms->GetEntriesFast() > 0)
    {
        Int_t nHits = fCalItemsAms->GetEntriesFast();
        // std::cout << "hit:"<<nHits << std::endl;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BAmsStripCalData* hit = dynamic_cast<R3BAmsStripCalData*>(fCalItemsAms->At(ihit));
            if (!hit)
                continue;
            fh_Ams_energy_allCalStrips[hit->GetDetId() * 2 + hit->GetSideId()]->Fill(hit->GetStripId(),
                                                                                     hit->GetEnergy());
        }
    }

    // Fill hit data
    if (fHitItemsAms && fHitItemsAms->GetEntriesFast() > 0)
    {
        Int_t mulhit[fNbDet][2];
        for (Int_t i = 0; i < fNbDet; i++)
            for (Int_t j = 0; j < 2; j++)
                mulhit[i][j] = 0;

        Float_t Emaxhit[fNbDet]; // just look for the max. energy per detector
        Float_t Thetamaxhit[fNbDet];
        Float_t Phimaxhit[fNbDet];
        for (Int_t i = 0; i < fNbDet; i++)
        {
            Emaxhit[i] = 0.;
            Thetamaxhit[i] = 0.; // at 0 degrees we have nothing!
            Phimaxhit[i] = 90.;  // at 90 degrees we have nothing!
        }

        Int_t nHits = fHitItemsAms->GetEntriesFast();
        Int_t DetId = -1;
        // std::cout << nHits << std::endl;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BAmsHitData* hit = dynamic_cast<R3BAmsHitData*>(fHitItemsAms->At(ihit));
            if (!hit)
                continue;
            DetId = hit->GetDetId();
            fh_Ams_hit_Pos[DetId]->Fill(hit->GetPos_S(), hit->GetPos_K());
            fh_Ams_hit_E[DetId]->Fill(hit->GetEnergyS(), hit->GetEnergyK());
            fh_Ams_hit_E_theta[DetId]->Fill(hit->GetTheta() * TMath::RadToDeg(), hit->GetEnergyS());
            mulhit[DetId][0] = hit->GetMulS();
            mulhit[DetId][1] = hit->GetMulK();

            // Only for experiment s455
            if (fNbDet == 6)
            {
                if (DetId == 0 || DetId == 3) // inner detectors, layer 1
                    fh2_ams_theta_phi[0]->Fill(hit->GetTheta() * TMath::RadToDeg(), hit->GetPhi() * TMath::RadToDeg());
                else // outer detectors, layer 2
                    fh2_ams_theta_phi[1]->Fill(hit->GetTheta() * TMath::RadToDeg(), hit->GetPhi() * TMath::RadToDeg());
                // look for the max. energy per AMS detector
                if (Emaxhit[DetId] < hit->GetEnergyS())
                {
                    Emaxhit[DetId] = hit->GetEnergyS();
                    Thetamaxhit[DetId] = hit->GetTheta() * TMath::RadToDeg();
                    if (DetId > 2)
                        Phimaxhit[DetId] = hit->GetPhi() * TMath::RadToDeg();
                    else
                    {
                        if (hit->GetPhi() < 0.)
                            Phimaxhit[DetId] = (TMath::Pi() + hit->GetPhi()) * TMath::RadToDeg() + 180.;
                        else
                            Phimaxhit[DetId] = hit->GetPhi() * TMath::RadToDeg();
                    }
                }
            }
        }
        for (Int_t i = 0; i < fNbDet; i++)
        {
            fh_Ams_hit_Mul[i][0]->Fill(mulhit[i][0]);
            fh_Ams_hit_Mul[i][1]->Fill(mulhit[i][1]);
        }

        // Only for experiment s455
        if (fNbDet == 6)
        {
            // right up
            if (Emaxhit[0] > 0. && Emaxhit[1] > 0.)
            {
                fh2_ams_e1_e2[0]->Fill(Emaxhit[0], TMath::Max(Emaxhit[1], Emaxhit[2]));
                fh2_ams_theta1_theta2[0]->Fill(Thetamaxhit[0], Thetamaxhit[1]);
                fh2_ams_phi1_phi2[0]->Fill(Phimaxhit[0], Phimaxhit[1]);
            }
            // right down
            if (Emaxhit[0] > 0. && Emaxhit[2] > 0.)
            {
                fh2_ams_e1_e2[0]->Fill(Emaxhit[0], TMath::Max(Emaxhit[1], Emaxhit[2]));
                fh2_ams_theta1_theta2[1]->Fill(Thetamaxhit[0], Thetamaxhit[2]);
                fh2_ams_phi1_phi2[0]->Fill(Phimaxhit[0], Phimaxhit[2]);
            }
            // left up
            if (Emaxhit[3] > 0. && Emaxhit[4] > 0.)
            {
                fh2_ams_e1_e2[1]->Fill(Emaxhit[3], TMath::Max(Emaxhit[4], Emaxhit[5]));
                fh2_ams_theta1_theta2[2]->Fill(Thetamaxhit[3], Thetamaxhit[4]);
                fh2_ams_phi1_phi2[1]->Fill(Phimaxhit[3], Phimaxhit[4]);
            }
            // left down
            if (Emaxhit[3] > 0. && Emaxhit[5] > 0.)
            {
                fh2_ams_e1_e2[1]->Fill(Emaxhit[3], TMath::Max(Emaxhit[4], Emaxhit[5]));
                fh2_ams_theta1_theta2[3]->Fill(Thetamaxhit[3], Thetamaxhit[5]);
                fh2_ams_phi1_phi2[1]->Fill(Phimaxhit[3], Phimaxhit[5]);
            }
        }
    }

    fNEvents += 1;
}

void R3BAmsOnlineSpectra::FinishEvent()
{

    if (fMappedItemsAms)
    {
        fMappedItemsAms->Clear();
    }
    if (fCalItemsAms)
    {
        fCalItemsAms->Clear();
    }
    if (fHitItemsAms)
    {
        fHitItemsAms->Clear();
    }
}

void R3BAmsOnlineSpectra::FinishTask()
{

    if (fMappedItemsAms)
    {
        cMap->Write();
    }
    if (fCalItemsAms)
    {
        cCalR->Write();
        cCalL->Write();
    }
    if (fHitItemsAms)
    {
        for (Int_t i = 0; i < fNbDet; i++)
            cHit[i]->Write();
        // Only for experiment s455
        if (fNbDet == 6)
        {
            cHitAngles->Write();
            cHitThetaCor->Write();
            cHitPhiCor->Write();
            cHitEnergyCor->Write();
        }
    }
}

ClassImp(R3BAmsOnlineSpectra)
