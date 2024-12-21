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

// --------------------------------------------------------------
// -----       R3BAmsCalifaCorrelatedOnlineSpectra             --
// -----    Created 17/02/19  by J.L. Rodriguez-Sanchez        --
// ----- Fill AMS and CALIFA correlations in online histograms --
// --------------------------------------------------------------

/*
 *  This taks reads hit data from AMS and califa detectors and plots
 *  online histograms
 */

#include "R3BAmsCalifaCorrelatedOnlineSpectra.h"
#include "R3BAmsHitData.h"
#include "R3BCalifaClusterData.h"
#include "R3BEventHeader.h"
#include "R3BLosCalData.h"
#include "THttpServer.h"
#include "TRandom.h"
#include "TVector3.h"

#include "R3BTCalEngine.h"

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
#include <algorithm>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#define IS_NAN(x) TMath::IsNaN(x)

using namespace std;

R3BAmsCalifaCorrelatedOnlineSpectra::R3BAmsCalifaCorrelatedOnlineSpectra()
    : FairTask("AmsCalifaOnlineSpectra", 1)
    , fHitItemsAms(NULL)
    , fHitItemsCalifa(NULL)
    , fCalItemsLos(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(4)
{
}

R3BAmsCalifaCorrelatedOnlineSpectra::R3BAmsCalifaCorrelatedOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHitItemsAms(NULL)
    , fHitItemsCalifa(NULL)
    , fCalItemsLos(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fTrigger(-1)
    , fNEvents(0)
    , fNbDet(4)
{
}

R3BAmsCalifaCorrelatedOnlineSpectra::~R3BAmsCalifaCorrelatedOnlineSpectra() {}

InitStatus R3BAmsCalifaCorrelatedOnlineSpectra::Init()
{

    LOG(info) << "R3BAmsCalifaCorrelatedOnlineSpectra::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BAmsCalifaCorrelatedOnlineSpectra::Init FairRootManager not found";
    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // create histograms of all detectors

    // get access to Hit data
    fHitItemsAms = dynamic_cast<TClonesArray*>(mgr->GetObject("AmsHitData"));
    if (!fHitItemsAms)
    {
        LOG(info) << "R3BAmsCalifaCorrelatedOnlineSpectra::Init AmsHitData not found";
    }

    // get access to Hit data
    fHitItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaClusterData"));
    if (!fHitItemsCalifa)
    {
        LOG(info) << "R3BAmsCalifaCorrelatedOnlineSpectra::Init CalifaClusterData not found";
    }

    // get access to Hit data
    fCalItemsLos = dynamic_cast<TClonesArray*>(mgr->GetObject("LosCal"));
    if (!fCalItemsLos)
    {
        LOG(info) << "R3BAmsCalifaCorrelatedOnlineSpectra::Init LosCalData not found";
    }

    // CANVAS 7
    cAMSangles = new TCanvas("AMS_Theta_vs_Phi", "Theta vs Phi", 10, 10, 500, 500);
    fh_ams_theta_phi = new TH2F("Ams_theta_vs_phi", "AMS theta vs phi", 90, 0, 90, 180, -180, 180);
    fh_ams_theta_phi->GetXaxis()->SetTitle("Theta (degrees)");
    fh_ams_theta_phi->GetYaxis()->SetTitle("Phi (degrees)");
    fh_ams_theta_phi->GetXaxis()->SetTitleOffset(1.2);
    fh_ams_theta_phi->GetYaxis()->SetTitleOffset(1.2);
    fh_ams_theta_phi->GetXaxis()->CenterTitle(true);
    fh_ams_theta_phi->GetYaxis()->CenterTitle(true);
    fh_ams_theta_phi->Draw("COLZ");

    char Name1[255];
    char Name2[255];

    // CANVAS 8
    cHitone = new TCanvas("AMS_theta_top_bottom", "theta between opposite detectors", 10, 10, 500, 500);
    cHitone->Divide(1, 2);

    fh_ams_thetatop_thetabottom =
        new TH2F("Ams_theta_top_vs_theta_bottom", "AMS theta top vs theta bottom", 90, 0, 90, 90, 0, 90);
    fh_ams_thetatop_thetabottom->GetXaxis()->SetTitle("Theta top (degrees)");
    fh_ams_thetatop_thetabottom->GetYaxis()->SetTitle("Theta bottom (degrees)");
    fh_ams_thetatop_thetabottom->GetXaxis()->SetTitleOffset(1.2);
    fh_ams_thetatop_thetabottom->GetYaxis()->SetTitleOffset(1.2);
    fh_ams_thetatop_thetabottom->GetXaxis()->CenterTitle(true);
    fh_ams_thetatop_thetabottom->GetYaxis()->CenterTitle(true);
    cHitone->cd(1);
    fh_ams_thetatop_thetabottom->Draw("COLZ");

    sprintf(Name1, "fh_Ams_opening_angle1");
    sprintf(Name2, "AMS Opening Angle (beam hits target at (0,0,0))");
    fh_Ams_openangle1 = new TH1F(Name1, Name2, 140, 0, 140);
    fh_Ams_openangle1->GetXaxis()->SetTitle("Opening angle [degrees]");
    fh_Ams_openangle1->GetYaxis()->SetTitle("Counts");
    fh_Ams_openangle1->GetXaxis()->CenterTitle(true);
    fh_Ams_openangle1->GetYaxis()->CenterTitle(true);
    fh_Ams_openangle1->GetYaxis()->SetTitleOffset(1.2);
    fh_Ams_openangle1->GetXaxis()->SetTitleOffset(1.2);
    cHitone->cd(2);
    fh_Ams_openangle1->Draw();

    //  CANVAS 9  -------------------------------
    cHittwo = new TCanvas("AMS_theta_left_right", "theta between opposite detectors", 10, 10, 500, 500);
    cHittwo->Divide(1, 2);

    fh_ams_thetaright_thetaleft =
        new TH2F("Ams_theta_left_vs_theta_right", "AMS theta left vs theta right", 90, 0, 90, 90, 0, 90);
    fh_ams_thetaright_thetaleft->GetXaxis()->SetTitle("Theta left (degrees)");
    fh_ams_thetaright_thetaleft->GetYaxis()->SetTitle("Theta right (degrees)");
    fh_ams_thetaright_thetaleft->GetXaxis()->SetTitleOffset(1.2);
    fh_ams_thetaright_thetaleft->GetYaxis()->SetTitleOffset(1.2);
    fh_ams_thetaright_thetaleft->GetXaxis()->CenterTitle(true);
    fh_ams_thetaright_thetaleft->GetYaxis()->CenterTitle(true);
    cHittwo->cd(1);
    fh_ams_thetaright_thetaleft->Draw("COLZ");

    sprintf(Name1, "fh_Ams_opening_angle2");
    sprintf(Name2, "AMS Opening Angle (beam hits target at (0,0,0))");
    fh_Ams_openangle2 = new TH1F(Name1, Name2, 140, 0, 140);
    fh_Ams_openangle2->GetXaxis()->SetTitle("Opening angle [degrees]");
    fh_Ams_openangle2->GetYaxis()->SetTitle("Counts");
    fh_Ams_openangle2->GetXaxis()->CenterTitle(true);
    fh_Ams_openangle2->GetYaxis()->CenterTitle(true);
    fh_Ams_openangle2->GetYaxis()->SetTitleOffset(1.2);
    fh_Ams_openangle2->GetXaxis()->SetTitleOffset(1.2);
    cHittwo->cd(2);
    fh_Ams_openangle2->Draw();

    //  CANVAS 10  ------------------------------
    cHitPhis = new TCanvas("AMS_Phi_Angles", "Phi between opposite detectors", 10, 10, 500, 500);
    fh_ams_phis = new TH2F("Ams_phis", "AMS Phi angles between opposite detectors", 90, 0, 360, 90, 0, 360);
    fh_ams_phis->GetXaxis()->SetTitle("Phi1 (degrees)");
    fh_ams_phis->GetYaxis()->SetTitle("Phi2 (degrees)");
    fh_ams_phis->GetXaxis()->SetTitleOffset(1.2);
    fh_ams_phis->GetYaxis()->SetTitleOffset(1.2);
    fh_ams_phis->GetXaxis()->CenterTitle(true);
    fh_ams_phis->GetYaxis()->CenterTitle(true);
    fh_ams_phis->Draw("COLZ");

    //  CANVAS 1  -------------------------------
    cMap1 = new TCanvas("AMS_top_right_Petal_1", "Correlation AMS_top_right and Califa Petals 1", 10, 10, 500, 500);
    // cMap->Divide(2,2);

    fh_Califa_coinc_petal1 =
        new TH2F("Ams_top_right_vs_petal_1", "Correlation AMS_top_right and Califa Petal 1", 90, 0, 90, 30, 0, 90);
    fh_Califa_coinc_petal1->GetXaxis()->SetTitle("Theta AMS (degrees)");
    fh_Califa_coinc_petal1->GetYaxis()->SetTitle("Theta petal 1 (degrees)");
    fh_Califa_coinc_petal1->GetYaxis()->SetTitleOffset(1.2);
    fh_Califa_coinc_petal1->GetXaxis()->CenterTitle(true);
    fh_Califa_coinc_petal1->GetYaxis()->CenterTitle(true);
    fh_Califa_coinc_petal1->Draw("COLZ");

    //  CANVAS 1  -------------------------------
    cMap1e = new TCanvas("AMS_top_right_Petal_1_E", "Correlation AMS_top_right and Califa Petals 1", 10, 10, 500, 500);
    // cMap->Divide(2,2);

    fh_Califa_coinc_petal1e = new TH2F(
        "Ams_top_right_vs_petal_1_E", "Correlation AMS_top_right and Califa Petal 1", 1000, 0, 5000, 250, 0, 600);
    fh_Califa_coinc_petal1e->GetXaxis()->SetTitle("Energy AMS (arb. units)");
    fh_Califa_coinc_petal1e->GetYaxis()->SetTitle("Energy petal 1 (MeV)");
    fh_Califa_coinc_petal1e->GetYaxis()->SetTitleOffset(1.3);
    fh_Califa_coinc_petal1e->GetXaxis()->CenterTitle(true);
    fh_Califa_coinc_petal1e->GetYaxis()->CenterTitle(true);
    fh_Califa_coinc_petal1e->Draw("COLZ");

    //  CANVAS 1  -------------------------------
    cMap2 = new TCanvas(
        "AMS_left_bottom_Petal_3_4", "Correlation AMS_left_bottom and Califa Petals 3 and 4", 10, 10, 500, 500);
    // cMap->Divide(2,2);

    fh_Califa_coinc_petal2 = new TH2F("Ams_left_bottom_vs_petal_3and4",
                                      "Correlation AMS_left_bottom and Califa Petals 3 and 4",
                                      90,
                                      0,
                                      90,
                                      30,
                                      0,
                                      90);
    fh_Califa_coinc_petal2->GetXaxis()->SetTitle("Theta AMS (degrees)");
    fh_Califa_coinc_petal2->GetYaxis()->SetTitle("Theta petal 3 and 4 (degrees)");
    fh_Califa_coinc_petal2->GetYaxis()->SetTitleOffset(1.2);
    fh_Califa_coinc_petal2->GetXaxis()->CenterTitle(true);
    fh_Califa_coinc_petal2->GetYaxis()->CenterTitle(true);
    fh_Califa_coinc_petal2->Draw("COLZ");

    cMap2e = new TCanvas(
        "AMS_left_bottom_Petal_3_4_E", "Correlation AMS_left_bottom and Califa Petals 3 and 4", 10, 10, 500, 500);
    // cMap->Divide(2,2);

    fh_Califa_coinc_petal2e = new TH2F("Ams_left_bottom_vs_petal_3and4_E",
                                       "Correlation AMS_left_bottom and Califa Petals 3 and 4",
                                       1000,
                                       0,
                                       5000,
                                       250,
                                       0,
                                       600);
    fh_Califa_coinc_petal2e->GetXaxis()->SetTitle("Energy AMS (arb. units)");
    fh_Califa_coinc_petal2e->GetYaxis()->SetTitle("Energy petal 3 and 4 (MeV)");
    fh_Califa_coinc_petal2e->GetYaxis()->SetTitleOffset(1.3);
    fh_Califa_coinc_petal2e->GetXaxis()->CenterTitle(true);
    fh_Califa_coinc_petal2e->GetYaxis()->CenterTitle(true);
    fh_Califa_coinc_petal2e->Draw("COLZ");

    TCanvas* cLos = new TCanvas("Los", "LOS", 10, 10, 750, 850);
    cLos->Divide(1, 2);
    fh_los_pos_MCFD = new TH2F("los_pos_MCFD", "LOS MCFD Position emp. cal.", 1000, -10., 10., 1000, -10., 10.);
    fh_los_pos_MCFD->GetXaxis()->SetTitle("X position / cm");
    fh_los_pos_MCFD->GetYaxis()->SetTitle("Y position / cm");
    fh_los_pos_TAMEX = new TH2F("los_pos_TAMEX", "LOS TAMEX Position emp. cal.", 1000, -10., 10., 1000, -10., 10.);
    fh_los_pos_TAMEX->GetXaxis()->SetTitle("X position / cm");
    fh_los_pos_TAMEX->GetYaxis()->SetTitle("Y position / cm");

    cLos->cd(1);
    gPad->SetLogz();
    fh_los_pos_MCFD->Draw("colz");
    cLos->cd(2);
    gPad->SetLogz();
    fh_los_pos_TAMEX->Draw("colz");

    // MAIN FOLDER-AMS-CALIFA
    TFolder* mainfolAmsCal = new TFolder("LOS_AMS_CALIFA", "LOS, AMS and CALIFA correlations");
    mainfolAmsCal->Add(cAMSangles);
    mainfolAmsCal->Add(cHitone);
    mainfolAmsCal->Add(cHittwo);
    mainfolAmsCal->Add(cHitPhis);
    mainfolAmsCal->Add(cMap1);
    mainfolAmsCal->Add(cMap2);
    mainfolAmsCal->Add(cMap1e);
    mainfolAmsCal->Add(cMap2e);
    mainfolAmsCal->Add(cLos);
    run->AddObject(mainfolAmsCal);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_LOS_AMS_CALIFA",
                                          Form("/Objects/%s/->Reset_AMS_CALIFA_Histo()", GetName()));

    return kSUCCESS;
}

void R3BAmsCalifaCorrelatedOnlineSpectra::Reset_AMS_CALIFA_Histo()
{
    LOG(info) << "R3BAmsCalifaCorrelatedOnlineSpectra::Reset_LOS_AMS_CALIFA_Histo";

    fh_Califa_coinc_petal1->Reset();
    fh_Califa_coinc_petal2->Reset();

    fh_Califa_coinc_petal1e->Reset();
    fh_Califa_coinc_petal2e->Reset();

    fh_ams_theta_phi->Reset();
    fh_ams_thetatop_thetabottom->Reset();
    fh_ams_thetaright_thetaleft->Reset();
    fh_ams_phis->Reset();
    fh_Ams_openangle1->Reset();
    fh_Ams_openangle2->Reset();
    fh_los_pos_MCFD->Reset();
    fh_los_pos_TAMEX->Reset();

    // Hit data
    for (Int_t i = 0; i < fNbDet; i++)
    {
        // fh_Ams_hit_Mul[i]->Reset();
        // fh_Ams_hit_Pos[i]->Reset();
        // fh_Ams_hit_E[i]->Reset();
    }
}

void R3BAmsCalifaCorrelatedOnlineSpectra::Exec(Option_t* option)
{

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BAmsCalifaCorrelatedOnlineSpectra::Exec FairRootManager not found";

    // if(header->GetTrigger()!=1){cout << header->GetTrigger()<<endl;}
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != 1))
        return;

    // los

    Double_t time_V[10][8] = { { 0.0 / 0.0 } }; // [multihit][pm]
    Double_t time_L[10][8] = { { 0.0 / 0.0 } };
    Double_t time_T[10][8] = { { 0.0 / 0.0 } };
    Double_t timeLosM[10] = { 0.0 };
    Double_t LosTresM[10] = { 0.0 / 0.0 };
    Double_t timeLosT[10] = { 0.0 };
    Double_t LosTresT[10] = { 0.0 / 0.0 };
    Double_t timeLos[10] = { 0.0 };
    Double_t totsum[10] = { 0.0 };
    Double_t tot[10][8] = { { 0.0 / 0.0 } };
    Double_t xT_cm[10] = { 0.0 / 0.0 };
    Double_t yT_cm[10] = { 0.0 / 0.0 };
    Double_t xV_cm[10] = { 0.0 / 0.0 };
    Double_t yV_cm[10] = { 0.0 / 0.0 };
    TVector3 PosBeam(0., 0., 0.);

    Int_t nPart;
    if (fCalItemsLos && fCalItemsLos->GetEntriesFast())
    {
        nPart = fCalItemsLos->GetEntriesFast();

        // if (nPart>0) fh_los_multihit->Fill(nPart);

        Int_t iDet = 0;
        Int_t nPart_VFTX[8] = { 0 };
        Int_t nPart_LEAD[8] = { 0 };
        Int_t nPart_TRAI[8] = { 0 };

        for (Int_t iPart = 0; iPart < nPart; iPart++)
        {

            /*
             * nPart is the number of particle passing through LOS detector in one event
             */
            R3BLosCalData* calData = dynamic_cast<R3BLosCalData*>(fCalItemsLos->At(iPart));
            iDet = calData->GetDetector();

            // lt=0, l=1,lb=2,b=3,rb=4,r=5,rt=6,t=7

            for (Int_t iCha = 0; iCha < 8; iCha++)
            {

                time_V[iPart][iCha] = 0. / 0.;
                if (!(IS_NAN(calData->GetTimeV_ns(iCha))))
                { // VFTX
                    nPart_VFTX[iCha] += 1;
                    time_V[iPart][iCha] = calData->GetTimeV_ns(iCha);
                }

                time_L[iPart][iCha] = 0. / 0.;
                if (!(IS_NAN(calData->GetTimeL_ns(iCha))))
                { // TAMEX leading
                    nPart_LEAD[iCha] += 1;
                    time_L[iPart][iCha] = calData->GetTimeL_ns(iCha);
                }
                time_T[iPart][iCha] = 0. / 0.;
                if (!(IS_NAN(calData->GetTimeT_ns(iCha))))
                { // TAMEX trailing
                    nPart_TRAI[iCha] += 1;
                    time_T[iPart][iCha] = calData->GetTimeT_ns(iCha);
                }
            }

            // Sorting VFTX data:

            std::qsort(time_V,
                       nPart,
                       sizeof(*time_V),
                       [](const void* arg1, const void* arg2) -> int
                       {
                           double const* lhs = static_cast<double const*>(arg1);
                           double const* rhs = static_cast<double const*>(arg2);

                           return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
                       });
            // End sorting

            /*
                  if(iPart > 0 && Multip%8 == 0)
                  {
                for(int k=0; k<8; k++)
                {
                  if(time_V[iPart][k] > 0. && time_V[iPart-1][k] > 0. && !(IS_NAN(time_V[iPart][k])) &&
               !(IS_NAN(time_V[iPart-1][k])))
                  {
                    fh_los_dt_hits->Fill(time_V[iPart][k]-time_V[iPart-1][k]);
                  }
                  if(time_L[iPart][k] > 0. && time_L[iPart-1][k] > 0. && !(IS_NAN(time_L[iPart][k])) &&
               !(IS_NAN(time_L[iPart-1][k])))
                  {
                    fh_los_dt_hits_l->Fill(time_L[iPart][k]-time_L[iPart-1][k]);
                  }
                  if(time_T[iPart][k] > 0. && time_T[iPart-1][k] > 0. && !(IS_NAN(time_T[iPart][k])) &&
               !(IS_NAN(time_T[iPart-1][k])))
                  {
                    fh_los_dt_hits_t->Fill(time_T[iPart][k]-time_T[iPart-1][k]);
                  }
                }
                  }	  */

            if (iDet == 1)
            {

                // lt=0, l=1,lb=2,b=3,rb=4,r=5,rt=6,t=7

                // if(iPart >= 0)
                if (1 == 1)
                {
                    int nPMT = 0;
                    int nPMV = 0;

                    for (int ipm = 0; ipm < 8; ipm++)
                    {
                        if (time_T[iPart][ipm] > 0. && time_L[iPart][ipm] > 0. && !(IS_NAN(time_T[iPart][ipm])) &&
                            !(IS_NAN(time_L[iPart][ipm])))
                        {
                            while (time_T[iPart][ipm] - time_L[iPart][ipm] < 0.)
                            {
                                time_T[iPart][ipm] = time_T[iPart][ipm] + 2048. * fClockFreq;
                            }

                            nPMT = nPMT + 1;
                            tot[iPart][ipm] = time_T[iPart][ipm] - time_L[iPart][ipm];
                        }

                        totsum[iPart] += tot[iPart][ipm];

                        // if(tot[iPart][ipm] != 0. && !(IS_NAN(tot[iPart][ipm])))
                        // fh_los_tot->Fill(ipm+1,tot[iPart][ipm]);

                        if (time_L[iPart][ipm] > 0. && !(IS_NAN(time_L[iPart][ipm])))
                            timeLosT[iPart] += time_L[iPart][ipm];

                        if (time_V[iPart][ipm] > 0. && !(IS_NAN(time_V[iPart][ipm])))
                        {
                            timeLosM[iPart] += time_V[iPart][ipm];
                            nPMV = nPMV + 1;
                        }
                    }

                    totsum[iPart] = totsum[iPart] / nPMT;
                    /*
                       if(totsum[iPart] < 88.) cout<<fNEvents<<"; "<<nPart<<"; "<<iPart<<",
                       "<<totsum[iPart]<<tot[iPart][0]<<
                       ", "<<tot[iPart][1]<<", " <<tot[iPart][2]<<", "<<tot[iPart][3]<<", "<<tot[iPart][4]<<", "
                       <<tot[iPart][5]<<", "<<tot[iPart][6]<<", "<<tot[iPart][7]<<endl;
                     */

                    timeLosM[iPart] = timeLosM[iPart] / nPMV;

                    timeLosT[iPart] = timeLosT[iPart] / nPMT;

                    timeLos[iPart] = timeLosM[iPart];

                    //    if(!(timeLosM[iPart] > 0.) && IS_NAN(timeLosM[iPart]))
                    //        cout<<"R3BOnline warn!! LOS VFTX time < 0 or nan! If nan we take TAMEX time, if <0 no
                    //        ToF info for this event!! "<<timeLosM[iPart]<<endl;
                    // Line below only after VFTX and TAMEX clocks are synhronized!
                    //  if(IS_NAN(timeLosM[iPart])) timeLos[iPart] = timeLosT[iPart];

                    if (nPMV == 8)
                        LosTresM[iPart] =
                            ((time_V[iPart][0] + time_V[iPart][2] + time_V[iPart][4] + time_V[iPart][6]) -
                             (time_V[iPart][1] + time_V[iPart][3] + time_V[iPart][5] + time_V[iPart][7])) /
                            4.;
                    if (nPMT == 8)
                        LosTresT[iPart] =
                            ((time_L[iPart][0] + time_L[iPart][2] + time_L[iPart][4] + time_L[iPart][6]) -
                             (time_L[iPart][1] + time_L[iPart][3] + time_L[iPart][5] + time_L[iPart][7])) /
                            4.;

                    // if(nPMV == 8) fh_los_tres_MCFD->Fill(LosTresM[iPart]);
                    // if(nPMT == 8) fh_los_tres_TAMEX->Fill(LosTresT[iPart]);
                    // if(nPMT == 8) fh_los_tot_mean->Fill(totsum[iPart]);
                    // Position from tamex:
                    xT_cm[iPart] =
                        (time_L[iPart][5] + time_L[iPart][6]) / 2. - (time_L[iPart][1] + time_L[iPart][2]) / 2.;
                    yT_cm[iPart] =
                        (time_L[iPart][7] + time_L[iPart][0]) / 2. - (time_L[iPart][3] + time_L[iPart][4]) / 2.;
                    // Empirical calibration (to be properly done):
                    xT_cm[iPart] = (xT_cm[iPart] - 0.) * 1.;
                    yT_cm[iPart] = (yT_cm[iPart] - 0.) * 1.;

                    // Position from VFTX:
                    xV_cm[iPart] =
                        (time_V[iPart][5] + time_V[iPart][6]) / 2. - (time_V[iPart][1] + time_V[iPart][2]) / 2.;
                    yV_cm[iPart] =
                        (time_V[iPart][7] + time_V[iPart][0]) / 2. - (time_V[iPart][3] + time_V[iPart][4]) / 2.;
                    // Empirical calibration (to be properly done):
                    xV_cm[iPart] = (xV_cm[iPart] - 0.) * 3.;
                    yV_cm[iPart] = (yV_cm[iPart] - 0.) * 3.;

                    if (nPMV == 8)
                        fh_los_pos_MCFD->Fill(xV_cm[iPart], yV_cm[iPart]);
                    if (nPMT == 8)
                    {
                        fh_los_pos_TAMEX->Fill(xT_cm[iPart], yT_cm[iPart]);
                        PosBeam.SetXYZ(xT_cm[iPart], yT_cm[iPart], 0.);
                    }
                }
            }
            else
            {
                cout << "Wrong detector ID for LOS!" << endl;
            }
        }
        /*
            for(int ik=0; ik<8; ik++)
            {
              fh_los_multihitVFTX->Fill(ik+1,nPart_VFTX[ik]);
              fh_los_multihitLEAD->Fill(ik+1,nPart_LEAD[ik]);
              fh_los_multihitTRAI->Fill(ik+1,nPart_TRAI[ik]);
            }   */
    }

    // califa and ams

    double hitams[100][2];
    double hitamse[100][2];
    for (Int_t i = 0; i < 100; i++)
        for (Int_t j = 0; j < 2; j++)
        {
            hitams[i][j] = 0.;
            hitamse[i][j] = 0.;
        }
    Int_t nbhitams[2];
    for (Int_t j = 0; j < 2; j++)
        nbhitams[j] = 0;
    TVector3 master[1000][fNbDet];
    Int_t mulhit[fNbDet];
    for (Int_t i = 0; i < fNbDet; i++)
        mulhit[i] = 0;

    // Fill hit data
    if (fHitItemsAms && fHitItemsAms->GetEntriesFast())
    {
        Int_t nHits = fHitItemsAms->GetEntriesFast();

        Double_t theta = 0., phi = 0.;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BAmsHitData* hit = dynamic_cast<R3BAmsHitData*>(fHitItemsAms->At(ihit));
            if (!hit)
                continue;
            if (hit->GetEnergyS() < 80 || hit->GetEnergyK() < 80 || hit->GetEnergyS() > 8500 ||
                hit->GetEnergyK() > 8500)
                continue;

            if (mulhit[hit->GetDetId()] < 1000)
            {
                TVector3 ams = hit->GetPosLab() - PosBeam;
                // master[mulhit[hit->GetDetId()]][hit->GetDetId()].SetMagThetaPhi(1.,hit->GetTheta(),hit->GetPhi());
                master[mulhit[hit->GetDetId()]][hit->GetDetId()] = ams;

                theta = ams.Theta() / TMath::Pi() * 180.;
                phi = ams.Phi() / TMath::Pi() * 180.;
                fh_ams_theta_phi->Fill(theta, phi);

                // std::cout << hit->GetDetId() << std::endl;
            }
            mulhit[hit->GetDetId()]++;

            if (hit->GetDetId() == 2 || hit->GetDetId() == 3)
            {

                TVector3 ams = hit->GetPosLab() - PosBeam;

                hitams[nbhitams[1]][1] = ams.Theta() / TMath::Pi() * 180.;
                hitamse[nbhitams[1]][1] = hit->GetEnergyS() + hit->GetEnergyK();
                nbhitams[1]++;
            }
            else
            {

                TVector3 ams = hit->GetPosLab() - PosBeam;

                hitams[nbhitams[0]][0] = ams.Theta() / TMath::Pi() * 180.;
                hitamse[nbhitams[0]][0] = hit->GetEnergyS() + hit->GetEnergyK();
                nbhitams[0]++;
            }
        }
    }

    // CALIFA
    Double_t hitcal[100][2];
    Double_t hitcale[100][2];
    for (Int_t i = 0; i < 100; i++)
        for (Int_t j = 0; j < 2; j++)
        {
            hitcal[i][j] = 0.;
            hitcale[i][j] = 0.;
        }
    Int_t nbhitcal[2];
    for (Int_t j = 0; j < 2; j++)
        nbhitcal[j] = 0;
    if (fHitItemsCalifa && fHitItemsCalifa->GetEntriesFast())
    {
        Int_t nHits = fHitItemsCalifa->GetEntriesFast();

        Double_t theta = 0., phi = 0.;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BCalifaClusterData* hit = dynamic_cast<R3BCalifaClusterData*>(fHitItemsCalifa->At(ihit));
            if (!hit)
                continue;
            theta = hit->GetTheta() / TMath::Pi() * 180.;
            phi = hit->GetPhi() / TMath::Pi() * 180.;
            if (phi > 120 && theta > 10)
            { // petal 1
                hitcal[nbhitcal[0]][0] = theta;
                hitcale[nbhitcal[0]][0] = hit->GetEnergy() / 1000.;
                nbhitcal[0]++;
            }
            else if (phi > -100 && phi < -10 && theta > 10)
            { // petal 3 and 4
                hitcal[nbhitcal[0]][1] = theta;
                hitcale[nbhitcal[0]][1] = hit->GetEnergy() / 1000.;
                nbhitcal[1]++;
            }
        }
    }

    // comparison AMS and Califa clusters petal 1
    Double_t window = 14.3; // degrees
    if (nbhitams[0] > 0 && nbhitcal[0] > 0)
    {
        for (Int_t i = 0; i < nbhitams[0]; i++)
        {
            for (Int_t j = 0; j < nbhitcal[0]; j++)
            {
                if (std::abs(hitcal[j][0] - hitams[i][0]) < window && hitcal[j][0] > 10.)
                {
                    fh_Califa_coinc_petal1->Fill(hitams[i][0], hitcal[j][0]);
                    fh_Califa_coinc_petal1e->Fill(hitamse[i][0], hitcale[j][0]);
                }
            }
        }
    }

    // comparison AMS and Califa clusters petal 3 and 4
    if (nbhitams[1] > 0 && nbhitcal[1] > 0)
    {
        for (Int_t i = 0; i < nbhitams[1]; i++)
        {
            for (Int_t j = 0; j < nbhitcal[1]; j++)
            {
                if (std::abs(hitcal[j][1] - hitams[i][1]) < window && hitcal[j][1] > 10.)
                {
                    fh_Califa_coinc_petal2->Fill(hitams[i][1], hitcal[j][1]);
                    fh_Califa_coinc_petal2e->Fill(hitamse[i][1], hitcale[j][1]);
                }
            }
        }
    }

    // AMS
    Double_t thetaopen = 0., phi1 = 0., phi2 = 0.;

    for (Int_t i = 0; i < mulhit[0]; i++)
    {
        for (Int_t j = i + 1; j < mulhit[1]; j++)
        {
            phi1 = master[i][0].Phi();
            phi2 = master[j][1].Phi();
            if (phi1 < 0)
                phi1 = phi1 + 2. * TMath::Pi();
            if (phi2 < 0)
                phi2 = phi2 + 2. * TMath::Pi();
            if (0.5 < gRandom->Uniform(0., 1.))
                fh_ams_phis->Fill(phi1 / TMath::Pi() * 180., phi2 / TMath::Pi() * 180.);
            else
                fh_ams_phis->Fill(phi2 / TMath::Pi() * 180., phi1 / TMath::Pi() * 180.);
        }
    }

    for (Int_t i = 0; i < mulhit[0]; i++)
    {
        for (Int_t j = i + 1; j < mulhit[2]; j++)
        {
            // std::cout<< master[i].Angle(master[j])/TMath::Pi()*180. <<std::endl;
            thetaopen = master[i][0].Angle(master[j][2]) / TMath::Pi() * 180.;
            if (thetaopen > 0)
            {
                fh_ams_thetatop_thetabottom->Fill(master[i][0].Theta() / TMath::Pi() * 180.,
                                                  master[j][2].Theta() / TMath::Pi() * 180.);
                fh_Ams_openangle1->Fill(thetaopen);
                phi1 = master[i][0].Phi();
                phi2 = master[j][2].Phi();
                // std::cout<< phi1 << " " <<phi2 <<std::endl;
                if (phi1 < 0)
                    phi1 = phi1 + 2. * TMath::Pi();
                if (phi2 < 0)
                    phi2 = phi2 + 2. * TMath::Pi();
                if (0.5 < gRandom->Uniform(0., 1.))
                    fh_ams_phis->Fill(phi1 / TMath::Pi() * 180., phi2 / TMath::Pi() * 180.);
                else
                    fh_ams_phis->Fill(phi2 / TMath::Pi() * 180., phi1 / TMath::Pi() * 180.);
            }
            else
            {
                std::cout << "Problem with opening angles" << std::endl;
            }
        }
    }

    for (Int_t i = 0; i < mulhit[0]; i++)
    {
        for (Int_t j = i + 1; j < mulhit[3]; j++)
        {
            phi1 = master[i][0].Phi();
            phi2 = master[j][3].Phi();
            if (phi1 < 0)
                phi1 = phi1 + 2. * TMath::Pi();
            if (phi2 < 0)
                phi2 = phi2 + 2. * TMath::Pi();
            if (0.5 < gRandom->Uniform(0., 1.))
                fh_ams_phis->Fill(phi1 / TMath::Pi() * 180., phi2 / TMath::Pi() * 180.);
            else
                fh_ams_phis->Fill(phi2 / TMath::Pi() * 180., phi1 / TMath::Pi() * 180.);
        }
    }

    for (Int_t i = 0; i < mulhit[1]; i++)
    {
        for (Int_t j = i + 1; j < mulhit[2]; j++)
        {
            phi1 = master[i][1].Phi();
            phi2 = master[j][2].Phi();
            if (phi1 < 0)
                phi1 = phi1 + 2. * TMath::Pi();
            if (phi2 < 0)
                phi2 = phi2 + 2. * TMath::Pi();
            if (0.5 < gRandom->Uniform(0., 1.))
                fh_ams_phis->Fill(phi1 / TMath::Pi() * 180., phi2 / TMath::Pi() * 180.);
            else
                fh_ams_phis->Fill(phi2 / TMath::Pi() * 180., phi1 / TMath::Pi() * 180.);
        }
    }

    for (Int_t i = 0; i < mulhit[1]; i++)
    {
        for (Int_t j = i + 1; j < mulhit[3]; j++)
        {

            // std::cout<< master[i].Angle(master[j])/TMath::Pi()*180. <<std::endl;
            thetaopen = master[i][1].Angle(master[j][3]) / TMath::Pi() * 180.;
            if (thetaopen > 0)
            {
                fh_ams_thetaright_thetaleft->Fill(master[i][1].Theta() / TMath::Pi() * 180.,
                                                  master[j][3].Theta() / TMath::Pi() * 180.);
                fh_Ams_openangle2->Fill(thetaopen);
                phi1 = master[i][1].Phi();
                phi2 = master[j][3].Phi();
                if (phi1 < 0)
                    phi1 = phi1 + 2. * TMath::Pi();
                if (phi2 < 0)
                    phi2 = phi2 + 2. * TMath::Pi();
                if (0.5 < gRandom->Uniform(0., 1.))
                    fh_ams_phis->Fill(phi1 / TMath::Pi() * 180., phi2 / TMath::Pi() * 180.);
                else
                    fh_ams_phis->Fill(phi2 / TMath::Pi() * 180., phi1 / TMath::Pi() * 180.);
            }
            else
            {
                std::cout << "Problem with opening angles" << std::endl;
            }
        }
    }

    for (Int_t i = 0; i < mulhit[2]; i++)
    {
        for (Int_t j = i + 1; j < mulhit[3]; j++)
        {
            phi1 = master[i][2].Phi();
            phi2 = master[j][3].Phi();
            if (phi1 < 0)
                phi1 = phi1 + 2. * TMath::Pi();
            if (phi2 < 0)
                phi2 = phi2 + 2. * TMath::Pi();
            if (0.5 < gRandom->Uniform(0., 1.))
                fh_ams_phis->Fill(phi1 / TMath::Pi() * 180., phi2 / TMath::Pi() * 180.);
            else
                fh_ams_phis->Fill(phi2 / TMath::Pi() * 180., phi1 / TMath::Pi() * 180.);
        }
    }

    fNEvents += 1;
}

void R3BAmsCalifaCorrelatedOnlineSpectra::FinishEvent()
{

    if (fHitItemsAms)
    {
        fHitItemsAms->Clear();
    }
    if (fHitItemsCalifa)
    {
        fHitItemsCalifa->Clear();
    }
    if (fCalItemsLos)
    {
        fCalItemsLos->Clear();
    }
}

void R3BAmsCalifaCorrelatedOnlineSpectra::FinishTask()
{
    if (fHitItemsAms)
    {
        cAMSangles->Write();
        cHitone->Write();
        cHittwo->Write();
        cHitPhis->Write();
    }
    if (fHitItemsAms && fHitItemsCalifa)
    {
        cMap1->Write();
        cMap2->Write();
        cMap1e->Write();
        cMap2e->Write();
    }
}

ClassImp(R3BAmsCalifaCorrelatedOnlineSpectra)
