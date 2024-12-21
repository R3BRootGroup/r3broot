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
// -----                  R3BOnlineSpectra                -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BOnlineSpectraLosVsSci2.h"
#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"
#include "R3BMusicHitData.h"
#include "R3BSci2TcalData.h"

#include "R3BEventHeader.h"
#include "R3BSamplerMappedData.h"
#include "R3BTCalEngine.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"

#include "TClonesArray.h"
#include "TFolder.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <algorithm>
#include <array>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

#define SPEED_OF_LIGHT_MNS 0.299792458

R3BOnlineSpectraLosVsSci2::R3BOnlineSpectraLosVsSci2()
    : R3BOnlineSpectraLosVsSci2("R3BOnlineSpectraLosVsSci2", 1)
{
}

R3BOnlineSpectraLosVsSci2::R3BOnlineSpectraLosVsSci2(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fClockFreq(5.)
    , fNEvents(0)
    , fTcalSci2(NULL)
    , fHitItemsMus(NULL)
    //, fFrsDataCA(NULL)
    , fToFoffset(0)
    , fToFmin(-5000)
    , fToFmax(5000)
    , fTof2InvV_p0(-7.8)
    , fTof2InvV_p1(0.0073)
    , fL2(137)
    , fPos_p0(-11)
    , fPos_p1(54.7)
    , fDispersionS2(7000.)
    , fBrho0_S2toCC(9.458)
    , fP0(-2.12371e7)
    , fP1(4.9473e7)
    , fP2(-2.87635e7)
    , fZprimary(50.)
    , fZoffset(-1.3)
    , ZMUSIC_cut(50)
    , ZMUSIC_wcut(0.5)
    , fOnline(kFALSE)
{
}

R3BOnlineSpectraLosVsSci2::~R3BOnlineSpectraLosVsSci2()
{
    if (fTcalSci2)
        delete fTcalSci2;
    if (fHitItemsMus)
        delete fHitItemsMus;
    /* if (fFrsDataCA)
     {
         delete fFrsDataCA;
     }*/
}

InitStatus R3BOnlineSpectraLosVsSci2::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BOnlineSpectraLosVsSci2::Init()";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    FairRunOnline* run = FairRunOnline::Instance();

    run->GetHttpServer()->Register("", this);

    // Get objects for detectors on all levels

    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    for (int det = 0; det < DET_MAX; det++)
    {
        fMappedItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sMapped", fDetectorNames[det]))));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        fCalItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sCal", fDetectorNames[det]))));
    }

    // --- Get access to Sci2 data at Tcal level --- //
    fTcalSci2 = dynamic_cast<TClonesArray*>(mgr->GetObject("Sci2Tcal"));
    if (!fTcalSci2)
    {
        LOG(info) << "R3BOnlineSpectraLosVsSci2::Init()  Could not find ci2Tcal";
    }

    // get access to hit data of the MUSIC
    fHitItemsMus = dynamic_cast<TClonesArray*>(mgr->GetObject("MusicHitData"));
    if (!fHitItemsMus)
        LOG(warn) << "R3BOnlineSpectraLosVsSci2: MusicHitData not found";

    // OUTPUT DATA
    /*fFrsDataCA = new TClonesArray("R3BFrsData", 5);
    if (!fOnline)
    {
        mgr->Register("FrsData", "Analysis FRS", fFrsDataCA, kTRUE);
    }
    else
    {
        mgr->Register("FrsData", "Analysis FRS", fFrsDataCA, kFALSE);
    }*/

    // MAIN FOLDER-ID
    TFolder* mainfolId = new TFolder("FRS_ID", "FRS ID info");

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------
    cPos = new TCanvas("PosSci2_m1", "PosSci2_m1", 10, 10, 800, 700);
    cPos->Divide(2, 2);
    fh1_RawPos_m1 = new TH1F("RawPosS2_mult1", "RawPosS2 mult=1", 100000, -50, 50);
    cPos->cd(1);
    fh1_RawPos_m1->Draw();
    fh1_CalPos_m1 = new TH1F("CalPosS2_mult1", "CalPosS2 mm mult=1", 30000, -150, 150);
    cPos->cd(2);
    fh1_CalPos_m1->Draw();
    // run->AddObject(cPos);
    fh1_RawPos = new TH1F("RawPosS2", "RawPosS2 ", 100000, -50, 50);
    cPos->cd(3);
    fh1_RawPos->Draw();
    fh1_CalPos = new TH1F("CalPosS2", "CalPosS2 mm ", 30000, -150, 150);
    cPos->cd(4);
    fh1_CalPos->Draw();
    mainfolId->Add(cPos);
    // run->AddObject(cPos);

    cTofFromS2 = new TCanvas("Tof_Sci2_Los_m1", "Tof_Sci2_Los_m1", 10, 10, 800, 700);
    cTofFromS2->Divide(1, 2);
    fh1_RawTofFromS2_TcalMult1 =
        new TH1D("RawTofNs_Tcal_m1_Sci2_Los", "RawTofNs_Tcal_m1_Sci2_Los", 10000000, -50000, 50000);
    fh1_RawTofFromS2_TcalMult1->GetXaxis()->SetTitle("Raw Tof [ns]");
    fh1_RawTofFromS2_TcalMult1->GetYaxis()->SetTitle("Counts");
    cTofFromS2->cd(1);
    fh1_RawTofFromS2_TcalMult1->Draw();
    fh1_RawTofFromS2 = new TH1D("RawTofNs_Tcal_Sci2_Los", "RawTofNs_Tcal_Sci2_Los", 10000000, -50000, 50000);
    fh1_RawTofFromS2->GetXaxis()->SetTitle("Raw Tof [ns]");
    fh1_RawTofFromS2->GetYaxis()->SetTitle("Counts");
    cTofFromS2->cd(2);
    fh1_RawTofFromS2->Draw();
    mainfolId->Add(cTofFromS2);
    // run->AddObject(cTofFromS2);

    cTofFromS2vsZ = new TCanvas("ZvsBeta_m1", "Z vs Beta", 10, 10, 800, 700);
    fh2_ZvsBeta_m1 = new TH2F("ZvsBeta_m1", "Z Music vs Beta", 3000, 0.75, 0.9, 2000, 42, 62);
    fh2_ZvsBeta_m1->GetXaxis()->SetTitle("Beta");
    fh2_ZvsBeta_m1->GetYaxis()->SetTitle("Z-Music");
    cTofFromS2vsZ->cd();
    fh2_ZvsBeta_m1->Draw("col");

    TCanvas* cTofFromS2vsZcor = new TCanvas("ZcorvsBeta_m1", "Zcor vs Beta", 10, 10, 800, 700);
    fh2_ZcorvsBeta_m1 = new TH2F("ZcorvsBeta_m1", "Zcor Music vs Beta", 3000, 0.75, 0.9, 2000, 42, 62);
    fh2_ZcorvsBeta_m1->GetXaxis()->SetTitle("Beta");
    fh2_ZcorvsBeta_m1->GetYaxis()->SetTitle("Z-Music corrBeta");
    cTofFromS2vsZcor->cd();
    fh2_ZcorvsBeta_m1->Draw("col");

    TCanvas* cZcorvsang = new TCanvas("ZcorvsAngle_m1", "Zcor vs Angle", 10, 10, 800, 700);
    fh2_MusZcorvsAng = new TH2F("ZcorvsAngle_m1", "Zcor Music vs Angle-Music", 1500, -20., 20., 2000, 42, 62);
    fh2_MusZcorvsAng->GetXaxis()->SetTitle("Angle [mrad]");
    fh2_MusZcorvsAng->GetYaxis()->SetTitle("Z-Music corrBeta");
    cZcorvsang->cd();
    fh2_MusZcorvsAng->Draw("col");

    TCanvas* cBetavsMusE = new TCanvas("BetavsMusE_m1", "MusicE vs Beta", 10, 10, 800, 700);
    fh2_EvsBeta_m1 = new TH2F("MusEvsBeta_m1", "E Music vs Beta", 3000, 0.75, 0.9, 5000, 2000, 62000);
    fh2_EvsBeta_m1->GetXaxis()->SetTitle("Beta");
    fh2_EvsBeta_m1->GetYaxis()->SetTitle("E-Music");
    cBetavsMusE->cd();
    fh2_EvsBeta_m1->Draw("col");

    // R3B-Music
    TCanvas* cMus_Z = new TCanvas("R3BMus_charge_z", "Mus: Charge Z", 10, 10, 800, 700);
    fh1_Mushit_z = new TH1F("fh1_R3BMus_charge_z", "Music: Charge Z", 60 * 20., 0, 60);
    fh1_Mushit_z->GetXaxis()->SetTitle("Charge (Z)");
    fh1_Mushit_z->GetYaxis()->SetTitle("Counts");
    fh1_Mushit_z->GetYaxis()->SetTitleOffset(1.1);
    fh1_Mushit_z->GetXaxis()->CenterTitle(true);
    fh1_Mushit_z->GetYaxis()->CenterTitle(true);
    fh1_Mushit_z->GetXaxis()->SetLabelSize(0.045);
    fh1_Mushit_z->GetXaxis()->SetTitleSize(0.045);
    fh1_Mushit_z->GetYaxis()->SetLabelSize(0.045);
    fh1_Mushit_z->GetYaxis()->SetTitleSize(0.045);
    fh1_Mushit_z->SetFillColor(28);
    fh1_Mushit_z->SetLineColor(1);
    cMus_Z->cd();
    fh1_Mushit_z->Draw("");

    cBeta = new TCanvas("Beta_Sci2_Los", "Beta_Sci2_Los", 10, 10, 800, 700);
    cBeta->Divide(1, 2);
    fh1_Beta_m1 = new TH1F("fh1_betaS2_m1", "BetaS2 to Cave-C with mult==1", 6000, 0.6, 0.9);
    fh1_Beta_m1->GetXaxis()->SetTitle("Beta");
    fh1_Beta_m1->GetYaxis()->SetTitle("Counts");
    fh1_Beta_m1->GetYaxis()->SetTitleOffset(1.15);
    fh1_Beta_m1->GetXaxis()->CenterTitle(true);
    fh1_Beta_m1->GetYaxis()->CenterTitle(true);
    fh1_Beta_m1->GetXaxis()->SetLabelSize(0.045);
    fh1_Beta_m1->GetXaxis()->SetTitleSize(0.045);
    fh1_Beta_m1->GetYaxis()->SetLabelSize(0.045);
    fh1_Beta_m1->GetYaxis()->SetTitleSize(0.045);
    fh1_Beta_m1->SetFillColor(2);
    fh1_Beta_m1->SetLineColor(1);
    cBeta->cd(1);
    fh1_Beta_m1->Draw("");
    fh1_Beta = new TH1F("fh1_betaS2", "BetaS2 to Cave-C", 6000, 0.6, 0.9);
    fh1_Beta->GetXaxis()->SetTitle("Beta");
    fh1_Beta->GetYaxis()->SetTitle("Counts");
    fh1_Beta->GetYaxis()->SetTitleOffset(1.15);
    fh1_Beta->GetXaxis()->CenterTitle(true);
    fh1_Beta->GetYaxis()->CenterTitle(true);
    fh1_Beta->GetXaxis()->SetLabelSize(0.045);
    fh1_Beta->GetXaxis()->SetTitleSize(0.045);
    fh1_Beta->GetYaxis()->SetLabelSize(0.045);
    fh1_Beta->GetYaxis()->SetTitleSize(0.045);
    fh1_Beta->SetFillColor(2);
    fh1_Beta->SetLineColor(1);
    cBeta->cd(2);
    fh1_Beta->Draw("");
    mainfolId->Add(cBeta);
    // run->AddObject(cBeta);

    cZvsAoQ = new TCanvas("ZvsAoQ", "ZvsAoQ", 10, 10, 800, 700);
    fh2_ZvsAoQ_m1 = new TH2F("fh2_ZvsAoQ", "Z-Music vs AoQ with mult==1", 1600, 2.35, 2.75, 2000, 42, 62);
    fh2_ZvsAoQ_m1->GetXaxis()->SetTitle("AoQ");
    fh2_ZvsAoQ_m1->GetYaxis()->SetTitle("Z [Charge units]");
    fh2_ZvsAoQ_m1->GetYaxis()->SetTitleOffset(1.1);
    fh2_ZvsAoQ_m1->GetXaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1->GetYaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1->GetXaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1->GetXaxis()->SetTitleSize(0.045);
    fh2_ZvsAoQ_m1->GetYaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1->GetYaxis()->SetTitleSize(0.045);
    fh2_ZvsAoQ_m1->Draw("colz");

    TCanvas* cZvsAoQcor = new TCanvas("ZvsAoQ_cor", "ZvsAoQ", 10, 10, 800, 700);
    fh2_ZvsAoQ_m1_cor =
        new TH2F("fh2_ZvsAoQ_cor", "Z-Music_betaCor vs AoQ with mult==1", 1600, 2.35, 2.75, 2000, 42, 62);
    fh2_ZvsAoQ_m1_cor->GetXaxis()->SetTitle("AoQ");
    fh2_ZvsAoQ_m1_cor->GetYaxis()->SetTitle("Z [Charge units]");
    fh2_ZvsAoQ_m1_cor->GetYaxis()->SetTitleOffset(1.1);
    fh2_ZvsAoQ_m1_cor->GetXaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1_cor->GetYaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1_cor->GetXaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1_cor->GetXaxis()->SetTitleSize(0.045);
    fh2_ZvsAoQ_m1_cor->GetYaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1_cor->GetYaxis()->SetTitleSize(0.045);
    fh2_ZvsAoQ_m1_cor->Draw("colz");

    TCanvas* cZvsAoQcor_tpat = new TCanvas("ZvsAoQ_cor_tpat", "ZvsAoQtpat", 10, 10, 800, 700);

    fh2_ZvsAoQ_m1_cor_tpat_tofd = new TH2F("fh2_ZvsAoQ_cor_tpat_tofd",
                                           "Z-Music_betaCor vs AoQ with mult==1 and (Tpat&2)=2",
                                           1600,
                                           2.35,
                                           2.75,
                                           2000,
                                           42,
                                           62);
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetXaxis()->SetTitle("AoQ");
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetYaxis()->SetTitle("Z [Charge units]");
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetYaxis()->SetTitleOffset(1.1);
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetXaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetYaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetXaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetXaxis()->SetTitleSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetYaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_tofd->GetYaxis()->SetTitleSize(0.045);

    fh2_ZvsAoQ_m1_cor_tpat_islos = new TH2F(
        "fh2_ZvsAoQ_cor_tpat_islos", "Z-Music_betaCor vs AoQ with mult==1 and Tpat==1", 1600, 2.35, 2.75, 2000, 42, 62);
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetXaxis()->SetTitle("AoQ");
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetYaxis()->SetTitle("Z [Charge units]");
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetYaxis()->SetTitleOffset(1.1);
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetXaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetYaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetXaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetXaxis()->SetTitleSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetYaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_islos->GetYaxis()->SetTitleSize(0.045);

    fh2_ZvsAoQ_m1_cor_tpat_los = new TH2F("fh2_ZvsAoQ_cor_tpat_los",
                                          "Z-Music_betaCor vs AoQ with mult==1 and (Tpat&1)=1",
                                          1600,
                                          2.35,
                                          2.75,
                                          2000,
                                          42,
                                          62);
    fh2_ZvsAoQ_m1_cor_tpat_los->GetXaxis()->SetTitle("AoQ");
    fh2_ZvsAoQ_m1_cor_tpat_los->GetYaxis()->SetTitle("Z [Charge units]");
    fh2_ZvsAoQ_m1_cor_tpat_los->GetYaxis()->SetTitleOffset(1.1);
    fh2_ZvsAoQ_m1_cor_tpat_los->GetXaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1_cor_tpat_los->GetYaxis()->CenterTitle(true);
    fh2_ZvsAoQ_m1_cor_tpat_los->GetXaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_los->GetXaxis()->SetTitleSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_los->GetYaxis()->SetLabelSize(0.045);
    fh2_ZvsAoQ_m1_cor_tpat_los->GetYaxis()->SetTitleSize(0.045);

    fh2_AoQ_cut = new TH1F("fh2_AoQ_cut", "fh2_AoQ_cut", 800, 2.45, 2.75);

    cZvsAoQcor_tpat->Divide(2, 2);
    cZvsAoQcor_tpat->cd(1);
    fh2_ZvsAoQ_m1_cor_tpat_los->Draw("colz");
    cZvsAoQcor_tpat->cd(2);
    fh2_ZvsAoQ_m1_cor_tpat_tofd->Draw("colz");
    cZvsAoQcor_tpat->cd(3);
    fh2_ZvsAoQ_m1_cor_tpat_islos->Draw("colz");
    cZvsAoQcor_tpat->cd(4);
    fh2_AoQ_cut->Draw();

    TCanvas* cAoQvsPosS2 = new TCanvas("AoQvsPosS2", "Pos-S2 vs AoQ", 10, 10, 800, 700);
    fh2_Pos2vsAoQ_m1 = new TH2F("fh2_PosS2vsAoQ", "Pos-S2 vs AoQ with mult==1", 2000, -100, 100, 1600, 2.35, 2.75);
    fh2_Pos2vsAoQ_m1->GetXaxis()->SetTitle("Pos-S2 [mm]");
    fh2_Pos2vsAoQ_m1->GetYaxis()->SetTitle("AoQ");
    fh2_Pos2vsAoQ_m1->GetYaxis()->SetTitleOffset(1.1);
    fh2_Pos2vsAoQ_m1->GetXaxis()->CenterTitle(true);
    fh2_Pos2vsAoQ_m1->GetYaxis()->CenterTitle(true);
    fh2_Pos2vsAoQ_m1->GetXaxis()->SetLabelSize(0.045);
    fh2_Pos2vsAoQ_m1->GetXaxis()->SetTitleSize(0.045);
    fh2_Pos2vsAoQ_m1->GetYaxis()->SetLabelSize(0.045);
    fh2_Pos2vsAoQ_m1->GetYaxis()->SetTitleSize(0.045);
    fh2_Pos2vsAoQ_m1->Draw("colz");

    if (fHitItemsMus)
    {
        mainfolId->Add(cMus_Z);
        mainfolId->Add(cTofFromS2vsZ);
        mainfolId->Add(cTofFromS2vsZcor);
        mainfolId->Add(cBetavsMusE);
        mainfolId->Add(cZcorvsang);
        mainfolId->Add(cZvsAoQ);
        mainfolId->Add(cZvsAoQcor);
        mainfolId->Add(cAoQvsPosS2);
    }
    run->AddObject(mainfolId);

    // Trigger and Tpat
    fhTpat = new TH1F("Tpat", "Tpat", 20, 0, 20);
    fhTpat->GetXaxis()->SetTitle("Tpat value");

    fhTrigger = new TH1F("Trigger", "Trigger all", 20, 0, 20);
    fhTrigger->GetXaxis()->SetTitle("Trigger value");

    TCanvas* cTrigg = new TCanvas("Trigg", "Triggers", 10, 10, 650, 350);
    cTrigg->Divide(2, 1);
    cTrigg->cd(1);
    gPad->SetLogy();
    fhTrigger->Draw();
    cTrigg->cd(2);
    gPad->SetLogy();
    fhTpat->Draw();
    cTrigg->cd(0);

    //------------------------------------------------------------------------
    // Los detector

    fNofLosDetectors = 1;

    // MAIN FOLDER-LOS
    TFolder* mainfol = new TFolder("LOS", "LOS info");

    TCanvas* cLos[fNofLosDetectors];

    if (fMappedItems.at(DET_LOS))
    {

        for (Int_t iloscount = 0; iloscount < fNofLosDetectors; iloscount++)
        {

            char detName[255];
            sprintf(detName, "LOS%d", iloscount + 1);

            cLos[iloscount] = new TCanvas(detName, detName, 10, 10, 1010, 810);

            fh_los_channels[iloscount] =
                new TH1F(Form("%s_channels", detName), Form("%s channels", detName), 10, 0., 10.);
            fh_los_channels[iloscount]->GetXaxis()->SetTitle("Channel number");

            fh_los_multihit[iloscount] =
                new TH1F(Form("%s_multihit", detName), Form("%s multihit && all 8 PMs", detName), 30, 0., 30.);
            fh_los_multihit[iloscount]->GetXaxis()->SetTitle("Multihit");

            fh_los_pos_MCFD[iloscount] =
                new TH2F(Form("%s_pos_MCFD", detName), Form("%s MCFD Position ", detName), 500, -5., 5., 500, -5., 5.);
            fh_los_pos_MCFD[iloscount]->GetXaxis()->SetTitle("X position / cm");
            fh_los_pos_MCFD[iloscount]->GetYaxis()->SetTitle("Y position / cm");

            fh_los_pos_TAMEX[iloscount] = new TH2F(
                Form("%s_pos_TAMEX", detName), Form("%s TAMEX Position ", detName), 500, -5., 5., 500, -5., 5.);
            fh_los_pos_TAMEX[iloscount]->GetXaxis()->SetTitle("X position / cm");
            fh_los_pos_TAMEX[iloscount]->GetYaxis()->SetTitle("Y position / cm");

            fh_los_pos_ToT[iloscount] =
                new TH2F(Form("%s_pos_ToT", detName), Form("%s ToT Position ", detName), 500, -5., 5., 500, -5., 5.);
            fh_los_pos_ToT[iloscount]->GetXaxis()->SetTitle("X position / cm");
            fh_los_pos_ToT[iloscount]->GetYaxis()->SetTitle("Y position / cm");

            fh_los_dt_hits_ToT[iloscount] =
                new TH2F(Form("%s_dt_ToT", detName), Form("%s ToT dt ", detName), 4000, -4., 4., 300, 0, 300.);
            fh_los_dt_hits_ToT[iloscount]->GetXaxis()->SetTitle("Time MCFD / ns"); // dt between two hits / ns
            fh_los_dt_hits_ToT[iloscount]->GetYaxis()->SetTitle("ToT / ns");       // ToT / ns

            fh_los_tres_MCFD[iloscount] = new TH1F(Form("%s_dt_4vs4_MCFD", detName),
                                                   Form("%s MCFD Time resolution - 4pmts vs 4pmts", detName),
                                                   4000,
                                                   -4.,
                                                   4.);
            fh_los_tres_MCFD[iloscount]->GetXaxis()->SetTitle("Time MCFD / ns");

            fh_los_tres_TAMEX[iloscount] = new TH1F(Form("%s_dt_4vs4_TAMEX", detName),
                                                    Form("%s TAMEX Time resolution - 4pmts vs 4pmts ", detName),
                                                    4000,
                                                    -4.,
                                                    4.);
            fh_los_tres_TAMEX[iloscount]->GetXaxis()->SetTitle("Time TAMEX / ns");

            fh_los_tot[iloscount] =
                new TH2F(Form("%s_tot", detName), Form("%s ToT vs PMT", detName), 10, 0, 10, 1500, 0., 300.);
            fh_los_tot[iloscount]->GetXaxis()->SetTitle("PMT number");
            fh_los_tot[iloscount]->GetYaxis()->SetTitle("ToT / ns");

            fh_los_tot_mean[iloscount] =
                new TH1F(Form("%s_tot_mean", detName), Form("%s mean ToT", detName), 1500, 0., 300.);
            fh_los_tot_mean[iloscount]->GetYaxis()->SetTitle("Counts");
            fh_los_tot_mean[iloscount]->GetXaxis()->SetTitle("ToT / ns");

            fh_los_ihit_ToT[iloscount] =
                new TH2F(Form("%s_tot_ihit", detName), Form("%s ToT vs ihit", detName), 10, 0, 10, 600, 0., 300.);
            // fh_los_ihit_ToT[iloscount]= new TH2F(Form("%s_tot_ihit",detName),Form("%s ToT vs
            // ihit",detName),15,0,15,300,0.,300.);
            fh_los_ihit_ToT[iloscount]->GetXaxis()->SetTitle("iHit");
            fh_los_ihit_ToT[iloscount]->GetYaxis()->SetTitle("ToT / ns");

            fh_losToT_vs_Events[iloscount] = new TH2F(Form("%s_tot_vs_event", detName),
                                                      Form("%s ToT vs. Event #", detName),
                                                      10000,
                                                      0,
                                                      10000000,
                                                      300,
                                                      0.,
                                                      300.);
            fh_losToT_vs_Events[iloscount]->GetYaxis()->SetTitle("ToT / ns");
            fh_losToT_vs_Events[iloscount]->GetXaxis()->SetTitle("Event number");

            fh_losTAMEX_vs_Events[iloscount] = new TH2F(Form("%s_dtTAMEX_vs_event", detName),
                                                        Form("%s dtTAMEX vs. Event #", detName),
                                                        10000,
                                                        0,
                                                        10000000,
                                                        1000,
                                                        -4.,
                                                        4.);
            fh_losTAMEX_vs_Events[iloscount]->GetYaxis()->SetTitle("Time TAMEX  / ns");
            fh_losTAMEX_vs_Events[iloscount]->GetXaxis()->SetTitle("Event number");

            fh_losMCFD_vs_Events[iloscount] = new TH2F(Form("%s_dtMCFD_vs_event", detName),
                                                       Form("%s dtMCFD vs. Event #", detName),
                                                       10000,
                                                       0,
                                                       10000000,
                                                       1000,
                                                       -4.,
                                                       4.);
            fh_losMCFD_vs_Events[iloscount]->GetYaxis()->SetTitle("Time MCFD  / ns");
            fh_losMCFD_vs_Events[iloscount]->GetXaxis()->SetTitle("Event number");

            fh_los_vftx_tamex[iloscount] =
                new TH1F(Form("%s_vftx_tamex", detName), Form("%s vftx_tamex", detName), 100000, -5000, 5000.);
            fh_los_vftx_tamex[iloscount]->GetYaxis()->SetTitle("Counts");
            fh_los_vftx_tamex[iloscount]->GetXaxis()->SetTitle("Ttamex-Tvftx / ns");

            cLos[iloscount]->Divide(3, 3);
            cLos[iloscount]->cd(1);
            fh_los_channels[iloscount]->Draw();
            cLos[iloscount]->cd(2);
            gPad->SetLogy();
            fh_los_multihit[iloscount]->Draw();
            cLos[iloscount]->cd(3);
            gPad->SetLogz();
            fh_los_tot[iloscount]->Draw("colz");
            cLos[iloscount]->cd(4); // gPad->SetLogy();
            fh_los_tot_mean[iloscount]->Draw();
            cLos[iloscount]->cd(5); // gPad->SetLogy();
            fh_los_tres_MCFD[iloscount]->Draw();
            cLos[iloscount]->cd(6);
            gPad->SetLogz();
            fh_los_dt_hits_ToT[iloscount]->Draw("colz");
            cLos[iloscount]->cd(7);
            gPad->SetLogz();
            fh_los_pos_ToT[iloscount]->Draw("colz");
            cLos[iloscount]->cd(8);
            gPad->SetLogz();
            fh_los_pos_MCFD[iloscount]->Draw("colz");
            cLos[iloscount]->cd(9);
            gPad->SetLogz();
            fh_los_ihit_ToT[iloscount]->Draw("colz");
            cLos[iloscount]->cd(0);
            mainfol->Add(cLos[iloscount]);
        }

        run->AddObject(mainfol);
        run->GetHttpServer()->RegisterCommand("Reset_LosVsSci2",
                                              Form("/Objects/%s/->Reset_LosVsSci2_Histo()", GetName()));
    }

    // -------------------------------------------------------------------------

    return kSUCCESS;
}

void R3BOnlineSpectraLosVsSci2::Reset_LosVsSci2_Histo()
{

    if (fMappedItems.at(DET_LOS))
    {
        for (Int_t iloscount = 0; iloscount < fNofLosDetectors; iloscount++)
        {
            fh_los_channels[iloscount]->Reset();
            fh_los_tres_MCFD[iloscount]->Reset();
            fh_los_tres_TAMEX[iloscount]->Reset();
            fh_los_pos_ToT[iloscount]->Reset();
            fh_los_tot[iloscount]->Reset();
            fh_los_tot_mean[iloscount]->Reset();
            fh_los_pos_MCFD[iloscount]->Reset();
            fh_los_pos_TAMEX[iloscount]->Reset();
            fh_los_multihit[iloscount]->Reset();
            fh_los_ihit_ToT[iloscount]->Reset();
            fh_losTAMEX_vs_Events[iloscount]->Reset();
            fh_losToT_vs_Events[iloscount]->Reset();
            fh_losMCFD_vs_Events[iloscount]->Reset();
            fh_los_dt_hits_ToT[iloscount]->Reset();
        }
        if (fTcalSci2)
        {
            fh1_RawTofFromS2_TcalMult1->Reset();
            fh1_Beta_m1->Reset();
            fh1_RawPos_m1->Reset();
            fh1_CalPos_m1->Reset();
            fh1_RawTofFromS2->Reset();
            fh1_Beta->Reset();
            fh1_RawPos->Reset();
            fh1_CalPos->Reset();
        }
    }

    if (fHitItemsMus)
    {
        fh1_Mushit_z->Reset();
        fh2_ZvsAoQ_m1->Reset();
        fh2_ZvsAoQ_m1_cor->Reset();
        fh2_ZvsAoQ_m1_cor_tpat_tofd->Reset();
        fh2_ZvsAoQ_m1_cor_tpat_islos->Reset();
        fh2_ZvsAoQ_m1_cor_tpat_los->Reset();
        fh2_AoQ_cut->Reset();
        fh2_EvsBeta_m1->Reset();
        fh2_MusZcorvsAng->Reset();
        fh2_ZcorvsBeta_m1->Reset();
        fh2_ZvsBeta_m1->Reset();
        fh2_Pos2vsAoQ_m1->Reset();
    }
}

void R3BOnlineSpectraLosVsSci2::Exec(Option_t* option)
{
    double Zmusic = 0., Music_ang = 0.;
    if (fHitItemsMus && fHitItemsMus->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitItemsMus->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusicHitData* hit = dynamic_cast<R3BMusicHitData*>(fHitItemsMus->At(ihit));
            if (!hit)
                continue;
            Zmusic = hit->GetZcharge();
            Music_ang = hit->GetTheta() * 1000.;
            // fh1_Mushit_z->Fill(Zmusic);
        }
    }

    // --- --------------- --- //
    // --- local variables --- //
    // --- --------------- --- //
    Int_t multTcal[3];
    Double_t iRawTimeNs[3][64];
    UInt_t nHits, iCh, cpt = 0;
    Double_t ToFraw = 0., PosRaw = -10000., ToFcal = 0., PosCal = -10000.;
    Double_t Velo = 0., Beta = 0., Gamma = 0., Brho = 0., AoQ = 0.;
    Double_t ToFraw_m1 = 0., PosRaw_m1 = -10000., ToFcal_m1 = 0., PosCal_m1 = -10000.;
    Double_t Velo_m1 = 0., Beta_m1 = 0., Gamma_m1 = 0., Brho_m1 = 0., AoQ_m1 = 0.;
    for (UShort_t pmt = 0; pmt < 3; pmt++)
    {
        multTcal[pmt] = 0;
        for (UShort_t m = 0; m < 64; m++)
            iRawTimeNs[pmt][m] = 0.;
    }

    // --- -------------- --- //
    // --- read tcal data --- //
    // --- -------------- --- //
    if (fTcalSci2 && fTcalSci2->GetEntriesFast())
    {
        nHits = fTcalSci2->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSci2TcalData* hittcal = dynamic_cast<R3BSci2TcalData*>(fTcalSci2->At(ihit));
            if (!hittcal)
                continue;
            iCh = hittcal->GetChannel() - 1;
            iRawTimeNs[iCh][multTcal[iCh]] = hittcal->GetRawTimeNs();
            multTcal[iCh]++;
        } // --- end of loop over Tcal data --- //
    }

    time = header->GetTimeStamp(); // / 1.6; // divided by 1.6 for stand alone daq with Vulom generating time stamps

    // time = 0;

    if (time_start == 0 && time > 0)
    {
        time_start = time;
        fNEvents_start = fNEvents;
    }

    if (header->GetTrigger() == 12)
        time_spill_start = time; // header->GetTimeStamp();    // spill start in nsec
    if (header->GetTrigger() == 13)
        time_spill_end = time; // header->GetTimeStamp();    // spill end  in nsec

    fhTrigger->Fill(header->GetTrigger());

    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    if (header->GetTrigger() == 12)
        cout << "Spill start: " << double(time_spill_start - time_start) / 1.e9 << " sec" << endl;
    if (header->GetTrigger() == 13)
        cout << "Spill stop: " << double(time_spill_end - time_start) / 1.e9 << " sec" << endl;

    Int_t tpatbin;
    for (int i = 0; i < 16; i++)
    {
        tpatbin = (header->GetTpat() & (1 << i));
        if (tpatbin != 0)
        {
            // if (i != 11)
            //    cout << "Tpat = " << i + 1 << endl;
            fhTpat->Fill(i + 1);
        }
    }

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    Int_t itpat;
    Int_t tpatvalue;
    if (fTpat_bit >= 0)
    {
        itpat = header->GetTpat();
        tpatvalue = (itpat && (1 < fTpat_bit)) >> fTpat_bit;
        if ((tpatvalue == 0))
            return;
    }

    if (fMappedItems.at(DET_SAMPLER))
    {

        auto det = fMappedItems.at(DET_SAMPLER);
        Int_t nHitsSamp = det->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHitsSamp; ihit++)
        {
            auto hit = dynamic_cast<R3BSamplerMappedData*>(det->At(ihit));
            // time is in steps of 10 ns
            // is is a 34 bit number, so max 1073741823
            samplerCurr = hit->GetTime();
            //		cout<<"Sampler vs Tstemp: " << samplerCurr<<", "<<time<<", "<<fNEvents<<endl;

            /*
                        if(spill_on_sampler && ihit == 0) {
                            samplerSpill = samplerCurr;
                            spill_on_sampler=false;
                        }


                        long samp = (samplerCurr - samplerSpill); //time in 10 ns
                        if(samp<0) samp += 1073741823;
                        if((double)samp/1e8 > fSpillLength && spillCounter>0) {
                            cout <<"Missed spill end!" << endl;
                            std::cout << samplerCurr << "  " << samplerPrev << "  " << samplerSpill<< "  " <<
               (double)samp/1e8 << '\n';
                        }
            */
            Double_t dt = ((double)(samplerCurr - samplerPrev)) / 100.;

            samplerPrev = samplerCurr;
        }
    }

    //----------------------------------------------------------------------
    // LOS detector
    //----------------------------------------------------------------------
    // Time nomenclature for LOS:
    // lt=0, l=1,lb=2,b=3,rb=4,r=5,rt=6,t=7
    //   t / \lt
    //   rt| |l
    //   r | |lb
    //   rb\ /b

    Double_t timeTofd = 0;

    Double_t timeLosV[fNofLosDetectors][32]; // mean time from VFTX
    Double_t LosTresV[fNofLosDetectors][32];
    Double_t timeLosT[fNofLosDetectors][32];
    Double_t LosTresT[fNofLosDetectors][32];
    Double_t timeLos[fNofLosDetectors][32];
    Double_t totsum[fNofLosDetectors][32];
    Double_t xT_cm[fNofLosDetectors][32];
    Double_t yT_cm[fNofLosDetectors][32];
    Double_t xToT_cm[fNofLosDetectors][32];
    Double_t yToT_cm[fNofLosDetectors][32];
    Double_t xV_cm[fNofLosDetectors][32];
    Double_t yV_cm[fNofLosDetectors][32];

    Double_t time_V[fNofLosDetectors][32][8]; // [det][multihit][pm]
    Double_t time_L[fNofLosDetectors][32][8];
    Double_t time_T[fNofLosDetectors][32][8];
    Double_t tot[fNofLosDetectors][32][8];
    Double_t time_MTDC[32][8] = { { 0. } };
    Double_t LosTresMTDC[32];

    for (Int_t idet = 0; idet < fNofLosDetectors; idet++)
    {
        for (Int_t imult = 0; imult < 32; imult++)
        {
            timeLosV[idet][imult] = 0.0;
            LosTresV[idet][imult] = 0.0 / 0.0;
            timeLosT[idet][imult] = 0.0;
            LosTresT[idet][imult] = 0.0 / 0.0;
            timeLos[idet][imult] = 0.0;
            totsum[idet][imult] = 0.0;
            xT_cm[idet][imult] = 0.0 / 0.0;
            yT_cm[idet][imult] = 0.0 / 0.0;
            xToT_cm[idet][imult] = -100000.;
            yToT_cm[idet][imult] = -100000.;
            xV_cm[idet][imult] = 0.0 / 0.0;
            yV_cm[idet][imult] = 0.0 / 0.0;
            for (Int_t icha = 0; icha < 8; icha++)
            {
                time_V[idet][imult][icha] = 0.0 / 0.0; // [det][multihit][pm]
                time_L[idet][imult][icha] = 0.0 / 0.0;
                time_T[idet][imult][icha] = 0.0 / 0.0;
                tot[idet][imult][icha] = 0.0 / 0.0;
            }
        }
    }
    Int_t nPartLOS = 0;
    Int_t nPartc[fNofLosDetectors];
    for (Int_t d = 0; d < fNofLosDetectors; d++)
        nPartc[d] = 0;

    Bool_t iLOSType[fNofLosDetectors][32];
    Bool_t iLOSPileUp[fNofLosDetectors][32];
    for (Int_t idet = 0; idet < fNofLosDetectors; idet++)
    {
        for (Int_t imult = 0; imult < 32; imult++)
        {
            iLOSType[idet][imult] = false;
            iLOSPileUp[idet][imult] = false;
        }
    }

    Int_t Multip;

    if (fMappedItems.at(DET_LOS))
    {
        auto det = fMappedItems.at(DET_LOS);
        nHits = det->GetEntriesFast();

        //       cout<<"nHits in LOS: "<<nHits<<endl;
        if (nHits > 0)
            nLosEvents += 1;
        Multip = nHits;

        //    cout<<"LOS: nHits"<<nHits<<endl;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BLosMappedData* hit = dynamic_cast<R3BLosMappedData*>(det->At(ihit));
            if (!hit)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hit->GetDetector(); // 1..
            Int_t iCha = hit->GetChannel();  // 1..
            Int_t iTyp = hit->GetType();     // 0,1,2,3

            fh_los_channels[iDet - 1]->Fill(iCha); // exclude MTDC data
        }
    }

    if (fCalItems.at(DET_LOS))
    {
        auto det = fCalItems.at(DET_LOS);
        nPartLOS = det->GetEntriesFast();

        Int_t iDet = 0;
        Double_t time_V_LOS1[32][8] = { { 0. } };
        Double_t time_V_LOS2[32][8] = { { 0. } };

        for (Int_t iPart = 0; iPart < nPartLOS; iPart++)
        {
            /*
             * nPart is the number of particle passing through LOS detector in one event
             */
            R3BLosCalData* calData = dynamic_cast<R3BLosCalData*>(det->At(iPart));
            iDet = calData->GetDetector();

            Double_t sumvtemp = 0, sumltemp = 0, sumttemp = 0;
            for (Int_t iCha = 0; iCha < 8; iCha++)
            {
                sumvtemp += calData->GetTimeV_ns(iCha);
                sumltemp += calData->GetTimeL_ns(iCha);
                sumttemp += calData->GetTimeT_ns(iCha);
            }
            if (!(IS_NAN(sumvtemp)) && !(IS_NAN(sumltemp)) && !(IS_NAN(sumltemp)))
            {
                nPartc[iDet - 1]++;
            }
            else
            {
                continue;
            }

            for (Int_t iCha = 0; iCha < 8; iCha++)
            {
                if (iDet == 1)
                {
                    time_V_LOS1[nPartc[iDet - 1] - 1][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeV_ns(iCha))))
                    { // VFTX
                        time_V_LOS1[nPartc[iDet - 1] - 1][iCha] = calData->GetTimeV_ns(iCha);
                    }
                    if (!(IS_NAN(calData->GetTimeL_ns(iCha))))
                    { // TAMEX leading
                        time_L[iDet - 1][nPartc[iDet - 1] - 1][iCha] = calData->GetTimeL_ns(iCha);
                    }

                    if (!(IS_NAN(calData->GetTimeT_ns(iCha))))
                    { // TAMEX trailing
                        time_T[iDet - 1][nPartc[iDet - 1] - 1][iCha] = calData->GetTimeT_ns(iCha);
                    }
                    time_MTDC[nPartc[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeM_ns(iCha))))
                    { // MTDC
                        time_MTDC[nPartc[iDet - 1] - 1][iCha] = calData->GetTimeM_ns(iCha);
                    }
                }
                if (iDet == 2)
                {
                    time_V_LOS2[nPartc[iDet - 1] - 1][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeV_ns(iCha))))
                    { // VFTX
                        time_V_LOS2[nPartc[iDet - 1] - 1][iCha] = calData->GetTimeV_ns(iCha);
                    }
                    if (!(IS_NAN(calData->GetTimeL_ns(iCha))))
                    { // TAMEX leading
                        time_L[iDet - 1][nPartc[iDet - 1] - 1][iCha] = calData->GetTimeL_ns(iCha);
                    }
                    if (!(IS_NAN(calData->GetTimeT_ns(iCha))))
                    { // TAMEX trailing
                        time_T[iDet - 1][nPartc[iDet - 1] - 1][iCha] = calData->GetTimeT_ns(iCha);
                    }
                }
            }

            if (!calData)
            {
                cout << "LOS !calData" << endl;
                continue; // can this happen?
            }
        }

        // Sorting VFTX data:

        // detector 1
        if (nPartc[0] > 0)
        {
            std::qsort(time_V_LOS1,
                       nPartc[0],
                       sizeof(*time_V_LOS1),
                       [](const void* arg1, const void* arg2) -> int
                       {
                           double const* lhs = static_cast<double const*>(arg1);
                           double const* rhs = static_cast<double const*>(arg2);

                           return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
                       });
            for (Int_t iPart = 0; iPart < nPartc[0]; iPart++)
            {
                for (int ipm = 0; ipm < 8; ipm++)
                {
                    time_V[0][iPart][ipm] = time_V_LOS1[iPart][ipm];
                }
            }
        }

        // detector 2
        if (fNofLosDetectors > 1 && nPartc[1] > 0)
        {
            std::qsort(time_V_LOS2,
                       nPartc[1],
                       sizeof(*time_V_LOS2),
                       [](const void* arg1, const void* arg2) -> int
                       {
                           double const* lhs = static_cast<double const*>(arg1);
                           double const* rhs = static_cast<double const*>(arg2);

                           return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
                       });
            for (Int_t iPart = 0; iPart < nPartc[1]; iPart++)
            {
                for (int ipm = 0; ipm < 8; ipm++)
                {
                    time_V[1][iPart][ipm] = time_V_LOS2[iPart][ipm];
                }
            }
        }

        // End sorting

        std::vector<double_t> time_first(fNofLosDetectors, -1.);
        std::vector<double_t> time0(fNofLosDetectors, -1.);
        std::vector<double_t> time1(fNofLosDetectors, -1.);
        std::vector<double_t> time_abs(fNofLosDetectors, -1.);

        for (iDet = 1; iDet <= fNofLosDetectors; iDet++)
        {
            cpt = 0;
            for (Int_t iPart = 0; iPart < nPartc[iDet - 1]; iPart++)
            {
                Bool_t iLOSTypeMCFD = false;
                Bool_t iLOSTypeTAMEX = false;

                if (time_V[iDet - 1][iPart][0] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][0])) &&
                    time_V[iDet - 1][iPart][1] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][1])) &&
                    time_V[iDet - 1][iPart][2] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][2])) &&
                    time_V[iDet - 1][iPart][3] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][3])) &&
                    time_V[iDet - 1][iPart][4] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][4])) &&
                    time_V[iDet - 1][iPart][5] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][5])) &&
                    time_V[iDet - 1][iPart][6] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][6])) &&
                    time_V[iDet - 1][iPart][7] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][7])))
                {
                    iLOSTypeMCFD = true; // all 8 MCFD times
                }

                if (time_L[iDet - 1][iPart][0] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][0])) &&
                    time_L[iDet - 1][iPart][1] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][1])) &&
                    time_L[iDet - 1][iPart][2] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][2])) &&
                    time_L[iDet - 1][iPart][3] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][3])) &&
                    time_L[iDet - 1][iPart][4] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][4])) &&
                    time_L[iDet - 1][iPart][5] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][5])) &&
                    time_L[iDet - 1][iPart][6] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][6])) &&
                    time_L[iDet - 1][iPart][7] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][7])) &&

                    time_T[iDet - 1][iPart][0] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][0])) &&
                    time_T[iDet - 1][iPart][1] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][1])) &&
                    time_T[iDet - 1][iPart][2] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][2])) &&
                    time_T[iDet - 1][iPart][3] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][3])) &&
                    time_T[iDet - 1][iPart][4] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][4])) &&
                    time_T[iDet - 1][iPart][5] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][5])) &&
                    time_T[iDet - 1][iPart][6] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][6])) &&
                    time_T[iDet - 1][iPart][7] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][7])))
                {
                    iLOSTypeTAMEX = true; // all 8 leading and trailing times
                }

                // We will consider only events in which booth MCFD and TAMEX see same number of channels:
                if (iLOSTypeTAMEX && iLOSTypeMCFD)
                    iLOSType[iDet - 1][iPart] = true;

                // if(iDet < fNofLosDetectors+1)
                // {

                if (iLOSType[iDet - 1][iPart])
                {
                    int nPMT = 0;
                    int nPMV = 0;

                    for (int ipm = 0; ipm < 8; ipm++)
                    {

                        if (time_T[iDet - 1][iPart][ipm] > 0. && time_L[iDet - 1][iPart][ipm] > 0. &&
                            !(IS_NAN(time_T[iDet - 1][iPart][ipm])) && !(IS_NAN(time_L[iDet - 1][iPart][ipm])))
                        {
                            while (time_T[iDet - 1][iPart][ipm] - time_L[iDet - 1][iPart][ipm] <= 0.)
                            {
                                time_T[iDet - 1][iPart][ipm] = time_T[iDet - 1][iPart][ipm] + 2048. * fClockFreq;
                            }

                            nPMT = nPMT + 1;
                            tot[iDet - 1][iPart][ipm] = time_T[iDet - 1][iPart][ipm] - time_L[iDet - 1][iPart][ipm];

                            // pileup rejection
                            if (tot[iDet - 1][iPart][ipm] > fEpileup)
                                iLOSPileUp[iDet - 1][iPart] = true;
                        }

                        if (tot[iDet - 1][iPart][ipm] != 0. && !(IS_NAN(tot[iDet - 1][iPart][ipm])))
                            totsum[iDet - 1][iPart] += tot[iDet - 1][iPart][ipm];

                        if (time_L[iDet - 1][iPart][ipm] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][ipm])))
                            timeLosT[iDet - 1][iPart] += time_L[iDet - 1][iPart][ipm];

                        // Calculate detector time
                        if (time_V[iDet - 1][iPart][ipm] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][ipm])))
                        {
                            timeLosV[iDet - 1][iPart] += time_V[iDet - 1][iPart][ipm];
                            nPMV = nPMV + 1;
                        }

                        /*  if (fNEvents < 5000000)
                              myFile << setprecision(10) << iDet << " " << iPart << " " << ipm << " "
                                     << time_V[iDet - 1][iPart][ipm] << " " << time_L[iDet - 1][iPart][ipm] << " "
                                     << tot[iDet - 1][iPart][ipm] << endl;*/
                    }

                    totsum[iDet - 1][iPart] = totsum[iDet - 1][iPart] / nPMT;

                    timeLosV[iDet - 1][iPart] = timeLosV[iDet - 1][iPart] / nPMV;

                    timeLosT[iDet - 1][iPart] = timeLosT[iDet - 1][iPart] / nPMT;

                    timeLos[iDet - 1][iPart] = timeLosV[iDet - 1][iPart];

                    // cout<<"LOS: "<<iDet<<", "<<iPart<<"; "<<timeLos[iDet-1][iPart]<<endl;

                    LosTresV[iDet - 1][iPart] = ((time_V[iDet - 1][iPart][0] + time_V[iDet - 1][iPart][2] +
                                                  time_V[iDet - 1][iPart][4] + time_V[iDet - 1][iPart][6]) -
                                                 (time_V[iDet - 1][iPart][1] + time_V[iDet - 1][iPart][3] +
                                                  time_V[iDet - 1][iPart][5] + time_V[iDet - 1][iPart][7])) /
                                                4.;

                    LosTresT[iDet - 1][iPart] = ((time_L[iDet - 1][iPart][0] + time_L[iDet - 1][iPart][2] +
                                                  time_L[iDet - 1][iPart][4] + time_L[iDet - 1][iPart][6]) -
                                                 (time_L[iDet - 1][iPart][1] + time_L[iDet - 1][iPart][3] +
                                                  time_L[iDet - 1][iPart][5] + time_L[iDet - 1][iPart][7])) /
                                                4.;

                    LosTresT[iDet - 1][iPart] = ((time_V[iDet - 1][iPart][0] * 0. + time_V[iDet - 1][iPart][2] +
                                                  time_V[iDet - 1][iPart][4] + time_V[iDet - 1][iPart][6]) /
                                                     3. -
                                                 (time_V[iDet - 1][iPart][1] + time_V[iDet - 1][iPart][3] +
                                                  time_V[iDet - 1][iPart][5] + time_V[iDet - 1][iPart][7]) /
                                                     4.);

                    // right koord.-system, Z-axis beam direction:
                    // Position from tamex:
                    xT_cm[iDet - 1][iPart] = (time_L[iDet - 1][iPart][1] + time_L[iDet - 1][iPart][2]) / 2. -
                                             (time_L[iDet - 1][iPart][5] + time_L[iDet - 1][iPart][6]) / 2.;
                    yT_cm[iDet - 1][iPart] = (time_L[iDet - 1][iPart][3] + time_L[iDet - 1][iPart][4]) / 2. -
                                             (time_L[iDet - 1][iPart][7] + time_L[iDet - 1][iPart][0]) / 2.;
                    xT_cm[iDet - 1][iPart] = (xT_cm[iDet - 1][iPart] - flosOffsetXT[iDet - 1]) * flosVeffXT[iDet - 1];
                    yT_cm[iDet - 1][iPart] = (yT_cm[iDet - 1][iPart] - flosOffsetYT[iDet - 1]) * flosVeffYT[iDet - 1];

                    // Position from VFTX:
                    xV_cm[iDet - 1][iPart] = (time_V[iDet - 1][iPart][1] + time_V[iDet - 1][iPart][2]) / 2. -
                                             (time_V[iDet - 1][iPart][5] + time_V[iDet - 1][iPart][6]) / 2.;
                    yV_cm[iDet - 1][iPart] = (time_V[iDet - 1][iPart][3] + time_V[iDet - 1][iPart][4]) / 2. -
                                             (time_V[iDet - 1][iPart][7] + time_V[iDet - 1][iPart][0]) / 2.;
                    xV_cm[iDet - 1][iPart] = (xV_cm[iDet - 1][iPart] - flosOffsetXV[iDet - 1]) * flosVeffXV[iDet - 1];
                    yV_cm[iDet - 1][iPart] = (yV_cm[iDet - 1][iPart] - flosOffsetYV[iDet - 1]) * flosVeffYV[iDet - 1];

                    // Position from ToT:
                    if (tot[iDet - 1][iPart][1] > 0. && tot[iDet - 1][iPart][2] > 0. && tot[iDet - 1][iPart][5] > 0. &&
                        tot[iDet - 1][iPart][6] > 0. && tot[iDet - 1][iPart][0] > 0. && tot[iDet - 1][iPart][3] > 0. &&
                        tot[iDet - 1][iPart][4] > 0. && tot[iDet - 1][iPart][7] > 0.)
                    {
                        xToT_cm[iDet - 1][iPart] = (((tot[iDet - 1][iPart][5] + tot[iDet - 1][iPart][6]) / 2. -
                                                     (tot[iDet - 1][iPart][1] + tot[iDet - 1][iPart][2]) / 2.) /
                                                    ((tot[iDet - 1][iPart][1] + tot[iDet - 1][iPart][2] +
                                                      tot[iDet - 1][iPart][5] + tot[iDet - 1][iPart][6]) /
                                                     4.));

                        yToT_cm[iDet - 1][iPart] = (((tot[iDet - 1][iPart][0] + tot[iDet - 1][iPart][7]) / 2. -
                                                     (tot[iDet - 1][iPart][3] + tot[iDet - 1][iPart][4]) / 2.) /
                                                    ((tot[iDet - 1][iPart][7] + tot[iDet - 1][iPart][0] +
                                                      tot[iDet - 1][iPart][3] + tot[iDet - 1][iPart][4]) /
                                                     4.));
                    }

                    xToT_cm[iDet - 1][iPart] =
                        (xToT_cm[iDet - 1][iPart] - flosOffsetXQ[iDet - 1]) * flosVeffXQ[iDet - 1];
                    yToT_cm[iDet - 1][iPart] =
                        (yToT_cm[iDet - 1][iPart] - flosOffsetYQ[iDet - 1]) * flosVeffYQ[iDet - 1];

                    // time difference between two hits
                    /*
                      if (time0[iDet - 1] < 0)
                          time0[iDet - 1] = timeLosV[iDet - 1][iPart];
                      time1[iDet - 1] = timeLosV[iDet - 1][iPart];
                      time_abs[iDet - 1] = time1[iDet - 1] - time0[iDet - 1]; // + (double)(time - time_prev[iDet - 1]);
                      if (time_abs[iDet - 1] > 0 && time_abs[iDet - 1] < 5.E8)
                      {
                          fh_los_dt_hits_ToT[iDet - 1]->Fill(time_abs[iDet - 1] / 1000., totsum[iDet - 1][iPart]);
                      }

                      time_prev[iDet - 1] = time;
                      time0[iDet - 1] = time1[iDet - 1];
                    */

                    if (timeLosV[iDet - 1][iPart] > 0. && timeLosV[iDet - 1][iPart] < 8192. * 5. &&
                        !(IS_NAN(timeLosV[iDet - 1][iPart])))
                    {
                        while (timeLosT[iDet - 1][iPart] - timeLosV[iDet - 1][iPart] < 2048. * 5. / 2.)
                        {
                            timeLosT[iDet - 1][iPart] = timeLosT[iDet - 1][iPart] + 2048. * 5.;
                        }
                        while (timeLosT[iDet - 1][iPart] - timeLosV[iDet - 1][iPart] > 2048. * 5. / 2.)
                        {
                            timeLosT[iDet - 1][iPart] = timeLosT[iDet - 1][iPart] - 2048. * 5.;
                        }

                        fh_los_vftx_tamex[iDet - 1]->Fill(timeLosT[iDet - 1][iPart] - timeLosV[iDet - 1][iPart]);
                    }

                    fh_losMCFD_vs_Events[iDet - 1]->Fill(fNEvents, LosTresV[iDet - 1][iPart]);
                    fh_losTAMEX_vs_Events[iDet - 1]->Fill(fNEvents, LosTresT[iDet - 1][iPart]);
                    fh_losToT_vs_Events[iDet - 1]->Fill(fNEvents, totsum[iDet - 1][iPart]);

                    if (!(iLOSPileUp[iDet - 1][iPart])) // desregard pile-up events
                    {
                        for (int ipm = 0; ipm < 8; ipm++)
                        {
                            fh_los_tot[iDet - 1]->Fill(ipm + 1, tot[iDet - 1][iPart][ipm]);
                        }
                        fh_los_dt_hits_ToT[iDet - 1]->Fill(LosTresV[iDet - 1][iPart], totsum[iDet - 1][iPart]);
                        fh_los_tot_mean[iDet - 1]->Fill(totsum[iDet - 1][iPart]);
                        fh_los_tres_MCFD[iDet - 1]->Fill(LosTresV[iDet - 1][iPart]);
                        fh_los_tres_TAMEX[iDet - 1]->Fill(LosTresT[iDet - 1][iPart]);
                        fh_los_pos_MCFD[iDet - 1]->Fill(xV_cm[iDet - 1][iPart], yV_cm[iDet - 1][iPart]);
                        fh_los_pos_TAMEX[iDet - 1]->Fill(xT_cm[iDet - 1][iPart], yT_cm[iDet - 1][iPart]);
                        fh_los_pos_ToT[iDet - 1]->Fill(xToT_cm[iDet - 1][iPart], yToT_cm[iDet - 1][iPart]);
                        fh_los_ihit_ToT[iDet - 1]->Fill(iPart + 1, totsum[iDet - 1][iPart]);
                        fh_los_multihit[iDet - 1]->Fill(iPart + 1);
                    }
                } // if iLosType

                // }
                // else
                // {
                // cout<<"Wrong detector ID for LOS!"<<endl;
                // }
                if (iDet == 1)
                {
                    for (Int_t multR = 0; multR < 64; multR++)
                    {
                        for (Int_t multL = 0; multL < 64; multL++)
                        {
                            ToFraw = timeLosV[0][iPart] - 0.5 * (iRawTimeNs[0][multR] + iRawTimeNs[1][multL]);
                            if (fToFmin <= ToFraw && ToFraw <= fToFmax)
                            {
                                cpt++;
                                Velo = 1. / (fTof2InvV_p0 + fTof2InvV_p1 * (fToFoffset + ToFraw)); // [m/ns]
                                Beta = Velo / 0.299792458;
                                Gamma = 1. / (TMath::Sqrt(1. - TMath::Power(Beta, 2)));
                                PosRaw = iRawTimeNs[0][multR] - iRawTimeNs[1][multL]; // [ns]
                                PosCal = fPos_p0 + fPos_p1 * PosRaw;                  // [mm] at S2
                                Brho = fBrho0_S2toCC * (1. - PosCal / fDispersionS2);
                                AoQ = Brho / (3.10716 * Beta * Gamma);
                            }
                        }
                    }
                }
            } // for iPart
            if (cpt == 1)
            {
                fh1_RawTofFromS2->Fill(ToFraw);
                fh1_RawPos->Fill(PosRaw);
                fh1_CalPos->Fill(PosCal);
                fh1_Beta->Fill(Beta);
            }

            // select multiplicity == 1 at Cave C and S2
            if (iDet == 1 && nPartc[0] == 1 && multTcal[0] == 1 && multTcal[1] == 1)
            {

                // --- -----------------------------
                // --- secondary beam identification
                // --- -----------------------------
                // if X is increasing from left to right:
                //    Brho = fBhro0 * (1 - xMwpc0/fDCC + xS2/fDS2)
                // in R3BRoot, X is increasing from right to left
                //    Bro = fBrho0 * (1 + xMwpc0/fDCC - xS2/fDS2)

                ToFraw_m1 = timeLosV[0][0] - 0.5 * (iRawTimeNs[0][0] + iRawTimeNs[1][0]);
                fh1_RawTofFromS2_TcalMult1->Fill(ToFraw_m1);
                Velo_m1 = 1. / (fTof2InvV_p0 + fTof2InvV_p1 * (fToFoffset + ToFraw_m1)); // [m/ns]
                Beta_m1 = Velo_m1 / 0.299792458;
                Gamma_m1 = 1. / (TMath::Sqrt(1. - TMath::Power(Beta_m1, 2)));
                PosRaw_m1 = iRawTimeNs[0][0] - iRawTimeNs[1][0]; // [ns]
                PosCal_m1 = fPos_p0 + fPos_p1 * PosRaw_m1;       // [mm] at S2
                Brho_m1 = fBrho0_S2toCC * (1. - PosCal_m1 / fDispersionS2);
                AoQ_m1 = Brho_m1 / (3.10716 * Beta_m1 * Gamma_m1);
                fh1_RawPos_m1->Fill(PosRaw_m1);
                fh1_CalPos_m1->Fill(PosCal_m1);
                fh1_Beta_m1->Fill(Beta_m1);
                if (Zmusic > 0.)
                {
                    // fh2_ZvsAoQ_m1->Fill(AoQ_m1, Zmusic);
                    fh2_ZvsBeta_m1->Fill(Beta_m1, Zmusic);
                    double Emus = ((Zmusic + 4.7) / 0.28) * ((Zmusic + 4.7) / 0.28);
                    // double zcor= fZprimary*sqrt( Emus/(fP0+fP1*Beta_m1+fP2*Beta_m1*Beta_m1) )+fZoffset;
                    double zcor = sqrt(Emus * Beta_m1) * 0.277;
                    fh2_EvsBeta_m1->Fill(Beta_m1, Emus);
                    fh1_Mushit_z->Fill(zcor);
                    fh2_ZcorvsBeta_m1->Fill(Beta_m1, zcor);
                    fh2_ZvsAoQ_m1->Fill(AoQ_m1, Zmusic);
                    double x0_point = 0.;
                    double y0_point = 50.39;
                    double rot_ang = 0.0375;
                    double zcorang =
                        y0_point + (Music_ang - x0_point) * sin(rot_ang) + (zcor - y0_point) * cos(rot_ang) + 0.2;
                    fh2_MusZcorvsAng->Fill(Music_ang, zcorang);
                    if (PosCal_m1 > fPos2min && PosCal_m1 < fPos2max)
                        fh2_ZvsAoQ_m1_cor->Fill(AoQ_m1, zcorang);
                    if (PosCal_m1 > fPos2min && PosCal_m1 < fPos2max &&
                        ((header->GetTpat() & 2) == 2 || (header->GetTpat() & 64) == 64))
                        fh2_ZvsAoQ_m1_cor_tpat_tofd->Fill(AoQ_m1, zcorang);
                    if (PosCal_m1 > fPos2min && PosCal_m1 < fPos2max &&
                        ((header->GetTpat() == 1) || (header->GetTpat() == 32)))
                        fh2_ZvsAoQ_m1_cor_tpat_islos->Fill(AoQ_m1, zcorang);
                    if (PosCal_m1 > fPos2min && PosCal_m1 < fPos2max &&
                        ((header->GetTpat() & 1) == 1 || ((header->GetTpat() & 32) == 32)))
                        fh2_ZvsAoQ_m1_cor_tpat_los->Fill(AoQ_m1, zcorang);
                    if (PosCal_m1 > fPos2min && PosCal_m1 < fPos2max && abs(zcorang - ZMUSIC_cut) < ZMUSIC_wcut)
                        fh2_AoQ_cut->Fill(AoQ_m1);
                    fh2_Pos2vsAoQ_m1->Fill(PosCal_m1, AoQ_m1);

                    // AddData(1, 2, zcorang, AoQ_m1, Beta_m1, Brho_m1, PosCal_m1, 0.);
                }
            }
        } // for iDet

    } // if fCallItems

    fNEvents += 1;
}

void R3BOnlineSpectraLosVsSci2::FinishEvent()
{

    for (Int_t det = 0; det < DET_MAX; det++)
    {
        if (fMappedItems.at(det))
        {
            fMappedItems.at(det)->Clear();
        }
        if (fCalItems.at(det))
        {
            fCalItems.at(det)->Clear();
        }
    }
    if (fTcalSci2)
        fTcalSci2->Clear();
    if (fHitItemsMus)
    {
        fHitItemsMus->Clear();
    }
}

void R3BOnlineSpectraLosVsSci2::FinishTask()
{

    if (fMappedItems.at(DET_LOS))
    {
        for (Int_t iloscount = 0; iloscount < fNofLosDetectors; iloscount++)
        {
            fh_los_channels[iloscount]->Write();
            fh_los_tot[iloscount]->Write();
            fh_los_dt_hits_ToT[iloscount]->Write();
            fh_los_ihit_ToT[iloscount]->Write();
            fh_los_tot_mean[iloscount]->Write();
            fh_los_tres_MCFD[iloscount]->Write();
            fh_los_tres_TAMEX[iloscount]->Write();
            fh_los_pos_MCFD[iloscount]->Write();
            fh_los_pos_TAMEX[iloscount]->Write();
            fh_los_pos_ToT[iloscount]->Write();
        }
    }
    fhTpat->Write();
    fhTrigger->Write();

    fh1_RawTofFromS2->Write();
    fh1_RawPos->Write();
    fh1_CalPos->Write();
    fh1_Beta->Write();
    fh1_RawTofFromS2_TcalMult1->Write();
    fh2_EvsBeta_m1->Write();
    fh1_RawPos_m1->Write();
    fh1_CalPos_m1->Write();
    fh1_Beta_m1->Write();
    fh1_Mushit_z->Write();
    fh2_MusZcorvsAng->Write();
    fh2_ZvsAoQ_m1->Write();
    fh2_ZvsAoQ_m1_cor->Write();
    fh2_ZvsAoQ_m1_cor_tpat_los->Write();
    fh2_ZvsBeta_m1->Write();
    fh2_ZcorvsBeta_m1->Write();
    fh2_Pos2vsAoQ_m1->Write();

    cout << "FinishTask: All events: " << fNEvents << ", LOS events: " << nLosEvents << endl;
}
/*
// -----   Private method AddData  --------------------------------------------
R3BFrsData* R3BOnlineSpectraLosVsSci2::AddData(Int_t StaId,
                                               Int_t StoId,
                                               Double_t z,
                                               Double_t aq,
                                               Double_t beta,
                                               Double_t brho,
                                               Double_t xs2,
                                               Double_t xc)
{
    // It fills the R3BSofFrsData
    TClonesArray& clref = *fFrsDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BFrsData(StaId, StoId, z, aq, beta, brho, xs2, xc);
}
*/
ClassImp(R3BOnlineSpectraLosVsSci2);
