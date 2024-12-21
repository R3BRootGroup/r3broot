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

// -----------------------------------------------------------
// -----             R3BRpcOnlineSpectra                 -----
// -----   Created 23/01/22 by J.L. Rodriguez-Sanchez    -----
// -----      Modified 24/01/22 by G.García-Jiménez
// -----         Fill RPC online histograms              -----
// -----------------------------------------------------------

#include "R3BRpcOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BRpcCalData.h"
#include "R3BRpcHitData.h"
#include "R3BRpcMappedData.h"
#include "R3BRpcPreCalData.h"
#include "R3BShared.h"

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
#include "THttpServer.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

R3BRpcOnlineSpectra::R3BRpcOnlineSpectra()
    : R3BRpcOnlineSpectra("RpcOnlineSpectra", 1)
{
}

R3BRpcOnlineSpectra::R3BRpcOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fEventHeader(nullptr)
    , fMappedDataItems(NULL)
    , fCalDataItems(NULL)
    , fHitDataItems(NULL)
    , fcounter(0)
    , fTrigger(-1)
    , fNEvents(0)
    , fLeftTofLim(-500)
    , fRightTofLim(500)
    , fRightRpcTimeLim(1000)
    , fLeftRpcTimeLim(-1000)
    , fRpcTimeBins(500)
    , fTofBins(500)
    , fRpcToTBins(400)
    , fLeftRpcToTLim(-400)
    , fRightRpcToTLim(400)
    , fTimeStart(0)
    , fTimeEnd(0)
    , fBarHits(0)
    , fStrip21Hits(0)
    , fSpill(1)
    , fFirstTPat(1)
    , fLastTPat(12)
{
}

R3BRpcOnlineSpectra::~R3BRpcOnlineSpectra()
{
    R3BLOG(debug1, "");

    if (fMappedDataItems)
        delete fMappedDataItems;

    if (fPreCalDataItems)
        delete fPreCalDataItems;

    if (fHitDataItems)
        delete fHitDataItems;
}

InitStatus R3BRpcOnlineSpectra::Init()
{
    R3BLOG(info, "");

    // Looking for FairRootManager
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    // Look for the R3BEventHeader
    fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!fEventHeader)
    {
        R3BLOG(warn, "EventHeader. not found");
        fEventHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    }
    else
        R3BLOG(info, "EventHeader. found");

    // Get access to Mapped data

    fMappedDataItems = dynamic_cast<TClonesArray*>(mgr->GetObject("R3BRpcMappedData"));
    if (!fMappedDataItems)
    {
        R3BLOG(fatal, "RpcMappedData not found");
        return kFATAL;
    }

    fPreCalDataItems = dynamic_cast<TClonesArray*>(mgr->GetObject("R3BRpcPreCalData"));
    if (!fPreCalDataItems)
    {
        R3BLOG(fatal, "R3BRpcPreCalData not found");
        return kFATAL;
    }

    fCalDataItems = dynamic_cast<TClonesArray*>(mgr->GetObject("R3BRpcCalData"));
    if (!fCalDataItems)
    {
        R3BLOG(fatal, "R3BRpcCalData not found");
        return kFATAL;
    }

    fHitDataItems = dynamic_cast<TClonesArray*>(mgr->GetObject("R3BRpcHitData"));
    if (!fHitDataItems)
    {
        R3BLOG(fatal, "RpcHitData not found");
        return kFATAL;
    }

    char name[100];
    /*
      6 Pmt's
      41 Strips per side
      2 Sides
      4 Ref's
    */

    /* ------------- FOLDERS --------------*/

    TFolder* mainfol = new TFolder("Rpc", "RPC_info");

    // Folder for mapped data
    TFolder* mapfol = new TFolder("Map", "Map_RPC_info");
    mainfol->Add(mapfol);

    // Folder for PreCal Data
    TFolder* preCalFolder = new TFolder("Pre Cal", "Pre_Cal_RPC_info");
    mainfol->Add(preCalFolder);

    // Folder for Cal Data
    TFolder* calFolder = new TFolder("Cal", "Cal_RPC_info");
    mainfol->Add(calFolder);

    // Folder for Hit Data
    TFolder* hitFolder = new TFolder("Hit", "Hit_RPC_info");
    mainfol->Add(hitFolder);

    TFolder* stripBarCorrFolder = new TFolder("Strip_Bar_Correlations", "Strip_Bar_Correlations");
    mainfol->Add(stripBarCorrFolder);

    TFolder* stripLeftFolder = new TFolder("Strip_Map:Left", "Strip_Map:Left");
    mapfol->Add(stripLeftFolder);

    TFolder* stripRightFolder = new TFolder("Strip_Map:Right", "Strip_Map:Right");
    mapfol->Add(stripRightFolder);

    TFolder* pmtFolder = new TFolder("Pmt_Map", "Pmt_Map");
    mapfol->Add(pmtFolder);

    TFolder* refFolder = new TFolder("Ref_Map", "Ref_Map");
    mapfol->Add(refFolder);

    TFolder* stripLeftPreCalFolder = new TFolder("Strip_Pre_Cal:Left", "Strip_Pre_Cal:Left");
    preCalFolder->Add(stripLeftPreCalFolder);

    TFolder* pmtPreCalFolder = new TFolder("Pmt_Pre_Cal", "Pmt_Pre_Cal");
    preCalFolder->Add(pmtPreCalFolder);

    TFolder* stripRightPreCalFolder = new TFolder("Strip_Pre_Cal:Right", "Strip_Pre_Cal:Right");
    preCalFolder->Add(stripRightPreCalFolder);

    TFolder* stripCalFolder = new TFolder("Strip_Cal", "Strip_Cal");
    calFolder->Add(stripCalFolder);

    TFolder* pmtCalFolder = new TFolder("Pmt_Cal", "Pmt_Cal");
    calFolder->Add(pmtCalFolder);

    TFolder* stripVsTofFolder = new TFolder("Channel_Vs_ToF", "Channel_Vs_ToF");
    hitFolder->Add(stripVsTofFolder);

    TFolder* efficiencyFolder = new TFolder("Efficiencies", "Efficiencies");
    mainfol->Add(efficiencyFolder);

    /* ----------- CANVASES ------------- */

    rightStripCanvasCoarse = new TCanvas("Right_Coarse_Canvas", "rightStripCanvasCoarse");
    leftStripCanvasCoarse = new TCanvas("Left_Coarse_Canvas", "leftStripCanvasCoarse");

    rightStripCanvasFine = new TCanvas("Right_Fine_Canvas", "rightStripCanvasFine");
    leftStripCanvasFine = new TCanvas("Left_Fine_Canvas", "leftStripCanvasFine");

    rightStripCanvasFine->Divide(9, 5);
    leftStripCanvasFine->Divide(9, 5);

    rightStripCanvasCoarse->Divide(9, 5);
    leftStripCanvasCoarse->Divide(9, 5);

    pmtCoarseCanvas = new TCanvas("PmtCoarseCanvas", "pmtCoarseCanvas");
    pmtCoarseCanvas->Divide(4, 2);

    pmtFineCanvas = new TCanvas("PmtFineCanvas", "pmtFineCanvas");
    pmtFineCanvas->Divide(4, 2);

    refFineCanvas = new TCanvas("RefFineCanvas", "refFineCanvas");
    refFineCanvas->Divide(3, 4);

    refCoarseCanvas = new TCanvas("RefCoarseCanvas", "refCoarseCanvas");
    refCoarseCanvas->Divide(3, 4);

    stripCoarseLeftCorrCanvas = new TCanvas("Left:Strip_Vs_Coarse", "strip_Coarse_Left_Corr_Canvas");
    stripCoarseRightCorrCanvas = new TCanvas("Right:Strip_Vs_Coarse", "strip_Coarse_Right_Corr_Canvas");

    stripFineLeftCorrCanvas = new TCanvas("Left:Strip_Vs_Fine", "strip_Fine_Left_Corr_Canvas");
    stripFineRightCorrCanvas = new TCanvas("Right:Strip_Vs_Fine", "strip_Fine_Right_Corr_Canvas");

    /* ----- Pre Cal Canvases ----- */
    stripLeftTotCorrCanvas = new TCanvas("Left:Strip_Vs_Tot", "strip_Vs_Tot");
    stripRightTotCorrCanvas = new TCanvas("Right:Strip_Vs_Tot", "strip_Vs_Tot");
    stripLeftTimeCorrCanvas = new TCanvas("Left:Strip_Vs_Time", "strip_Vs_Time");
    stripRightTimeCorrCanvas = new TCanvas("Right:Strip_Vs_Time", "strip_Vs_Time");

    pmtPreCalTimeCanvas = new TCanvas("Pmt_PreCal_Time", "pmt_Pre_Cal_Time");
    pmtPreCalTimeCanvas->Divide(4, 2);

    pmtPreCalTotCanvas = new TCanvas("Pmt_PreCal_ToT_Time", "pmt_Pre_Cal_ToT");
    pmtPreCalTotCanvas->Divide(4, 2);

    /* ----- Cal Canvases ----- */
    stripCalTimeCorrCanvas = new TCanvas("Left_Time_Vs_Right_Time", "left_Time_Vs_Right_Time");
    stripCalTotCorrCanvas = new TCanvas("Left_ToT_Vs_Right_ToT", "left_ToT_Vs_Right_ToT");

    /* ------ Hit Canvases ------*/
    hitMapCanvas = new TCanvas("Hit_Map", "hit_Map_Canvas");
    hitMapCanvas->Divide(2, 1);

    /* ------ ToF Correlations -------*/
    tofCorrCanvas = new TCanvas("tof_Corr_Canvas", "ToF_Correlations");
    stripTofCanvas = new TCanvas("strip_Tof_Canvas", "Strip_ToF_Correlations");
    stripTofCanvas->Divide(5, 9);

    tofPosCanvas = new TCanvas("tof_Pos_Canvas", "ToF_Pos_Correlations");

    /* ------ Strip Bar Correlations ------------*/
    timeDiffStripPmtCanvas = new TCanvas("time_Diff_Strip_Pmt_Canvas", "Time_Differences");
    timeDiffStripPmtCanvas->Divide(5, 9);

    timeDiffStripPmtCorrCanvas = new TCanvas("time_Diff_Strip_Pmt_Corr_Canvas", "Time_Differences");
    posRpcCorrCanvas = new TCanvas("pos_Rpc_Corr_Canvas", "Pos_Vs_Rpc");
    posRpcCorrCanvas->Divide(2, 1);

    stripIdVsHitsCanvas = new TCanvas("strip_Id_Vs_Hits_Canvas", "Hits_Vs_Strip_Id(H_Bar_AND_Good_Beam)");

    efficiencyCanvas = new TCanvas("efficiency_Canvas", "Efficiency_Canvas");
    efficiencyCanvas->Divide(3, 1);

    /* ----- Map Histograms -----*/

    // stripCoarseRightHisto = R3B::root_owned<TH1F>*[41];
    // stripFineRightHisto = R3B::root_owned<TH1F>*[41];

    // stripCoarseLeftHisto = R3B::root_owned<TH1F>*[41];
    // stripFineLeftHisto = R3B::root_owned<TH1F>*[41];

    // pmtCoarseHistoTop = R3B::root_owned<TH1F>*[4];
    // pmtFineHistoTop = R3B::root_owned<TH1F>*[4];

    // pmtCoarseHistoBottom = R3B::root_owned<TH1F>*[4];
    // pmtFineHistoBottom = R3B::root_owned<TH1F>*[4];

    // refCoarseHisto = R3B::root_owned<TH1F>*[9];
    // refFineHisto = R3B::root_owned<TH1F>*[9];

    stripCoarseLeftCorr =
        R3B::root_owned<TH2F>("strip_Coarse_Left", " Strip_Vs_Coarse_Time:Left", 41, 0.5, 41.5, 400, 0, 2500);
    stripCoarseRightCorr =
        R3B::root_owned<TH2F>("strip_Coarse_Right", " Strip_Vs_Coarse_Time:Right", 41, 0.5, 41.5, 400, 0, 2500);
    stripFineLeftCorr =
        R3B::root_owned<TH2F>("strip_Fine_Left", " Strip_Vs_Fine_Time:Left", 41, 0.5, 41.5, 400, 0, 600);
    stripFineRightCorr =
        R3B::root_owned<TH2F>("strip_Fine_Right", " Strip_Vs_Fine_Time:Right", 41, 0.5, 41.5, 400, 0, 600);

    /* ----- Pre Cal Histograms ----- */
    //    pmtPreCalTimeHistoTop = R3B::root_owned<TH1F>*[4];
    //    pmtPreCalTotHistoTop = R3B::root_owned<TH1F>*[4];
    //
    //    pmtPreCalTimeHistoBottom = R3B::root_owned<TH1F>*[4];
    //    pmtPreCalTotHistoBottom = R3B::root_owned<TH1F>*[4];

    stripLeftTotCorr = R3B::root_owned<TH2F>(
        "strip_Left_Tot_Corr", "Strip_Vs_Tot:Left", 41, 0.5, 41.5, fRpcToTBins, fLeftRpcToTLim, fRightRpcToTLim);
    stripRightTotCorr = R3B::root_owned<TH2F>(
        "strip_Right_Tot_Corr", "Strip_Vs_Tot:Right", 41, 0.5, 41.5, fRpcToTBins, fLeftRpcToTLim, fRightRpcToTLim);
    stripLeftTimeCorr = R3B::root_owned<TH2F>(
        "strip_Left_Time_Corr", "Strip_Vs_Time:Left", 41, 0.5, 41.5, fRpcTimeBins, fLeftRpcTimeLim, fRightRpcTimeLim);
    stripRightTimeCorr = R3B::root_owned<TH2F>(
        "strip_Right_Time_Corr", "Strip_Vs_Time:Right", 41, 0.5, 41.5, fRpcTimeBins, fLeftRpcTimeLim, fRightRpcTimeLim);

    /* ----- Cal Histograms ----- */
    stripCalTimeCorr =
        R3B::root_owned<TH2F>("strip_Cal_Time_Corr", "Strip:Time_Left_Vs_Time_Right", 800, -680, -580, 800, -800, -600);
    stripCalToTCorr =
        R3B::root_owned<TH2F>("strip_Cal_ToT_Corr", "Strip:ToT_Left_Vs_ToT_Right", 400, -100, -350, 400, -100, 350);

    /* ------------- HIT Histograms ------------ */
    stripPosHitCorr = R3B::root_owned<TH2F>("strip_Pos_Hit_Corr", "Strip_Vs_Position", 150, 0, 1500, 41, 0.5, 41.5);
    totalChargeHist = R3B::root_owned<TH1F>("total_Charge_Hist", "Charge", 1000, -100, 400);
    meanChargeCorr = R3B::root_owned<TH2F>("mean_Charge_Corr", "Heat_Map:Mean_Charge", 50, 0, 1500, 41, 0.5, 41.5);

    tofVsPosCorr =
        R3B::root_owned<TH2F>("tof_Vs_Pos_Corr", "ToF_Vs_Position", 600, 0, 1500, fTofBins, fLeftTofLim, fRightTofLim);
    tofCorr = R3B::root_owned<TH2F>("tof_Corr", "ToF_Vs_Strip", fTofBins, fLeftTofLim, fRightTofLim, 41, 0.5, 41.5);

    // stripTofHisto = R3B::root_owned<TH1F>*[41];

    // timeDiffStripPmtHisto = R3B::root_owned<TH1F>*[41];

    timeDiffStripPmtCorr = R3B::root_owned<TH2F>(
        "time_Diff_Strip_Pmt_Corr", "Time_Differences_Strip_and_Bar ", 1000, -1000, -1000, 41, 0.5, 41.5);
    deltaTVsPosCorr =
        R3B::root_owned<TH2F>("deltaT_Vs_Pos_Corr", "Time_Difference_Vs_Position", 300, 0, 1500, 1000, -1000, 1000);
    posRpcVsPosBar =
        R3B::root_owned<TH2F>("pos_Rpc_Vs_Pos_Bar", "Rpc_Position_Vs_Bar_Position", 300, 0, 1500, 500, -1000, 3000);

    stripIdVsNHitsCorr =
        R3B::root_owned<TH2F>("strip_Id_Vs_N_Hits_Corr", "Hits_Vs_Strip_ID(Good_Beam)", 21, -0.5, 20.5, 41, 0.5, 41.5);

    stripMultHisto = R3B::root_owned<TH1F>("strip_Mult_Histo", "Strips:H_Bar_AND_Start", 41, 0.5, 41.5);

    hStripEffHisto_H = R3B::root_owned<TH1F>("h_Strip_Eff_Histo_H", "Strip_Efficiency(H Bar)", 350, -1000, 2500);
    hBarEffHisto_H = R3B::root_owned<TH1F>("h_Bar_Eff_Histo_H", "Bar_Efficiency", 350, -1000, 2500);
    hEfficiencyHisto_H = R3B::root_owned<TH1F>("h_Efficiency_Histo_H", "Strip_Efficiency(H Bar)", 350, -1000, 2500);

    hStripEffHisto_V1 = R3B::root_owned<TH1F>("h_Strip_Eff_Histo_V1", "Strip_Efficiency(VR_Bar)", 350, -1000, 2500);
    hBarEffHisto_V1 = R3B::root_owned<TH1F>("h_Bar_Eff_Histo_V1", "Bar_Efficiency", 350, -1000, 2500);
    hEfficiencyHisto_V1 = R3B::root_owned<TH1F>("h_Efficiency_Histo_V1", "Strip_Efficiency(VR_Bar)", 350, -1000, 2500);

    hStripEffHisto_V2 = R3B::root_owned<TH1F>("h_Strip_Eff_Histo_V2", "Strip_Efficiency(VL_Bar)", 350, 0, 3500);
    hBarEffHisto_V2 = R3B::root_owned<TH1F>("h_Bar_Eff_Histo_V2", "Bar_Efficiency", 350, 0, 3500);
    hEfficiencyHisto_V2 = R3B::root_owned<TH1F>("h_Efficiency_Histo_V2", "Strip_Efficiency(VL_Bar)", 350, 0, 3500);

    for (Int_t i = 0; i < 41; i++)
    {

        sprintf(name, "Left:Coarse_Time_Strip_%i", i + 1);
        stripCoarseLeftHisto[i] = R3B::root_owned<TH1F>(name, name, 200, 0, 2200);

        sprintf(name, "Left:Fine_Time_Strip_%i", i + 1);
        stripFineLeftHisto[i] = R3B::root_owned<TH1F>(name, name, 200, 0, 600);

        sprintf(name, "Right:Coarse_Time_Strip_%i", i + 1);
        stripCoarseRightHisto[i] = R3B::root_owned<TH1F>(name, name, 200, 0, 2200);

        sprintf(name, "Right:Fine_Time_Strip_%i", i + 1);
        stripFineRightHisto[i] = R3B::root_owned<TH1F>(name, name, 200, 0, 600);

        sprintf(name, "ToF:Strip_%i", i + 1);
        stripTofHisto[i] = R3B::root_owned<TH1F>(name, name, fTofBins, fLeftTofLim, fRightTofLim);

        sprintf(name, "Time_Strip-Time_Bar_2:Strip_%i", i + 1);
        timeDiffStripPmtHisto[i] = R3B::root_owned<TH1F>(name, name, 200, -70, -20);
    }

    for (Int_t i = 0; i < 8; i++)
    {
        sprintf(name, "Coarse_Time_Pmt_%i_TOP", i + 1);
        pmtCoarseHistoTop[i] = R3B::root_owned<TH1F>(name, name, 1000, 0, 3000);

        sprintf(name, "Fine_Time_Pmt_%i_TOP", i + 1);
        pmtFineHistoTop[i] = R3B::root_owned<TH1F>(name, name, 1000, 0, 3000);

        sprintf(name, "Time_Pmt_%i_TOP", i + 1);
        pmtPreCalTimeHistoTop[i] = R3B::root_owned<TH1F>(name, name, 1000, -550, 400);

        sprintf(name, "ToT_Pmt_%i_TOP", i + 1);
        pmtPreCalTotHistoTop[i] = R3B::root_owned<TH1F>(name, name, 1000, -550, 400);
    }

    for (Int_t i = 0; i < 8; i++)
    {

        sprintf(name, "Coarse_Time_Pmt_%i_BOTTOM", i + 1);
        pmtCoarseHistoBottom[i] = R3B::root_owned<TH1F>(name, name, 1000, 0, 3000);

        sprintf(name, "Fine_Time_Pmt_%i_BOTTOM", i + 1);
        pmtFineHistoBottom[i] = R3B::root_owned<TH1F>(name, name, 1000, 0, 3000);

        sprintf(name, "Time_Pmt_%i_BOTTOM", i + 1);
        pmtPreCalTimeHistoBottom[i] = R3B::root_owned<TH1F>(name, name, 1000, -550, 400);

        sprintf(name, "ToT_Pmt_%i_BOTTOM", i + 1);
        pmtPreCalTotHistoBottom[i] = R3B::root_owned<TH1F>(name, name, 1000, -550, 400);
    }

    for (Int_t i = 0; i < 9; i++)
    {

        sprintf(name, "Coarse_Time_Ref_%i", i + 1);
        refCoarseHisto[i] = R3B::root_owned<TH1F>(name, name, 1000, 0, 3000);

        sprintf(name, "Fine_Time_Ref_%i", i + 1);
        refFineHisto[i] = R3B::root_owned<TH1F>(name, name, 1000, 0, 3000);
    }

    for (Int_t i = 0; i < 41; i++)
    {

        stripCoarseLeftHisto[i]->GetXaxis()->SetTitle("Coarse Time");
        stripCoarseLeftHisto[i]->GetYaxis()->SetTitle("Counts");
        stripCoarseLeftHisto[i]->GetXaxis()->SetTitleOffset(1.1);
        stripCoarseLeftHisto[i]->GetYaxis()->SetTitleOffset(1.1);
        stripCoarseLeftHisto[i]->GetXaxis()->CenterTitle(true);
        stripCoarseLeftHisto[i]->GetYaxis()->CenterTitle(true);

        leftStripCanvasCoarse->cd(i + 1);
        stripCoarseLeftHisto[i]->Draw();

        stripFineRightHisto[i]->GetXaxis()->SetTitle("Fine Time");
        stripFineRightHisto[i]->GetYaxis()->SetTitle("Counts");
        stripFineRightHisto[i]->GetXaxis()->SetTitleOffset(1.1);
        stripFineRightHisto[i]->GetYaxis()->SetTitleOffset(1.1);
        stripFineRightHisto[i]->GetXaxis()->CenterTitle(true);
        stripFineRightHisto[i]->GetYaxis()->CenterTitle(true);

        rightStripCanvasFine->cd(i + 1);
        stripFineRightHisto[i]->Draw();

        stripFineLeftHisto[i]->GetXaxis()->SetTitle("Fine Time");
        stripFineLeftHisto[i]->GetYaxis()->SetTitle("Counts");
        stripFineLeftHisto[i]->GetXaxis()->SetTitleOffset(1.1);
        stripFineLeftHisto[i]->GetYaxis()->SetTitleOffset(1.1);
        stripFineLeftHisto[i]->GetXaxis()->CenterTitle(true);
        stripFineLeftHisto[i]->GetYaxis()->CenterTitle(true);

        leftStripCanvasFine->cd(i + 1);
        stripFineLeftHisto[i]->Draw();

        stripCoarseRightHisto[i]->GetXaxis()->SetTitle("Coarse Time");
        stripCoarseRightHisto[i]->GetYaxis()->SetTitle("Counts");
        stripCoarseRightHisto[i]->GetXaxis()->SetTitleOffset(1.1);
        stripCoarseRightHisto[i]->GetYaxis()->SetTitleOffset(1.1);
        stripCoarseRightHisto[i]->GetXaxis()->CenterTitle(true);
        stripCoarseRightHisto[i]->GetYaxis()->CenterTitle(true);

        rightStripCanvasCoarse->cd(i + 1);
        stripCoarseRightHisto[i]->Draw();

        stripTofHisto[i]->GetXaxis()->SetTitle("ToF (ns)");
        stripTofHisto[i]->GetYaxis()->SetTitle("Counts");
        stripTofCanvas->cd(i + 1);
        stripTofHisto[i]->Draw();

        timeDiffStripPmtHisto[i]->GetXaxis()->SetTitle("Time Difference (ns)");
        timeDiffStripPmtHisto[i]->GetYaxis()->SetTitle("Counts");
        timeDiffStripPmtCanvas->cd(i + 1);
        timeDiffStripPmtHisto[i]->Draw();
    }

    stripLeftFolder->Add(leftStripCanvasCoarse);
    stripLeftFolder->Add(leftStripCanvasFine);

    stripRightFolder->Add(rightStripCanvasCoarse);
    stripRightFolder->Add(rightStripCanvasFine);

    stripVsTofFolder->Add(stripTofCanvas);

    stripLeftTotCorr->GetXaxis()->SetTitle("Strip Number");
    stripLeftTotCorr->GetYaxis()->SetTitle("ToT");
    stripLeftTotCorrCanvas->cd();
    stripLeftTotCorr->Draw("COLZ");
    stripLeftPreCalFolder->Add(stripLeftTotCorrCanvas);

    stripRightTotCorr->GetXaxis()->SetTitle("Strip Number");
    stripRightTotCorr->GetYaxis()->SetTitle("ToT");
    stripRightTotCorrCanvas->cd();
    stripRightTotCorr->Draw("COLZ");
    stripRightPreCalFolder->Add(stripRightTotCorrCanvas);

    stripLeftTimeCorr->GetXaxis()->SetTitle("Strip Number");
    stripLeftTimeCorr->GetYaxis()->SetTitle("Time");
    stripLeftTimeCorrCanvas->cd();
    stripLeftTimeCorr->Draw("COLZ");
    stripLeftPreCalFolder->Add(stripLeftTimeCorrCanvas);

    stripRightTimeCorr->GetXaxis()->SetTitle("Strip Number");
    stripRightTimeCorr->GetYaxis()->SetTitle("Time");
    stripRightTimeCorrCanvas->cd();
    stripRightTimeCorr->Draw("COLZ");
    stripRightPreCalFolder->Add(stripRightTimeCorrCanvas);

    stripCalTimeCorr->GetXaxis()->SetTitle("Left Time");
    stripCalTimeCorr->GetYaxis()->SetTitle("Right Time");
    stripCalTimeCorrCanvas->cd();
    stripCalTimeCorr->Draw("COLZ");
    stripCalFolder->Add(stripCalTimeCorrCanvas);

    stripCalToTCorr->GetXaxis()->SetTitle("Left ToT");
    stripCalToTCorr->GetYaxis()->SetTitle("Right ToT");
    stripCalTotCorrCanvas->cd();
    stripCalToTCorr->Draw("COLZ");
    stripCalFolder->Add(stripCalTotCorrCanvas);

    stripPosHitCorr->GetXaxis()->SetTitle("Position");
    stripPosHitCorr->GetYaxis()->SetTitle("Strip");
    hitMapCanvas->cd(1);
    stripPosHitCorr->Draw("COLZ");
    hitFolder->Add(hitMapCanvas);

    totalChargeHist->GetXaxis()->SetTitle("Charge");
    totalChargeHist->GetYaxis()->SetTitle("Counts");
    hitFolder->Add(totalChargeHist);

    meanChargeCorr->GetXaxis()->SetTitle("Position");
    meanChargeCorr->GetYaxis()->SetTitle("Strip");
    hitMapCanvas->cd(2);
    meanChargeCorr->Draw("COLZ");

    tofCorr->GetXaxis()->SetTitle("ToF (ps)");
    tofCorr->GetYaxis()->SetTitle("Strip Number");
    tofCorrCanvas->cd();
    tofCorr->Draw("COLZ");
    stripVsTofFolder->Add(tofCorrCanvas);

    tofVsPosCorr->GetXaxis()->SetTitle("Pos (mm)");
    tofVsPosCorr->GetYaxis()->SetTitle("ToF (ps)");
    tofPosCanvas->cd();
    tofVsPosCorr->Draw("COLZ");

    stripVsTofFolder->Add(tofPosCanvas);

    timeDiffStripPmtCorr->GetXaxis()->SetTitle("Time Difference (ns)");
    timeDiffStripPmtCorr->GetYaxis()->SetTitle("Strip Number");
    timeDiffStripPmtCorrCanvas->cd();
    timeDiffStripPmtCorr->Draw("COLZ");
    stripBarCorrFolder->Add(timeDiffStripPmtCorrCanvas);
    stripBarCorrFolder->Add(timeDiffStripPmtCanvas);

    deltaTVsPosCorr->GetXaxis()->SetTitle("Position (mm)");
    deltaTVsPosCorr->GetYaxis()->SetTitle("Time Difference (ns)");
    posRpcCorrCanvas->cd(1);
    deltaTVsPosCorr->Draw("COLZ");

    posRpcVsPosBar->GetXaxis()->SetTitle("RPC Position (mm)");
    posRpcVsPosBar->GetYaxis()->SetTitle("Sci Position (mm)");
    posRpcCorrCanvas->cd(2);
    posRpcVsPosBar->Draw("COLZ");

    stripIdVsNHitsCorr->GetXaxis()->SetTitle("Number of Hits");
    stripIdVsNHitsCorr->GetYaxis()->SetTitle("Strip Number");
    stripIdVsHitsCanvas->cd();
    stripIdVsNHitsCorr->Draw("COLZ");
    stripBarCorrFolder->Add(stripIdVsHitsCanvas);

    hEfficiencyHisto_H->GetXaxis()->SetTitle("Position (mm)");
    hEfficiencyHisto_H->GetYaxis()->SetTitle("Efficiency");
    efficiencyCanvas->cd(1);
    hEfficiencyHisto_H->Draw("");

    hEfficiencyHisto_V1->GetXaxis()->SetTitle("Position (mm)");
    hEfficiencyHisto_V1->GetYaxis()->SetTitle("Efficiency");
    efficiencyCanvas->cd(2);
    hEfficiencyHisto_V1->Draw("");

    hEfficiencyHisto_V2->GetXaxis()->SetTitle("Position (mm)");
    hEfficiencyHisto_V2->GetYaxis()->SetTitle("Efficiency");
    efficiencyCanvas->cd(3);
    hEfficiencyHisto_V2->Draw("");

    efficiencyFolder->Add(efficiencyCanvas);

    /* ---------- Efficiency Canvases ---------- */
    stripMultCanvas = new TCanvas("stripMultCanvas", "Strips:Start_AND_H_Bar");
    // stripMultCanvas->Divide(2,1);

    stripMultCanvas->cd();

    stripMultHisto->GetXaxis()->SetTitle("Strip ID");
    stripMultHisto->GetYaxis()->SetTitle("Counts");
    stripMultHisto->Draw();

    tex1 = new TLatex(2, 80000, "Elapsed time (since last restart): ");
    tex1->SetTextFont(42);
    tex1->SetTextSize(0.04);
    tex1->SetLineWidth(1);

    tex2 = new TLatex(2, 80000, "Hits in Horizontal Bar : ");
    tex2->SetTextFont(42);
    tex2->SetTextSize(0.04);
    tex2->SetLineWidth(1);

    // tex1->Draw();

    efficiencyFolder->Add(stripMultCanvas);

    stripBarCorrFolder->Add(posRpcCorrCanvas);

    for (Int_t i = 0; i < 8; i++)
    {

        pmtCoarseHistoTop[i]->GetXaxis()->SetTitle("Coarse Time");
        pmtCoarseHistoTop[i]->GetYaxis()->SetTitle("Counts");
        pmtCoarseCanvas->cd(i + 1);
        pmtCoarseHistoTop[i]->Draw();

        pmtFineHistoTop[i]->GetXaxis()->SetTitle("Fine Time");
        pmtFineHistoTop[i]->GetYaxis()->SetTitle("Counts");
        pmtFineCanvas->cd(i + 1);
        pmtFineHistoTop[i]->Draw();

        pmtPreCalTimeHistoTop[i]->GetXaxis()->SetTitle("Time");
        pmtPreCalTimeHistoTop[i]->GetYaxis()->SetTitle("Counts");
        pmtPreCalTimeCanvas->cd(i + 1);
        pmtPreCalTimeHistoTop[i]->Draw();

        pmtPreCalTotHistoTop[i]->GetXaxis()->SetTitle("ToT");
        pmtPreCalTotHistoTop[i]->GetYaxis()->SetTitle("Counts");
        pmtPreCalTotCanvas->cd(i + 1);
        pmtPreCalTotHistoTop[i]->Draw();

        pmtCoarseHistoBottom[i]->GetXaxis()->SetTitle("Coarse Time");
        pmtCoarseHistoBottom[i]->GetYaxis()->SetTitle("Counts");
        pmtCoarseCanvas->cd(i + 1 + 4);
        pmtCoarseHistoBottom[i]->Draw();

        pmtFineHistoBottom[i]->GetXaxis()->SetTitle("Fine Time");
        pmtFineHistoBottom[i]->GetYaxis()->SetTitle("Counts");
        pmtFineCanvas->cd(i + 1 + 4);
        pmtFineHistoBottom[i]->Draw();

        pmtPreCalTimeHistoBottom[i]->GetXaxis()->SetTitle("Time");
        pmtPreCalTimeHistoBottom[i]->GetYaxis()->SetTitle("Counts");
        pmtPreCalTimeCanvas->cd(i + 1 + 4);
        pmtPreCalTimeHistoBottom[i]->Draw();

        pmtPreCalTotHistoBottom[i]->GetXaxis()->SetTitle("ToT");
        pmtPreCalTotHistoBottom[i]->GetYaxis()->SetTitle("Counts");
        pmtPreCalTotCanvas->cd(i + 1 + 4);
        pmtPreCalTotHistoBottom[i]->Draw();
    }

    for (Int_t i = 0; i < 9; i++)
    {

        refFineHisto[i]->GetXaxis()->SetTitle("Fine Time");
        refFineHisto[i]->GetYaxis()->SetTitle("Counts");
        refFineCanvas->cd(i + 1);
        refFineHisto[i]->Draw();

        refCoarseHisto[i]->GetXaxis()->SetTitle("Coarse Time");
        refCoarseHisto[i]->GetYaxis()->SetTitle("Counts");
        refCoarseCanvas->cd(i + 1);
        refCoarseHisto[i]->Draw();
    }

    refFolder->Add(refFineCanvas);
    refFolder->Add(refCoarseCanvas);
    pmtFolder->Add(pmtCoarseCanvas);
    pmtFolder->Add(pmtFineCanvas);

    pmtPreCalFolder->Add(pmtPreCalTimeCanvas);
    pmtPreCalFolder->Add(pmtPreCalTotCanvas);

    stripCoarseRightCorr->GetXaxis()->SetTitle("Strip Number");
    stripCoarseRightCorr->GetYaxis()->SetTitle("Coarse Time");
    stripCoarseRightCorrCanvas->cd();
    stripCoarseRightCorr->Draw("COLZ");
    stripRightFolder->Add(stripCoarseRightCorrCanvas);

    stripCoarseLeftCorr->GetXaxis()->SetTitle("Strip Number");
    stripCoarseLeftCorr->GetYaxis()->SetTitle("Coarse Time");
    stripCoarseLeftCorrCanvas->cd();
    stripCoarseLeftCorr->Draw("COLZ");
    stripLeftFolder->Add(stripCoarseLeftCorrCanvas);

    stripFineRightCorr->GetXaxis()->SetTitle("Strip Number");
    stripFineRightCorr->GetYaxis()->SetTitle("Fine Time");
    stripFineRightCorrCanvas->cd();
    stripFineRightCorr->Draw("COLZ");
    stripRightFolder->Add(stripFineRightCorrCanvas);

    stripFineLeftCorr->GetXaxis()->SetTitle("Strip Number");
    stripFineLeftCorr->GetYaxis()->SetTitle("Fine Time");
    stripFineLeftCorrCanvas->cd();
    stripFineLeftCorr->Draw("COLZ");
    stripLeftFolder->Add(stripFineLeftCorrCanvas);

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);
    run->AddObject(mainfol);

    run->GetHttpServer()->RegisterCommand("Reset_RPC", Form("/Objects/%s/->Reset_RPC_Histo()", GetName()));
    run->GetHttpServer()->RegisterCommand("Reset_Efficiencies", Form("/Objects/%s/->Reset_Efficiencies()", GetName()));
    run->GetHttpServer()->RegisterCommand("GO_ONSPILL", Form("/Objects/%s/->GO_ONSPILL()", GetName()));
    run->GetHttpServer()->RegisterCommand("GO_OFFSPILL", Form("/Objects/%s/->GO_OFFSPILL()", GetName()));

    return kSUCCESS;
}

void R3BRpcOnlineSpectra::Reset_RPC_Histo()
{

    for (Int_t i = 0; i < 41; i++)
    {
        stripCoarseLeftHisto[i]->Reset();
        stripFineLeftHisto[i]->Reset();
        stripCoarseRightHisto[i]->Reset();
        stripCoarseLeftHisto[i]->Reset();
        stripTofHisto[i]->Reset();
        timeDiffStripPmtHisto[i]->Reset();
        stripTofHisto[i]->Reset();
    }

    for (Int_t i = 0; i < 9; i++)
    {
        refFineHisto[i]->Reset();
        refCoarseHisto[i]->Reset();
    }

    for (Int_t i = 0; i < 8; i++)
    {

        pmtFineHistoTop[i]->Reset();
        pmtCoarseHistoTop[i]->Reset();
        pmtFineHistoBottom[i]->Reset();
        pmtCoarseHistoBottom[i]->Reset();
    }

    stripCoarseLeftCorr->Reset();
    stripCoarseRightCorr->Reset();
    stripFineLeftCorr->Reset();
    stripFineRightCorr->Reset();
    stripCoarseLeftCorr->Reset();
    stripCoarseRightCorr->Reset();
    stripFineLeftCorr->Reset();
    stripFineRightCorr->Reset();
    stripCalTimeCorr->Reset();
    stripCalToTCorr->Reset();
    totalChargeHist->Reset();
    stripPosHitCorr->Reset();
    meanChargeCorr->Reset();
    tofCorr->Reset();
    timeDiffStripPmtCorr->Reset();
    deltaTVsPosCorr->Reset();
    posRpcVsPosBar->Reset();
    tofCorr->Reset();
    tofVsPosCorr->Reset();
    stripIdVsNHitsCorr->Reset();

    for (Int_t i = 0; i < 41; i++)
    {
        for (Int_t j = 0; j < 50; j++)
        {
            meanCharges[i][j] = 0.0;
        }
    }

    return;
}

void R3BRpcOnlineSpectra::Reset_Efficiencies()
{

    stripMultHisto->Reset();
    fBarHits = 0.0;
    fTimeStart = 0;

    hEfficiencyHisto_H->Reset();
    hStripEffHisto_H->Reset();
    hBarEffHisto_H->Reset();

    hEfficiencyHisto_V1->Reset();
    hStripEffHisto_V1->Reset();
    hBarEffHisto_V1->Reset();

    hEfficiencyHisto_V2->Reset();
    hStripEffHisto_V2->Reset();
    hBarEffHisto_V2->Reset();

    return;
}

void R3BRpcOnlineSpectra::GO_ONSPILL()
{
    R3BLOG(info, "");
    fSpill = 1;
    return;
}

void R3BRpcOnlineSpectra::GO_OFFSPILL()
{
    R3BLOG(info, "");
    fSpill = 0;
    return;
}

void R3BRpcOnlineSpectra::Exec(Option_t* option)
{
    fcounter++;
    Float_t elapsedTime;
    if (fTimeStart == 0)
        fTimeStart = fEventHeader->GetTimeStamp();

    if (fEventHeader->GetTrigger() != fTrigger && fTrigger > -1)
    {
        return;
    }

    if (fEventHeader->GetTpat() > 0)
    {
        for (Int_t i = 0; i < 16; i++)
        {
            tpatbin = (fEventHeader->GetTpat() & (1 << i));
            if (tpatbin != 0)
            {
                fTPat = i + 1;
            }
        }
    }

    auto nMappedHits = fMappedDataItems->GetEntriesFast();

    if (fSpill == 1)
        execBool = (fTPat >= fFirstTPat && fTPat <= fLastTPat);

    if (fSpill == 0)
        execBool = !(fTPat >= fFirstTPat && fTPat <= fLastTPat);

    if (execBool)
    {

        /* ------------------- Map EventLoop ------------------*/
        for (Int_t ihit = 0; ihit < nMappedHits; ihit++)
        {

            R3BRpcMappedData* hit = dynamic_cast<R3BRpcMappedData*>(fMappedDataItems->At(ihit));

            if (!hit)
                continue;

            UInt_t idetector = hit->GetDetId();

            if (idetector == 0)
            {

                Int_t side = hit->GetSide();

                if (side == 0)
                {

                    stripFineLeftHisto[hit->GetChannelId() - 1]->Fill(hit->GetFineTime());
                    stripCoarseLeftHisto[hit->GetChannelId() - 1]->Fill(hit->GetCoarseTime());
                    stripCoarseLeftCorr->Fill(hit->GetChannelId(), hit->GetCoarseTime());
                    stripFineLeftCorr->Fill(hit->GetChannelId(), hit->GetFineTime());
                }

                if (side == 1)
                {

                    stripFineRightHisto[hit->GetChannelId() - 1]->Fill(hit->GetFineTime());
                    stripCoarseRightHisto[hit->GetChannelId() - 1]->Fill(hit->GetCoarseTime());
                    stripCoarseRightCorr->Fill(hit->GetChannelId(), hit->GetCoarseTime());
                    stripFineRightCorr->Fill(hit->GetChannelId(), hit->GetFineTime());
                }
            }

            if (idetector == 2)
            {

                refFineHisto[hit->GetChannelId() - 1]->Fill(hit->GetFineTime());
                refCoarseHisto[hit->GetChannelId() - 1]->Fill(hit->GetCoarseTime());
            }

            if (idetector == 1)
            {

                Int_t pmtSide = hit->GetSide();

                if (pmtSide == 0)
                {

                    pmtFineHistoTop[hit->GetChannelId() - 1]->Fill(hit->GetFineTime());
                    pmtCoarseHistoTop[hit->GetChannelId() - 1]->Fill(hit->GetCoarseTime());
                }

                if (pmtSide == 1)
                {

                    pmtFineHistoBottom[hit->GetChannelId() - 1]->Fill(hit->GetFineTime());
                    pmtCoarseHistoBottom[hit->GetChannelId() - 1]->Fill(hit->GetCoarseTime());
                }
            }
        }

        /* ----------------- PreCal EventLoop ------------------*/
        auto nPreCalHits = fPreCalDataItems->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nPreCalHits; ihit++)
        {

            R3BRpcPreCalData* hit = dynamic_cast<R3BRpcPreCalData*>(fPreCalDataItems->At(ihit));

            Int_t side = hit->GetSide();

            if (hit->GetDetId() == 0)
            {

                if (side == 0)
                {
                    stripLeftTotCorr->Fill(hit->GetChannelId(), hit->GetTot());
                    stripLeftTimeCorr->Fill(hit->GetChannelId(), hit->GetTime());
                }
                if (side == 1)
                {

                    stripRightTotCorr->Fill(hit->GetChannelId(), hit->GetTot());
                    stripRightTimeCorr->Fill(hit->GetChannelId(), hit->GetTime());
                }
            }

            if (hit->GetDetId() == 1)
            {

                if (side == 0)
                {
                    pmtPreCalTimeHistoTop[hit->GetChannelId() - 1]->Fill(hit->GetTime());
                    pmtPreCalTotHistoTop[hit->GetChannelId() - 1]->Fill(hit->GetTot());
                }

                if (side == 1)
                {
                    pmtPreCalTimeHistoBottom[hit->GetChannelId() - 1]->Fill(hit->GetTime());
                    pmtPreCalTotHistoBottom[hit->GetChannelId() - 1]->Fill(hit->GetTot());
                }
            }
        }

        /* -------------------- Cal EventLoop ---------------------*/
        auto nStripCalHits = fCalDataItems->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nStripCalHits; ihit++)
        {

            R3BRpcCalData* hit = dynamic_cast<R3BRpcCalData*>(fCalDataItems->At(ihit));

            if (hit->GetDetId() == 0)
            {

                stripCalTimeCorr->Fill(hit->GetTimeL_T(), hit->GetTimeR_B());
                stripCalToTCorr->Fill(hit->GetTotL_T(), hit->GetTotR_B());
            }
        }

        /*---------------- HIT EventLoop ---------------*/
        auto nStripHits = fHitDataItems->GetEntriesFast();

        Int_t channelId;
        Int_t detId;
        Float_t pos, charge;
        Int_t bin;

        Float_t barTime = 0.0;
        Float_t barID = 0.0;
        Float_t sciPos = 0.0;
        Int_t StripCounter = 0;

        Bool_t goodBar = kFALSE;

        for (Int_t ihit = 0; ihit < nStripHits; ihit++)
        {

            R3BRpcHitData* hit = dynamic_cast<R3BRpcHitData*>(fHitDataItems->At(ihit));

            detId = hit->GetDetId();
            channelId = hit->GetChannelId();
            pos = hit->GetPos();
            charge = hit->GetCharge();

            if (detId == 1 && channelId == 2)
            {

                goodBar = kTRUE;
                barTime = hit->GetTime();
                sciPos = hit->GetPos();
                fBarHits++;
                hBarEffHisto_H->Fill(pos);
            }

            if (detId == 1 && channelId == 1)
            {
                hBarEffHisto_V1->Fill(pos);
            }

            if (detId == 1 && channelId == 3)
            {
                hBarEffHisto_V2->Fill(pos);
            }
        }

        Int_t stripMultiplicityArray[41] = { 0 };

        for (Int_t ihit = 0; ihit < nStripHits; ihit++)
        {

            R3BRpcHitData* hit = dynamic_cast<R3BRpcHitData*>(fHitDataItems->At(ihit));
            detId = hit->GetDetId();
            channelId = hit->GetChannelId();
            pos = hit->GetPos();
            charge = hit->GetCharge();

            if (detId == 1)
                continue;

            StripCounter++;

            stripMultiplicityArray[channelId - 1] += 1;

            Double_t rpc_tof = hit->GetTof();

            stripTofHisto[channelId - 1]->Fill(rpc_tof);

            tofCorr->Fill(rpc_tof, channelId);

            tofVsPosCorr->Fill(pos, rpc_tof);

            if (goodBar)
            {

                timeDiffStripPmtHisto[channelId - 1]->Fill(hit->GetTime() - barTime);
                timeDiffStripPmtCorr->Fill(hit->GetTime() - barTime, channelId);

                posRpcVsPosBar->Fill(pos, sciPos);
                deltaTVsPosCorr->Fill(pos, hit->GetTime() - barTime);
                stripMultHisto->Fill(channelId);
            }

            stripPosHitCorr->Fill(pos, channelId);
            totalChargeHist->Fill(hit->GetCharge());
            counts++;

            bin = Int_t(50 - (1500 - pos) / 30.0 + 1);

            if (charge >= 0.0 && bin <= 50 && bin > 0)
            {

                meanCharges[channelId - 1][bin - 1] = meanCharges[channelId - 1][bin - 1] +
                                                      (1.0 / counts) * (charge - meanCharges[channelId - 1][bin - 1]);

                meanChargeCorr->SetBinContent(bin, channelId, meanCharges[channelId - 1][bin - 1]);
            }
        }

        for (Int_t ihit = 0; ihit < nStripHits; ihit++)
        {

            R3BRpcHitData* hit = dynamic_cast<R3BRpcHitData*>(fHitDataItems->At(ihit));
            detId = hit->GetDetId();
            channelId = hit->GetChannelId();
            pos = hit->GetPos();
            charge = hit->GetCharge();

            if (detId == 1 && channelId == 2 && StripCounter != 0)
                hStripEffHisto_H->Fill(pos);

            if (detId == 1 && channelId == 1 && StripCounter != 0)
                hStripEffHisto_V1->Fill(pos);

            if (detId == 1 && channelId == 3 && StripCounter != 0)
                hStripEffHisto_V2->Fill(pos);
        }

        for (int i = 0; i < 41; i++)
            stripIdVsNHitsCorr->Fill(stripMultiplicityArray[i], i + 1);

        if (fNEvents % 1000 == 0)
        {

            efficiencyCanvas->cd(1);
            hEfficiencyHisto_H->Divide(hStripEffHisto_H, hBarEffHisto_H, 1, 1, "");
            efficiencyCanvas->Modified();

            efficiencyCanvas->cd(2);
            hEfficiencyHisto_V1->Divide(hStripEffHisto_V1, hBarEffHisto_V1, 1, 1, "");
            efficiencyCanvas->Modified();

            efficiencyCanvas->cd(3);
            hEfficiencyHisto_V2->Divide(hStripEffHisto_V2, hBarEffHisto_V2, 1, 1, "");
            efficiencyCanvas->Modified();

            fTimeEnd = fEventHeader->GetTimeStamp();
            elapsedTime = (fTimeEnd - fTimeStart) / 1E9;

            Int_t y = (int)stripMultHisto->GetMaximum();

            char name1[255];

            sprintf(name1, "Elapsed time (since last restart) : %d  s", (int)elapsedTime);

            char name2[255];

            sprintf(name2, "Hits in Horizontal Bar  : %d  ", (int)fBarHits);

            stripMultCanvas->Clear();
            stripMultHisto->Draw();
            tex1->DrawLatex(2, y, name1);
            tex1->Draw();
            tex2->DrawLatex(2, 0.95 * y, name2);
            tex2->Draw();

            stripMultCanvas->Modified();
            stripMultCanvas->cd();
        }
    }

    fNEvents += 1;
    return;
}

void R3BRpcOnlineSpectra::FinishEvent()
{

    if (fMappedDataItems)
    {
        fMappedDataItems->Clear();
    }

    if (fPreCalDataItems)
    {
        fPreCalDataItems->Clear();
    }

    if (fCalDataItems)
    {
        fCalDataItems->Clear();
    }

    if (fHitDataItems)
    {
        fHitDataItems->Clear();
    }
}

void R3BRpcOnlineSpectra::FinishTask()
{
    if (fMappedDataItems)
    {
        for (Int_t i = 0; i < 7; i++)
        {
            refFineHisto[i]->Write();
            refCoarseHisto[i]->Write();
        }

        for (Int_t i = 0; i < 4; i++)
        {
            pmtFineHistoTop[i]->Write();
            pmtCoarseHistoTop[i]->Write();
            pmtFineHistoBottom[i]->Write();
            pmtCoarseHistoBottom[i]->Write();
        }

        for (Int_t i = 0; i < 41; i++)
        {

            stripFineLeftHisto[i]->Write();
            stripFineRightHisto[i]->Write();

            stripCoarseLeftHisto[i]->Write();
            stripCoarseRightHisto[i]->Write();
            stripTofHisto[i]->Write();
        }

        leftStripCanvasCoarse->Write();
        rightStripCanvasCoarse->Write();

        leftStripCanvasFine->Write();
        rightStripCanvasFine->Write();

        stripCoarseLeftCorrCanvas->Write();
        stripCoarseRightCorrCanvas->Write();

        stripFineLeftCorrCanvas->Write();
        stripFineRightCorrCanvas->Write();

        pmtCoarseCanvas->Write();
        pmtFineCanvas->Write();
        refFineCanvas->Write();
        refCoarseCanvas->Write();

        stripLeftTotCorrCanvas->Write();
        stripRightTotCorrCanvas->Write();
        stripLeftTimeCorrCanvas->Write();
        stripRightTimeCorrCanvas->Write();
    }
}

ClassImp(R3BRpcOnlineSpectra);
