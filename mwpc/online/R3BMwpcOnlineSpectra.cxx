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
// -----                  R3BMwpcOnlineSpectra            -----
// -----    Created 29/09/19 by J.L. Rodriguez-Sanchez    -----
// -----           Fill MWPC online histograms            -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with MWPC online data
 */

#include "R3BMwpcOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BMwpcCalData.h"
#include "R3BMwpcHitData.h"
#include "R3BMwpcMappedData.h"
#include "R3BShared.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFolder.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdlib.h>

using namespace std;

R3BMwpcOnlineSpectra::R3BMwpcOnlineSpectra()
    : R3BMwpcOnlineSpectra("R3BMwpcOnlineSpectra", 1, "Mwpc0")
{
}

R3BMwpcOnlineSpectra::R3BMwpcOnlineSpectra(const TString& name, Int_t iVerbose, const TString& namedet)
    : FairTask(name, iVerbose)
    , fNameDet(namedet)
{
}

InitStatus R3BMwpcOnlineSpectra::Init()
{
    R3BLOG(info, "for " << fNameDet);

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3B" + fNameDet + "OnlineSpectra::Init FairRootManager not found";

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // get access to mapped data of mwpcs
    fMapItemsMwpc = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet + "MappedData"));
    R3BLOG_IF(fatal, fMapItemsMwpc == nullptr, fNameDet + "MappedData not found");

    // get access to cal data of mwpcs
    fCalItemsMwpc = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet + "CalData"));
    R3BLOG_IF(fatal, fCalItemsMwpc == nullptr, fNameDet + "CalData not found");

    // get access to hit data of mwpcs
    fHitItemsMwpc = dynamic_cast<TClonesArray*>(mgr->GetObject(fNameDet + "HitData"));
    R3BLOG_IF(warn, fHitItemsMwpc == nullptr, fNameDet + "HitData not found");

    // Create histograms for detectors
    TString Name1;
    TString Name2;

    // MWPC: Mapped data
    cMwpc_nPads = new TCanvas(fNameDet + "_nPads", fNameDet + "_nPads", 10, 10, 800, 700);
    cMwpcMap = new TCanvas(fNameDet + "_Q", fNameDet + "_Q", 10, 10, 800, 700);

    if (fNameDet == "Mwpc0")
    {
        cMwpc_nPads->Divide(2, 1);
        cMwpcMap->Divide(2, 1);

        Name1 = fNameDet + "_Plane1_nPadsPerEvent";
        fh1_mwpc_map_nPads[0] = R3B::root_owned<TH1I>(Name1, Name1, 67, -1.5, 65.5);
        fh1_mwpc_map_nPads[0]->SetFillColor(kBlue + 3);
        cMwpc_nPads->cd(1);
        fh1_mwpc_map_nPads[0]->Draw();

        Name1 = fNameDet + "_Plane3_nPadsPerEvent";
        fh1_mwpc_map_nPads[2] = R3B::root_owned<TH1I>(Name1, Name1, 67, -1.5, 65.5);
        fh1_mwpc_map_nPads[2]->SetFillColor(kRed + 3);
        cMwpc_nPads->cd(2);
        fh1_mwpc_map_nPads[2]->Draw();

        Name1 = fNameDet + "_Plane1Q";
        fh2_mwpc_map_q[0] = R3B::root_owned<TH2I>(Name1, Name1, 67, -1.5, 65.5, 1000, 0, 5000);
        fh2_mwpc_map_q[0]->GetXaxis()->SetTitle("X Pad number");
        fh2_mwpc_map_q[0]->GetYaxis()->SetTitle("Charge [channels]");
        cMwpcMap->cd(1);
        fh2_mwpc_map_q[0]->Draw("col");

        Name1 = fNameDet + "_Plane3Q";
        fh2_mwpc_map_q[2] = R3B::root_owned<TH2I>(Name1, Name1, 67, -1.5, 65.5, 1000, 0, 5000);
        fh2_mwpc_map_q[2]->SetFillColor(kRed + 3);
        fh2_mwpc_map_q[2]->GetXaxis()->SetTitle("Y Pad number");
        fh2_mwpc_map_q[2]->GetYaxis()->SetTitle("Charge [channels]");
        cMwpcMap->cd(2);
        fh2_mwpc_map_q[2]->Draw("col");
    }
    else if (fNameDet == "Mwpc1" || fNameDet == "Mwpc2")
    {
        cMwpc_nPads->Divide(3, 1);
        cMwpcMap->Divide(3, 1);

        Name1 = fNameDet + "_Plane1_nPadsPerEvent";
        fh1_mwpc_map_nPads[0] = R3B::root_owned<TH1I>(Name1, Name1, 67, -1.5, 65.5);
        fh1_mwpc_map_nPads[0]->SetFillColor(kBlue + 3);
        cMwpc_nPads->cd(1);
        fh1_mwpc_map_nPads[0]->Draw();

        Name1 = fNameDet + "_Plane2_nPadsPerEvent";
        fh1_mwpc_map_nPads[1] = R3B::root_owned<TH1I>(Name1, Name1, 67, -1.5, 65.5);
        fh1_mwpc_map_nPads[1]->SetFillColor(kBlue - 3);
        cMwpc_nPads->cd(2);
        fh1_mwpc_map_nPads[1]->Draw();

        Name1 = fNameDet + "_Plane3_nPadsPerEvent";
        fh1_mwpc_map_nPads[2] = R3B::root_owned<TH1I>(Name1, Name1, 43, -1.5, 41.5);
        fh1_mwpc_map_nPads[2]->SetFillColor(kRed + 3);
        cMwpc_nPads->cd(3);
        fh1_mwpc_map_nPads[2]->Draw();

        Name1 = fNameDet + "_Plane1Q";
        fh2_mwpc_map_q[0] = R3B::root_owned<TH2I>(Name1, Name1, 67, -1.5, 65.5, 1000, 0, 5000);
        fh2_mwpc_map_q[0]->GetXaxis()->SetTitle("X Pad number");
        fh2_mwpc_map_q[0]->GetYaxis()->SetTitle("Charge [channels]");
        cMwpcMap->cd(1);
        fh2_mwpc_map_q[0]->Draw("col");

        Name1 = fNameDet + "_Plane2Q";
        fh2_mwpc_map_q[1] = R3B::root_owned<TH2I>(Name1, Name1, 67, -1.5, 65.5, 1000, 0, 5000);
        fh2_mwpc_map_q[1]->GetXaxis()->SetTitle("X Pad number");
        fh2_mwpc_map_q[1]->GetYaxis()->SetTitle("Charge [channels]");
        cMwpcMap->cd(2);
        fh2_mwpc_map_q[1]->Draw("col");

        Name1 = fNameDet + "_Plane3Q";
        fh2_mwpc_map_q[2] = R3B::root_owned<TH2I>(Name1, Name1, 43, -1.5, 41.5, 1000, 0, 5000);
        fh2_mwpc_map_q[2]->GetXaxis()->SetTitle("Y Pad number");
        fh2_mwpc_map_q[2]->GetYaxis()->SetTitle("Charge [channels]");
        cMwpcMap->cd(3);
        fh2_mwpc_map_q[2]->Draw("col");
    }
    else
    {
        cMwpc_nPads->Divide(2, 1);
        cMwpcMap->Divide(2, 1);

        Name1 = fNameDet + "_Plane1_nPadsPerEvent";
        fh1_mwpc_map_nPads[0] = R3B::root_owned<TH1I>(Name1, Name1, 243, -1.5, 241.5);
        fh1_mwpc_map_nPads[0]->SetFillColor(kBlue + 3);
        cMwpc_nPads->cd(1);
        fh1_mwpc_map_nPads[0]->Draw();

        Name1 = fNameDet + "_Plane3_nPadsPerEvent";
        fh1_mwpc_map_nPads[2] = R3B::root_owned<TH1I>(Name1, Name1, 403, -1.5, 401.5);
        fh1_mwpc_map_nPads[2]->SetFillColor(kRed + 3);
        cMwpc_nPads->cd(2);
        fh1_mwpc_map_nPads[2]->Draw();

        Name1 = fNameDet + "_Plane1Q";
        fh2_mwpc_map_q[0] = R3B::root_owned<TH2I>(Name1, Name1, 402, -1.5, 400.5, 1000, 0, 5000);
        fh2_mwpc_map_q[0]->GetXaxis()->SetTitle("X Pad number");
        fh2_mwpc_map_q[0]->GetYaxis()->SetTitle("Charge [channels]");
        cMwpcMap->cd(1);
        fh2_mwpc_map_q[0]->Draw("col");

        Name1 = fNameDet + "_Plane3Q";
        fh2_mwpc_map_q[2] = R3B::root_owned<TH2I>(Name2, Name1, 152, -1.5, 150.5, 1000, 0, 5000);
        fh2_mwpc_map_q[2]->GetXaxis()->SetTitle("Y Pad number");
        fh2_mwpc_map_q[2]->GetYaxis()->SetTitle("Charge [channels]");
        cMwpcMap->cd(2);
        fh2_mwpc_map_q[2]->Draw("col");
    }

    // MWPC: Cal data
    cMWPCCal = new TCanvas(fNameDet + "_cal", fNameDet + " cal info", 10, 10, 800, 700);
    cMWPCCal->Divide(2, 1);
    Name1 = "fh1_" + fNameDet + "_posx_cal";
    Name2 = fNameDet + ": Position X";
    if (fNameDet != "Mwpc3")
        fh1_mwpc_cal[0] = R3B::root_owned<TH1F>(Name1, Name2, 64, 0.5, 64.5);
    else
        fh1_mwpc_cal[0] = R3B::root_owned<TH1F>(Name1, Name2, 400, 0.5, 400.5);
    fh1_mwpc_cal[0]->GetXaxis()->SetTitle("Position X [pads]");
    fh1_mwpc_cal[0]->GetYaxis()->SetTitle("Counts");
    fh1_mwpc_cal[0]->GetYaxis()->SetTitleOffset(1.15);
    fh1_mwpc_cal[0]->GetXaxis()->CenterTitle(true);
    fh1_mwpc_cal[0]->GetYaxis()->CenterTitle(true);
    fh1_mwpc_cal[0]->GetXaxis()->SetLabelSize(0.045);
    fh1_mwpc_cal[0]->GetXaxis()->SetTitleSize(0.045);
    fh1_mwpc_cal[0]->GetYaxis()->SetLabelSize(0.045);
    fh1_mwpc_cal[0]->GetYaxis()->SetTitleSize(0.045);
    fh1_mwpc_cal[0]->SetFillColor(8);
    fh1_mwpc_cal[0]->SetLineColor(1);
    cMWPCCal->cd(1);
    fh1_mwpc_cal[0]->Draw("");

    Name1 = "fh1_" + fNameDet + "_posy_cal";
    Name2 = fNameDet + ": Position Y";

    if (fNameDet == "Mwpc0")
        fh1_mwpc_cal[1] = R3B::root_owned<TH1F>(Name1, Name2, 64, 0.5, 64.5);
    else if (fNameDet == "Mwpc1" || fNameDet == "Mwpc2")
        fh1_mwpc_cal[1] = R3B::root_owned<TH1F>(Name1, Name2, 40, 0.5, 40.5);
    else
        fh1_mwpc_cal[1] = R3B::root_owned<TH1F>(Name1, Name2, 120, 0.5, 120.5);
    fh1_mwpc_cal[1]->GetXaxis()->SetTitle("Position Y [pads]");
    fh1_mwpc_cal[1]->GetYaxis()->SetTitle("Counts");
    fh1_mwpc_cal[1]->GetYaxis()->SetTitleOffset(1.15);
    fh1_mwpc_cal[1]->GetXaxis()->CenterTitle(true);
    fh1_mwpc_cal[1]->GetYaxis()->CenterTitle(true);
    fh1_mwpc_cal[1]->GetXaxis()->SetLabelSize(0.045);
    fh1_mwpc_cal[1]->GetXaxis()->SetTitleSize(0.045);
    fh1_mwpc_cal[1]->GetYaxis()->SetLabelSize(0.045);
    fh1_mwpc_cal[1]->GetYaxis()->SetTitleSize(0.045);
    fh1_mwpc_cal[1]->SetFillColor(8);
    fh1_mwpc_cal[1]->SetLineColor(1);
    cMWPCCal->cd(2);
    fh1_mwpc_cal[1]->Draw("");

    cMWPCCal2D = new TCanvas(fNameDet + "_cal2D", fNameDet + " cal 2D info", 10, 10, 800, 700);

    Name1 = "fh2_" + fNameDet + "_posXy_cal";
    Name2 = fNameDet + ": Position X vs Y";
    if (fNameDet == "Mwpc0")
        fh2_mwpc_cal = R3B::root_owned<TH2F>(Name1, Name2, 64, 0.5, 64.5, 64, 0.5, 64.5);
    else if (fNameDet == "Mwpc1" || fNameDet == "Mwpc2")
        fh2_mwpc_cal = R3B::root_owned<TH2F>(Name1, Name2, 64, 0.5, 64.5, 40, 0.5, 40.5);
    else
        fh2_mwpc_cal = R3B::root_owned<TH2F>(Name1, Name2, 288, 0.5, 288.5, 120, 0.5, 120.5);
    fh2_mwpc_cal->GetXaxis()->SetTitle("Position X [pads]");
    fh2_mwpc_cal->GetYaxis()->SetTitle("Position Y [pads]");
    fh2_mwpc_cal->GetYaxis()->SetTitleOffset(1.1);
    fh2_mwpc_cal->GetXaxis()->CenterTitle(true);
    fh2_mwpc_cal->GetYaxis()->CenterTitle(true);
    fh2_mwpc_cal->GetXaxis()->SetLabelSize(0.045);
    fh2_mwpc_cal->GetXaxis()->SetTitleSize(0.045);
    fh2_mwpc_cal->GetYaxis()->SetLabelSize(0.045);
    fh2_mwpc_cal->GetYaxis()->SetTitleSize(0.045);
    fh2_mwpc_cal->Draw("col");

    TCanvas* cx = new TCanvas(fNameDet + "_XQ", "", 10, 10, 800, 700);
    Name1 = "fh2_" + fNameDet + "_XQ";
    Name2 = fNameDet + ": Pad vs XQ";

    if (fNameDet == "Mwpc0" || fNameDet == "Mwpc1" || fNameDet == "Mwpc2")
        fh2_mwpc_xq = R3B::root_owned<TH2F>(Name1, Name2, 64, 0.5, 64.5, 4200, 0, 4200);
    else
        fh2_mwpc_xq = R3B::root_owned<TH2F>(Name1, Name2, 400, 0.5, 400.5, 4200, 0, 4200);

    fh2_mwpc_xq->GetXaxis()->SetTitle("Pad number");
    fh2_mwpc_xq->GetYaxis()->SetTitle("Charge [channels]");
    fh2_mwpc_xq->GetYaxis()->SetTitleOffset(1.1);
    fh2_mwpc_xq->GetXaxis()->CenterTitle(true);
    fh2_mwpc_xq->GetYaxis()->CenterTitle(true);
    fh2_mwpc_xq->GetXaxis()->SetLabelSize(0.045);
    fh2_mwpc_xq->GetXaxis()->SetTitleSize(0.045);
    fh2_mwpc_xq->GetYaxis()->SetLabelSize(0.045);
    fh2_mwpc_xq->GetYaxis()->SetTitleSize(0.045);
    fh2_mwpc_xq->Draw("col");

    TCanvas* cy = new TCanvas(fNameDet + "_YQ", "", 10, 10, 800, 700);
    Name1 = "fh2_" + fNameDet + "_YQ";
    Name2 = fNameDet + ": Pad vs YQ";
    if (fNameDet == "Mwpc0")
        fh2_mwpc_yq = R3B::root_owned<TH2F>(Name1, Name2, 64, 0.5, 64.5, 4200, 0, 4200);
    else if (fNameDet == "Mwpc1" || fNameDet == "Mwpc2")
        fh2_mwpc_yq = R3B::root_owned<TH2F>(Name1, Name2, 40, 0.5, 40.5, 4200, 0, 4200);
    else
        fh2_mwpc_yq = R3B::root_owned<TH2F>(Name1, Name2, 120, 0.5, 120.5, 4200, 0, 4200);

    fh2_mwpc_yq->GetXaxis()->SetTitle("Pad number");
    fh2_mwpc_yq->GetYaxis()->SetTitle("Charge [channels]");
    fh2_mwpc_yq->GetYaxis()->SetTitleOffset(1.1);
    fh2_mwpc_yq->GetXaxis()->CenterTitle(true);
    fh2_mwpc_yq->GetYaxis()->CenterTitle(true);
    fh2_mwpc_yq->GetXaxis()->SetLabelSize(0.045);
    fh2_mwpc_yq->GetXaxis()->SetTitleSize(0.045);
    fh2_mwpc_yq->GetYaxis()->SetLabelSize(0.045);
    fh2_mwpc_yq->GetYaxis()->SetTitleSize(0.045);
    fh2_mwpc_yq->Draw("col");

    // Hit data
    chitx = new TCanvas(fNameDet + "_Xpos", "", 10, 10, 800, 700);
    Name1 = "fh1_" + fNameDet + "_Xpos";
    Name2 = fNameDet + ": X (mm)";
    if (fNameDet == "Mwpc3")
        fh1_Xpos = R3B::root_owned<TH1F>(Name1, Name2, 1800, -450, 450);
    else
        fh1_Xpos = R3B::root_owned<TH1F>(Name1, Name2, 400, -100, 100);
    fh1_Xpos->GetXaxis()->SetTitle("(Wixhausen)<---  X (mm)  ---> (Messel)");
    fh1_Xpos->GetYaxis()->SetTitle("counts per bin");
    fh1_Xpos->GetYaxis()->SetTitleOffset(1.1);
    fh1_Xpos->GetXaxis()->CenterTitle(true);
    fh1_Xpos->GetYaxis()->CenterTitle(true);
    fh1_Xpos->GetXaxis()->SetLabelSize(0.045);
    fh1_Xpos->GetXaxis()->SetTitleSize(0.045);
    fh1_Xpos->GetYaxis()->SetLabelSize(0.045);
    fh1_Xpos->GetYaxis()->SetTitleSize(0.045);
    fh1_Xpos->SetFillColor(8);
    fh1_Xpos->SetLineColor(1);
    fh1_Xpos->Draw("col");

    // Hit data
    chity = new TCanvas(fNameDet + "_Ypos", "", 10, 10, 800, 700);
    Name1 = "fh1_" + fNameDet + "_Ypos";
    Name2 = fNameDet + ": Y (mm)";
    if (fNameDet == "Mwpc3")
        fh1_Ypos = R3B::root_owned<TH1F>(Name1, Name2, 1200, -300, 300);
    else
        fh1_Ypos = R3B::root_owned<TH1F>(Name1, Name2, 400, -100, 100);
    fh1_Ypos->GetXaxis()->SetTitle("Y (mm)");
    fh1_Ypos->GetYaxis()->SetTitle("counts per bin");
    fh1_Ypos->GetYaxis()->SetTitleOffset(1.1);
    fh1_Ypos->GetXaxis()->CenterTitle(true);
    fh1_Ypos->GetYaxis()->CenterTitle(true);
    fh1_Ypos->GetXaxis()->SetLabelSize(0.045);
    fh1_Ypos->GetXaxis()->SetTitleSize(0.045);
    fh1_Ypos->GetYaxis()->SetLabelSize(0.045);
    fh1_Ypos->GetYaxis()->SetTitleSize(0.045);
    fh1_Ypos->SetFillColor(8);
    fh1_Ypos->SetLineColor(1);
    fh1_Ypos->Draw("col");

    // Hit data
    chitxy = new TCanvas(fNameDet + "_XYpos", "", 10, 10, 800, 700);
    Name1 = "fh2_" + fNameDet + "_XYpos";
    Name2 = fNameDet + ": Y vs X (mm)";
    if (fNameDet == "Mwpc3")
        fh2_XYpos = R3B::root_owned<TH2F>(Name1, Name2, 1800, -450, 450, 1200, -300, 300);
    else
        fh2_XYpos = R3B::root_owned<TH2F>(Name1, Name2, 400, -100, 100, 400, -100, 100);
    fh2_XYpos->GetXaxis()->SetTitle("(Wixhausen)<---  X (mm)  ---> (Messel)");
    fh2_XYpos->GetYaxis()->SetTitle("Y (mm)");
    fh2_XYpos->GetYaxis()->SetTitleOffset(1.1);
    fh2_XYpos->GetXaxis()->CenterTitle(true);
    fh2_XYpos->GetYaxis()->CenterTitle(true);
    fh2_XYpos->GetXaxis()->SetLabelSize(0.045);
    fh2_XYpos->GetXaxis()->SetTitleSize(0.045);
    fh2_XYpos->GetYaxis()->SetLabelSize(0.045);
    fh2_XYpos->GetYaxis()->SetTitleSize(0.045);
    fh2_XYpos->Draw("col");

    // MAIN FOLDER-MWPC
    TFolder* mainfolMW = new TFolder(fNameDet, fNameDet + " info");
    mainfolMW->Add(cMwpc_nPads);
    mainfolMW->Add(cMwpcMap);
    mainfolMW->Add(cMWPCCal);
    mainfolMW->Add(cMWPCCal2D);
    mainfolMW->Add(cx);
    mainfolMW->Add(cy);
    if (fHitItemsMwpc)
    {
        mainfolMW->Add(chitx);
        mainfolMW->Add(chity);
        mainfolMW->Add(chitxy);
    }
    run->AddObject(mainfolMW);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_" + fNameDet + "_HIST",
                                          Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}

void R3BMwpcOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "for " << fNameDet);
    // Mapped data
    if (fMapItemsMwpc)
    {
        if (fNameDet == "Mwpc1" || fNameDet == "Mwpc2")
            for (Int_t i = 0; i < 3; i++)
            {
                fh1_mwpc_map_nPads[i]->Reset();
                fh2_mwpc_map_q[i]->Reset();
            }
        if (fNameDet == "Mwpc0" || fNameDet == "Mwpc3")
        {
            fh1_mwpc_map_nPads[0]->Reset();
            fh1_mwpc_map_nPads[2]->Reset();
            fh2_mwpc_map_q[0]->Reset();
            fh2_mwpc_map_q[2]->Reset();
        }
    }

    // Cal data
    if (fCalItemsMwpc)
    {
        for (Int_t i = 0; i < 2; i++)
            fh1_mwpc_cal[i]->Reset();
        fh2_mwpc_cal->Reset();
        fh2_mwpc_xq->Reset();
        fh2_mwpc_yq->Reset();
    }
    // Hit data
    if (fHitItemsMwpc)
    {
        fh1_Xpos->Reset();
        fh1_Ypos->Reset();
        fh2_XYpos->Reset();
    }
}

void R3BMwpcOnlineSpectra::Exec(Option_t* /*option*/)
{

    if (fMapItemsMwpc && fMapItemsMwpc->GetEntriesFast() > 0)
    {
        Int_t nPadsPerEvent[3];
        for (Int_t plane = 0; plane < 3; plane++)
        {
            nPadsPerEvent[plane] = 0;
        }
        Int_t nHits = fMapItemsMwpc->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMwpcMappedData* hit = dynamic_cast<R3BMwpcMappedData*>(fMapItemsMwpc->At(ihit));
            if (!hit)
                continue;
            nPadsPerEvent[hit->GetPlane() - 1]++;
            fh2_mwpc_map_q[hit->GetPlane() - 1]->Fill(hit->GetPad(), hit->GetQ());
        }
        if (fNameDet == "Mwpc1" || fNameDet == "Mwpc2")
        {
            for (int plane = 0; plane < 3; plane++)
                fh1_mwpc_map_nPads[plane]->Fill(nPadsPerEvent[plane]);
        }
        else
        {
            fh1_mwpc_map_nPads[0]->Fill(nPadsPerEvent[0]);
            fh1_mwpc_map_nPads[2]->Fill(nPadsPerEvent[2]);
        }
    }

    // Fill Cal data
    if (fCalItemsMwpc && fCalItemsMwpc->GetEntriesFast() > 0)
    {
        Int_t nHits = fCalItemsMwpc->GetEntriesFast();
        Int_t maxpadx = -1, maxpady = -1, maxqx = 0, maxqy = 0;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMwpcCalData* hit = dynamic_cast<R3BMwpcCalData*>(fCalItemsMwpc->At(ihit));
            if (!hit)
                continue;
            if (hit->GetPlane() == 1 || hit->GetPlane() == 2)
            // if (hit->GetPlane()==2)
            {
                fh1_mwpc_cal[0]->Fill(hit->GetPad());
                fh2_mwpc_xq->Fill(hit->GetPad() + gRandom->Uniform(-0.5, 0.5), hit->GetQ());
                if (hit->GetQ() > maxqx)
                {
                    maxpadx = hit->GetPad();
                    maxqx = hit->GetQ();
                }
            }
            if (hit->GetPlane() == 3)
            {
                fh1_mwpc_cal[1]->Fill(hit->GetPad());
                fh2_mwpc_yq->Fill(hit->GetPad() + gRandom->Uniform(-0.5, 0.5), hit->GetQ());
                if (hit->GetQ() > maxqy)
                {
                    maxpady = hit->GetPad();
                    maxqy = hit->GetQ();
                }
            }
        }
        if (maxpadx > -1 && maxpady > -1)
            fh2_mwpc_cal->Fill(maxpadx, maxpady);
    }

    // Fill Hit data
    if (fHitItemsMwpc && fHitItemsMwpc->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitItemsMwpc->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMwpcHitData* hit = dynamic_cast<R3BMwpcHitData*>(fHitItemsMwpc->At(ihit));
            if (!hit)
                continue;
            fh1_Xpos->Fill(hit->GetX());
            fh1_Ypos->Fill(hit->GetY());
            fh2_XYpos->Fill(hit->GetX(), hit->GetY());
        }
    }

    fNEvents += 1;
    return;
}

void R3BMwpcOnlineSpectra::FinishEvent()
{
    if (fMapItemsMwpc)
    {
        fMapItemsMwpc->Clear();
    }
    if (fCalItemsMwpc)
    {
        fCalItemsMwpc->Clear();
    }
    if (fHitItemsMwpc)
    {
        fHitItemsMwpc->Clear();
    }
}

void R3BMwpcOnlineSpectra::FinishTask()
{
    if (fMapItemsMwpc)
    {
        cMwpc_nPads->Write();
        cMwpcMap->Write();
    }

    if (fCalItemsMwpc)
    {
        cMWPCCal->Write();
        cMWPCCal2D->Write();
        fh2_mwpc_xq->Write();
        fh2_mwpc_yq->Write();
    }
    if (fHitItemsMwpc)
    {
        chitx->Write();
        chity->Write();
        chitxy->Write();
    }
}

ClassImp(R3BMwpcOnlineSpectra)
