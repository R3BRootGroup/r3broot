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
// -----               R3BTofDOnlineSpectra               -----
// -----    Created 18/03/22 by J.L. Rodriguez-Sanchez    -----
// -----    Updated 18/03/22 by M. Feijoo                 -----
// -----           Fill TofD online histograms            -----
// ------------------------------------------------------------

// ROOT headers
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFolder.h>
#include <TGaxis.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TMath.h>
#include <TVector3.h>

// FAIR headers
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>

// R3B headers
#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BShared.h"
#include "R3BTCalEngine.h"
#include "R3BTofDMappingPar.h"
#include "R3BTofDOnlineSpectra.h"
#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

namespace
{
    uint64_t n1, n2;
} // namespace

// R3BTofDOnlineSpectra::Default Constructor --------------------------
R3BTofDOnlineSpectra::R3BTofDOnlineSpectra()
    : R3BTofDOnlineSpectra("TofDOnlineSpectra", 1)
{
}

// R3BTofDOnlineSpectra::Standard Constructor --------------------------
R3BTofDOnlineSpectra::R3BTofDOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

// Virtual R3BTofDOnlineSpectra::Destructor
R3BTofDOnlineSpectra::~R3BTofDOnlineSpectra()
{
    R3BLOG(debug1, "Destructor");
    if (fMappedItems)
        delete fMappedItems;
    if (fCalItems)
        delete fCalItems;
    if (fHitItems)
        delete fHitItems;
}

void R3BTofDOnlineSpectra::SetParContainers()
{
    fMapPar = dynamic_cast<R3BTofDMappingPar*>(FairRuntimeDb::instance()->getContainer("tofdMappingPar"));
    R3BLOG_IF(warn, !fMapPar, "Could not get access to tofdMappingPar container");
    return;
}

void R3BTofDOnlineSpectra::SetParameter()
{
    R3BLOG_IF(info, fMapPar, "Nb of planes " << fMapPar->GetNbPlanes() << " and paddles " << fMapPar->GetNbPaddles());
    if (fMapPar)
    {
        fNofPlanes = fMapPar->GetNbPlanes();
        fPaddlesPerPlane = fMapPar->GetNbPaddles();
    }

    return;
}

InitStatus R3BTofDOnlineSpectra::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, mgr == nullptr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(error, header == nullptr, "EventHeader. not found");

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdTriggerCal"));
    R3BLOG_IF(warn, fCalTriggerItems == nullptr, "TofdTriggerCal not found");

    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdMapped"));
    R3BLOG_IF(fatal, fMappedItems == nullptr, "TofdMapped not found");

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdCal"));
    R3BLOG_IF(fatal, fCalItems == nullptr, "TofdCal not found");

    fHitItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdHit"));
    R3BLOG_IF(warn, fHitItems == nullptr, "TofdHit not found");

    SetParameter();

    // MAIN FOLDER-Twim-Foot
    auto* maintofd = new TFolder("TofD", "TofD info");

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------
    // TofD detector

    if (fMappedItems && fCalItems)
    {
        auto* cTofd_planes = new TCanvas("TofD_planes_Cal", "TOFD planes CAL data", 10, 10, 1100, 1000);
        cTofd_planes->Divide(6, fNofPlanes);

        fh_tofd_channels.resize(fNofPlanes);
        fh_tofd_TotPm.resize(fNofPlanes);
        fh_tofd_TotPm_coinc.resize(fNofPlanes);
        fh_num_side.resize(fNofPlanes);
        fh_tofd_multihit.resize(fNofPlanes);
        fh_tofd_multihit_coinc.resize(fNofPlanes);
        fh_tofd_dt.resize(fNofPlanes - 1);

        // Canvas to display the Y position as a function of paddle and per plane  ---------------
        auto* cTofd_Y_per_planes =
            new TCanvas("TofD_Ypos_planes_Cal", "TOFD: Y-pos per plane with CAL data", 10, 10, 1100, 1000);
        cTofd_Y_per_planes->Divide(2, fNofPlanes);

        fh2_tofd_ypos_cal.resize(fNofPlanes);
        fh2_tofd_timedif_cal.resize(fNofPlanes);

        for (Int_t j = 0; j < fNofPlanes; j++)
        {
            char strName1[255];
            sprintf(strName1, "tofd_channels_plane_%d", j + 1);
            char strName2[255];
            sprintf(strName2, "Tofd channels plane %d", j + 1);
            fh_tofd_channels[j] = R3B::root_owned<TH1F>(strName1, strName2, 90, -45., 45.);
            fh_tofd_channels[j]->GetXaxis()->SetTitle("Channel");
            fh_tofd_channels[j]->GetYaxis()->SetTitle("Counts");
            fh_tofd_channels[j]->GetYaxis()->SetTitleOffset(1.);
            fh_tofd_channels[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_channels[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_channels[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_channels[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_channels[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_channels[j]->GetYaxis()->SetTitleSize(0.045);
            fh_tofd_channels[j]->SetFillColor(31);

            char strName3[255];
            sprintf(strName3, "tofd_ToT_plane_%d", j + 1);
            char strName4[255];
            sprintf(strName4, "Tofd ToT plane %d", j + 1);
            fh_tofd_TotPm[j] = R3B::root_owned<TH2F>(strName3, strName4, 90, -45, 45, 1500, 0., 500.);
            fh_tofd_TotPm[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_TotPm[j]->GetYaxis()->SetTitle("ToT / ns");
            fh_tofd_TotPm[j]->GetYaxis()->SetTitleOffset(1.1);
            fh_tofd_TotPm[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_TotPm[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_TotPm[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_TotPm[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_TotPm[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_TotPm[j]->GetYaxis()->SetTitleSize(0.045);

            char strName5[255];
            sprintf(strName5, "tofd_ToT_coinc_plane_%d", j + 1);
            char strName6[255];
            sprintf(strName6, "Tofd ToT coinc plane %d", j + 1);
            fh_tofd_TotPm_coinc[j] = R3B::root_owned<TH2F>(strName5, strName6, 90, -45, 45, 1500, 0., 500.);
            fh_tofd_TotPm_coinc[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_TotPm_coinc[j]->GetYaxis()->SetTitle("ToT / ns");
            fh_tofd_TotPm_coinc[j]->GetYaxis()->SetTitleOffset(1.1);
            fh_tofd_TotPm_coinc[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_TotPm_coinc[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_TotPm_coinc[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_TotPm_coinc[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_TotPm_coinc[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_TotPm_coinc[j]->GetYaxis()->SetTitleSize(0.045);

            char strName7[255];
            sprintf(strName7, "tofd_multihit_plane_%d", j + 1);
            char strName8[255];
            sprintf(strName8, "Tofd multihit plane %d", j + 1);
            fh_tofd_multihit[j] = R3B::root_owned<TH2F>(strName7, strName8, 90, -45., 45., 30, 0, 30);
            fh_tofd_multihit[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_multihit[j]->GetYaxis()->SetTitle("Multihit");
            fh_tofd_multihit[j]->GetYaxis()->SetTitleOffset(1.);
            fh_tofd_multihit[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_multihit[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_multihit[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_multihit[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_multihit[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_multihit[j]->GetYaxis()->SetTitleSize(0.045);

            char strName9[255];
            sprintf(strName9, "tofd_multihit_coinc_plane_%d", j + 1);
            char strName10[255];
            sprintf(strName10, "Tofd multihit coinc plane %d", j + 1);
            fh_tofd_multihit_coinc[j] = R3B::root_owned<TH2F>(strName9, strName10, 45, 0., 45., 30, 0, 30);
            fh_tofd_multihit_coinc[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_multihit_coinc[j]->GetYaxis()->SetTitle("Multihit");
            fh_tofd_multihit_coinc[j]->GetYaxis()->SetTitleOffset(1.);
            fh_tofd_multihit_coinc[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_multihit_coinc[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_multihit_coinc[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_multihit_coinc[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_multihit_coinc[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_multihit_coinc[j]->GetYaxis()->SetTitleSize(0.045);

            if (j < fNofPlanes - 1)
            {
                int jk = j + 1;
                char strName11[255];
                sprintf(strName11, "tofd_ToF_plane_%d_%d", jk, jk + 1);
                char strName12[255];
                sprintf(strName12, "Tofd ToF plane %d and %d ", jk, jk + 1);
                fh_tofd_dt[j] = R3B::root_owned<TH2F>(strName11, strName12, 50, 0, 50, 4000, -2000., 2000);
                fh_tofd_dt[j]->GetXaxis()->SetTitle("Bar number");
                fh_tofd_dt[j]->GetYaxis()->SetTitle("dt / ns");
                fh_tofd_dt[j]->GetYaxis()->SetTitleOffset(1.);
                fh_tofd_dt[j]->GetXaxis()->CenterTitle(true);
                fh_tofd_dt[j]->GetYaxis()->CenterTitle(true);
                fh_tofd_dt[j]->GetXaxis()->SetLabelSize(0.045);
                fh_tofd_dt[j]->GetXaxis()->SetTitleSize(0.045);
                fh_tofd_dt[j]->GetYaxis()->SetLabelSize(0.045);
                fh_tofd_dt[j]->GetYaxis()->SetTitleSize(0.045);
            }

            char strName13[255];
            sprintf(strName13, "tofd_numHits_top_vs_bottom_%d", j + 1);
            char strName14[255];
            sprintf(strName14, "Tofd numHitsMapped top vs bottom %d", j + 1);
            fh_num_side[j] = R3B::root_owned<TH2F>(strName13, strName14, 45, 0., 45., 45, 0, 45);
            fh_num_side[j]->GetXaxis()->SetTitle("Num hits up");
            fh_num_side[j]->GetYaxis()->SetTitle("Num hits bottom");
            fh_num_side[j]->GetYaxis()->SetTitleOffset(1.);
            fh_num_side[j]->GetXaxis()->CenterTitle(true);
            fh_num_side[j]->GetYaxis()->CenterTitle(true);
            fh_num_side[j]->GetXaxis()->SetLabelSize(0.045);
            fh_num_side[j]->GetXaxis()->SetTitleSize(0.045);
            fh_num_side[j]->GetYaxis()->SetLabelSize(0.045);
            fh_num_side[j]->GetYaxis()->SetTitleSize(0.045);

            sprintf(strName13, "tofd_Ypos_plane_%d", j + 1);
            sprintf(strName14, "Tofd Ypos for plane %d", j + 1);
            fh2_tofd_ypos_cal[j] = R3B::root_owned<TH2F>(strName13, strName14, 45, 0, 45, 2000., -50., 50.);
            fh2_tofd_ypos_cal[j]->GetXaxis()->SetTitle("Bar number");
            fh2_tofd_ypos_cal[j]->GetYaxis()->SetTitle("Y-position [ns]");
            fh2_tofd_ypos_cal[j]->GetYaxis()->SetTitleOffset(1.05);
            fh2_tofd_ypos_cal[j]->GetXaxis()->CenterTitle(true);
            fh2_tofd_ypos_cal[j]->GetYaxis()->CenterTitle(true);
            fh2_tofd_ypos_cal[j]->GetXaxis()->SetLabelSize(0.045);
            fh2_tofd_ypos_cal[j]->GetXaxis()->SetTitleSize(0.045);
            fh2_tofd_ypos_cal[j]->GetYaxis()->SetLabelSize(0.045);
            fh2_tofd_ypos_cal[j]->GetYaxis()->SetTitleSize(0.045);
            cTofd_Y_per_planes->cd(j * 2 + 1);
            fh2_tofd_ypos_cal[j]->Draw("colz");

            sprintf(strName13, "tofd_timediff_plane_%d", j + 1);
            sprintf(strName14, "Tofd time diff. per PMT for plane %d", j + 1);
            fh2_tofd_timedif_cal[j] = R3B::root_owned<TH2F>(strName13, strName14, 90, -45, 45, 2000., -5000., 1000.);
            fh2_tofd_timedif_cal[j]->GetXaxis()->SetTitle("Bar number");
            fh2_tofd_timedif_cal[j]->GetYaxis()->SetTitle("PMT_time - Trig_time [ns]");
            fh2_tofd_timedif_cal[j]->GetYaxis()->SetTitleOffset(1.05);
            fh2_tofd_timedif_cal[j]->GetXaxis()->CenterTitle(true);
            fh2_tofd_timedif_cal[j]->GetYaxis()->CenterTitle(true);
            fh2_tofd_timedif_cal[j]->GetXaxis()->SetLabelSize(0.045);
            fh2_tofd_timedif_cal[j]->GetXaxis()->SetTitleSize(0.045);
            fh2_tofd_timedif_cal[j]->GetYaxis()->SetLabelSize(0.045);
            fh2_tofd_timedif_cal[j]->GetYaxis()->SetTitleSize(0.045);
            cTofd_Y_per_planes->cd(j * 2 + 2);
            fh2_tofd_timedif_cal[j]->Draw("colz");
        }

        cTofd_planes->cd(1);
        fh_tofd_channels[0]->Draw();
        cTofd_planes->cd(2);
        gPad->SetLogz();
        fh_tofd_TotPm[0]->Draw("colz");
        cTofd_planes->cd(3);
        gPad->SetLogz();
        fh_tofd_multihit[0]->SetAxisRange(1, 30, "Y");
        fh_tofd_multihit[0]->Draw("colz");
        cTofd_planes->cd(4);
        gPad->SetLogz();
        fh_tofd_TotPm_coinc[0]->Draw("colz");
        cTofd_planes->cd(5);
        gPad->SetLogz();
        fh_tofd_multihit_coinc[0]->SetAxisRange(1, 30, "Y");
        fh_tofd_multihit_coinc[0]->Draw("colz");

        cTofd_planes->cd(7);
        fh_tofd_channels[1]->Draw();
        cTofd_planes->cd(8);
        gPad->SetLogz();
        fh_tofd_TotPm[1]->Draw("colz");
        cTofd_planes->cd(9);
        gPad->SetLogz();
        fh_tofd_multihit[1]->SetAxisRange(1, 30, "Y");
        fh_tofd_multihit[1]->Draw("colz");
        cTofd_planes->cd(10);
        gPad->SetLogz();
        fh_tofd_TotPm_coinc[1]->Draw("colz");
        cTofd_planes->cd(11);
        gPad->SetLogz();
        fh_tofd_multihit_coinc[1]->SetAxisRange(1, 30, "Y");
        fh_tofd_multihit_coinc[1]->Draw("colz");
        cTofd_planes->cd(12);
        gPad->SetLogz();
        fh_tofd_dt[0]->Draw("colz");

        cTofd_planes->cd(13);
        fh_tofd_channels[2]->Draw();
        cTofd_planes->cd(14);
        gPad->SetLogz();
        fh_tofd_TotPm[2]->Draw("colz");
        cTofd_planes->cd(15);
        gPad->SetLogz();
        fh_tofd_multihit[2]->SetAxisRange(1, 30, "Y");
        fh_tofd_multihit[2]->Draw("colz");
        cTofd_planes->cd(16);
        gPad->SetLogz();
        fh_tofd_TotPm_coinc[2]->Draw("colz");
        cTofd_planes->cd(17);
        gPad->SetLogz();
        fh_tofd_multihit_coinc[2]->SetAxisRange(1, 30, "Y");
        fh_tofd_multihit_coinc[2]->Draw("colz");
        cTofd_planes->cd(18);
        gPad->SetLogz();
        fh_tofd_dt[1]->Draw("colz");

        cTofd_planes->cd(19);
        fh_tofd_channels[3]->Draw();
        cTofd_planes->cd(20);
        gPad->SetLogz();
        fh_tofd_TotPm[3]->Draw("colz");
        cTofd_planes->cd(21);
        gPad->SetLogz();
        fh_tofd_multihit[3]->SetAxisRange(1, 30, "Y");
        fh_tofd_multihit[3]->Draw("colz");
        cTofd_planes->cd(22);
        gPad->SetLogz();
        fh_tofd_TotPm_coinc[3]->Draw("colz");
        cTofd_planes->cd(23);
        gPad->SetLogz();
        fh_tofd_multihit_coinc[3]->SetAxisRange(1, 30, "Y");
        fh_tofd_multihit_coinc[3]->Draw("colz");
        cTofd_planes->cd(24);
        gPad->SetLogz();
        fh_tofd_dt[2]->Draw("colz");

        // Adding this canvas to the main folder
        maintofd->Add(cTofd_planes);
        maintofd->Add(cTofd_Y_per_planes);
    }

    if (fHitItems)
    {
        TCanvas* cTofd_planes_hit = new TCanvas("TofD_planes_Hit", "TOFD planes HIT data", 20, 20, 1120, 1020);
        cTofd_planes_hit->Divide(5, fNofPlanes);

        fh_tofd_Tot_hit.resize(fNofPlanes);
        fh_tofd_dt_hit.resize(fNofPlanes - 1);
        fh_tofd_multihit_hit.resize(fNofPlanes);
        fh_tofd_bars.resize(fNofPlanes);
        fh_tofd_time_hit.resize(fNofPlanes);

        for (Int_t j = 0; j < fNofPlanes; j++)
        {
            char strName3[255];
            sprintf(strName3, "tofd_hit_Q_plane_%d", j + 1);
            char strName4[255];
            sprintf(strName4, "Tofd hit Charge plane %d", j + 1);
            fh_tofd_Tot_hit[j] = R3B::root_owned<TH2F>(strName3, strName4, 45, 0, 45, 1000, 0., 12.);
            fh_tofd_Tot_hit[j]->GetXaxis()->SetTitle("BarId");
            fh_tofd_Tot_hit[j]->GetYaxis()->SetTitle("Charge");
            fh_tofd_Tot_hit[j]->GetYaxis()->SetTitleOffset(1.);
            fh_tofd_Tot_hit[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_Tot_hit[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_Tot_hit[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_Tot_hit[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_Tot_hit[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_Tot_hit[j]->GetYaxis()->SetTitleSize(0.045);

            char strName7[255];
            sprintf(strName7, "tofd_hit_multihit_plane_%d", j + 1);
            char strName8[255];
            sprintf(strName8, "Tofd hit multihit plane %d", j + 1);
            fh_tofd_multihit_hit[j] = R3B::root_owned<TH1F>(strName7, strName8, fMaxmul, 0, fMaxmul);
            fh_tofd_multihit_hit[j]->GetYaxis()->SetTitle("Counts");
            fh_tofd_multihit_hit[j]->GetXaxis()->SetTitle("Multihit");
            fh_tofd_multihit_hit[j]->GetYaxis()->SetTitleOffset(1.);
            fh_tofd_multihit_hit[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_multihit_hit[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_multihit_hit[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_multihit_hit[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_multihit_hit[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_multihit_hit[j]->GetYaxis()->SetTitleSize(0.045);
            fh_tofd_multihit_hit[j]->SetFillColor(31);

            char strName21[255];
            sprintf(strName21, "tofd_bars_plane_%d", j + 1);
            char strName22[255];
            sprintf(strName22, "Tofd bars plane %d", j + 1);
            fh_tofd_bars[j] = R3B::root_owned<TH1F>(strName21, strName22, 45, 0., 45.);
            fh_tofd_bars[j]->GetXaxis()->SetTitle("Wix<--- BarId --->Mes");
            fh_tofd_bars[j]->GetYaxis()->SetTitle("Counts");
            fh_tofd_bars[j]->GetYaxis()->SetTitleOffset(1.);
            fh_tofd_bars[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_bars[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_bars[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_bars[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_bars[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_bars[j]->GetYaxis()->SetTitleSize(0.045);
            fh_tofd_bars[j]->SetFillColor(31);

            char strName23[255];
            sprintf(strName23, "tofd_hit_time_plane_%d", j + 1);
            char strName24[255];
            sprintf(strName24, "Tofd hit time plane %d", j + 1);
            fh_tofd_time_hit[j] = R3B::root_owned<TH2F>(strName23, strName24, 45, 0, 45, 4000, -200., 200.);
            fh_tofd_time_hit[j]->GetXaxis()->SetTitle("BarId");
            fh_tofd_time_hit[j]->GetYaxis()->SetTitle("time / ns");
            fh_tofd_time_hit[j]->GetYaxis()->SetTitleOffset(1.1);
            fh_tofd_time_hit[j]->GetXaxis()->CenterTitle(true);
            fh_tofd_time_hit[j]->GetYaxis()->CenterTitle(true);
            fh_tofd_time_hit[j]->GetXaxis()->SetLabelSize(0.045);
            fh_tofd_time_hit[j]->GetXaxis()->SetTitleSize(0.045);
            fh_tofd_time_hit[j]->GetYaxis()->SetLabelSize(0.045);
            fh_tofd_time_hit[j]->GetYaxis()->SetTitleSize(0.045);

            if (j < fNofPlanes - 1)
            {
                int jk = j + 1;
                char strName11[255];
                sprintf(strName11, "tofd_hit_ToF_dt_plane_%d_%d", jk, jk + 1);
                char strName12[255];
                sprintf(strName12, "Tofd hit ToF dt plane %d and %d ", jk, jk + 1);
                fh_tofd_dt_hit[j] = R3B::root_owned<TH2F>(strName11, strName12, 45, 0, 45, 4000, -20., 20);
                char strName9[255];
                sprintf(strName9, "BarId plane %d ", j + 1);
                fh_tofd_dt_hit[j]->GetXaxis()->SetTitle(strName9);
                fh_tofd_dt_hit[j]->GetYaxis()->SetTitle("dt / ns");
                fh_tofd_dt_hit[j]->GetYaxis()->SetTitleOffset(1.);
                fh_tofd_dt_hit[j]->GetXaxis()->CenterTitle(true);
                fh_tofd_dt_hit[j]->GetYaxis()->CenterTitle(true);
                fh_tofd_dt_hit[j]->GetXaxis()->SetLabelSize(0.045);
                fh_tofd_dt_hit[j]->GetXaxis()->SetTitleSize(0.045);
                fh_tofd_dt_hit[j]->GetYaxis()->SetLabelSize(0.045);
                fh_tofd_dt_hit[j]->GetYaxis()->SetTitleSize(0.045);
            }
        }

        cTofd_planes_hit->cd(1);
        gPad->SetLogy();
        fh_tofd_multihit_hit[0]->Draw();
        cTofd_planes_hit->cd(3);
        gPad->SetLogz();
        fh_tofd_Tot_hit[0]->Draw("colz");
        cTofd_planes_hit->cd(2);
        gPad->SetLogy();
        fh_tofd_bars[0]->Draw();
        cTofd_planes_hit->cd(4);
        gPad->SetLogz();
        fh_tofd_time_hit[0]->Draw("colz");

        cTofd_planes_hit->cd(6);
        gPad->SetLogy();
        fh_tofd_multihit_hit[1]->Draw();
        cTofd_planes_hit->cd(8);
        gPad->SetLogz();
        fh_tofd_Tot_hit[1]->Draw("colz");
        cTofd_planes_hit->cd(7);
        gPad->SetLogy();
        fh_tofd_bars[1]->Draw();
        cTofd_planes_hit->cd(9);
        gPad->SetLogz();
        fh_tofd_time_hit[1]->Draw("colz");
        cTofd_planes_hit->cd(10);
        gPad->SetLogz();
        fh_tofd_dt_hit[0]->Draw("colz");

        cTofd_planes_hit->cd(11);
        gPad->SetLogy();
        fh_tofd_multihit_hit[2]->Draw();
        cTofd_planes_hit->cd(13);
        gPad->SetLogz();
        fh_tofd_Tot_hit[2]->Draw("colz");
        cTofd_planes_hit->cd(12);
        gPad->SetLogy();
        fh_tofd_bars[2]->Draw();
        cTofd_planes_hit->cd(14);
        gPad->SetLogz();
        fh_tofd_time_hit[2]->Draw("colz");
        cTofd_planes_hit->cd(15);
        gPad->SetLogz();
        fh_tofd_dt_hit[1]->Draw("colz");

        cTofd_planes_hit->cd(16);
        gPad->SetLogy();
        fh_tofd_multihit_hit[3]->Draw();
        cTofd_planes_hit->cd(18);
        gPad->SetLogz();
        fh_tofd_Tot_hit[3]->Draw("colz");
        cTofd_planes_hit->cd(17);
        gPad->SetLogy();
        fh_tofd_bars[3]->Draw();
        cTofd_planes_hit->cd(19);
        gPad->SetLogz();
        fh_tofd_time_hit[3]->Draw("colz");
        cTofd_planes_hit->cd(20);
        gPad->SetLogz();
        fh_tofd_dt_hit[2]->Draw("colz");

        // Adding this canvas to the main folder
        maintofd->Add(cTofd_planes_hit);

        auto* cToFd_los_h2 = new TCanvas("ToFD_time_Los_time", "ToFD time - Los time", 20, 20, 1120, 1020);
        cToFd_los_h2->Divide(2, 2);
        fh_tofd_time_los_h2.resize(fNofPlanes);

        auto* cToFd_los_h2_wt =
            new TCanvas("ToFD_Los_time_without_trigger", "ToFD time - Los time without trigger", 20, 20, 1120, 1020);
        cToFd_los_h2_wt->Divide(2, 2);
        fh2_tofd_time_los_cal.resize(fNofPlanes);

        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            fh_tofd_time_los[j].resize(fNofPlanes);

        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            char strNameLos_c[455];
            sprintf(strNameLos_c, "tofd_los_timediff_plane_%d", i + 1);
            fh_tofd_time_los_h2[i] = R3B::root_owned<TH2F>(strNameLos_c, strNameLos_c, 44, 1, 45, 20000, 40, 190);
            fh_tofd_time_los_h2[i]->GetXaxis()->SetTitle("Bar");
            fh_tofd_time_los_h2[i]->GetYaxis()->SetTitle("ToF [ns]");
            fh_tofd_time_los_h2[i]->GetXaxis()->CenterTitle(true);
            fh_tofd_time_los_h2[i]->GetYaxis()->CenterTitle(true);
            cToFd_los_h2->cd(i + 1);
            gPad->SetLogz();
            fh_tofd_time_los_h2[i]->Draw("colz");

            char strNameLos_c2[255];
            snprintf(strNameLos_c2, sizeof(strNameLos_c2), "tofd_los_time_without_trigger_%d", i + 1);
            fh2_tofd_time_los_cal[i] = R3B::root_owned<TH2F>(strNameLos_c2, strNameLos_c2, 44, 1, 45, 10000, 40, 70);
            fh2_tofd_time_los_cal[i]->GetXaxis()->SetTitle("Bar");
            fh2_tofd_time_los_cal[i]->GetYaxis()->SetTitle("ToF [ns]");
            fh2_tofd_time_los_cal[i]->GetXaxis()->CenterTitle(true);
            fh2_tofd_time_los_cal[i]->GetYaxis()->CenterTitle(true);
            cToFd_los_h2_wt->cd(i + 1);
            gPad->SetLogz();
            fh2_tofd_time_los_cal[i]->Draw("colz");

            auto cToFd_los = new TCanvas(strNameLos_c, strNameLos_c, 20, 20, 1120, 1020);
            cToFd_los->Divide(5, 9);
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                char strNameLos[255];
                sprintf(strNameLos, "tofd_los_timediff_bar_%d_plane_%d", j + 1, i + 1);
                char strNameLos2[255];
                sprintf(strNameLos2, "Tofd_time - Los_time bar %d plane %d", j + 1, i + 1);
                fh_tofd_time_los[j][i] = R3B::root_owned<TH1F>(strNameLos, strNameLos2, 20000, 40, 190);
                fh_tofd_time_los[j][i]->GetXaxis()->SetTitle("ToF [ns]");
                fh_tofd_time_los[j][i]->GetYaxis()->SetTitle("counts");
                fh_tofd_time_los[j][i]->SetFillColor(31);
                cToFd_los->cd(j + 1);
                fh_tofd_time_los[j][i]->Draw("");
            }
            // Adding this canvas to the main folder
            maintofd->Add(cToFd_los);
        }
        maintofd->Add(cToFd_los_h2);
        maintofd->Add(cToFd_los_h2_wt);

        auto cToFd_time_charge = new TCanvas("tofd_time_vs_charge", "", 20, 20, 1120, 1020);
        fh2_tofd_time_vs_charge =
            R3B::root_owned<TH2F>("fh2_tofd_time_vs_charge", "Time vs Charge", 10000, 50, 150, 1000, 0, 12);
        fh2_tofd_time_vs_charge->GetXaxis()->SetTitle("ToF [ns]");
        fh2_tofd_time_vs_charge->GetYaxis()->SetTitle("Charge");
        fh2_tofd_time_vs_charge->GetXaxis()->CenterTitle(true);
        fh2_tofd_time_vs_charge->GetYaxis()->CenterTitle(true);
        gPad->SetLogz();
        fh2_tofd_time_vs_charge->Draw("colz");
        maintofd->Add(cToFd_time_charge);

        auto cToFd_time_wtrigger_charge = new TCanvas("tofd_time_withoutTrig_vs_charge", "", 20, 20, 1120, 1020);
        fh2_tofd_time_wouttrig_vs_charge =
            R3B::root_owned<TH2F>("fh2_tofd_time_vs_charge_withoutTrig", "Time vs Charge", 10000, 25, 70, 1000, 0, 12);
        fh2_tofd_time_wouttrig_vs_charge->GetXaxis()->SetTitle("ToF [ns]");
        fh2_tofd_time_wouttrig_vs_charge->GetYaxis()->SetTitle("Charge");
        fh2_tofd_time_wouttrig_vs_charge->GetXaxis()->CenterTitle(true);
        fh2_tofd_time_wouttrig_vs_charge->GetYaxis()->CenterTitle(true);
        gPad->SetLogz();
        fh2_tofd_time_wouttrig_vs_charge->Draw("colz");
        maintofd->Add(cToFd_time_wtrigger_charge);

        auto cToFd_charge_cor = new TCanvas("tofd_Charge_correlations", "", 20, 20, 1120, 1020);
        cToFd_charge_cor->Divide(2, 1);
        cToFd_charge_cor->cd(1);
        fh2_tofd_charges13 =
            R3B::root_owned<TH2F>("fh2_tofd_charge1_vs_charge3", "Charges: Plane1 vs Plane3", 1000, 0, 12, 1000, 0, 12);
        fh2_tofd_charges13->GetXaxis()->SetTitle("Charge-Z P1");
        fh2_tofd_charges13->GetYaxis()->SetTitle("Charge-Z P3");
        fh2_tofd_charges13->GetXaxis()->CenterTitle(true);
        fh2_tofd_charges13->GetYaxis()->CenterTitle(true);
        fh2_tofd_charges13->Draw("colz");
        fh2_tofd_charges13->SetStats(0);

        cToFd_charge_cor->cd(2);
        fh2_tofd_charges24 =
            R3B::root_owned<TH2F>("fh2_tofd_charge2_vs_charge4", "Charges: Plane2 vs Plane4", 1000, 0, 12, 1000, 0, 12);
        fh2_tofd_charges24->GetXaxis()->SetTitle("Charge-Z P2");
        fh2_tofd_charges24->GetYaxis()->SetTitle("Charge-Z P4");
        fh2_tofd_charges24->GetXaxis()->CenterTitle(true);
        fh2_tofd_charges24->GetYaxis()->CenterTitle(true);
        fh2_tofd_charges24->Draw("colz");
        fh2_tofd_charges24->SetStats(0);

        maintofd->Add(cToFd_charge_cor);

        auto cToFd_charge_cor_v2 = new TCanvas("tofd_Charge_correlations_all", "", 20, 20, 1120, 1020);
        cToFd_charge_cor_v2->Divide(3, 1);
        cToFd_charge_cor_v2->cd(1);
        fh2_tofd_charges12 =
            R3B::root_owned<TH2F>("fh2_tofd_charge1_vs_charge2", "Charges: Plane1 vs Plane2", 1000, 0, 12, 1000, 0, 12);
        fh2_tofd_charges12->GetXaxis()->SetTitle("Charge-Z P1");
        fh2_tofd_charges12->GetYaxis()->SetTitle("Charge-Z P2");
        fh2_tofd_charges12->GetXaxis()->CenterTitle(true);
        fh2_tofd_charges12->GetYaxis()->CenterTitle(true);
        fh2_tofd_charges12->Draw("colz");
        fh2_tofd_charges12->SetStats(0);

        cToFd_charge_cor_v2->cd(2);
        fh2_tofd_charges23 =
            R3B::root_owned<TH2F>("fh2_tofd_charge2_vs_charge3", "Charges: Plane2 vs Plane3", 1000, 0, 12, 1000, 0, 12);
        fh2_tofd_charges23->GetXaxis()->SetTitle("Charge-Z P2");
        fh2_tofd_charges23->GetYaxis()->SetTitle("Charge-Z P3");
        fh2_tofd_charges23->GetXaxis()->CenterTitle(true);
        fh2_tofd_charges23->GetYaxis()->CenterTitle(true);
        fh2_tofd_charges23->Draw("colz");
        fh2_tofd_charges23->SetStats(0);

        cToFd_charge_cor_v2->cd(3);
        fh2_tofd_charges34 =
            R3B::root_owned<TH2F>("fh2_tofd_charge3_vs_charge4", "Charges: Plane3 vs Plane4", 1000, 0, 12, 1000, 0, 12);
        fh2_tofd_charges34->GetXaxis()->SetTitle("Charge-Z P3");
        fh2_tofd_charges34->GetYaxis()->SetTitle("Charge-Z P4");
        fh2_tofd_charges34->GetXaxis()->CenterTitle(true);
        fh2_tofd_charges34->GetYaxis()->CenterTitle(true);
        fh2_tofd_charges34->Draw("colz");
        fh2_tofd_charges34->SetStats(0);

        maintofd->Add(cToFd_charge_cor_v2);

        auto cToFd_zcharge = new TCanvas("tofd_Charge_Z", "", 20, 20, 1120, 1020);
        fh1_Zcharge = R3B::root_owned<TH1F>("fh1_tofd_chargeZ", "Charges Z from Plane 1", 1000, 0, 12);
        fh1_Zcharge->GetXaxis()->SetTitle("Charge-Z");
        fh1_Zcharge->GetYaxis()->SetTitle("Counts");
        fh1_Zcharge->GetXaxis()->CenterTitle(true);
        fh1_Zcharge->GetYaxis()->CenterTitle(true);
        fh1_Zcharge->SetFillColor(31);
        fh1_Zcharge->Draw("");

        maintofd->Add(cToFd_zcharge);

        auto cToFd_zcharge_tpat = new TCanvas("Tofd_chargeZ_vs_tpat", "Charges Z vs tpat", 10, 10, 500, 500);
        fh2_zcharge_tpat =
            R3B::root_owned<TH2F>("fh2_Tofd_chargeZ_vs_tpat", "Charges Z vs tpat", 17, -0.5, 16.5, 1000, 0, 12);
        fh2_zcharge_tpat->GetXaxis()->SetTitle("Tpat");
        fh2_zcharge_tpat->GetYaxis()->SetTitle("Charge-Z");
        fh2_zcharge_tpat->GetXaxis()->CenterTitle(true);
        fh2_zcharge_tpat->GetYaxis()->CenterTitle(true);
        fh2_zcharge_tpat->GetYaxis()->SetTitleOffset(1.2);
        fh2_zcharge_tpat->GetXaxis()->SetTitleOffset(1.2);
        fh2_zcharge_tpat->Draw("colz");

        maintofd->Add(cToFd_zcharge_tpat);
    }

    run->AddObject(maintofd);
    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_TofD_HIST", Form("/Objects/%s/->Reset_Histo()", GetName()));

    // Definition of a time stich object to correlate times coming from different systems
    fTimeStitch = std::make_unique<R3BCoarseTimeStitch>();

    return kSUCCESS;
}

void R3BTofDOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "");
    for (int i = 0; i < fNofPlanes; i++)
    {
        fh_tofd_channels[i]->Reset();
        fh_tofd_multihit[i]->Reset();
        fh_tofd_TotPm[i]->Reset();
        fh_tofd_multihit_coinc[i]->Reset();
        fh_tofd_TotPm_coinc[i]->Reset();
        fh2_tofd_ypos_cal[i]->Reset();
        fh2_tofd_timedif_cal[i]->Reset();
    }
    fh_tofd_dt[0]->Reset();
    fh_tofd_dt[1]->Reset();
    fh_tofd_dt[2]->Reset();

    if (fHitItems)
    {
        for (int i = 0; i < fNofPlanes; i++)
        {
            fh_tofd_Tot_hit[i]->Reset();
            fh_tofd_time_hit[i]->Reset();
            fh_tofd_multihit_hit[i]->Reset();
            fh_tofd_bars[i]->Reset();
            fh_tofd_time_los_h2[i]->Reset();
            fh2_tofd_time_los_cal[i]->Reset();
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
                fh_tofd_time_los[j][i]->Reset();
        }
        for (int i = 0; i < fNofPlanes - 1; i++)
        {
            fh_tofd_dt_hit[i]->Reset();
        }
        fh2_tofd_time_vs_charge->Reset();
        fh2_tofd_time_wouttrig_vs_charge->Reset();
        fh2_tofd_charges13->Reset();
        fh2_tofd_charges24->Reset();
        fh2_tofd_charges12->Reset();
        fh2_tofd_charges23->Reset();
        fh2_tofd_charges34->Reset();
        fh1_Zcharge->Reset();
        fh2_zcharge_tpat->Reset();
    }
    return;
}

void R3BTofDOnlineSpectra::Exec(Option_t* option)
{
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;
    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit1 = fTpat1 - 1;
    Int_t fTpat_bit2 = fTpat2 - 1;
    Int_t tpatbin = 0;
    std::vector<int> tpatindex;
    if (header && fTpat1 >= 0 && fTpat2 >= 0)
    {
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));

            if (tpatbin != 0)
                tpatindex.push_back(i + 1);

            if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2))
            {
                return;
            }
        }
    }
    else if (header)
    {
        for (Int_t i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0)
                tpatindex.push_back(i + 1);
        }
    }

    Int_t NumPaddles[fNofPlanes];
    for (int i = 0; i < fNofPlanes; i++)
        NumPaddles[i] = 0;

    if (fMappedItems)
    {
        Int_t nMapped = fMappedItems->GetEntriesFast();
        Int_t iPlaneMem = 1, iBarMem = 0;

        Int_t nsum_top[fNofPlanes];
        Int_t nsum_bot[fNofPlanes];
        for (Int_t p = 0; p < fNofPlanes; p++)
        {
            nsum_top[p] = 0;
            nsum_bot[p] = 0;
        }

        for (Int_t imapped = 0; imapped < nMapped; imapped++)
        {
            auto mapped = dynamic_cast<R3BTofdMappedData*>(fMappedItems->At(imapped));
            if (!mapped)
                continue; // should not happen

            Int_t const iPlane = mapped->GetDetectorId(); // 1..n
            Int_t const iBar = mapped->GetBarId();        // 1..n
            Int_t const iSide = mapped->GetSideId();      // 1..n
            Int_t const iEdge = mapped->GetEdgeId();

            if (iSide == 1 && iEdge == 1 && iPlane < fNofPlanes)
                nsum_bot[iPlane - 1] += 1;
            if (iSide == 2 && iEdge == 1 && iPlane < fNofPlanes)
                nsum_top[iPlane - 1] += 1;

            if (iPlaneMem != iPlane)
                iBarMem = 0;
            if (iPlane == 1 && iBarMem != iBar && iSide == 1 && iEdge == 1)
                NumPaddles[0] += 1;
            if (iPlane == 2 && iBarMem != iBar && iSide == 1 && iEdge == 1)
                NumPaddles[1] += 1;
            if (iPlane == 3 && iBarMem != iBar && iSide == 1 && iEdge == 1)
                NumPaddles[2] += 1;
            if (iPlane == 3 && iBarMem != iBar && iSide == 1 && iEdge == 1)
                NumPaddles[3] += 1;

            iPlaneMem = iPlane;
            iBarMem = iBar;

            if (iPlane <= fNofPlanes)
            {
                if (iSide == 1) // bottom
                    fh_tofd_channels[iPlane - 1]->Fill(-iBar - 1);
                if (iSide == 2) // top
                    fh_tofd_channels[iPlane - 1]->Fill(iBar);
            }
        }
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            fh_num_side[i]->Fill(nsum_bot[i], nsum_top[i]);
        }
    }

    std::vector<double> fTof_without_trig(44, 0.);
    if (fCalItems)
    {
        UInt_t vmultihits[fNofPlanes + 1][fPaddlesPerPlane];
        UInt_t vmultihits_top[fNofPlanes + 1][fPaddlesPerPlane];
        UInt_t vmultihits_bot[fNofPlanes + 1][fPaddlesPerPlane];
        Double_t time_bar[fNofPlanes + 1][fPaddlesPerPlane][fMaxmul];
        for (Int_t i = 0; i < fNofPlanes + 1; i++)
        {
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                vmultihits[i][j] = 0;
                vmultihits_top[i][j] = 0;
                vmultihits_bot[i][j] = 0;
                for (Int_t l = 0; l < fMaxmul; l++)
                {
                    time_bar[i][j][l] = 0. / 0.;
                }
            }
        }

        //    std::cout<<"new event!*************************************\n";
        Int_t nHits = fCalItems->GetEntriesFast();

        Int_t nHitsEvent = 0;
        // Organize cals into bars.
        struct Entry
        {
            std::vector<R3BTofdCalData*> top;
            std::vector<R3BTofdCalData*> bot;
        };

        std::map<size_t, Entry> bar_map;
        //   puts("Event");
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto* hit = dynamic_cast<R3BTofdCalData*>(fCalItems->At(ihit));
            size_t idx = (hit->GetDetectorId() - 1) * fPaddlesPerPlane + hit->GetBarId() - 1;

            auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));

            // hit->GetSideId() gives 1 for bottom and 2 for top
            auto& vec = 1 == hit->GetSideId() ? ret.first->second.bot : ret.first->second.top;
            vec.push_back(hit);
        }

        // Build trigger map.
        std::vector<R3BTofdCalData*> trig_map;
        for (int i = 0; i < fCalTriggerItems->GetEntries(); ++i)
        {
            auto trig = dynamic_cast<R3BTofdCalData*>(fCalTriggerItems->At(i));
            if (trig_map.size() < trig->GetBarId())
            {
                trig_map.resize(trig->GetBarId());
            }
            trig_map.at(trig->GetBarId() - 1) = trig;
        }

        // Without coincidences
        bool s_was_trig_missing = false;

        for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
        {
            auto const& top_vec = it->second.top;
            size_t top_i = 0;
            auto const& bot_vec = it->second.bot;
            size_t bot_i = 0;
            for (; top_i < top_vec.size();)
            {
                auto top = top_vec.at(top_i);
                Int_t top_trig_i = 0;

                if (fMapPar)
                    top_trig_i = fMapPar->GetTrigMap(top->GetDetectorId(), top->GetBarId(), top->GetSideId());

                Double_t top_trig_ns = 0;
                if (top_trig_i < trig_map.size() && trig_map.at(top_trig_i))
                {
                    auto top_trig = trig_map.at(top_trig_i);
                    top_trig_ns = top_trig->GetTimeLeading_ns();
                }
                else
                {
                    if (!s_was_trig_missing)
                    {
                        R3BLOG(error, "Missing trigger information!");
                        s_was_trig_missing = true;
                    }
                }

                Int_t iPlane = top->GetDetectorId(); // 1..n
                Int_t iBar = top->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    R3BLOG(error, "More detectors than expected! Det: " << iPlane << " allowed are 1.." << fNofPlanes);
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    R3BLOG(error, "More bars then expected! Det: " << iBar << " allowed are 1.." << fPaddlesPerPlane);
                    continue;
                }

                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + fC_range_ns, fC_range_ns);
                fh_tofd_TotPm[iPlane - 1]->Fill(iBar, top_tot);
                vmultihits_top[iPlane - 1][iBar - 1] += 1;
                ++top_i;
            }

            for (; bot_i < bot_vec.size();)
            {
                auto bot = bot_vec.at(bot_i);
                Int_t bot_trig_i = 0;

                if (fMapPar)
                    bot_trig_i = fMapPar->GetTrigMap(bot->GetDetectorId(), bot->GetBarId(), bot->GetSideId());

                Double_t bot_trig_ns = 0;
                if (bot_trig_i < trig_map.size() && trig_map.at(bot_trig_i))
                {
                    auto bot_trig = trig_map.at(bot_trig_i);
                    bot_trig_ns = bot_trig->GetTimeLeading_ns();
                }
                else
                {
                    if (!s_was_trig_missing)
                    {
                        R3BLOG(error, "Missing trigger information!");
                        s_was_trig_missing = true;
                    }
                }

                // Shift the cyclic difference window by half a window-length and move it back,
                // this way the trigger time will be at 0.
                auto bot_ns =
                    fmod(bot->GetTimeLeading_ns() - bot_trig_ns + fC_range_ns + fC_range_ns / 2, fC_range_ns) -
                    fC_range_ns / 2;

                Int_t iPlane = bot->GetDetectorId(); // 1..n
                Int_t iBar = bot->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    R3BLOG(error, "More detectors than expected! Det: " << iPlane << " allowed are 1.." << fNofPlanes);
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    R3BLOG(error, "More bars then expected! Det: " << iBar << " allowed are 1.." << fPaddlesPerPlane);
                    continue;
                }

                auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + fC_range_ns, fC_range_ns);

                fh_tofd_TotPm[iPlane - 1]->Fill(-iBar - 1, bot_tot);

                // register multi hits
                vmultihits_bot[iPlane - 1][iBar - 1] += 1;

                ++bot_i;
            }
        }

        for (Int_t ipl = 0; ipl < fNofPlanes; ipl++)
        {
            for (Int_t ibr = 1; ibr < fPaddlesPerPlane + 1; ibr++)
            {
                fh_tofd_multihit[ipl]->Fill(-ibr - 1, vmultihits_bot[ipl][ibr - 1]);
                fh_tofd_multihit[ipl]->Fill(ibr, vmultihits_top[ipl][ibr - 1]);
            }
        }

        // With coincidences:
        bool s_was_trig_missingc = false;
        for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
        {
            auto const& topc_vec = it->second.top;
            auto const& botc_vec = it->second.bot;
            size_t topc_i = 0;
            size_t botc_i = 0;
            for (; topc_i < topc_vec.size() && botc_i < botc_vec.size();)
            {
                auto topc = topc_vec.at(topc_i);
                auto botc = botc_vec.at(botc_i);

                Int_t topc_trig_i = 0;
                Int_t botc_trig_i = 0;
                if (fMapPar)
                {
                    topc_trig_i = fMapPar->GetTrigMap(topc->GetDetectorId(), topc->GetBarId(), topc->GetSideId());
                    botc_trig_i = fMapPar->GetTrigMap(botc->GetDetectorId(), botc->GetBarId(), botc->GetSideId());
                }

                Double_t topc_trig_ns = 0, botc_trig_ns = 0;
                if (topc_trig_i < trig_map.size() && trig_map.at(topc_trig_i) && botc_trig_i < trig_map.size() &&
                    trig_map.at(botc_trig_i))
                {
                    auto topc_trig = trig_map.at(topc_trig_i);
                    auto botc_trig = trig_map.at(botc_trig_i);
                    topc_trig_ns = topc_trig->GetTimeLeading_ns();
                    botc_trig_ns = botc_trig->GetTimeLeading_ns();
                    ++n1;
                }
                else
                {
                    if (!s_was_trig_missingc)
                    {
                        R3BLOG(error, "Missing trigger information!");
                        s_was_trig_missingc = true;
                    }
                    ++n2;
                }

                // Shift the cyclic difference window by half a window-length and move it back,
                // this way the trigger time will be at 0.
                auto topc_ns = fTimeStitch->GetTime(topc->GetTimeLeading_ns() - topc_trig_ns);
                auto botc_ns = fTimeStitch->GetTime(botc->GetTimeLeading_ns() - botc_trig_ns);

                auto dt = topc_ns - botc_ns;
                // Handle wrap-around.
                auto dt_mod = fmod(dt + fC_range_ns, fC_range_ns);
                if (dt < 0)
                {
                    // We're only interested in the short time-differences, so we
                    // want to move the upper part of the coarse counter range close
                    // to the lower range, i.e. we cut the middle of the range and
                    // glue zero and the largest values together.
                    dt_mod -= fC_range_ns;
                }

                if (std::abs(dt_mod) < fC_bar_coincidence_ns * 10.)
                {
                    int iPlane = topc->GetDetectorId(); // 1..n
                    int iBar = topc->GetBarId();        // 1..n
                    // Histograms to display Y position
                    fh2_tofd_ypos_cal[iPlane - 1]->Fill(iBar, dt_mod);
                    fh2_tofd_timedif_cal[iPlane - 1]->Fill(iBar, topc_ns);
                    fh2_tofd_timedif_cal[iPlane - 1]->Fill(-1 * iBar, botc_ns);

                    if (topc->GetTimeLeading_ns() > 0 && header->GetTStartSimple() > 0)
                    {
                        auto tof_without_trig_top =
                            fTimeStitch->GetTime(topc->GetTimeLeading_ns() - header->GetTStartSimple());
                        auto tof_without_trig_bot =
                            fTimeStitch->GetTime(botc->GetTimeLeading_ns() - header->GetTStartSimple());
                        auto mean_tof_trig = (tof_without_trig_top + tof_without_trig_bot) / 2.;
                        if (iPlane == 1)
                        {
                            fTof_without_trig[iBar - 1] = mean_tof_trig + fTofcor[iBar - 1];
                            fh2_tofd_time_los_cal[topc->GetDetectorId() - 1]->Fill(topc->GetBarId(),
                                                                                   fTof_without_trig[iBar - 1]);
                        }
                        else
                            fh2_tofd_time_los_cal[topc->GetDetectorId() - 1]->Fill(
                                topc->GetBarId(), mean_tof_trig + fTofcor[44 * (topc->GetDetectorId() - 1) + iBar - 1]);
                    }
                }

                if (std::abs(dt_mod) < fC_bar_coincidence_ns)
                {
                    // Hit!
                    // std::cout << "Hit!\n";
                    Int_t iPlane = topc->GetDetectorId(); // 1..n
                    Int_t iBar = topc->GetBarId();        // 1..n
                    if (iPlane > fNofPlanes)              // this also errors for iDetector==0
                    {
                        R3BLOG(error,
                               "More detectors than expected! Det: " << iPlane << " allowed are 1.." << fNofPlanes);
                        continue;
                    }
                    if (iBar > fPaddlesPerPlane) // same here
                    {
                        R3BLOG(error,
                               "More bars then expected! Det: " << iBar << " allowed are 1.." << fPaddlesPerPlane);
                        continue;
                    }

                    auto topc_tot = fTimeStitch->GetTime(topc->GetTimeTrailing_ns() - topc->GetTimeLeading_ns());
                    auto botc_tot = fTimeStitch->GetTime(botc->GetTimeTrailing_ns() - botc->GetTimeLeading_ns());

                    fh_tofd_TotPm_coinc[iPlane - 1]->Fill(-iBar - 1, botc_tot);
                    fh_tofd_TotPm_coinc[iPlane - 1]->Fill(iBar, topc_tot);

                    // std::cout<<"ToT: "<<top_tot << " "<<bot_tot<<"\n";

                    // register multi hits
                    Int_t imlt = vmultihits[iPlane - 1][iBar - 1];
                    time_bar[iPlane - 1][iBar - 1][imlt] = (topc_ns + botc_ns) / 2.;
                    vmultihits[iPlane - 1][iBar - 1] += 1;

                    ++topc_i;
                    ++botc_i;
                }
                else if (dt < 0 && dt > -fC_range_ns / 2)
                {
                    ++topc_i;
                }
                else
                {
                    ++botc_i;
                }
            }
        }

        for (Int_t ipl = 0; ipl < fNofPlanes; ipl++)
        {
            for (Int_t ibr = 1; ibr < fPaddlesPerPlane + 1; ibr++)
            {
                fh_tofd_multihit_coinc[ipl]->Fill(ibr, vmultihits[ipl][ibr - 1]);
                if (ipl > 0)
                {
                    for (Int_t imult1 = 0; imult1 < vmultihits[ipl][ibr - 1]; imult1++)
                    {
                        for (Int_t imult2 = 0; imult2 < vmultihits[ipl - 1][ibr - 1]; imult2++)
                        {
                            Double_t tof_plane = 0. / 0.;
                            tof_plane = fTimeStitch->GetTime(time_bar[ipl][ibr - 1][imult1] -
                                                             time_bar[ipl - 1][ibr - 1][imult2]);
                            fh_tofd_dt[ipl - 1]->Fill(ibr, tof_plane);
                        }
                    }
                }
            }
        }

    } // endi if fCalItems

    if (fHitItems)
    {
        Int_t nHits = fHitItems->GetEntriesFast();
        if (nHits > fMaxmul)
            return;

        Double_t x[fNofPlanes][fMaxmul], y[fNofPlanes][fMaxmul], t[fNofPlanes][fMaxmul], q[fNofPlanes][fMaxmul],
            bar[fNofPlanes][fMaxmul];
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            for (Int_t k = 0; k < fMaxmul; k++)
            {
                x[i][k] = -1000.;
                y[i][k] = -1000.;
                q[i][k] = -1000.;
                t[i][k] = 0. / 0.;
                bar[i][k] = 0;
            }
        }

        Int_t nMulti[fNofPlanes], iCounts[fNofPlanes];
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            nMulti[i] = 0;
            iCounts[i] = 0;
        }

        double charges[4] = { 0 };
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hitTofd = dynamic_cast<R3BTofdHitData*>(fHitItems->At(ihit));
            if (IS_NAN(hitTofd->GetTime()))
                continue;
            Int_t iPlane = hitTofd->GetDetId();
            Double_t randx = (std::rand() / (float)RAND_MAX) - 0.5;
            Int_t ictemp = iCounts[iPlane - 1];
            x[iPlane - 1][ictemp] = hitTofd->GetX() + 2.7 * randx;
            y[iPlane - 1][ictemp] = hitTofd->GetY();
            t[iPlane - 1][ictemp] = hitTofd->GetTime();
            q[iPlane - 1][ictemp] = hitTofd->GetEloss();
            bar[iPlane - 1][ictemp] = hitTofd->GetBarId();
            fh_tofd_Tot_hit[iPlane - 1]->Fill(bar[iPlane - 1][ictemp], q[iPlane - 1][ictemp]);
            fh_tofd_time_hit[iPlane - 1]->Fill(bar[iPlane - 1][ictemp], t[iPlane - 1][ictemp]);
            fh_tofd_bars[iPlane - 1]->Fill(hitTofd->GetBarId());
            fh_tofd_time_los_h2[iPlane - 1]->Fill(hitTofd->GetBarId(), hitTofd->GetTof());
            fh_tofd_time_los[hitTofd->GetBarId() - 1][iPlane - 1]->Fill(hitTofd->GetTof());
            if (iPlane == 1)
            {
                fh2_tofd_time_vs_charge->Fill(hitTofd->GetTof(), hitTofd->GetEloss());
                fh2_tofd_time_wouttrig_vs_charge->Fill(fTof_without_trig[hitTofd->GetBarId() - 1], hitTofd->GetEloss());
            }

            if (hitTofd->GetEloss() > charges[iPlane - 1])
                charges[iPlane - 1] = hitTofd->GetEloss() - 0.2;

            iCounts[iPlane - 1] += 1;
            nMulti[iPlane - 1] += 1;
        }
        if (charges[0] > 0 && charges[2] > 0)
            fh2_tofd_charges13->Fill(charges[0], charges[2]);
        if (charges[1] > 0 && charges[3] > 0)
            fh2_tofd_charges24->Fill(charges[1], charges[3]);
        if (charges[0] > 0 && charges[1] > 0)
        {
            fh2_tofd_charges12->Fill(charges[0], charges[1]);
            if (abs(charges[0] - charges[1]) < 0.5)
            {
                auto ZCharge =
                    0.6213293 + 0.65953563445425 * TMath::Sqrt(charges[0] * charges[1]) +
                    0.039440822 * TMath::Sqrt(charges[0] * charges[1]) * TMath::Sqrt(charges[0] * charges[1]);
                fh1_Zcharge->Fill(ZCharge);
                for (const auto& itpat : tpatindex)
                    fh2_zcharge_tpat->Fill(itpat, ZCharge);
            }
        }
        if (charges[1] > 0 && charges[2] > 0)
            fh2_tofd_charges23->Fill(charges[1], charges[2]);
        if (charges[2] > 0 && charges[3] > 0)
            fh2_tofd_charges34->Fill(charges[2], charges[3]);

        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            fh_tofd_multihit_hit[i]->Fill(nMulti[i]);
            if (i > 0)
            {
                for (Int_t im1 = 0; im1 < iCounts[i]; im1++)
                {
                    for (Int_t im2 = 0; im2 < iCounts[i - 1]; im2++)
                    {
                        Double_t tdif = fTimeStitch->GetTime(t[i][im1] - t[i - 1][im2]);
                        fh_tofd_dt_hit[i - 1]->Fill(bar[i][im1], tdif);
                    }
                }
            }
        }
    }

    fNEvents += 1;
}

void R3BTofDOnlineSpectra::FinishEvent()
{
    R3BLOG(debug1, "Cleaning data structures");
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
}

void R3BTofDOnlineSpectra::FinishTask()
{
    if (fCalItems)
    {
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            fh_tofd_TotPm[i]->Write();
            fh_tofd_TotPm_coinc[i]->Write();
            fh_tofd_channels[i]->Write();
            fh_num_side[i]->Write();
            fh_tofd_multihit[i]->Write();
            fh_tofd_multihit_coinc[i]->Write();
            fh2_tofd_ypos_cal[i]->Write();
            fh2_tofd_timedif_cal[i]->Write();
        }
        for (Int_t i = 0; i < fNofPlanes - 1; i++)
        {
            fh_tofd_dt[i]->Write();
        }
    }
    if (fHitItems)
    {
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            fh_tofd_Tot_hit[i]->Write();
            fh_tofd_time_hit[i]->Write();
            fh_tofd_multihit_hit[i]->Write();
            fh_tofd_bars[i]->Write();
            fh_tofd_time_hit[i]->Write();
            fh_tofd_time_los_h2[i]->Write();
            fh2_tofd_time_los_cal[i]->Write();
        }
        for (Int_t i = 0; i < fNofPlanes - 1; i++)
        {
            fh_tofd_dt_hit[i]->Write();
        }
        fh2_tofd_time_vs_charge->Write();
        fh2_tofd_time_wouttrig_vs_charge->Write();
        fh2_tofd_charges13->Write();
        fh2_tofd_charges24->Write();
        fh2_tofd_charges12->Write();
        fh2_tofd_charges23->Write();
        fh2_tofd_charges34->Write();
        fh1_Zcharge->Write();
        fh2_zcharge_tpat->Write();
    }
}

ClassImp(R3BTofDOnlineSpectra)
