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
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BOnlineSpectraToFD_S494.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BLosCalData.h"

#include "R3BPaddleCalData.h"
#include "R3BPaddleTamexMappedData.h"

#include "R3BEventHeader.h"
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
#include "TMath.h"
#include "mapping_tofd_trig.hh"
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

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 15; // nanoseconds.
} // namespace

R3BOnlineSpectraToFD_S494::R3BOnlineSpectraToFD_S494()
    : R3BOnlineSpectraToFD_S494("OnlineSpectraToFD_S494", 1)
{
}

R3BOnlineSpectraToFD_S494::R3BOnlineSpectraToFD_S494(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fCalTriggerItems(NULL)
    , fNofPlanes(N_PLANE_MAX_TOFD_S494)
    , fPaddlesPerPlane(N_PADDLE_MAX_TOFD_S494)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvents(0)
    , fLosTriggerCalDataItems(NULL)
    , fLosCalDataItems(NULL)

{
}

R3BOnlineSpectraToFD_S494::~R3BOnlineSpectraToFD_S494() {}

InitStatus R3BOnlineSpectraToFD_S494::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BOnlineSpectraToFD_S494::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
    {
        LOG(warn) << "R3BOnlineSpectraToFD_S494::Init() EventHeader. not found";
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    }
    else
        LOG(info) << "R3BOnlineSpectraToFD_S494::Init() EventHeader. found";

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("/Tasks", this);

    fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdTriggerCal"));
    if (NULL == fCalTriggerItems)
        printf("Branch TofdTriggerCal not found.\n");

    fLosTriggerCalDataItems = dynamic_cast<TClonesArray*>(mgr->GetObject("LosTriggerCal"));
    if (!fLosTriggerCalDataItems)
    {
        // R3BLOG(fatal, "LOS Data not found");
        // return kFATAL;
        LOG(fatal) << "los trigger not found";
    }

    fLosCalDataItems = dynamic_cast<TClonesArray*>(mgr->GetObject("LosCal"));
    if (!fLosCalDataItems)
    {
        // R3BLOG(fatal, "LOS Hit Data not found");
        // return kFATAL;
        LOG(fatal) << "los data not found";
    }

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
        fHitItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sHit", fDetectorNames[det]))));
    }

    tofd_trig_map_setup();

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------
    // TofD detector

    if (fMappedItems.at(DET_TOFD))
    {
        TCanvas* cTofd_planes = new TCanvas("TOFD_planes_CalLevel", "TOFD planes CAL Level", 10, 10, 1100, 1000);
        cTofd_planes->Divide(6, 4);

        for (Int_t j = 0; j < 4; j++)
        {
            char strName1[255];
            sprintf(strName1, "tofd_channels_plane_%d", j + 1);
            char strName2[255];
            sprintf(strName2, "Tofd channels plane %d", j + 1);
            fh_tofd_channels[j] = new TH1F(strName1, strName2, 90, -45., 45.);
            fh_tofd_channels[j]->GetXaxis()->SetTitle("Channel");
            fh_tofd_channels[j]->GetYaxis()->SetTitle("Counts");

            char strName3[255];
            sprintf(strName3, "tofd_ToT_plane_%d", j + 1);
            char strName4[255];
            sprintf(strName4, "Tofd ToT plane %d", j + 1);
            fh_tofd_TotPm[j] = new TH2F(strName3, strName4, 90, -45, 45, 900, 0., 300.);
            fh_tofd_TotPm[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_TotPm[j]->GetYaxis()->SetTitle("ToT / ns");

            char strName5[255];
            sprintf(strName5, "tofd_ToT_coinc_plane_%d", j + 1);
            char strName6[255];
            sprintf(strName6, "Tofd ToT coinc plane %d", j + 1);
            fh_tofd_TotPm_coinc[j] = new TH2F(strName5, strName6, 90, -45, 45, 900, 0., 300.);
            fh_tofd_TotPm_coinc[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_TotPm_coinc[j]->GetYaxis()->SetTitle("ToT / ns");

            char strName7[255];
            sprintf(strName7, "tofd_multihit_plane_%d", j + 1);
            char strName8[255];
            sprintf(strName8, "Tofd multihit plane %d", j + 1);
            fh_tofd_multihit[j] = new TH2F(strName7, strName8, 90, -45., 45., 30, 0, 30);
            fh_tofd_multihit[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_multihit[j]->GetYaxis()->SetTitle("Multihit");

            char strName9[255];
            sprintf(strName9, "tofd_multihit_coinc_plane_%d", j + 1);
            char strName10[255];
            sprintf(strName10, "Tofd multihit coinc plane %d", j + 1);
            fh_tofd_multihit_coinc[j] = new TH2F(strName9, strName10, 45, 0., 45., 30, 0, 30);
            fh_tofd_multihit_coinc[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_multihit_coinc[j]->GetYaxis()->SetTitle("Multihit");

            if (j < 3)
            {
                int jk = j + 1;
                char strName11[255];
                sprintf(strName11, "tofd_ToF_plane_%d_%d", jk, jk + 1);
                char strName12[255];
                sprintf(strName12, "Tofd ToF plane %d and %d ", jk, jk + 1);
                fh_tofd_dt[j] = new TH2F(strName11, strName12, 50, 0, 50, 4000, -2000., 2000);
                fh_tofd_dt[j]->GetXaxis()->SetTitle("Bar number");
                fh_tofd_dt[j]->GetYaxis()->SetTitle("dt / ns");
            }

            char strName13[255];
            sprintf(strName13, "tofd_numHits_top_vs_bottom_%d", j + 1);
            char strName14[255];
            sprintf(strName14, "Tofd numHitsMapped top vs bottom %d", j + 1);
            fh_num_side[j] = new TH2F(strName13, strName14, 45, 0., 45., 45, 0, 45);
            fh_num_side[j]->GetXaxis()->SetTitle("Num hits up");
            fh_num_side[j]->GetYaxis()->SetTitle("Num hits bottom");
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

        cTofd_planes->cd(0);
        run->AddObject(cTofd_planes);

        run->GetHttpServer()->RegisterCommand("Reset_TOFD", Form("/Tasks/%s/->Reset_TOFD_Histo()", GetName()));
    }
    if (fHitItems.at(DET_TOFD))
    {
        TCanvas* cTofd_planes_hit = new TCanvas("TOFD_planes_HitLevel", "TOFD planes HIT Level", 20, 20, 1120, 1020);
        cTofd_planes_hit->Divide(5, 4);

        for (Int_t j = 0; j < 4; j++)
        {

            char strName3[255];
            sprintf(strName3, "tofd_hit_Q_plane_%d", j + 1);
            char strName4[255];
            sprintf(strName4, "Tofd hit Charge plane %d", j + 1);

            fh_tofd_Tot_hit[j] = new TH2F(strName3, strName4, 45, 0, 45, 200, 0., 20.);
            fh_tofd_Tot_hit[j]->GetXaxis()->SetTitle("BarId");
            fh_tofd_Tot_hit[j]->GetYaxis()->SetTitle("Charge");

            char strName7[255];
            sprintf(strName7, "tofd_hit_multihit_plane_%d", j + 1);
            char strName8[255];
            sprintf(strName8, "Tofd hit multihit plane %d", j + 1);
            fh_tofd_multihit_hit[j] = new TH1F(strName7, strName8, 100, 0, 100);
            fh_tofd_multihit_hit[j]->GetXaxis()->SetTitle("Multihit");

            char strName21[255];
            sprintf(strName21, "tofd_bars_plane_%d", j + 1);
            char strName22[255];
            sprintf(strName22, "Tofd bars plane %d", j + 1);
            fh_tofd_bars[j] = new TH1F(strName21, strName22, 45, 0., 45.);
            fh_tofd_bars[j]->GetXaxis()->SetTitle("BarId");
            fh_tofd_bars[j]->GetYaxis()->SetTitle("Counts");

            char strName23[255];
            sprintf(strName23, "tofd_hit_time_plane_%d", j + 1);
            char strName24[255];
            sprintf(strName24, "Tofd hit time plane %d", j + 1);
            fh_tofd_time_hit[j] = new TH2F(strName23, strName24, 45, 0, 45, 4000, -2000., 2000.);
            fh_tofd_time_hit[j]->GetXaxis()->SetTitle("BarId");
            fh_tofd_time_hit[j]->GetYaxis()->SetTitle("time / ns");

            if (j < 3)
            {
                int jk = j + 1;
                char strName11[255];
                sprintf(strName11, "tofd_hit_ToF_dt_plane_%d_%d", jk, jk + 1);
                char strName12[255];
                sprintf(strName12, "Tofd hit ToF dt plane %d and %d ", jk, jk + 1);
                fh_tofd_dt_hit[j] = new TH2F(strName11, strName12, 45, 0, 45, 4000, -2000., 2000);
                char strName9[255];
                sprintf(strName9, "BarId plane %d ", j + 1);
                fh_tofd_dt_hit[j]->GetXaxis()->SetTitle(strName9);
                fh_tofd_dt_hit[j]->GetYaxis()->SetTitle("dt / ns");
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

        cTofd_planes_hit->cd(0);
        run->AddObject(cTofd_planes_hit);

        run->GetHttpServer()->RegisterCommand("Reset_TOFD", Form("/Tasks/%s/->Reset_TOFD_Histo()", GetName()));
    }

    // -------------------------------------------------------------------------

    TCanvas* cToFd_los[N_PLANE_MAX_TOFD_S494];
    for (Int_t i = 0; i < 4; i++)
    {
        cToFd_los[i] = new TCanvas("TOFD_correlations_LOS", "TOFD correlations with LOS", 20, 20, 1120, 1020);
        cToFd_los[i]->Divide(5, 9);
        for (Int_t j = 0; j < 44; j++)
        {
            char strNameLos[255];
            sprintf(strNameLos, "tofd_hit_timeLos_bar_%d", j + 1);
            char strNameLos2[255];
            sprintf(strNameLos2, "Tofd_time - Los_time bar %d plane %d", j + 1, i + 1);
            fh_tofd_time_los[i][j] = new TH1F(strNameLos, strNameLos2, 50, 0, 1000);
            fh_tofd_time_los[i][j]->GetXaxis()->SetTitle("Time");
            fh_tofd_time_los[i][j]->GetYaxis()->SetTitle("counts");
            cToFd_los[i]->cd(j);
            fh_tofd_time_los[i][j]->Draw("");
        }
    }
    return kSUCCESS;
}

void R3BOnlineSpectraToFD_S494::Reset_TOFD_Histo()
{
    for (int i = 0; i < N_PLANE_MAX_TOFD_S494; i++)
    {
        fh_tofd_channels[i]->Reset();
        fh_tofd_multihit[i]->Reset();
        fh_tofd_TotPm[i]->Reset();
        fh_tofd_multihit_coinc[i]->Reset();
        fh_tofd_TotPm_coinc[i]->Reset();
    }
    fh_tofd_dt[0]->Reset();
    fh_tofd_dt[1]->Reset();
    fh_tofd_dt[2]->Reset();

    if (fHitItems.at(DET_TOFD))
    {
        for (int i = 0; i < N_PLANE_MAX_TOFD_S494; i++)
        {
            fh_tofd_Tot_hit[i]->Reset();
            fh_tofd_time_hit[i]->Reset();
            fh_tofd_multihit_hit[i]->Reset();
            fh_tofd_bars[i]->Reset();

            for (int j = 0; i < N_PADDLE_MAX_TOFD_S494; j++)
            {
                fh_tofd_time_los[i][j]->Reset();
            }
        }
        for (int i = 0; i < N_PLANE_MAX_TOFD_S494 - 1; i++)
        {
            fh_tofd_dt_hit[i]->Reset();
        }
    }
}

namespace
{
    uint64_t n1, n2;
};

void R3BOnlineSpectraToFD_S494::Exec(Option_t* option)
{
    //  cout << "fNEvents " << fNEvents << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        // FairLogger::GetLogger()->Fatal(MESSAGE_ORIGIN, "FairRootManager not found");
        LOG(error) << "FairRootManager not found";
        return;
    }

    time = header->GetTimeStamp();

    if (time_start == 0 && time > 0)
    {
        time_start = time;
        fNEvents_start = fNEvents;
    }
    Bool_t spectra_clear = false;
    Double_t xtime = double(time - time_start) / 1.e9;
    // for reseting spectra
    Int_t icount = iclear_count * reset_time; // reset after reset_time (sec)
    if (time_clear < 0. && int(xtime) % icount == 0 && xtime > 1.)
    {
        time_clear = xtime;
        spectra_clear = true;
    }

    if (header->GetTrigger() == 12)
        time_spill_start = header->GetTimeStamp(); // spill start in nsec
    if (header->GetTrigger() == 13)
        time_spill_end = header->GetTimeStamp(); // spill end  in nsec

    //   if (header->GetTrigger() == 12)
    //       cout << "Spill start: " << double(time_spill_start - time_start) / 1.e9 << " sec" << endl;
    //   if (header->GetTrigger() == 13)
    //       cout << "Spill stop: " << double(time_spill_end - time_start) / 1.e9 << " sec" << endl;

    //   check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit1 = fTpat1 - 1;
    Int_t fTpat_bit2 = fTpat2 - 1;
    Int_t tpatbin;
    for (int i = 0; i < 16; i++)
    {
        tpatbin = (header->GetTpat() & (1 << i));
        if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2))
        {
            return;
        }
    }
    /*
        // fTpat = 1-16; fTpat_bit = 0-15
        Int_t fTpat_bit = fTpat - 1;
        Int_t itpat;
        Int_t tpatvalue;
        if (fTpat_bit >= 0)
        {
            itpat = header->GetTpat();
            tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
            if ((tpatvalue == 0))
                return;
        }
    */
    //----------------------------------------------------------------------
    // TOFD
    //----------------------------------------------------------------------

    Int_t NumPaddles[4] = { 0 };
    if (fMappedItems.at(DET_TOFD))
    {
        auto det = fMappedItems.at(DET_TOFD);
        Int_t nMapped = det->GetEntriesFast();
        Int_t iPlaneMem = 1, iBarMem = 0;
        Int_t nsum_top[N_PLANE_MAX_TOFD_S494] = { 0 };
        Int_t nsum_bot[N_PLANE_MAX_TOFD_S494] = { 0 };
        for (Int_t imapped = 0; imapped < nMapped; imapped++)
        {
            auto mapped = dynamic_cast<R3BTofdMappedData const*>(det->At(imapped));
            if (!mapped)
                continue; // should not happen

            Int_t const iPlane = mapped->GetDetectorId(); // 1..n
            Int_t const iBar = mapped->GetBarId();        // 1..n
            Int_t const iSide = mapped->GetSideId();      // 1..n
            Int_t const iEdge = mapped->GetEdgeId();

            if (iSide == 1 && iEdge == 1 && iPlane < N_PLANE_MAX_TOFD_S494)
                nsum_bot[iPlane - 1] += 1;
            if (iSide == 2 && iEdge == 1 && iPlane < N_PLANE_MAX_TOFD_S494)
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

            //  cout<<imapped<<", "<<iPlane<<"; "<<iPlaneMem<<"; "<<iBar<<"; "<<iBarMem<<",
            //  "<<NumPaddles[iPlane-1]<<endl;

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
        for (Int_t i = 0; i < N_PLANE_MAX_TOFD_S494; i++)
        {
            fh_num_side[i]->Fill(nsum_bot[i], nsum_top[i]);
        }
    }

    Float_t losTime = 0.0;
    Float_t losTriggerTime = 0.0;

    R3BLosCalData* losHit = dynamic_cast<R3BLosCalData*>(fLosCalDataItems->At(0));
    Int_t losChannel = losHit->GetDetector();
    // std::cout<<"LOS Time : "<<losHit->GetTimeT_ns(losChannel)<<std::endl;
    losTime = losHit->GetTimeT_ns(losChannel);

    R3BLosCalData* losTriggerHit = dynamic_cast<R3BLosCalData*>(fLosTriggerCalDataItems->At(0));
    Int_t losChannelTrigger = losTriggerHit->GetDetector();
    // std::cout<<"LOS Time (Trigger) :
    // "<<losTriggerHit->Ge(losHit->GetTime()-losCalTriggerHits->GetTimeL_ns(channelLos)tTimeL_ns(0)<<"
    // "<<losChannelTrigger<<std::endl;
    losTriggerTime = losTriggerHit->GetTimeL_ns(0);

    if (fCalItems.at(DET_TOFD))
    {

        UInt_t vmultihits[N_PLANE_MAX_TOFD_S494 + 1][N_PADDLE_MAX_TOFD_S494];
        UInt_t vmultihits_top[N_PLANE_MAX_TOFD_S494 + 1][N_PADDLE_MAX_TOFD_S494];
        UInt_t vmultihits_bot[N_PLANE_MAX_TOFD_S494 + 1][N_PADDLE_MAX_TOFD_S494];
        Double_t time_bar[N_PLANE_MAX_TOFD_S494 + 1][N_PADDLE_MAX_TOFD_S494][100];
        for (Int_t i = 0; i < fNofPlanes + 1; i++)
        {
            for (Int_t j = 0; j < N_PADDLE_MAX_TOFD_S494; j++)
            {
                vmultihits[i][j] = 0;
                vmultihits_top[i][j] = 0;
                vmultihits_bot[i][j] = 0;
                for (Int_t l = 0; l < 100; l++)
                {
                    time_bar[i][j][l] = 0. / 0.;
                }
            }
        }

        //    std::cout<<"new event!*************************************\n";
        auto det = fCalItems.at(DET_TOFD);
        Int_t nHits = det->GetEntries();

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
            auto* hit = dynamic_cast<R3BTofdCalData*>(det->At(ihit));
            size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();

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
                auto top_trig_i = g_tofd_trig_map[top->GetDetectorId() - 1][top->GetSideId() - 1][top->GetBarId() - 1];

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
                        LOG(error) << "R3BOnlineSpectraToFD::Exec() : Missing trigger information!";
                        s_was_trig_missing = true;
                    }
                }

                Int_t iPlane = top->GetDetectorId(); // 1..n
                Int_t iBar = top->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    LOG(error) << "R3BTOnlineSpectraToFD::Exec() : more detectors than expected! Det: " << iPlane
                               << " allowed are 1.." << fNofPlanes;
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    LOG(error) << "R3BTOnlineSpectraToFD::Exec() : more bars then expected! Det: " << iBar
                               << " allowed are 1.." << fPaddlesPerPlane;
                    continue;
                }

                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                fh_tofd_TotPm[iPlane - 1]->Fill(iBar, top_tot);
                vmultihits_top[iPlane - 1][iBar - 1] += 1;

                ++top_i;
            }

            for (; bot_i < bot_vec.size();)
            {
                auto bot = bot_vec.at(bot_i);
                auto bot_trig_i = g_tofd_trig_map[bot->GetDetectorId() - 1][bot->GetSideId() - 1][bot->GetBarId() - 1];

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
                        LOG(error) << "R3BOnlineSpectraToFD::Exec() : Missing trigger information!";
                        s_was_trig_missing = true;
                    }
                }

                // Shift the cyclic difference window by half a window-length and move it back,
                // this way the trigger time will be at 0.
                auto bot_ns = fmod(bot->GetTimeLeading_ns() - bot_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                              c_range_ns / 2;

                Int_t iPlane = bot->GetDetectorId(); // 1..n
                Int_t iBar = bot->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    LOG(error) << "R3BTOnlineSpectraToFD::Exec() : more detectors than expected! Det: " << iPlane
                               << " allowed are 1.." << fNofPlanes;
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    LOG(error) << "R3BTOnlineSpectraToFD::Exec() : more bars then expected! Det: " << iBar
                               << " allowed are 1.." << fPaddlesPerPlane;
                    continue;
                }

                auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                fh_tofd_TotPm[iPlane - 1]->Fill(-iBar - 1, bot_tot);

                // register multi hits
                vmultihits_bot[iPlane - 1][iBar - 1] += 1;

                ++bot_i;
            }
        }

        for (Int_t ipl = 0; ipl < N_PLANE_MAX_TOFD_S494; ipl++)
        {
            for (Int_t ibr = 1; ibr < N_PADDLE_MAX_TOFD_S494 + 1; ibr++)
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
                auto topc_trig_i =
                    g_tofd_trig_map[topc->GetDetectorId() - 1][topc->GetSideId() - 1][topc->GetBarId() - 1];
                auto botc_trig_i =
                    g_tofd_trig_map[botc->GetDetectorId() - 1][botc->GetSideId() - 1][botc->GetBarId() - 1];

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
                        LOG(error) << "R3BOnlineSpectraToFD::Exec() : Missing trigger information!";
                        s_was_trig_missingc = true;
                    }
                    ++n2;
                }

                // Shift the cyclic difference window by half a window-length and move it back,
                // this way the trigger time will be at 0.
                auto topc_ns =
                    fmod(topc->GetTimeLeading_ns() - topc_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                    c_range_ns / 2;
                auto botc_ns =
                    fmod(botc->GetTimeLeading_ns() - botc_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                    c_range_ns / 2;

                auto dt = topc_ns - botc_ns;
                auto mean_time = (topc_ns + botc_ns) / 2;

                // Handle wrap-around.
                auto dt_mod = fmod(dt + c_range_ns, c_range_ns);
                if (dt < 0)
                {
                    // We're only interested in the short time-differences, so we
                    // want to move the upper part of the coarse counter range close
                    // to the lower range, i.e. we cut the middle of the range and
                    // glue zero and the largest values together.
                    dt_mod -= c_range_ns;
                }

                if (std::abs(dt_mod) < c_bar_coincidence_ns)
                {
                    // Hit!
                    // std::cout << "Hit!\n";
                    Int_t iPlane = topc->GetDetectorId(); // 1..n
                    Int_t iBar = topc->GetBarId();        // 1..n
                    if (iPlane > fNofPlanes)              // this also errors for iDetector==0
                    {
                        LOG(error) << "R3BTOnlineSpectraToFD::Exec() : more detectors than expected! Det: " << iPlane
                                   << " allowed are 1.." << fNofPlanes;
                        continue;
                    }
                    if (iBar > fPaddlesPerPlane) // same here
                    {
                        LOG(error) << "R3BTOnlineSpectraToFD::Exec() : more bars then expected! Det: " << iBar
                                   << " allowed are 1.." << fPaddlesPerPlane;
                        continue;
                    }

                    auto topc_tot =
                        fmod(topc->GetTimeTrailing_ns() - topc->GetTimeLeading_ns() + c_range_ns + c_range_ns / 2,
                             c_range_ns) -
                        c_range_ns / 2;
                    auto botc_tot =
                        fmod(botc->GetTimeTrailing_ns() - botc->GetTimeLeading_ns() + c_range_ns + c_range_ns / 2,
                             c_range_ns) -
                        c_range_ns / 2;

                    fh_tofd_TotPm_coinc[iPlane - 1]->Fill(-iBar - 1, botc_tot);
                    fh_tofd_TotPm_coinc[iPlane - 1]->Fill(iBar, topc_tot);
                    if (losTime > 0)
                    {
                        fh_tofd_time_los[iPlane - 1][iBar - 1]->Fill(mean_time - losTime);
                    }
                    // std::cout<<"ToT: "<<top_tot << " "<<bot_tot<<"\n";

                    // register multi hits
                    Int_t imlt = vmultihits[iPlane - 1][iBar - 1];
                    time_bar[iPlane - 1][iBar - 1][imlt] = (topc_ns + botc_ns) / 2.;
                    vmultihits[iPlane - 1][iBar - 1] += 1;

                    ++topc_i;
                    ++botc_i;
                }
                else if (dt < 0 && dt > -c_range_ns / 2)
                {
                    ++topc_i;
                }
                else
                {
                    ++botc_i;
                }
            }
        }

        for (Int_t ipl = 0; ipl < N_PLANE_MAX_TOFD_S494; ipl++)
        {
            for (Int_t ibr = 1; ibr < N_PADDLE_MAX_TOFD_S494 + 1; ibr++)
            {
                fh_tofd_multihit_coinc[ipl]->Fill(ibr, vmultihits[ipl][ibr - 1]);
                if (ipl > 0)
                {
                    for (Int_t imult1 = 0; imult1 < vmultihits[ipl][ibr - 1]; imult1++)
                    {
                        for (Int_t imult2 = 0; imult2 < vmultihits[ipl - 1][ibr - 1]; imult2++)
                        {
                            Double_t tof_plane = 0. / 0.;
                            tof_plane = fmod(time_bar[ipl][ibr - 1][imult1] - time_bar[ipl - 1][ibr - 1][imult2] +
                                                 c_range_ns + c_range_ns / 2,
                                             c_range_ns) -
                                        c_range_ns / 2;
                            fh_tofd_dt[ipl - 1]->Fill(ibr, tof_plane);
                        }
                    }
                }
            }
        }

    } // endi if fCalItems

    if (fHitItems.at(DET_TOFD))
    {
        auto detTofd = fHitItems.at(DET_TOFD);
        Int_t nHits = detTofd->GetEntriesFast();

        if (nHits > 100)
            return;

        Double_t x[N_PLANE_MAX_TOFD_S494][10], y[N_PLANE_MAX_TOFD_S494][10], t[N_PLANE_MAX_TOFD_S494][10],
            q[N_PLANE_MAX_TOFD_S494][10], bar[N_PLANE_MAX_TOFD_S494][10];
        for (Int_t i = 0; i < N_PLANE_MAX_TOFD_S494; i++)
        {
            for (Int_t k = 0; k < 10; k++)
            {
                x[i][k] = -1000.;
                y[i][k] = -1000.;
                q[i][k] = -1000.;
                t[i][k] = 0. / 0.;
                bar[i][k] = 0;
            }
        }
        Int_t nMulti[N_PLANE_MAX_TOFD_S494] = { 0 }, iCounts[N_PLANE_MAX_TOFD_S494] = { 0 };
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BTofdHitData* hitTofd = dynamic_cast<R3BTofdHitData*>(detTofd->At(ihit));

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
            fh_tofd_bars[iPlane - 1]->Fill(bar[iPlane - 1][ictemp]);
            iCounts[iPlane - 1] += 1;
            nMulti[iPlane - 1] += 1;
        }

        for (Int_t i = 0; i < N_PLANE_MAX_TOFD_S494; i++)
        {
            fh_tofd_multihit_hit[i]->Fill(nMulti[i]);
            if (i > 0)
            {
                for (Int_t im1 = 0; im1 < iCounts[i]; im1++)
                {
                    for (Int_t im2 = 0; im2 < iCounts[i - 1]; im2++)
                    {
                        Double_t tdif =
                            fmod(t[i][im1] - t[i - 1][im2] + c_range_ns + c_range_ns / 2, c_range_ns) - c_range_ns / 2;
                        fh_tofd_dt_hit[i - 1]->Fill(bar[i][im1], tdif);
                    }
                }
            }
        }
    }

    fNEvents += 1;
}

void R3BOnlineSpectraToFD_S494::FinishEvent()
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
        if (fHitItems.at(det))
        {
            fHitItems.at(det)->Clear();
        }
        fLosCalDataItems->Clear();
        fLosTriggerCalDataItems->Clear();
    }
}

void R3BOnlineSpectraToFD_S494::FinishTask()
{
    if (fCalItems.at(DET_TOFD))
    {
        for (Int_t i = 0; i < N_PLANE_MAX_TOFD_S494; i++)
        {
            fh_tofd_TotPm[i]->Write();
            fh_tofd_TotPm_coinc[i]->Write();
            fh_tofd_channels[i]->Write();
            fh_num_side[i]->Write();
            fh_tofd_multihit[i]->Write();
            fh_tofd_multihit_coinc[i]->Write();

            for (int j = 0; j < N_PADDLE_MAX_TOFD_S494; j++)
            {
                fh_tofd_time_los[i][j]->Write();
            }
        }
        for (Int_t i = 0; i < N_PLANE_MAX_TOFD_S494 - 1; i++)
        {
            fh_tofd_dt[i]->Write();
        }
    }
    if (fHitItems.at(DET_TOFD))
    {
        for (Int_t i = 0; i < N_PLANE_MAX_TOFD_S494; i++)
        {
            fh_tofd_Tot_hit[i]->Write();
            fh_tofd_time_hit[i]->Write();
            fh_tofd_multihit_hit[i]->Write();
            fh_tofd_bars[i]->Write();
            fh_tofd_time_hit[i]->Write();
        }
        for (Int_t i = 0; i < N_PLANE_MAX_TOFD_S494 - 1; i++)
        {
            fh_tofd_dt_hit[i]->Write();
        }
    }
}

ClassImp(R3BOnlineSpectraToFD_S494)
