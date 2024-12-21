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

#include "R3BOnlineSpectraToFI_S494.h"

#include "R3BTofiCalData.h"
#include "R3BTofiHitData.h"
#include "R3BTofiMappedData.h"

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

#include "mapping_tofi_trig.hh"

extern unsigned g_tofi_trig_map[2][24];
void tofi_trig_map_setup();

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
} // namespace

R3BOnlineSpectraToFI_S494::R3BOnlineSpectraToFI_S494()
    : R3BOnlineSpectraToFI_S494("OnlineSpectraToFI_S494", 1)
{
}

R3BOnlineSpectraToFI_S494::R3BOnlineSpectraToFI_S494(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fCalTriggerItems(NULL)
    , fNofPlanes(N_PLANE_MAX_TOFI)
    , fPaddlesPerPlane(N_PADDLE_MAX_TOFI)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvents(0)
{
}

R3BOnlineSpectraToFI_S494::~R3BOnlineSpectraToFI_S494() {}

InitStatus R3BOnlineSpectraToFI_S494::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BOnlineSpectraToFI_S494::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("/Tasks", this);

    fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofiTriggerCal"));
    if (NULL == fCalTriggerItems)
        printf("Branch TofiTriggerCal not found.\n");

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
        if (NULL == fCalItems.at(det))
        {
            printf("Could not find cal data for '%s'.\n", fDetectorNames[det]);
        }
        fHitItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sHit", fDetectorNames[det]))));
        if (NULL == fHitItems.at(det))
        {
            printf("Could not find hit data for '%s'.\n", fDetectorNames[det]);
        }
    }

    tofi_trig_map_setup();

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------
    // TofI detector

    if (fMappedItems.at(DET_TOFI))
    {
        TCanvas* ctofi_planes = new TCanvas("TOFI_CalLevel_planes", "TOFI CalLevel planes", 10, 10, 1100, 1000);

        ctofi_planes->Divide(4, 3);

        for (Int_t j = 0; j < N_PLANE_MAX_TOFI; j++)
        {
            char strName1[255];
            sprintf(strName1, "tofi_channels_plane_%d", j + 1);
            char strName2[255];
            sprintf(strName2, "Tofi Mapped channels plane %d", j + 1);
            fh_tofi_channels[j] = new TH1F(strName1, strName2, 60, -30, 30);
            fh_tofi_channels[j]->GetXaxis()->SetTitle("Channel");
            fh_tofi_channels[j]->GetYaxis()->SetTitle("Counts");

            char strName3[255];
            sprintf(strName3, "tofi_ToT_plane_%d", j + 1);
            char strName4[255];
            sprintf(strName4, "Tofi ToT plane %d", j + 1);
            fh_tofi_TotPm[j] = new TH2F(strName3, strName4, 60, -30, 30, 1500, 0., 300.);
            fh_tofi_TotPm[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofi_TotPm[j]->GetYaxis()->SetTitle("ToT / ns");

            char strName30[255];
            sprintf(strName30, "tofi_time_plane_%d", j + 1);
            char strName40[255];
            sprintf(strName40, "Tofi time plane %d", j + 1);
            fh_tofi_timePm[j] = new TH2F(strName30, strName40, 60, -30, 30, 3000, -2000., 1000);
            fh_tofi_timePm[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofi_timePm[j]->GetYaxis()->SetTitle("time / ns");

            char strName9[255];
            sprintf(strName9, "tofi_multihit_plane_%d", j + 1);
            char strName10[255];
            sprintf(strName10, "Tofi multihit plane %d", j + 1);
            fh_tofi_multihit[j] = new TH2F(strName9, strName10, 60, -30., 30., 10, 0, 10);
            fh_tofi_multihit[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofi_multihit[j]->GetYaxis()->SetTitle("Multihit");

            fh_tofi_TotPm_coinc[j] = new TH2F("TOFI_ToT_with_coinc", "ToFI ToT with coinc", 30, 0, 30, 6000, 0., 300.);
            fh_tofi_TotPm_coinc[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofi_TotPm_coinc[j]->GetYaxis()->SetTitle("sqrt(ToT1*ToT2) / ns");

            fh_tofi_multihit_coinc[j] =
                new TH2F("ToFI_multihit_with_coinc", "ToFI multihit with coinc", 30, 0., 30., 10, 0, 10);
            fh_tofi_multihit_coinc[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofi_multihit_coinc[j]->GetYaxis()->SetTitle("Multihit");
        }

        fh_tofi_dt = new TH2F("tofi_dt", "Tofi dt", 25, 0, 25, 2000, -100., 100);
        fh_tofi_dt->GetXaxis()->SetTitle("Bar number");
        fh_tofi_dt->GetYaxis()->SetTitle("timeTop-timeBottom / ns");

        fh_tofi_time = new TH2F("tofi_time", "Tofi not-sync time", 25, 0, 25, 4000, -2000., 2000);
        fh_tofi_time->GetXaxis()->SetTitle("Bar number");
        fh_tofi_time->GetYaxis()->SetTitle("time / ns");

        fh_num_bars = new TH1F("tofi_numBars", "Tofi Mapped num bars", 30, 0, 30);
        fh_num_bars->GetXaxis()->SetTitle("Num of bars with hit");
        fh_num_bars->GetYaxis()->SetTitle("Counts");

        fh_num_side = new TH2F("tofi_top_vs_bottom", "Tofi Mapped top vs bottom", 30, 0, 30, 30, 0, 30);
        fh_num_side->GetXaxis()->SetTitle("Num of hits with bottom");
        fh_num_side->GetYaxis()->SetTitle("Num of hits with top");

        fh_test = new TH2F("tofi_ToT_numBars_eq_1", "Tofi ToT for Multipl==1", 60, -30., 30., 600, 0., 300.);
        fh_test->GetXaxis()->SetTitle("Bar number");
        fh_test->GetYaxis()->SetTitle("ToT  / ns");

        fh_test1 = new TH2F("tofi_ToT_vs_ToT1", "Tofi ToT corel iBar vs iBar+1", 600, 0., 300., 600, 0., 300.);
        fh_test1->GetXaxis()->SetTitle("ToT iBar=5 / ns");
        fh_test1->GetYaxis()->SetTitle("ToT  iBar=4,6 / ns");

        fh_test2 = new TH2F("tofi_ToT_vs_ToT2", "Tofi ToT corel iBar vs iBar+2", 600, 0., 300., 600, 0., 300.);
        fh_test2->GetXaxis()->SetTitle("ToT iBar=5 / ns");
        fh_test2->GetYaxis()->SetTitle("ToT  iBar=3,7 / ns");

        fh_test3 = new TH2F("tofi_iBar_vs_iBar", "Tofi iBar vs iBar", 30, 0., 30., 30, 0., 30.);
        fh_test3->GetXaxis()->SetTitle("iBar");
        fh_test3->GetYaxis()->SetTitle("iBar");

        ctofi_planes->cd(1);
        fh_tofi_channels[0]->Draw();
        ctofi_planes->cd(2);
        gPad->SetLogy();
        fh_num_bars->Draw();
        ctofi_planes->cd(3);
        gPad->SetLogz();
        fh_num_side->Draw("colz");

        ctofi_planes->cd(6);
        gPad->SetLogz();
        fh_tofi_TotPm[0]->Draw("colz");
        ctofi_planes->cd(5);
        gPad->SetLogz();
        fh_tofi_multihit[0]->SetAxisRange(1, 10, "Y");
        fh_tofi_multihit[0]->Draw("colz");
        ctofi_planes->cd(7);
        gPad->SetLogz();
        fh_tofi_timePm[0]->Draw("colz");

        ctofi_planes->cd(10);
        gPad->SetLogz();
        fh_tofi_TotPm_coinc[0]->Draw("colz");
        ctofi_planes->cd(9);
        gPad->SetLogz();
        fh_tofi_multihit_coinc[0]->Draw("colz");
        ctofi_planes->cd(11);
        gPad->SetLogz();
        fh_tofi_dt->SetAxisRange(1, 25, "X");
        fh_tofi_dt->Draw("colz");
        ctofi_planes->cd(12);
        gPad->SetLogz();
        gPad->SetLogz();
        fh_tofi_time->Draw("colz");

        ctofi_planes->cd(0);
        run->AddObject(ctofi_planes);

        run->GetHttpServer()->RegisterCommand("Reset_TOFI", Form("/Tasks/%s/->Reset_TOFI_Histo()", GetName()));
    }
    if (fHitItems.at(DET_TOFI))
    {
        TCanvas* ctofi_planes_hit = new TCanvas("TOFI_planes_HitLevel", "TOFI planes HIT Level", 20, 20, 1120, 1020);
        ctofi_planes_hit->Divide(3, 2);

        for (Int_t j = 0; j < N_PLANE_MAX_TOFI; j++)
        {

            char strName11[255];
            sprintf(strName11, "tofi_hit_ToT_plane_%d", j + 1);
            char strName12[255];
            sprintf(strName12, "tofi hit ToT plane %d", j + 1);
            fh_tofi_Tot_hit[j] = new TH2F(strName11, strName12, 25, 0, 25, 200, 0., 20.);
            fh_tofi_Tot_hit[j]->GetXaxis()->SetTitle("BarId");
            fh_tofi_Tot_hit[j]->GetYaxis()->SetTitle("Nuclear charge");

            char strName13[255];
            sprintf(strName13, "tofi_hit_multihit_plane_%d", j + 1);
            char strName14[255];
            sprintf(strName14, "tofi hit multihit plane %d", j + 1);
            fh_tofi_multihit_hit[j] = new TH1F(strName13, strName14, 100, 0, 100);
            fh_tofi_multihit_hit[j]->GetXaxis()->SetTitle("Multihit");

            char strName21[255];
            sprintf(strName21, "tofi_bars_plane_%d", j + 1);
            char strName22[255];
            sprintf(strName22, "tofi bars plane %d", j + 1);
            fh_tofi_bars[j] = new TH1F(strName21, strName22, 25, 0., 25.);
            fh_tofi_bars[j]->GetXaxis()->SetTitle("BarId");
            fh_tofi_bars[j]->GetYaxis()->SetTitle("Counts");

            char strName23[255];
            sprintf(strName23, "tofi_hit_time_plane_%d", j + 1);
            char strName24[255];
            sprintf(strName24, "tofi hit time plane %d", j + 1);
            fh_tofi_time_hit[j] = new TH2F(strName23, strName24, 25, 0, 25, 4000, -2000., 2000.);
            fh_tofi_time_hit[j]->GetXaxis()->SetTitle("BarId");
            fh_tofi_time_hit[j]->GetYaxis()->SetTitle("time / ns");

            char strName25[255];
            sprintf(strName25, "tofi_hit_y_plane_%d", j + 1);
            char strName26[255];
            sprintf(strName26, "tofi hit y plane %d ", j + 1);
            fh_tofi_dt_hit[j] = new TH2F(strName25, strName26, 25, 0, 25, 4000, -2000., 2000);
            fh_tofi_dt_hit[j]->GetXaxis()->SetTitle("Bar Id");
            fh_tofi_dt_hit[j]->GetYaxis()->SetTitle("dt / ns");

            char strName43[255];
            sprintf(strName43, "tofi_Z_plane_%d", j + 1);
            char strName44[255];
            sprintf(strName44, "tofi Z plane %d", j + 1);
            fh_tofi_Z[j] = new TH1F(strName43, strName44, 500, 0, 20);
            fh_tofi_Z[j]->GetXaxis()->SetTitle("Nuclear charge");
        }

        ctofi_planes_hit->cd(1);
        gPad->SetLogy();
        fh_tofi_multihit_hit[0]->Draw();
        ctofi_planes_hit->cd(4);
        gPad->SetLogy();
        fh_tofi_bars[0]->Draw();
        ctofi_planes_hit->cd(5);
        gPad->SetLogy();
        fh_tofi_Z[0]->Draw();
        ctofi_planes_hit->cd(2);
        gPad->SetLogz();
        fh_tofi_Tot_hit[0]->Draw("colz");
        ctofi_planes_hit->cd(3);
        gPad->SetLogz();
        fh_tofi_dt_hit[0]->Draw("colz"); // fh_tofi_bars[0]->Draw();
        ctofi_planes_hit->cd(6);
        gPad->SetLogz();
        fh_tofi_time_hit[0]->Draw("colz");

        ctofi_planes_hit->cd(0);
        run->AddObject(ctofi_planes_hit);

        run->GetHttpServer()->RegisterCommand("Reset_TOFI", Form("/Tasks/%s/->Reset_TOFI_Histo()", GetName()));
    }

    // -------------------------------------------------------------------------

    return kSUCCESS;
}

void R3BOnlineSpectraToFI_S494::Reset_TOFI_Histo()
{
    for (int i = 0; i < N_PLANE_MAX_TOFI; i++)
    {
        fh_tofi_channels[i]->Reset();
        fh_tofi_multihit[i]->Reset();
        fh_tofi_TotPm[i]->Reset();
        fh_tofi_timePm[i]->Reset();
        fh_tofi_multihit_coinc[i]->Reset();
        fh_tofi_TotPm_coinc[i]->Reset();
    }
    fh_tofi_dt->Reset();
    fh_tofi_time->Reset();
    fh_num_bars->Reset();
    fh_test->Reset();
    fh_test1->Reset();
    fh_test2->Reset();
    fh_test3->Reset();
    fh_num_side->Reset();

    if (fHitItems.at(DET_TOFI))
    {
        for (int i = 0; i < N_PLANE_MAX_TOFI; i++)
        {
            fh_tofi_Tot_hit[i]->Reset();
            fh_tofi_time_hit[i]->Reset();
            fh_tofi_multihit_hit[i]->Reset();
            fh_tofi_bars[i]->Reset();
            fh_tofi_Z[i]->Reset();
        }
        for (int i = 0; i < N_PLANE_MAX_TOFI - 1; i++)
        {
            fh_tofi_dt_hit[i]->Reset();
        }
    }
}

namespace
{
    uint64_t n1, n2;
};

void R3BOnlineSpectraToFI_S494::Exec(Option_t* option)
{
    // cout << "fNEvents " << fNEvents << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        // FairLogger::GetLogger()->Fatal(MESSAGE_ORIGIN, "FairRootManager not found");
        LOG(error) << "FairRootManager not found";
        return;
    }
    Bool_t time_header = false;
    time = header->GetTimeStamp();

    if (time > 0)
    {
        time_header = true;
        time_mem = time;
    }

    time = time_mem;

    // cout<<"TIMES: "<<time<<", "<<time_mem<<endl;

    if (time_start == 0 && time > 0)
    {
        time_start = time;
        fNEvents_start = fNEvents;
    }
    Double_t xtime = double(time - time_start) / 1.e9;

    if (header->GetTrigger() == 12)
        time_spill_start = header->GetTimeStamp(); // spill start in nsec
    if (header->GetTrigger() == 13)
        time_spill_end = header->GetTimeStamp(); // spill end  in nsec

    if (header->GetTrigger() == 12)
        cout << "Spill start: " << double(time_spill_start - time_start) / 1.e9 << " sec " << endl;
    if (header->GetTrigger() == 13)
        cout << "Spill stop: " << double(time_spill_end - time_start) / 1.e9 << " sec " << endl;

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
    // TOFI
    //----------------------------------------------------------------------

    Int_t NumPaddles = 0;

    if (fMappedItems.at(DET_TOFI))
    {
        Int_t nsum_lead = 0;
        Int_t nsum_trail = 0;
        Int_t nsum_top = 0;
        Int_t nsum_bot = 0;
        auto det = fMappedItems.at(DET_TOFI);
        Int_t nMapped = det->GetEntriesFast();
        Int_t iPlaneMem = 1, iBarMem = 0;
        for (Int_t imapped = 0; imapped < nMapped; imapped++)
        {
            auto mapped = dynamic_cast<R3BTofiMappedData const*>(det->At(imapped));
            if (!mapped)
                continue; // should not happen

            Int_t const iPlane = mapped->GetDetectorId(); // 1..n
            Int_t const iBar = mapped->GetBarId();        // 1..n
            Int_t const iSide = mapped->GetSideId();      // 1..n
            Int_t const iEdge = mapped->GetEdgeId();

            if (iEdge == 1)
                nsum_lead += 1;
            if (iEdge == 2)
                nsum_trail += 1;
            if (iSide == 1 && iEdge == 1 && iPlane == 1)
                nsum_bot += 1;
            if (iSide == 2 && iEdge == 1 && iPlane == 1)
                nsum_top += 1;

            if (iBarMem != iBar && iSide == 1 && iEdge == 1)
                NumPaddles += 1;

            iBarMem = iBar;
            if (iPlane <= fNofPlanes)
            {
                if (iSide == 1 && iEdge == 1 && iPlane == 1) // only leading edges iEdge == 1
                    fh_tofi_channels[iPlane - 1]->Fill(-iBar - 1);
                if (iSide == 2 && iEdge == 1 && iPlane == 1)
                    fh_tofi_channels[iPlane - 1]->Fill(iBar);
            }
            Int_t maxentry = fh_tofi_channels[0]->GetBinContent(fh_tofi_channels[0]->GetMaximumBin());
        }
        fh_num_bars->Fill(NumPaddles);
        fh_num_side->Fill(nsum_bot, nsum_top);
    }

    if (fCalItems.at(DET_TOFI)) // without coincidances top/bottom
    {

        UInt_t vmultihits[N_PLANE_MAX_TOFI][N_PADDLE_MAX_TOFI];
        UInt_t vmultihits_top[N_PLANE_MAX_TOFI][N_PADDLE_MAX_TOFI];
        UInt_t vmultihits_bot[N_PLANE_MAX_TOFI][N_PADDLE_MAX_TOFI];
        Double_t time_bar[N_PLANE_MAX_TOFI][N_PADDLE_MAX_TOFI][32];
        Double_t tot_bar[N_PLANE_MAX_TOFI][N_PADDLE_MAX_TOFI][32];
        auto time0_top = -1;
        auto time0_bot = -1;
        auto time_abs = -1;
        auto time1 = -1;
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            for (Int_t j = 0; j < N_PADDLE_MAX_TOFI; j++)
            {
                vmultihits[i][j] = 0;
                vmultihits_top[i][j] = 0;
                vmultihits_bot[i][j] = 0;
                for (Int_t k = 0; k < 32; k++)
                {
                    tot_bar[i][j][k] = 0;
                    time_bar[i][j][k] = 0;
                }
            }
        }

        auto det = fCalItems.at(DET_TOFI);
        Int_t nHits = det->GetEntries();

        Int_t nHitsEvent = 0;

        // Organize cals into bars.
        struct Entry
        {
            std::vector<R3BTofiCalData*> top;
            std::vector<R3BTofiCalData*> bot;
        };
        std::map<size_t, Entry> bar_map;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {

            auto* hit = dynamic_cast<R3BTofiCalData*>(det->At(ihit));
            size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();

            auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
            auto& vec = 1 == hit->GetSideId() ? ret.first->second.bot : ret.first->second.top;
            vec.push_back(hit);

            if (hit->GetDetectorId() == 1 && hit->GetSideId() == 1 && hit->GetTimeLeading_ns() > 0)
                nHitsEvent += 1;
        }

        //     if(NumPaddles > 0)cout <<"nHits cal level: "<<nHitsEvent<<", "<<NumPaddles<<endl;

        // Build trigger map.
        std::vector<R3BTofiCalData*> trig_map;
        for (int i = 0; i < fCalTriggerItems->GetEntries(); ++i)
        {
            auto trig = dynamic_cast<R3BTofiCalData*>(fCalTriggerItems->At(i));
            if (trig_map.size() < trig->GetBarId())
            {
                trig_map.resize(trig->GetBarId());
            }
            trig_map.at(trig->GetBarId() - 1) = trig;
        }

        // ******************* Without coincidences ****************************************************************
        bool s_was_trig_missing = false;

        for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
        {
            auto const& top_vec = it->second.top;
            auto const& bot_vec = it->second.bot;
            size_t top_i = 0;
            size_t bot_i = 0;

            // TOP PMTS
            for (; top_i < top_vec.size();)
            {
                auto top = top_vec.at(top_i);
                auto top_trig_i = g_tofi_trig_map[top->GetSideId() - 1][top->GetBarId() - 1];
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
                        LOG(error) << "R3BOnlineSpectraToFI::Exec() : Missing trigger information!";
                        s_was_trig_missing = true;
                    }
                }

                // Shift the cyclic difference window by half a window-length and move it back,
                // this way the trigger time will be at 0.
                auto top_ns = fmod(top->GetTimeLeading_ns() - top_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                              c_range_ns / 2;

                auto top_raw = fmod(top->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                //      cout<<"Time: "<<top_ns<<"; "<<bot_ns<<endl;
                Int_t iPlane = top->GetDetectorId(); // 1..n
                Int_t iBar = top->GetBarId();        // 1..n

                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                fh_tofi_TotPm[iPlane - 1]->Fill(iBar, top_tot);
                fh_tofi_timePm[iPlane - 1]->Fill(iBar, top_ns);
                if (nHitsEvent == 1)
                    fh_test->Fill(iBar, top_tot);
                // register multi hits
                vmultihits_top[iPlane - 1][iBar - 1] += 1;

                ++top_i;
            }

            time_abs = -1;
            time1 = -1;

            // BOTTOM PMTS
            for (; bot_i < bot_vec.size();)
            {
                auto bot = bot_vec.at(bot_i);
                auto bot_trig_i = g_tofi_trig_map[bot->GetSideId() - 1][bot->GetBarId() - 1];
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
                        LOG(error) << "R3BOnlineSpectraToFI::Exec() : Missing trigger information!";
                        s_was_trig_missing = true;
                    }
                }

                // Shift the cyclic difference window by half a window-length and move it back,
                // this way the trigger time will be at 0.
                auto bot_ns = fmod(bot->GetTimeLeading_ns() - bot_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                              c_range_ns / 2;

                auto bot_raw = fmod(bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                //      cout<<"Time: "<<top_ns<<"; "<<bot_ns<<endl;
                Int_t iPlane = bot->GetDetectorId(); // 1..n
                Int_t iBar = bot->GetBarId();        // 1..n
                // auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                auto bot_tot = bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns();

                fh_tofi_TotPm[iPlane - 1]->Fill(-iBar - 1, bot_tot);
                fh_tofi_timePm[iPlane - 1]->Fill(-iBar - 1, bot_ns);
                if (nHitsEvent == 1)
                    fh_test->Fill(-iBar - 1, bot_tot);
                // register multi hits
                vmultihits_bot[iPlane - 1][iBar - 1] += 1;

                ++bot_i;
            }
        }
        for (Int_t ipl = 0; ipl < N_PLANE_MAX_TOFI; ipl++)
        {
            for (Int_t ibr = 1; ibr < N_PADDLE_MAX_TOFI + 1; ibr++)
            {
                if (vmultihits_bot[ipl][ibr - 1] > 0)
                    fh_tofi_multihit[ipl]->Fill(-ibr - 1, vmultihits_bot[ipl][ibr - 1]);
                if (vmultihits_top[ipl][ibr - 1] > 0)
                    fh_tofi_multihit[ipl]->Fill(ibr, vmultihits_top[ipl][ibr - 1]);
            }
        }

        // With coincidences ****************************************************************
        bool s_was_trig_missingc = false;
        Int_t imlt = 0;
        for (auto itc = bar_map.begin(); bar_map.end() != itc; ++itc)
        {
            auto const& topc_vec = itc->second.top;
            auto const& botc_vec = itc->second.bot;
            size_t topc_i = 0;
            size_t botc_i = 0;
            for (; topc_i < topc_vec.size() && botc_i < botc_vec.size();)
            {
                auto topc = topc_vec.at(topc_i);
                auto botc = botc_vec.at(botc_i);
                auto topc_trig_i = g_tofi_trig_map[topc->GetSideId() - 1][topc->GetBarId() - 1];
                auto botc_trig_i = g_tofi_trig_map[botc->GetSideId() - 1][botc->GetBarId() - 1];
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
                        LOG(error) << "R3BOnlineSpectraToFI::Exec() : Missing trigger information!";
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

                //      cout<<"Time: "<<top_ns<<"; "<<bot_ns<<endl;

                auto dt = topc_ns - botc_ns;
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
                fh_tofi_dt->Fill(topc->GetBarId(), dt_mod);

                //       cout<<"dt_mod: "<<dt_mod <<endl;

                // std::cout << top_i << ' ' << bot_i << ": " << top_ns << ' ' << bot_ns << " = " << dt << ' ' <<
                // std::abs(dt_mod) << '\n';
                if (std::abs(dt_mod) < c_bar_coincidence_ns)
                {
                    // Hit!
                    // std::cout << "Hit!\n";
                    Int_t iPlane = topc->GetDetectorId(); // 1..n
                    Int_t iBar = topc->GetBarId();        // 1..n
                    if (iPlane > fNofPlanes)              // this also errors for iDetector==0
                    {
                        LOG(error) << "R3BTOnlineSpectraToFI::Exec() : more detectors than expected! Det: " << iPlane
                                   << " allowed are 1.." << fNofPlanes;
                        continue;
                    }
                    if (iBar > fPaddlesPerPlane) // same here
                    {
                        LOG(error) << "R3BTOnlineSpectraToFI::Exec() : more bars then expected! Det: " << iBar
                                   << " allowed are 1.." << fPaddlesPerPlane;
                        continue;
                    }

                    auto topc_tot =
                        fmod(topc->GetTimeTrailing_ns() - topc->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                    auto botc_tot =
                        fmod(botc->GetTimeTrailing_ns() - botc->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                    // register multi hits
                    vmultihits[iPlane - 1][iBar - 1] += 1;
                    imlt = vmultihits[iPlane - 1][iBar - 1];

                    time_bar[iPlane - 1][iBar - 1][imlt - 1] = (topc_ns + botc_ns) / 2.;
                    fh_tofi_time->Fill(iBar, time_bar[iPlane - 1][iBar - 1][imlt - 1]);

                    tot_bar[iPlane - 1][iBar - 1][imlt - 1] = sqrt(topc_tot * botc_tot);
                    fh_tofi_TotPm_coinc[iPlane - 1]->Fill(iBar, tot_bar[iPlane - 1][iBar - 1][imlt - 1]);

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

        for (Int_t ipl = 0; ipl < N_PLANE_MAX_TOFI; ipl++)
        {
            for (Int_t ibr = 0; ibr < N_PADDLE_MAX_TOFI; ibr++)
            {
                if (vmultihits[ipl][ibr] > 0)
                    fh_tofi_multihit_coinc[ipl]->Fill(ibr + 1, vmultihits[ipl][ibr]);

                for (Int_t imult = 0; imult < vmultihits[ipl][ibr]; imult++)
                {
                    if (ibr == 4 && tot_bar[ipl][ibr][imult] > 0) // Paddle 5
                    {
                        if (ibr < N_PADDLE_MAX_TOFI - 1)
                        {
                            for (Int_t imult1 = 0; imult1 < vmultihits[ipl][ibr + 1]; imult1++)
                            {
                                if (tot_bar[ipl][ibr + 1][imult1] > 0)
                                    fh_test1->Fill(tot_bar[ipl][ibr][imult], tot_bar[ipl][ibr + 1][imult1]);
                            }
                        }

                        if (ibr < N_PADDLE_MAX_TOFI - 2)
                        {
                            for (Int_t imult1 = 0; imult1 < vmultihits[ipl][ibr + 2]; imult1++)
                            {
                                if (tot_bar[ipl][ibr + 2][imult1] > 0)
                                    fh_test2->Fill(tot_bar[ipl][ibr][imult1], tot_bar[ipl][ibr + 2][imult1]);
                            }
                        }

                        if (ibr > 0)
                        {
                            for (Int_t imult1 = 0; imult1 < vmultihits[ipl][ibr - 1]; imult1++)
                            {
                                if (tot_bar[ipl][ibr - 1][imult1] > 0)
                                    fh_test1->Fill(tot_bar[ipl][ibr][imult1], tot_bar[ipl][ibr - 1][imult1]);
                            }
                        }

                        if (ibr > 1)
                        {
                            for (Int_t imult1 = 0; imult1 < vmultihits[ipl][ibr - 2]; imult1++)
                            {
                                if (tot_bar[ipl][ibr - 2][imult1] > 0)
                                    fh_test2->Fill(tot_bar[ipl][ibr][imult1], tot_bar[ipl][ibr - 2][imult1]);
                            }
                        }
                    }
                    if (tot_bar[ipl][ibr][imult] > 0)
                    {
                        for (Int_t ibr1 = 0; ibr1 < N_PADDLE_MAX_TOFI; ibr1++)
                        {
                            for (Int_t imult1 = 0; imult1 < vmultihits[ipl][ibr1]; imult1++)
                            {
                                if (tot_bar[ipl][ibr1][imult1] > 0)
                                    fh_test3->Fill(ibr + 1, ibr1 + 1);
                            }
                        }
                    }
                }
            }
        }

    } // endi if fCalItems

    if (fHitItems.at(DET_TOFI))
    {
        auto detTofi = fHitItems.at(DET_TOFI);
        Int_t nHits = detTofi->GetEntriesFast();

        if (nHits > 100)
            return;

        Double_t x[N_PLANE_MAX_TOFI][10], y[N_PLANE_MAX_TOFI][10], t[N_PLANE_MAX_TOFI][10], q[N_PLANE_MAX_TOFI][10],
            bar[N_PLANE_MAX_TOFI][10];
        for (Int_t i = 0; i < N_PLANE_MAX_TOFI; i++)
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
        Int_t nMulti[N_PLANE_MAX_TOFI] = { 0 }, iCounts[N_PLANE_MAX_TOFI] = { 0 };
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BTofiHitData* hitTofi = dynamic_cast<R3BTofiHitData*>(detTofi->At(ihit));

            if (IS_NAN(hitTofi->GetTime()))
                continue;

            Int_t iPlane = hitTofi->GetDetId();
            Double_t randx = (std::rand() / (float)RAND_MAX) - 0.5;
            Int_t ictemp = iCounts[iPlane - 1];
            x[iPlane - 1][ictemp] = hitTofi->GetX() + 2.7 * randx;
            y[iPlane - 1][ictemp] = hitTofi->GetY();
            t[iPlane - 1][ictemp] = hitTofi->GetTime();
            q[iPlane - 1][ictemp] = hitTofi->GetEloss();
            bar[iPlane - 1][ictemp] = hitTofi->GetBarId();
            if (std::abs(t[iPlane - 1][ictemp]) < 30.)
            {
                fh_tofi_Tot_hit[iPlane - 1]->Fill(bar[iPlane - 1][ictemp], q[iPlane - 1][ictemp]);
                fh_tofi_Z[iPlane - 1]->Fill(q[iPlane - 1][ictemp]);
            }

            fh_tofi_dt_hit[iPlane - 1]->Fill(bar[iPlane - 1][ictemp], y[iPlane - 1][ictemp]);
            fh_tofi_time_hit[iPlane - 1]->Fill(bar[iPlane - 1][ictemp], t[iPlane - 1][ictemp]);
            fh_tofi_bars[iPlane - 1]->Fill(bar[iPlane - 1][ictemp]);
            iCounts[iPlane - 1] += 1;
            nMulti[iPlane - 1] += 1;
        }

        for (Int_t i = 0; i < N_PLANE_MAX_TOFI; i++)
        {
            fh_tofi_multihit_hit[i]->Fill(nMulti[i]);
            /* if (i > 0)
             {
                 for (Int_t im1 = 0; im1 < iCounts[i]; im1++)
                 {
                     for (Int_t im2 = 0; im2 < iCounts[i - 1]; im2++)
                     {
                         Double_t tdif =
                             fmod(t[i][im1] - t[i - 1][im2] + c_range_ns + c_range_ns / 2, c_range_ns) - c_range_ns / 2;
                         fh_tofi_dt_hit[i - 1]->Fill(bar[i][im1], tdif);
                     }
                 }
             }*/
        }
    }

    time_prev = time;
    fNEvents += 1;
}

void R3BOnlineSpectraToFI_S494::FinishEvent()
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
}

void R3BOnlineSpectraToFI_S494::FinishTask()
{
    if (fCalItems.at(DET_TOFI))
    {
        for (Int_t i = 0; i < N_PLANE_MAX_TOFI; i++)
        {
            fh_tofi_TotPm[i]->Write();
            fh_tofi_timePm[i]->Write();
            fh_tofi_TotPm_coinc[i]->Write();
            fh_tofi_channels[i]->Write();
            fh_tofi_multihit[i]->Write();
            fh_tofi_multihit_coinc[i]->Write();
        }
        fh_tofi_dt->Write();
        fh_tofi_time->Write();
        fh_num_bars->Write();
        fh_test->Write();
        fh_test1->Write();
        fh_test2->Write();
        fh_test3->Write();
        fh_num_side->Write();
    }
    if (fHitItems.at(DET_TOFI))
    {
        for (Int_t i = 0; i < N_PLANE_MAX_TOFI; i++)
        {
            fh_tofi_Tot_hit[i]->Write();
            fh_tofi_time_hit[i]->Write();
            fh_tofi_multihit_hit[i]->Write();
            fh_tofi_bars[i]->Write();
            fh_tofi_time_hit[i]->Write();
            fh_tofi_dt_hit[i]->Write();
            fh_tofi_Z[i]->Write();
        }
        /*   for (Int_t i = 0; i < N_PLANE_MAX_TOFI - 1; i++)
           {
               fh_tofi_dt_hit[i]->Write();
           }*/
    }
}

ClassImp(R3BOnlineSpectraToFI_S494)
