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
// -----                R3BOnlineSpillAnalysis            -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

#include "R3BOnlineSpillAnalysis.h"

#include "R3BEventHeader.h"

#include "R3BSamplerMappedData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"

#include "THttpServer.h"

#include "TVirtualFFT.h"

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

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BOnlineSpillAnalysis::R3BOnlineSpillAnalysis()
    : R3BOnlineSpillAnalysis("OnlineSpillAnalysis", 1)
{
}

R3BOnlineSpillAnalysis::R3BOnlineSpillAnalysis(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fSamplerMappedItems(NULL)
    , fTrigger(-1)
    , fTpat(-1)
    , fSpillLength(2.)
    , fNEvents(0)
{
}

R3BOnlineSpillAnalysis::~R3BOnlineSpillAnalysis()
{
    if (fSamplerMappedItems)
        delete fSamplerMappedItems;
}

InitStatus R3BOnlineSpillAnalysis::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BOnlineSpillAnalysis::Init() ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    // Look for the R3BEventHeader
    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
    {
        LOG(warn) << "R3BOnlineSpillAnalysis::Init() EventHeader. not found";
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    }
    else
        LOG(info) << "R3BOnlineSpillAnalysis::Init() EventHeader. found";

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // fServer->Register(“/Server”, fServer);
    // fServer->RegisterCommand(“ProcessRequests”, Form(“/Server/%s/->ProcessRequests()”, fServer->GetName()));

    // Get objects for detectors on all levels
    fSamplerMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("SamplerMapped"));
    if (NULL == fSamplerMappedItems)
    {
        cout << "Could not find mapped data for Sampler." << endl;
    }

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    fh_dt_hits = new TH1F("los_dt_hits", "los dt hits", 200000, 0, 1000.);
    fh_dt_hits->GetXaxis()->SetTitle("dt between two hits / micros");
    fh_dt_hits->GetYaxis()->SetTitle("Counts");

    // Spill
    fh_spill_times = new TH1F("CountsPerSpillLOSRaw",
                              "counts per spill RAW [100us bin]",
                              fSpillLength * 10000,
                              0.,
                              fSpillLength); // bin 100µs
    fh_spill_times_Fine = new TH1F("CountsPerSpillLOSfine",
                                   "counts per spill fine [10us bin] ",
                                   fSpillLength * 100000,
                                   0.,
                                   fSpillLength); // 10µs bin
    fh_spill_times_Coarse = new TH1F("CountsPerSpillLOScoarse",
                                     "counts per spill coarse [10ms bin] ",
                                     fSpillLength * 100,
                                     0.,
                                     fSpillLength); // 10ms bin

    fh_spill_times_Fine_adj =
        new TH1F("CountsPerSpillLOSfineadj", "counts per spill fine [10us bin] ", 1000000, 0., 10.);      // 10µs bin
    fh_spill_times_FFT = new TH1F("CountsPerSpillFFT", "counts per spill [10us bin] ", 1000000, 0., 10.); // 10µs bin

    fh_SAMP_tDiff = new TH1F("Time_difference_SAMP", "Time-difference Sampler", 20000, 0., 1.); // 1ms, 5ns bin
    fh_SAMP_tDiff_long =
        new TH1F("Time_difference_SAMP_long", "Time-difference Sampler", 20000, 0., 1.); // 1ms, 5ns bin
    fh_SAMP_freq = new TH1F("SAMP_Frequencies", "Sampler Frequencies", 6000., 0., 6000.);
    fh_SAMP_freq_long = new TH1F("SAMP_Frequencies_long", "Sampler Frequencies", 10000., 0., 1000.);

    fh_SAMP_tDiff_pois = new TH1F("Random_distributed_tDIFF", "Random_tDiff", 20000, 0., 1.); // 0-1 ms, 5ns bin

    fh_FFT = new TH1F("FFT", "FFT", 1000000, 0., 10.);
    fh_FFT_adj = new TH1F("adj_FFT", "FFT", 1000000, 0., 10.);
    fh_FFT_add = new TH1F("added_FFT", "FFT summed", 1000000, 0., 10.);

    fh_hans_sum = new TH1F("Spillquality_sum", "Spill quality HT", 600, -300., 300.);
    fh_hans_sum_byMax = new TH1F("Spillquality_sum_byMax", "Spill quality HT", 600, -300., 300.);
    fh_MtA_sum = new TH1F("MaxToAverage_sum", "Maximum to average", 1000, 0., 100.);

    fh_rate = new TH1F("rate", "accepted events relativ to poisson", 360, 0, 180);
    fh_rate->GetXaxis()->SetTitle("spill number");
    fh_rate->GetYaxis()->SetTitle("Percent");
    fh_rate->GetXaxis()->SetTitle("spill number");
    fh_rate->GetYaxis()->SetTitle("Percent");

    fh_spill_hans = new TH1F("Spill_quality", "Spill Quality", 360., 0., 360.);
    fh_spill_hans_byMax = new TH1F("Spill_quality_HT_byMax", "Spill Quality method 1", 100., 0., 100.);
    fh_DutyFactor = new TH1F(
        "DutyFactor", "Duty Factor", fSpillLength * 100., 0., fSpillLength /*, 2000., 0., 2.*/); // timebin in 10ms
    fh_DutyFactor_MtA = new TH1F("DutyFactor_MtA",
                                 "Duty Factor MtA",
                                 fSpillLength * 100.,
                                 0.,
                                 fSpillLength /*, 2000., 0., 2.*/); // timebin in 10ms
    fh_DutyFactor_pois = new TH1F(
        "DutyFactor_pois", "Duty Factor", fSpillLength * 100., 0., fSpillLength /*, 2000., 0., 2.*/); // timebin in 10ms
    fh_DutyFactor_MaxToAvg = new TH1F("MaxToAvg",
                                      "Duty Factor - Maximum to Average",
                                      fSpillLength * 100.,
                                      0.,
                                      fSpillLength * 50. /*, 2000., 0., 2.*/); // timebin in 10ms
    fh_DutyFactor_MaxRun = new TH1F("DFMaxRun",
                                    "Duty Factor - Maximum over run",
                                    fSpillLength * 100.,
                                    0.,
                                    fSpillLength * 50. /*, 2000., 0., 2.*/); // timebin in 10ms
    fh_DutyFactor_AvgRun = new TH1F("DFAvgRun",
                                    "Duty Factor - Average over run",
                                    fSpillLength * 100.,
                                    0.,
                                    fSpillLength * 50. /*, 2000., 0., 2.*/); // timebin in 10ms
    fh_DutyFactor_SAMP_clean = new TH1F("DutyFactor_SAMP_clean",
                                        "Duty Factor - rel. to poisson",
                                        fSpillLength * 100.,
                                        0.,
                                        fSpillLength * 50. /*, 2000., 0., 2.*/); // timebin in 10ms
    fh_rate_sum = new TH1F("rate_sum", "Rate accepted relative to poisson", 200, -100., 100.);

    fh_DutyFactor_Max = new TH1F("DutyFactor_Max", "DutyFactor Max", 100, 0., 1.);
    fh_DutyFactor_Avg = new TH1F("DutyFactor_Avg", "DutyFactor Avg", 100, 0., 1.);
    fh_DutyFactor_PLD = new TH1F("DutyFactor_PLD", "Spill Quality method2", 100, 0., 100.);

    cSpill = new TCanvas("Spill_SAMP", "Spill_SAMP", 10, 10, 910, 810);
    cSpill->Divide(3, 2);

    cSpill->cd(1);
    fh_spill_times->GetXaxis()->SetTitle("time / sec");
    fh_spill_times->GetYaxis()->SetTitle("counts");
    fh_spill_times->Draw();

    cSpill->cd(2);
    fh_SAMP_tDiff->GetXaxis()->SetTitle("dt/ms");
    fh_SAMP_tDiff->GetYaxis()->SetTitle("counts");
    gPad->SetLogx();
    gPad->SetLogy();
    fh_SAMP_tDiff->Draw();
    fh_SAMP_tDiff_pois->SetLineColor(kRed);
    fh_SAMP_tDiff_pois->Draw("same L");

    cSpill->cd(3);
    fh_spill_hans_byMax->GetYaxis()->SetTitle("Spill quality in percent");
    fh_spill_hans_byMax->GetXaxis()->SetTitle("Spill number");
    fh_spill_hans_byMax->SetFillColor(3);
    fh_spill_hans_byMax->GetYaxis()->SetRange(-100, 100);
    fh_spill_hans_byMax->Draw("column");

    cSpill->cd(4);
    fh_DutyFactor->GetXaxis()->SetTitle("time [s]");
    fh_DutyFactor->GetYaxis()->SetTitle("Duty Factor");
    fh_DutyFactor->SetFillColor(kGreen);
    fh_DutyFactor->Draw("HIST");

    cSpill->cd(5);
    fh_DutyFactor_SAMP_clean->GetXaxis()->SetTitle("time bins");
    fh_DutyFactor_SAMP_clean->GetYaxis()->SetTitle("Duty Factor");
    fh_DutyFactor_SAMP_clean->SetFillColor(kGreen);
    fh_DutyFactor_SAMP_clean->Draw("Hist");

    cSpill->cd(6);
    fh_DutyFactor_PLD->GetXaxis()->SetTitle("Spill number");
    fh_DutyFactor_PLD->GetYaxis()->SetTitle("Difference Avg to Pois-Limit");
    fh_DutyFactor_PLD->SetFillColor(kGreen);
    fh_DutyFactor_PLD->Draw("column");

    run->AddObject(cSpill);
    run->GetHttpServer()->RegisterCommand("Reset", Form("/Tasks/%s/->Reset_Histo()", GetName()));

    cFFT = new TCanvas("FFT", "FFT", 10, 10, 910, 810);
    cFFT->Divide(2, 2);

    cFFT->cd(1);
    gPad->SetLogy();
    fh_SAMP_freq_long->GetXaxis()->SetTitle("kHz");
    fh_SAMP_freq_long->GetYaxis()->SetTitle("counts");
    fh_SAMP_freq_long->Draw();

    cFFT->cd(2);
    fh_SAMP_tDiff_long->GetXaxis()->SetTitle("dt/ms");
    fh_SAMP_tDiff_long->GetYaxis()->SetTitle("counts");
    gPad->SetLogx();
    gPad->SetLogy();
    fh_SAMP_tDiff_long->Draw();

    cFFT->cd(3);
    gPad->SetLogy();
    fh_FFT_adj->GetXaxis()->SetTitle("Frequency [10kHz]");
    fh_FFT_adj->GetYaxis()->SetTitle("Magnitude arb. units");
    fh_FFT_adj->Draw();

    cFFT->cd(4);
    gPad->SetLogy();
    fh_FFT_add->GetXaxis()->SetTitle("Frequency [10kHz]");
    fh_FFT_add->GetYaxis()->SetTitle("Magnitude arb. units");
    fh_FFT_add->Draw();

    run->AddObject(cFFT);

    run->GetHttpServer()->RegisterCommand("Reset", Form("/Objects/%s/->Reset_Histo()", GetName()));

    // run->GetHttpServer()->RegisterCommand("Update", Form("/Tasks/%s/->Update_Histo()", GetName()));

    return kSUCCESS;
}

void R3BOnlineSpillAnalysis::Update_Histo()
{

    // fServer->Register(“/Server”, fServer);
    // fServer->RegisterCommand(“ProcessRequests”, Form(“/Server/%s/->ProcessRequests()”, fServer->GetName()));
}

void R3BOnlineSpillAnalysis::Reset_Histo()
{
    fh_spill_times->Reset();
    fh_spill_times_FFT->Reset();
    fh_spill_times_Fine->Reset();
    fh_spill_times_Coarse->Reset();
    fh_SAMP_tDiff->Reset();
    fh_SAMP_tDiff_long->Reset();
    fh_SAMP_freq->Reset();
    fh_SAMP_freq_long->Reset();
    fh_SAMP_tDiff_pois->Reset();
    fh_DutyFactor->Reset();
    fh_DutyFactor_pois->Reset();
    fh_DutyFactor_MtA->Reset();
    fh_DutyFactor_MaxToAvg->Reset();
    fh_DutyFactor_PLD->Reset();
    fh_DutyFactor_Max->Reset();
    fh_DutyFactor_Avg->Reset();
    fh_rate->Reset();
    fh_rate_sum->Reset();
    fh_spill_hans->Reset();
    fh_FFT->Reset();
    fh_FFT_adj->Reset();
    fh_FFT_add->Reset();
    fh_hans_sum->Reset();
    fh_hans_sum_byMax->Reset();
    fh_MtA_sum->Reset();
    fh_dt_hits->Reset();
    fh_DutyFactor_MaxRun->Reset();
    fh_DutyFactor_AvgRun->Reset();
}

void R3BOnlineSpillAnalysis::Exec(Option_t* option)
{

    if (0 == (fNEvents % 1000))
    {
        // cout << "Event: " << fNEvents << endl;
    }

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(error) << "FairRootManager not found";
        return;
    }

    Int_t trigger = header->GetTrigger();
    // cout << "Trigger: " << trigger << " requested trigger: " << fTrigger << endl;
    if (fTrigger >= 0 && (header) && trigger != fTrigger && trigger != 12 && trigger != 13)
    {
        // cout << "Wrong trigger!" << endl;
        return;
    }

    Int_t tpatbin;
    for (int i = 0; i < 16; i++)
    {
        tpatbin = (header->GetTpat() & (1 << i));
        // if (tpatbin != 0)
        // fh_Tpat->Fill(i + 1);
    }

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    Int_t itpat;
    Int_t tpatvalue;
    if (fTpat_bit >= 0)
    {
        itpat = header->GetTpat();
        tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if (tpatvalue == 0)
        {
            // counterWrongTpat++;
            return;
        }
    }

    if (header->GetTrigger() == 12)
    {
        cout << "Start of spill!" << endl;
        // time_spill_start = header->GetTimeStamp() / 1.6; // spill start in nsec
        // spill_on = true;
        spill_on_sampler = true;
        if (spillCounter > 99)
        {
            spillCounter = 0;
            fh_spill_hans_byMax->Reset();
            fh_DutyFactor_PLD->Reset();
        }

        spillCounter++;

        fh_spill_times->Reset();
        fh_spill_times_Fine->Reset();
        fh_spill_times_Coarse->Reset();
        fh_SAMP_tDiff->Reset();
        fh_dt_hits->Reset();

        fh_DutyFactor->Reset();
        fh_DutyFactor_SAMP_clean->Reset();
    }
    /*
        else
        {
            spill_on = false;
        }
    */
    if (spill_on && time > 0)
    {
        // time_spill_start = time; // spill start in nsec
        // cout << "Spill start: " << double(time_spill_start - time_start) / 1.e9 << " " << endl;
        // dt_LOS = 0.;
        // tt_prev = 0.;
        // spill_ctr_LOS = 0;
    }
    if (fSamplerMappedItems)
    {
        auto det = fSamplerMappedItems;
        Int_t nHitsSamp = det->GetEntriesFast();
        // cout << "Hits sampler: " << nHitsSamp << endl;
        for (Int_t ihit = 0; ihit < nHitsSamp; ihit++)
        {
            auto hit = dynamic_cast<R3BSamplerMappedData*>(det->At(ihit));
            // time is in steps of 10 ns
            // is is a 34 bit number, so max 1073741823
            samplerCurr = hit->GetTime();

            // cout << "TSampler: " << samplerCurr << endl;

            if (spill_on_sampler && ihit == 0)
            {
                samplerSpill = samplerCurr;
                spill_on_sampler = false;
                // cout << "new start time sample spill: " << samplerSpill << endl;
            }

            Int_t SE_ctr = 0;
            // spill_off_calc = false;
            Double_t SL_var = fSpillLength;
            long samp = (samplerCurr - samplerSpill); // time in 10 ns
            if (samp < 0)
                samp += 1073741823;
            if ((double)samp / (1e8 * 1.6) > fSpillLength && spill_ctr_LOS > 0)
            {
                SE_ctr = 1;
                MissedSpillEnd = true;
                SL_var = (double)samp / (1e8 * 1.6) - SL_var;
            }
            else if (MissedSpillEnd && SE_ctr == 0)
            {
                cout << "Missed spill end!" << endl;

                MissedSpillEnd = false;
            }

            Double_t dt = 0.;
            if (spill_ctr_LOS > 0)
            {
                dt = ((double)(samplerCurr - samplerPrev));
                if (dt < 0)
                    dt += 1073741823.;
            }

            if (spill_ctr_LOS == 0) // first trigger in the spill
            {
                dt_LOS_first = (double)samp;
                // cout << "Setting Los_first: " << dt_LOS_first << endl;
            }

            Double_t DIFF = (double)samp;

            fh_spill_times->Fill(DIFF / (1.e8)); // time in seconds
            fh_spill_times_Fine->Fill(DIFF / (1.e8));
            fh_spill_times_Coarse->Fill(DIFF / (1.e8));
            fh_spill_times_FFT->Fill(DIFF / (1.e8));

            // cout << "sampler time: " << samplerCurr << " previous: " <<
            // samplerPrev << " dt: " << dt << " Tdiff: " << DIFF << endl;

            if (spill_ctr_LOS > 0)
            {
                if (dt > 0.)
                {
                    fh_SAMP_tDiff->Fill(dt / 100000.);         /// samp in 10 ns => dt in 1ms w\ 10µs bin
                    fh_SAMP_tDiff_long->Fill(dt / 100000.);    /// samp in 10 ns => dt in 1ms w\ 10µs bin
                    fh_SAMP_freq->Fill(100. / dt * 1.E3);      // kHz
                    fh_SAMP_freq_long->Fill(100. / dt * 1.E3); // kHz
                    dt_prev = dt;
                }
            }
            spill_ctr_LOS++;

            samplerPrev = samplerCurr;
        }
    }

    if (header->GetTrigger() == 13)
    {
        cout << "End of spill!" << endl;
        // time_spill_end = header->GetTimeStamp() / 1.6; // spill end  in nsec
        spill_off = true;
    }
    else if (spill_off_calc)
        spill_off = true;

    else
        spill_off = false;

    if (spill_off == true)
    {
        fh_DutyFactor_pois->Reset();
        fh_DutyFactor_MtA->Reset();
        fh_SAMP_tDiff_pois->Reset();
        fh_spill_times_Fine_adj->Reset();
        // time_spill_end = time; // spill end  in nsec
        // cout << "Spill stop: " << double(time_spill_end - time_start) / 1.e9 << " " << endl;

        Int_t pps_LOS = fh_spill_times->Integral(1, fh_spill_times->GetSize() - 2); /// pps_fib = parts per spill?

        /* Get Duty-Factor. Get the mean value over 10ms of 1000 10 µs bin of spill times (Rahul Singh)*/
        for (int k1 = 0; k1 < 2; k1++)
        { // Loop over "only sampler"(k1=0) & "sampler & fiber1" (k1=1)
            int Xmax = 0;
            if (k1 == 0)
                Xmax = fh_spill_times_Fine->FindLastBinAbove(0, 1);

            Int_t N_count = 0;
            Int_t N_MtA = 0;
            Int_t Np = 0;
            Int_t Nps = 0;
            Int_t NCTR = 0;
            Double_t FD = 0.;
            Double_t DutyAvg = 0.;
            Double_t DutyAvg_clean = 0.;
            Double_t PoisAvg = 0.;
            Int_t AvgCtr = 0;
            /*Get DutyFactor like RS */

            for (int i = 0; i <= Xmax; i++)
            {

                if (k1 == 0)
                    N_count = fh_spill_times_Fine->GetBinContent(i);
                if (N_MtA < N_count)
                    N_MtA = N_count;
                Np += N_count;
                Nps += N_count * N_count;
                NCTR++;
                if ((((i + 1) % 1000 == 0) || i == Xmax) && (Nps > 0))
                {
                    if (i < Xmax)
                    {
                        F_duty = (double)(Np * Np) / (double)Nps * 1. / 1000.; /* F_Duty = MEAN(Np)² / MEAN(Nps) */
                        FD = Np / 1000.;
                    }
                    else
                    {
                        F_duty = (double)(Np * Np) / (double)Nps * 1. / (double)((i) % 1000);
                        FD = Np / (double)((i) % 1000);
                    }

                    Np = 0;
                    Nps = 0;
                    Double_t Filler =
                        ((double)(i)) / 1000.; /// Not exactly sure why, but need to reduce i by 1, otherwise there
                                               /// occurs a problem in the binning and some bins arent getting filled,
                                               /// while the ones before are filled twice.

                    if (Filler < 0.)
                        Filler += fSpillLength;
                    if (F_duty > -1000. && F_duty < 1000.)
                    {
                        if (k1 == 0)
                        {
                            fh_DutyFactor->SetBinContent(Filler, F_duty);
                        }
                    }
                    if (FD / (FD + 1.) > -1000. && FD / (FD + 1.) < 1000.)
                    {
                        if (k1 == 0)
                            fh_DutyFactor_pois->SetBinContent(
                                Filler,
                                FD / (FD +
                                      1.)); /// Get the Poisson-Limit for a given N (=FD). This is the highest possible
                                            /// value. If DutyFactor is lower, then there is a problem with a device.
                    }
                    Double_t fd = (F_duty - (FD / (FD + 1.))) / (FD / (FD + 1.));
                    if (!IS_NAN(fd))
                    {
                        if (k1 == 0)
                            fh_DutyFactor_SAMP_clean->SetBinContent(Filler, fd * 100.);

                        DutyAvg_clean += (F_duty - FD / (FD + 1.)) / (FD / (FD + 1.));
                        DutyAvg += F_duty;
                        PoisAvg += FD / (FD + 1.);
                        AvgCtr++;
                        // cout << "DutyAvg_clean: " << DutyAvg_clean << "  DutyAvg: " << DutyAvg << "  PoisAvg: " <<
                        // PoisAvg << endl;
                    }
                    // cout<<"N_MtA: "<<N_MtA<<" pps: "<<pps_LOS<<endl;
                    if (k1 == 0)
                        fh_DutyFactor_MtA->SetBinContent(Filler, (double)N_MtA / (double)pps_LOS);
                    N_MtA = 0;
                }
            }
            if (k1 == 0)
            {
                DutyAvg = DutyAvg / (double)AvgCtr;
                PoisAvg = PoisAvg / (double)AvgCtr;

                // cout << "outside for   DutyAvg: " << DutyAvg << "  PoisAvg: " << PoisAvg << endl;

                Double_t DutyMax = fh_DutyFactor->GetBinContent(fh_DutyFactor->GetMaximumBin());
                Double_t DutyMax_clean = fh_DutyFactor->GetBinContent(fh_DutyFactor->GetMaximumBin());
                if (DutyMax / DutyAvg > -1000. && DutyMax / DutyAvg < 1000. && !IS_NAN(DutyAvg) && !IS_NAN(PoisAvg))
                {
                    // cout << "Test: " << spillCounter << "  " << DutyMax << "  " << DutyAvg << endl;

                    fh_DutyFactor_Avg->Fill(DutyAvg);
                    fh_DutyFactor_Max->Fill(DutyMax);
                    fh_DutyFactor_MaxToAvg->SetBinContent(spillCounter, DutyMax / DutyAvg);
                    if (DutyMax > -1000. && DutyMax < 1000.)
                        fh_DutyFactor_MaxRun->SetBinContent(spillCounter, 100. * DutyMax);
                    if (DutyAvg > -1000. && DutyAvg < 1000.)
                        fh_DutyFactor_AvgRun->SetBinContent(spillCounter, 100. * DutyAvg);
                    if (DutyAvg < 1000. && DutyAvg > -1000. && PoisAvg < 1000. && PoisAvg > -1000.)
                        fh_DutyFactor_PLD->SetBinContent(spillCounter, (DutyAvg - PoisAvg) / PoisAvg * 100.);
                    fh_MtA_sum->Fill(DutyMax / DutyAvg);
                }

                // cout<<"1Max. DUTYFACTOR: "<<DutyMax<<" Avg Duty: "<<DutyAvg<<" PoisAvg: "<<PoisAvg<<endl;
            }
            if (k1 == 1)
            {
                DutyAvg = DutyAvg / (double)AvgCtr;
                DutyAvg_clean = DutyAvg_clean / (double)AvgCtr;
                PoisAvg = PoisAvg / (double)AvgCtr;
                // cout << "Avg Duty: "<<DutyAvg<<" PoisAvg: "<<PoisAvg<<endl;
            }
        }

#define DO_FFT
#ifdef DO_FFT
        for (int ft = 0; ft < 2; ft++)
        {
            /*Get FFT of fh_spill_times_Fine */
            // if (spillCounter > 1)
            // fh_FFT->Reset();
            // if (ft == 0)
            // fh_spill_times_Fine_adj->Reset();
            // if (ft == 0)
            // fh_FFT_adj->Reset();

            Double_t mean_of_spill = 0;
            Int_t Last_of_spill = 0;
            Int_t First_of_spill = 0;
            Int_t Int_mean = 0;
            if (ft == 0)
                Last_of_spill = fh_spill_times_FFT->FindLastBinAbove(0, 1);
            if (ft == 0)
                First_of_spill = fh_spill_times_FFT->FindFirstBinAbove(0, 1);
            if (ft == 0)
                Int_mean = fh_spill_times_FFT->Integral(First_of_spill, Last_of_spill);
            Int_t mean_count = 0;
            for (int k1 = First_of_spill; k1 < Last_of_spill + 1; k1++)
            {
                if (ft == 0)
                    mean_of_spill += fh_spill_times_FFT->GetBinContent(k1);
                mean_count++;
            }

            /// Substract the mean (here calculated via integral) from the spill to eliminate the lowest freq. (0th bin
            /// of FFT)
            for (int i = First_of_spill; i < Last_of_spill + 1; i++)
            {
                if (ft == 0)
                    fh_spill_times_Fine_adj->SetBinContent(
                        i, fh_spill_times_FFT->GetBinContent(i) - (double)mean_of_spill / (double)mean_count);
            }

            TH1* fft = 0;
            // actual FFT
            TVirtualFFT::SetTransform(0);
            if (ft == 0)
                fh_FFT_adj = fh_spill_times_Fine_adj->FFT(fh_FFT_adj, "MAG");
            /// Add FFTs
            if (ft == 0)
                fh_FFT_add->Add(fh_FFT_adj);
        }
#endif

        /*Get the TimeDifferences of the randomized poisson Spill */
        Double_t lambda = 1. / (fh_SAMP_tDiff->GetMean());
        // cout << "lambda: " << lambda << endl;
        Double_t nPPS = 0.;
        /*
                if (pps_LOS < 100000.)
                    nPPS = pps_LOS * 1000.;
                if (pps_LOS >= 100000.)
                    nPPS = pps_LOS * 10.;
        */

        nPPS = pps_LOS * 10;
        for (int i = 0; i < nPPS; i++)
        {
            /* Make a poisson distribution with same lambda value as timedifferences*/
            Double_t num = ran_expo(lambda);
            fh_SAMP_tDiff_pois->Fill(num, 0.1);
        }

        // Calculate amount of low timedifferences relative to poisson (Hans Törnqvist)
        //        for (int k2 = 0; k2 < 2; k2++)
        //        { // Loop over "only sampler"(k2=0) & "sampler & fiber1" (k2=1)
        // Now get the Integral of the real, and the poisson dt. Get the maximum of pois. Ratio of Area  until the
        // Max is quality Int_t max_pois = fh_SAMP_tDiff_pois->GetMaximumBin(); Int_t max_real =
        // fh_SAMP_tDiff->GetMaximumBin();
        Int_t max_pois = (int)(fh_SAMP_tDiff_pois->GetMean() * 1.e5 / 5.);
        Int_t max_real = (int)(fh_SAMP_tDiff->GetMean() * 1.e5 / 5.); // in ns
        cout << "max_real: " << max_real << " max_pois: " << max_pois << endl;
        Double_t int_pois = 0.;
        Double_t int_real = 0.;
        Double_t int_pois_byMax = 0.;
        Double_t int_real_byMax = 0.;
        int_pois = fh_SAMP_tDiff_pois->Integral(1, 200);            // Get Integral of tDiffs equal and below 1µs
        int_real = fh_SAMP_tDiff->Integral(1, 200);                 // with 5ns bin from 0th to 200th bin
        int_pois_byMax = fh_SAMP_tDiff_pois->Integral(1, max_pois); // Get Integral of tDiffs equal and below 1µs
        int_real_byMax = fh_SAMP_tDiff->Integral(1, max_real);      // with 5ns bin from 0th to 200th bin
        Double_t int_real_LOS = fh_dt_hits->Integral(1, 200);       // with 5ns bin from 0th to 200th bin
        /**One can either use 1µs tDiff as Integration limit or the mean -> both show the same course, but the
         * absolute values are a bit different. Both seem valid and independent of particle rate, because we compare
         * it to a poisson distribution with the same mean value => the same particle rate. */

        cout << "int_pois: " << int_pois << " int_real: " << int_real << endl;
        cout << "int_pois max: " << int_pois_byMax << " int_real max: " << int_real_byMax << endl;
        cout << "Spill #: " << spillCounter << endl;
        Double_t frac_int_real = 0;
        frac_int_real = int_real / fh_SAMP_tDiff->Integral(1, fh_SAMP_tDiff->GetSize() - 2);
        Double_t frac_int_real_LOS = int_real_LOS / fh_dt_hits->Integral(1, 100000);
        Double_t frac_int_pois = int_pois / fh_SAMP_tDiff_pois->Integral(1, fh_SAMP_tDiff_pois->GetSize() - 2);
        Double_t frac_int_real_byMax = 0;
        frac_int_real_byMax = int_real_byMax / fh_SAMP_tDiff->Integral(1, fh_SAMP_tDiff->GetSize() - 2);
        Double_t frac_int_pois_byMax =
            int_pois_byMax / fh_SAMP_tDiff_pois->Integral(1, fh_SAMP_tDiff_pois->GetSize() - 2);

        cout << "Integral real: " << fh_SAMP_tDiff->Integral(1, fh_SAMP_tDiff->GetSize() - 2)
             << "Integral pois: " << fh_SAMP_tDiff_pois->Integral(1, fh_SAMP_tDiff_pois->GetSize() - 2) << endl;
        // cout<<"frac_int_real: "<<frac_int_real<<endl;
        // cout<<"frac_int_pois: "<<frac_int_pois<<endl;
        // cout<<"frac_int_real_byMax: "<<frac_int_real_byMax<<endl;
        // cout<<"frac_int_pois_byMax: "<<frac_int_pois_byMax<<endl;

        /* Since the exponential distribution has many entries (high statistics) compare the fractions of dt's under
         * 1µs */
        Double_t Hans2 = (frac_int_pois - frac_int_real) / frac_int_pois; // >0 = good; <0 = bad
        Double_t Hans = (frac_int_pois_byMax - frac_int_real_byMax) / frac_int_pois_byMax;

        if (!IS_NAN(Hans))
        {
            Hans_mean += Hans;
            hans_ctr++;
        }
        if (!IS_NAN(Hans2))
        {
            fh_spill_hans->SetBinContent(spillCounter, 100. * Hans2);
            fh_hans_sum->Fill(100. * Hans2);
        }
        if (!IS_NAN(Hans))
        {
            fh_spill_hans_byMax->SetBinContent(spillCounter, 100. * Hans);
            fh_hans_sum_byMax->Fill(100. * Hans2);
        }
        //        }

        fh_spill_times_Fine->Reset();
        fh_spill_times_Coarse->Reset();

        // cout << "updating" << endl;
        for (int i = 0; i < 6; i++)
        {
            TVirtualPad* pad = cSpill->cd(i + 1);
            pad->Modified();
            pad->Update();
        }
        for (int i = 0; i < 4; i++)
        {
            TVirtualPad* pad = cFFT->cd(i + 1);
            pad->Modified();
            pad->Update();
        }

        FairRunOnline::Instance()->GetHttpServer()->ProcessRequests();

    } // end spill_off

    fNEvents += 1;
}

void R3BOnlineSpillAnalysis::FinishEvent()
{
    if (fSamplerMappedItems)
        fSamplerMappedItems->Clear();
}

void R3BOnlineSpillAnalysis::FinishTask()
{

    fh_dt_hits->Write();
    fh_spill_times->Write();
    fh_spill_times_Fine->Write();
    fh_spill_times_FFT->Write();
    fh_spill_times_Coarse->Write();
    fh_SAMP_tDiff->Write();
    fh_SAMP_tDiff_long->Write();
    fh_SAMP_freq->Write();
    fh_SAMP_freq_long->Write();
    fh_SAMP_tDiff_pois->Write();
    fh_DutyFactor->Write();
    fh_DutyFactor_MtA->Write();
    fh_DutyFactor_pois->Write();
    fh_DutyFactor_MaxToAvg->Write();
    fh_DutyFactor_Avg->Write();
    fh_DutyFactor_Max->Write();
    fh_DutyFactor_PLD->Write();
    fh_spill_hans->Write();
    fh_FFT->Write();
    fh_FFT_adj->Write();
    fh_FFT_add->Write();
    fh_hans_sum->Write();
    fh_MtA_sum->Write();
    fh_rate_sum->Write();
    fh_rate->Write();
    fh_dt_hits->Write();
    fh_DutyFactor_MaxRun->Write();
    fh_DutyFactor_AvgRun->Write();
}

ClassImp(R3BOnlineSpillAnalysis);
