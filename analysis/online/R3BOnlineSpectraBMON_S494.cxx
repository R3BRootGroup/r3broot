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
// -----                  R3BOnlineSpectraBMON_S494       -----
// -----          Created May 2021 by A. Kelic-Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */
#include "R3BOnlineSpectraBMON_S494.h"

#include "R3BRoluCalData.h"
#include "R3BRoluHitData.h"
#include "R3BRoluMappedData.h"

#include "R3BBeamMonitorMappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdMappedData.h"

#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TDatime.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TMath.h"
#include <vector>
#define IS_NAN(x) TMath::IsNaN(x)

#include <boost/multi_array.hpp>
#include <cassert>
#include <cmath> // NAN

template <class T>
static void init_array(T& array, double init)
{
    // boost::multi_array prefers cascaded for loops over index sets
    // (e.g. python's dict().keys())
    // thus we can not really use a range based for loop
    // and instead wrangle pointers like it is 1980 again.
    for (double* e = array.origin(); e < array.origin() + array.num_elements(); ++e)
        *e = init;
}

using namespace std;
namespace
{
    double c_period = 2048 * 5;
} // namespace

R3BOnlineSpectraBMON_S494::R3BOnlineSpectraBMON_S494()
    : R3BOnlineSpectraBMON_S494("OnlineSpectra", 1)
{
}

R3BOnlineSpectraBMON_S494::R3BOnlineSpectraBMON_S494(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvents(0)
{
}

R3BOnlineSpectraBMON_S494::~R3BOnlineSpectraBMON_S494()
{
    //	delete fh_Tpat ;
    //	delete fh_Trigger;
    //	delete fh_SEETRAM;
}

InitStatus R3BOnlineSpectraBMON_S494::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BOnlineSpectraBMON_S494::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
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
        fHitItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sHit", fDetectorNames[det]))));
    }
    maxevent = mgr->CheckMaxEventNo();

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    // Trigger and Tpat
    TCanvas* cTrigg = new TCanvas("Trigg", "Triggers", 10, 10, 650, 350);
    fh_Tpat = new TH1F("Tpat", "Tpat", 20, 0, 20);
    fh_Tpat->GetXaxis()->SetTitle("Tpat value");

    fh_Trigger = new TH1F("Trigger", "Trigger all", 20, 0, 20);
    fh_Trigger->GetXaxis()->SetTitle("Trigger value");

    cTrigg->Divide(2, 1);
    cTrigg->cd(1);
    gPad->SetLogy();
    fh_Trigger->Draw();
    cTrigg->cd(2);
    gPad->SetLogy();
    fh_Tpat->Draw();
    cTrigg->cd(0);

    //------------------------------------------------------------------------
    // Rolu detector

    if (fMappedItems.at(DET_ROLU))
    {
        TCanvas* cROLU = new TCanvas("ROLU", "ROLU", 10, 10, 650, 350);
        fh_rolu_tot = new TH2F("Rolu_tot", "ROLU ToT", 10, 0, 10, 600, 0, 300);
        fh_rolu_tot->GetXaxis()->SetTitle("Channel number");
        fh_rolu_tot->GetYaxis()->SetTitle("ToT / ns");
        fh_rolu_tot->GetXaxis()->CenterTitle(true);
        fh_rolu_tot->GetYaxis()->CenterTitle(true);

        fh_rolu_tof = new TH2F("Rolu_tof", "ROLU-TOFD ToF", 9, 0, 9, 4000, -2000, 2000);
        fh_rolu_tof->GetXaxis()->SetTitle("Channel number");
        fh_rolu_tof->GetYaxis()->SetTitle("Rolu time / ns");

        fh_rolu_channels = new TH1F("Rolu_channels", "ROLU channels", 6, 0, 6);
        fh_rolu_channels->GetXaxis()->SetTitle("Channel number");
        fh_rolu_channels->GetYaxis()->SetTitle("Counts");
        fh_rolu_channels->SetFillColor(31);
        fh_rolu_channels->GetXaxis()->CenterTitle(true);
        fh_rolu_channels->GetYaxis()->CenterTitle(true);

        cROLU->Divide(1, 2);
        cROLU->cd(1);
        fh_rolu_channels->Draw();
        cROLU->cd(2);
        fh_rolu_tot->Draw("colz");

        auto mainfolRolu = new TFolder("ROLU", "ROLU info");
        mainfolRolu->Add(cROLU);

        run->AddObject(mainfolRolu);
        run->GetHttpServer()->RegisterCommand("Reset_ROLU", Form("/Objects/%s/->Reset_ROLU_Histo()", GetName()));
    }

    if (fMappedItems.at(DET_BMON))
    {

        // get the theoretical calib factors for SEETRAM
        Double_t fexp = float(fsens_SEE + 9);
        Double_t fpow = float(pow(10., fexp));
        calib_SEE = 104457.9 * fpow;
        cout << "SEETRAM: " << fsens_SEE << ", " << fexp << ", " << fpow << ", " << calib_SEE << endl;
        // get the theoretical calib factors for IC
        Double_t fexp_ic = float(fsens_IC + 7);
        Double_t fpow_ic = float(pow(10., fexp_ic));
        calib_IC = 1; // 5551. * fpow_ic;
        cout << "IC     : " << fsens_IC << ", " << fexp_ic << ", " << fpow_ic << ", " << calib_IC << endl;

        TCanvas* cbmon = new TCanvas("Beam_Monitor", "Beam Monitors", 820, 10, 900, 900);
        Int_t Nbin_bmon = reset_time / read_time;

        fh_spill_length = new TH1F("spill_length", "Spill ", 1000, 0, 1000);
        fh_spill_length->GetXaxis()->SetTitle("time / sec");

        fh_IC = new TH1F("IC", "IC ", Nbin_bmon, 0, reset_time);
        fh_IC->GetXaxis()->SetTitle("spill number");
        fh_IC->GetYaxis()->SetTitle("IC counts");

        fh_SEE = new TH1F("SEETRAM", "SEETRAM ", 1000, 0, 1000);
        fh_SEE->GetXaxis()->SetTitle("spill number");
        fh_SEE->GetYaxis()->SetTitle("SEETRAM counts");

        fh_SEEPS = new TH1F("SEETRAM_PS", "SEETRAM per spill", 1000, 0, 1000);
        fh_SEEPS->GetXaxis()->SetTitle("spill number");
        fh_SEEPS->GetYaxis()->SetTitle("SEETRAM counts");

        fh_SEER = new TH1F("SEETRAM_R", "SEETRAM rate per second", 1000, 0, 1000);
        fh_SEER->GetXaxis()->SetTitle("spill number");
        fh_SEER->GetYaxis()->SetTitle("SEETRAM <rate>/s");

        fh_TOFDOR = new TH1F("TOFDOR", "TOFDOR ", Nbin_bmon, 0, reset_time);
        fh_TOFDOR->GetXaxis()->SetTitle("time / sec");
        fh_TOFDOR->GetYaxis()->SetTitle("TOFDOR counts");

        fh_SROLU1 = new TH1F("SROLU1", "SROLU1 ", Nbin_bmon, 0, reset_time);
        fh_SROLU1->GetXaxis()->SetTitle("time / sec");
        fh_SROLU1->GetYaxis()->SetTitle("SROLU1 counts");

        fh_SROLU2 = new TH1F("SROLU2", "SROLU2 ", Nbin_bmon, 0, reset_time);
        fh_SROLU2->GetXaxis()->SetTitle("time / sec");
        fh_SROLU2->GetYaxis()->SetTitle("SROLU2 counts");

        fh_IC_spill = new TH1F("IC_spill", "IC particles rate in kHz ", Nbin_bmon, 0, reset_time);
        fh_IC_spill->GetXaxis()->SetTitle("time / sec");
        fh_IC_spill->GetYaxis()->SetTitle("Particles / kHz");

        fh_SEE_spill = new TH1F("SEE_spill", "SEE particles rate in kHz", Nbin_bmon, 0, reset_time);
        fh_SEE_spill->GetXaxis()->SetTitle("time / sec");
        fh_SEE_spill->GetYaxis()->SetTitle("Particles / kHz");

        fh_SEE_spill_long = new TH1F("SEE_spill_long", "SEE particles rate in kHz", 3600, 0, 3600);
        fh_SEE_spill_long->GetXaxis()->SetTitle("time / min");
        fh_SEE_spill_long->GetYaxis()->SetTitle("Particles / MHz");

        fh_SEE_spill_raw = new TH1F("SEE_spill_raw", "SEE count rate in kHz", Nbin_bmon, 0, reset_time);
        fh_SEE_spill_raw->GetXaxis()->SetTitle("time / sec");
        fh_SEE_spill_raw->GetYaxis()->SetTitle("Counts / kHz");

        fh_TOFDOR_spill = new TH1F("TOFDOR_spill", "TOFDOR rate in kHz ", Nbin_bmon, 0, reset_time);
        fh_TOFDOR_spill->GetXaxis()->SetTitle("time / sec");
        fh_TOFDOR_spill->GetYaxis()->SetTitle("TOFDOR rate / kHz");

        fh_SROLU1_spill = new TH1F("SROLU1_spill", "SROLU1 rate in kHz ", Nbin_bmon, 0, reset_time);
        fh_SROLU1_spill->GetXaxis()->SetTitle("time / sec");
        fh_SROLU1_spill->GetYaxis()->SetTitle("SROLU1 rate / kHz");

        fh_SROLU2_spill = new TH1F("SROLU2_spill", "SROLU2 rate in kHz ", Nbin_bmon, 0, reset_time);
        fh_SROLU2_spill->GetXaxis()->SetTitle("time / sec");
        fh_SROLU2_spill->GetYaxis()->SetTitle("SROLU2 rate / kHz");

        fh_SEE_TOFDOR = new TH1F("SEE_vs_TOFDOR", "SEETRAM vs TOFDOR ", 20000, 0, 20000);
        fh_SEE_TOFDOR->GetYaxis()->SetTitle("SEETRAM counts/spill");
        fh_SEE_TOFDOR->GetXaxis()->SetTitle("TOFDOR counts/spill");

        fh_IC_TOFDOR = new TH1F("IC_vs_TOFDOR", "IC vs TOFDOR ", 20000, 0, 20000);
        fh_IC_TOFDOR->GetYaxis()->SetTitle("IC counts/spill");
        fh_IC_TOFDOR->GetXaxis()->SetTitle("TOFDOR counts/spill");

        fh_IC_SEE = new TH1F("IC_vs_SEE", "IC vs SEE ", 10000, 0, 10000);
        fh_IC_SEE->GetYaxis()->SetTitle("SEE particles/spill");
        fh_IC_SEE->GetXaxis()->SetTitle("IC counts /spill");

        cbmon->Divide(2, 2);
        cbmon->cd(1);
        fh_IC->Draw();
        cbmon->cd(2);
        fh_SEE->Draw();
        cbmon->cd(3);
        fh_SEEPS->Draw();
        cbmon->cd(4);
        fh_SEER->Draw();

        cbmon->cd(0);

        run->AddObject(cbmon);

        run->GetHttpServer()->RegisterCommand("Reset_BMON", Form("/Objects/%s/->Reset_BMON_Histo()", GetName()));
    }

    return kSUCCESS;
}

void R3BOnlineSpectraBMON_S494::Reset_ROLU_Histo()
{
    fh_rolu_channels->Reset();
    fh_rolu_tot->Reset();
    fh_Trigger->Reset();
    fh_Tpat->Reset();
    if (fHitItems.at(DET_TOFD))
        fh_rolu_tof->Reset();
}

void R3BOnlineSpectraBMON_S494::Reset_BMON_Histo()
{

    fh_spill_length->Reset();
    fh_IC->Reset();
    fh_IC_spill->Reset();
    fh_SEE->Reset();
    fh_SEE_spill->Reset();
    fh_SEE_spill_long->Reset();
    fh_SEE_spill_raw->Reset();
    fh_TOFDOR->Reset();
    fh_TOFDOR_spill->Reset();
    fh_SROLU1->Reset();
    fh_SROLU1_spill->Reset();
    fh_SROLU2->Reset();
    fh_SROLU2_spill->Reset();
    fh_SEE_TOFDOR->Reset();
    fh_IC_TOFDOR->Reset();
    fh_IC_SEE->Reset();
    time_start = -1;
}

void R3BOnlineSpectraBMON_S494::Exec(Option_t* option)
{
    fNEvents += 1;
    /*
        if (fNEvents / 10000. == (int)fNEvents / 10000)
            std::cout << "\rEvents: " << fNEvents << " / " << maxevent << " (" << (int)(fNEvents * 100. / maxevent)
                      << " %) " << std::flush;
    */
    Bool_t debug = false;

    if (header)
    {
        time = header->GetTimeStamp();
        //		if (time > 0) cout << "header time: " << time << endl;
        if (time_begin == -1 && time > 0 && time_end < 1)
        {
            time_begin = time;
        }

        if (time_start == -1 && time > 0)
        {
            time_start = time;
            fNEvents_start = fNEvents;
            // cout << "Start event number " << fNEvents_start << endl;
        }

        if (header->GetTrigger() == 12)
        {
            // spill start in nsec
            // cout << "spill start" << endl;
            num_spills++;
            newSpill = true;
            // cout << "Spill num: " << num_spills << endl;

            time_spill_start = time; // header->GetTimeStamp();    // spill start in nsec
            // cout << "Spill start: " << double(time_spill_start - time_begin) / 1.e9 << " sec " << endl;
            // reset counters:
            see_spill = 0;
            ic_spill = 0;
            tofdor_spill = 0;
            spill_on = true;
            in_spill_off = 0;
            fNSpills += 1;
        }
        if (header->GetTrigger() == 13)
        {
            // spill end  in nsec
            // cout << "spill stop" << endl;
            time_spill_end = time; // header->GetTimeStamp();    // spill end  in nsec
            // cout << "Spill stop: " << double(time_spill_end - time_begin) / 1.e9 << " sec " << endl;
            spill_on = false;
        }

        fh_Trigger->Fill(header->GetTrigger());
        //   check for requested trigger (Todo: should be done globablly / somewhere else)
        if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        {
            counterWrongTrigger++;
            return;
        }

        Int_t tpatbin;
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0)
                fh_Tpat->Fill(i + 1);
        }

        // fTpat = 1-16; fTpat_bit = 0-15
        Int_t fTpat_bit1 = fTpat1 - 1;
        Int_t fTpat_bit2 = fTpat2 - 1;
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2))
            {
                counterWrongTpat++;
                return;
            }
        }
    }

    // **************   ROLU *************************************
    if (fMappedItems.at(DET_ROLU))
    {
        auto det = fMappedItems.at(DET_ROLU);
        Int_t nHits = det->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BRoluMappedData* hit = dynamic_cast<R3BRoluMappedData*>(det->At(ihit));
            if (!hit)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hit->GetDetector(); // 1..
            Int_t iCha = hit->GetChannel();  // 1..

            if (iDet < 2)
                fh_rolu_channels->Fill(iCha); // ROLU 1
            if (iDet > 1)
                fh_rolu_channels->Fill(iCha + 4); // ROLU 2
        }
    }

    if (fCalItems.at(DET_ROLU))
    {
        Int_t nParts = 0;
        auto det = fCalItems.at(DET_ROLU);
        nParts = det->GetEntriesFast();

        if (nParts > 0)
        {
            fNEventsRolu += 1;
            Int_t iDet = 0;
            /*
             * Note: double x[nParts][2][4]={NAN};
             * will initialize everything other than x[0][0][0] to 0.0.
             * (and also not compile with our ancient gcc 4.8.5.)
             */
            using A = boost::multi_array<double, 3>;
            auto dims = boost::extents[nParts][2][4];
            A timeRolu_L(dims);
            init_array(timeRolu_L, NAN);
            A timeRolu_T(dims);
            init_array(timeRolu_T, NAN);
            A totRolu(dims);
            init_array(totRolu, NAN);

            for (Int_t iPart = 0; iPart < nParts; iPart++)
            {
                /*
                 * nParts is the number of particle passing through detector in one event
                 */
                R3BRoluCalData* calData = dynamic_cast<R3BRoluCalData*>(det->At(iPart));
                iDet = calData->GetDetector();

                for (Int_t iCha = 0; iCha < 4; iCha++)
                {

                    if (!(IS_NAN(calData->GetTimeL_ns(iCha))))
                    { // TAMEX leading
                        timeRolu_L[iPart][iDet - 1][iCha] = calData->GetTimeL_ns(iCha);
                    }
                    if (!(IS_NAN(calData->GetTimeT_ns(iCha))))
                    { // TAMEX trailing
                        timeRolu_T[iPart][iDet - 1][iCha] = calData->GetTimeT_ns(iCha);
                    }

                    if (timeRolu_T[iPart][iDet - 1][iCha] > 0. && timeRolu_L[iPart][iDet - 1][iCha] > 0. &&
                        !(IS_NAN(timeRolu_T[iPart][iDet - 1][iCha])) && !(IS_NAN(timeRolu_L[iPart][iDet - 1][iCha])))
                    {
                        while (timeRolu_T[iPart][iDet - 1][iCha] - timeRolu_L[iPart][iDet - 1][iCha] <= 0.)
                        {
                            timeRolu_T[iPart][iDet - 1][iCha] = timeRolu_T[iPart][iDet - 1][iCha] + 2048. * fClockFreq;
                        }

                        totRolu[iPart][iDet - 1][iCha] =
                            timeRolu_T[iPart][iDet - 1][iCha] - timeRolu_L[iPart][iDet - 1][iCha];
                    }

                    if (iDet < 2)
                        fh_rolu_tot->Fill(iCha + 1, totRolu[iPart][iDet - 1][iCha]);
                    if (iDet > 1)
                        fh_rolu_tot->Fill(iCha + 5, totRolu[iPart][iDet - 1][iCha]);
                }

                if (!calData)
                {
                    cout << "Rolu !calData" << endl;
                    continue; // can this happen?
                }
            }
        }
    } // end if fCalItems(ROLU)

    if (fMappedItems.at(DET_BMON))
    {
        unsigned long IC;
        unsigned long SEETRAM_raw;
        Double_t SEETRAM;
        unsigned long TOFDOR;

        auto detBmon = fMappedItems.at(DET_BMON);
        Int_t nHitsbm = detBmon->GetEntriesFast();
        // cout<<"Bmon hits: "<<nHitsbm<<endl;

        for (Int_t ihit = 0; ihit < nHitsbm; ihit++)
        {

            R3BBeamMonitorMappedData* hit = dynamic_cast<R3BBeamMonitorMappedData*>(detBmon->At(ihit));
            if (!hit)
                continue;

            IC = hit->GetIC(); // negative values if offset not high enough
            counts_IC += (double)IC;

            SEETRAM_raw = hit->GetSEETRAM();           // raw counts
            SEETRAM = (double)SEETRAM_raw * calib_SEE; // calibrated SEETRAM counts
            // cout<<SEETRAM_raw<<" "<<calib_SEE<<" "<<SEETRAM<<"\n";
            counts_SEE += SEETRAM;

            TOFDOR = hit->GetTOFDOR(); // only positive values possible
            counts_TofD += TOFDOR;

            if (fNEvents == fNEvents_start)
            {
                see_start = SEETRAM;
                ic_start = IC;
                tofdor_start = TOFDOR;
                cout << "SEE_start " << see_start << "\n";
                cout << "IC_start " << ic_start << "\n";
                cout << "TOFDOR_start " << tofdor_start << "\n";
            }
            // IC:
            unsigned long yIC = IC - ic_start;

            // SEETRAM:
            ySEE = SEETRAM - see_start;

            if (num_spills == 0)
            {
                nBeamParticle += ySEE - ySEE_mem;
                ySEE_mem = ySEE;
                seeLastSpill = 0;
            }

            if (newSpill == true)
            {
                spilltime = (time / 1e9) - time_lastSpill;
                TDatime t;
                t.Set(time / 1e9);

                nBeamParticle += (ySEE - seeLastSpill);

                if (debug)
                {
                    cout << "*******************************************" << endl;
                    cout << "Time WR " << time << " -> " << t.GetDay() << "." << t.GetMonth() << "." << t.GetYear()
                         << " " << t.GetHour() << ":" << t.GetMinute() << ":" << t.GetSecond() << endl;
                    cout << "last spill length " << (time / 1e9) - time_lastSpill << " s; SEE last spill "
                         << seeLastSpill << "\n";
                    cout << "Spill: " << num_spills << endl;
                    cout << "SEE  " << SEETRAM_raw << "  " << SEETRAM << "  " << counts_SEE << "  " << ySEE << endl;
                    cout << "incr SEE " << ySEE - seeLastSpill << "\n";
                    cout << "nBeamPart: " << nBeamParticle << endl;
                }

                seeLastSpill = ySEE;
                time_lastSpill = time / 1e9;
                newSpill = false;
            }

            // TOFDOR:
            unsigned long yTOFDOR = TOFDOR - tofdor_start;
            if (1 == 1) // spill_on == false)
            {
                fh_IC->SetBinContent(num_spills, yIC);
                fh_SEE->SetBinContent(num_spills, ySEE);
                fh_SEEPS->SetBinContent(num_spills, ySEE - seeLastSpill);
                fh_SEER->SetBinContent(num_spills, (ySEE - seeLastSpill) / spilltime);
                fh_TOFDOR->SetBinContent(num_spills, yTOFDOR);
            }

            if (num_spills != num_spills_prev)
            {
                // cout << "IC   " << IC << "  " << counts_IC << "  " << yIC << endl;
                // cout << "TofD " << TOFDOR << "  " << counts_TofD << "  " << yTOFDOR << endl;
            }
            num_spills_prev = num_spills;

            if (fNEvents == maxevent && time / 1.e9 > time_lastSpill && newSpill == false)
            {
                cout << "******** LAST EVENT *****" << endl;
                cout << "time till spill start " << (time / 1e9) - time_lastSpill << " s; SEE last spill "
                     << seeLastSpill << "\n";
                cout << "Spill: " << num_spills << endl;
                cout << "SEE  " << SEETRAM_raw << "  " << SEETRAM << "  " << counts_SEE << "  " << ySEE << endl;
                cout << "incr SEE " << ySEE - seeLastSpill << "\n";
                nBeamParticle += (ySEE - seeLastSpill);
                cout << "nBeamPart: " << nBeamParticle << endl;
                cout << "******* END LAST EVENT **" << endl;
            }
        }
    }

    time_end = header->GetTimeStamp();
}

void R3BOnlineSpectraBMON_S494::FinishEvent()
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
    }
}

void R3BOnlineSpectraBMON_S494::FinishTask()
{
    R3BLOG(info,
           "\n"
               << "nEvents total " << fNEvents << "\n"
               << "nEvents Rolu " << fNEventsRolu << "\n"
               << "Time_start      : " << time_begin << "\n"
               << "Time end        : " << time_end << "\n"
               << "Time duration   : " << (double)(time_end - time_begin) / 1.e9 << " sec \n"
               << "nSpill          : " << fNSpills);

    if (fMappedItems.at(DET_ROLU))
    {
        fh_rolu_channels->Write();
        fh_rolu_tot->Write();
        if (fHitItems.at(DET_TOFD))
            fh_rolu_tof->Write();
    }

    if (fMappedItems.at(DET_BMON))
    {
        fh_Tpat->Write();
        fh_spill_length->Write();
        fh_Trigger->Write();
        fh_IC->Write();
        fh_SEE->Write();
        fh_SEEPS->Write();
        fh_SEER->Write();
        fh_TOFDOR->Write();
        fh_IC_spill->Write();
        fh_SEE_spill_long->Write();
        fh_SEE_spill->Write();
        fh_SEE_spill_raw->Write();
        fh_TOFDOR_spill->Write();
        fh_SEE_TOFDOR->Write();
        fh_IC_TOFDOR->Write();
        fh_IC_SEE->Write();
    }
}

ClassImp(R3BOnlineSpectraBMON_S494);
