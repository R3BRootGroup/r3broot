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
// -----                  R3BOnlineSpectraDec2019         -----
// -----          Created April 13th 2016 by M.Heil       -----
// ------------------------------------------------------------

/*
 * This task should fill histograms with detector variables which allow
 * to test the detectors online
 *
 */

#include "R3BOnlineSpectraDec2019.h"
#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"

#include "R3BRoluCalData.h"
#include "R3BRoluMappedData.h"

#include "R3BSci8CalData.h"
#include "R3BSci8MappedData.h"

#include "R3BBeamMonitorMappedData.h"

#include "R3BTofdCalData.h"
#include "R3BTofdMappedData.h"

#include "R3BPaddleCalData.h"
#include "R3BPaddleTamexMappedData.h"

#include "R3BPspxCalData.h"
#include "R3BPspxHitData.h"
#include "R3BPspxMappedData.h"
#include "R3BPspxPrecalData.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberMappedData.h"

#include "R3BSamplerMappedData.h"

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
#include "TLegend.h"

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

#include <assert.h>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

R3BOnlineSpectraDec2019::R3BOnlineSpectraDec2019()
    : R3BOnlineSpectraDec2019("OnlineSpectra", 1)
{
}

R3BOnlineSpectraDec2019::R3BOnlineSpectraDec2019(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fSamp(-1)
    , fNofPlanes(N_PLANE_MAX_TOFD)
    , fPaddlesPerPlane(N_PADDLE_MAX_TOFD)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fNEvents(0)
{
}

R3BOnlineSpectraDec2019::~R3BOnlineSpectraDec2019()
{
    //	delete fhTpat ;
    //	delete fhTrigger;
    //	delete fh_SEETRAM;

    for (int i = 0; i < NOF_FIB_DET; i++)
    {
        delete fh_channels_Fib[i];
        delete fh_fibers_Fib[i];
        delete fh_mult_Fib[i];
        delete fh_time_Fib[i];
        delete fh_multihit_m_Fib[i];
        delete fh_multihit_s_Fib[i];
        delete fh_ToT_m_Fib[i];
        delete fh_ToT_s_Fib[i];
        delete fh_ToT_single_Fib[i];
        delete fh_channels_single_Fib[i];
    }
}

InitStatus R3BOnlineSpectraDec2019::Init()
{

    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BOnlineSpectraDec2019::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    FairRunOnline* run = FairRunOnline::Instance();

    run->GetHttpServer()->Register("/Tasks", this);

    // Get objects for detectors on all levels
    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    printf("Have %d fiber detectors.\n", NOF_FIB_DET);
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

    //------------------------------------------------------------------------
    // create histograms of all detectors
    //------------------------------------------------------------------------

    // Trigger and Tpat
    TCanvas* cTrigg = new TCanvas("Trigg", "Triggers", 10, 10, 650, 350);
    fhTpat = new TH1F("Tpat", "Tpat", 20, 0, 20);
    fhTpat->GetXaxis()->SetTitle("Tpat value");

    fhTrigger = new TH1F("Trigger", "Trigger all", 20, 0, 20);
    fhTrigger->GetXaxis()->SetTitle("Trigger value");

    cTrigg->Divide(2, 1);
    cTrigg->cd(1);
    gPad->SetLogy();
    fhTrigger->Draw();
    cTrigg->cd(2);
    gPad->SetLogy();
    fhTpat->Draw();
    cTrigg->cd(0);

    //------------------------------------------------------------------------
    // Rolu detector

    if (fMappedItems.at(DET_ROLU))
    {
        TCanvas* cRolu = new TCanvas("Rolu", "Rolu", 10, 10, 800, 400);

        fh_rolu_tot = new TH2F("Rolu_tot", "ROLU ToT", 10, 0, 10, 3000, 0, 300);
        fh_rolu_tot->GetXaxis()->SetTitle("Channel number");
        fh_rolu_tot->GetYaxis()->SetTitle("ToT / ns");

        fh_rolu_channels = new TH1F("Rolu_channels", "ROLU channels", 10, 0, 10);
        fh_rolu_channels->GetXaxis()->SetTitle("Channel number");
        fh_rolu_channels->GetYaxis()->SetTitle("Counts");

        cRolu->Divide(2, 1);
        cRolu->cd(1);
        gPad->SetLogy();
        fh_rolu_channels->Draw();
        cRolu->cd(2);
        gPad->SetLogz();
        fh_rolu_tot->Draw("colz");
        cRolu->cd(0);
        run->AddObject(cRolu);

        run->GetHttpServer()->RegisterCommand("Reset_ROLU", Form("/Tasks/%s/->Reset_ROLU_Histo()", GetName()));
    }

    //------------------------------------------------------------------------
    // Sci8 detector

    if (fMappedItems.at(DET_SCI8))
    {
        TCanvas* cSci8 = new TCanvas("Sci8", "SCI8", 10, 10, 850, 850);

        fh_sci8_channels = new TH1F("sci8_channels", "SCI8 channels", 4, 0., 4.);
        fh_sci8_channels->GetXaxis()->SetTitle("Channel number");
        fh_sci8_multihit = new TH1F("sci8_multihit", "SCI8 multihit", 20, 0., 20.);
        fh_sci8_multihit->GetXaxis()->SetTitle("Multihit");
        fh_sci8_multihitVFTX = new TH2F("sci8_vftx_multihit", "SCI8 VFTX multihit", 4, 0., 4., 20, 0., 20.);
        fh_sci8_multihitVFTX->GetYaxis()->SetTitle("Multihit ");
        fh_sci8_multihitVFTX->GetXaxis()->SetTitle("Channel number");
        fh_sci8_multihitLEAD = new TH2F("sci8_lead_multihit", "SCI8 TAMEX Leading multihit", 4, 0., 4., 20, 0., 20.);
        fh_sci8_multihitLEAD->GetYaxis()->SetTitle("Multihit ");
        fh_sci8_multihitLEAD->GetXaxis()->SetTitle("Channel number");
        fh_sci8_multihitTRAI = new TH2F("sci8_trai_multihit", "SCI8 TAMEX Trailing multihit", 4, 0., 4., 20, 0., 20.);
        fh_sci8_multihitTRAI->GetYaxis()->SetTitle("Multihit ");
        fh_sci8_multihitTRAI->GetXaxis()->SetTitle("Channel number");
        fh_sci8_dt_hits = new TH1F("sci8_dt_hits", "SCI8 dt multihits", 4000, -2000, 2000);
        fh_sci8_dt_hits->GetXaxis()->SetTitle("dt_VFTX between two hits / ns");
        fh_sci8_dt_hits_l = new TH1F("sci8_dt_hits_l", "SCI8 dt multihits leading", 4000, -2000, 2000);
        fh_sci8_dt_hits_l->GetXaxis()->SetTitle("dt_LEAD between two hits / ns");
        fh_sci8_dt_hits_t = new TH1F("sci8_dt_hits_t", "SCI8 dt multihits trailing", 4000, -2000, 2000);
        fh_sci8_dt_hits_t->GetXaxis()->SetTitle("dt_TRAIL between two hits / ns");
        fh_sci8_tres_MCFD = new TH1F("sci8_time_res_MCFD", "SCI8 MCFD Time resolution - raw", 2000, -10., 10.);
        fh_sci8_tres_MCFD->GetXaxis()->SetTitle("Time MCFD / ns");
        fh_sci8_tres_TAMEX = new TH1F("sci8_time_res_TAMEX", "SCI8 TAMEX Time resolution -raw ", 8000, -4., 4.);
        fh_sci8_tres_TAMEX->GetXaxis()->SetTitle("Time TAMEX / ns");
        fh_sci8_tot = new TH2F("sci8_tot", "SCI8 ToT vs PM", 4, 0, 4, 3000, 0., 300.);
        fh_sci8_tot->GetXaxis()->SetTitle("PMT number");
        fh_sci8_tot->GetYaxis()->SetTitle("ToT / ns");
        fh_sci8_tot_mean = new TH1F("sci8_tot_mean", "SCI8 mean ToT", 3000, 0., 300.);
        fh_sci8_tot_mean->GetYaxis()->SetTitle("Counts");
        fh_sci8_tot_mean->GetXaxis()->SetTitle("ToT / ns");
        fh_tof_sci8 = new TH1F("sci8_tof", "SCI8-LOS ToF", 20000, -5000, 15000);
        fh_tof_sci8->GetXaxis()->SetTitle("ToF / ns");

        cSci8->Divide(2, 2);
        cSci8->cd(1);
        fh_sci8_channels->Draw();
        cSci8->cd(2);
        gPad->SetLogy();
        fh_sci8_multihit->Draw();
        // cSci8->cd(3);gPad->SetLogz();
        // fh_sci8_tot->Draw("colz");
        cSci8->cd(3);
        fh_sci8_tres_MCFD->Draw();
        // cSci8->cd(5);
        // fh_sci8_tres_TAMEX->Draw();
        cSci8->cd(4);
        fh_tof_sci8->Draw();
        // cSci8->cd(7);
        // fh_sci8_tot_mean->Draw();

        cSci8->cd(0);
        run->AddObject(cSci8);

        run->GetHttpServer()->RegisterCommand("Reset_SCI8", Form("/Tasks/%s/->Reset_SCI8_Histo()", GetName()));
    }

    //------------------------------------------------------------------------
    // Los detector
    TCanvas* cLos[fNofLosDetectors];

    if (fMappedItems.at(DET_LOS))
    {

        if (fNofLosDetectors == 2)
        {
            fh_los_ToT_vs_ToT = new TH2F("los_tot_vs_tot", "LOS ToT vs ToT", 3000, 0., 300., 3000, 0., 300.);
            fh_los_ToT_vs_ToT->GetXaxis()->SetTitle("ToT LOS1 / ns");
            fh_los_ToT_vs_ToT->GetYaxis()->SetTitle("ToT LOS2 / ns");
        }

        for (Int_t iloscount = 0; iloscount < fNofLosDetectors; iloscount++)
        {

            char detName[255];
            sprintf(detName, "LOS%d", iloscount + 1);

            cLos[iloscount] = new TCanvas(detName, detName, 10, 10, 910, 810);

            fh_los_channels[iloscount] =
                new TH1F(Form("%s_channels", detName), Form("%s channels", detName), 10, 0., 10.);
            fh_los_channels[iloscount]->GetXaxis()->SetTitle("Channel number");

            fh_los_multihit[iloscount] =
                new TH1F(Form("%s_multihit", detName), Form("%s multihit && all 8 PMs", detName), 30, 0., 30.);
            fh_los_multihit[iloscount]->GetXaxis()->SetTitle("Multihit");

            fh_los_pos_MCFD[iloscount] = new TH2F(
                Form("%s_pos_MCFD", detName), Form("%s MCFD Position ", detName), 1000, -5., 5., 1000, -5., 5.);
            fh_los_pos_MCFD[iloscount]->GetXaxis()->SetTitle("X position / cm");
            fh_los_pos_MCFD[iloscount]->GetYaxis()->SetTitle("Y position / cm");

            fh_los_pos_TAMEX[iloscount] = new TH2F(
                Form("%s_pos_TAMEX", detName), Form("%s TAMEX Position ", detName), 1000, -5., 5., 1000, -5., 5.);
            fh_los_pos_TAMEX[iloscount]->GetXaxis()->SetTitle("X position / cm");
            fh_los_pos_TAMEX[iloscount]->GetYaxis()->SetTitle("Y position / cm");

            fh_los_pos_ToT[iloscount] =
                new TH2F(Form("%s_pos_ToT", detName), Form("%s ToT Position ", detName), 1000, -5., 5., 1000, -5., 5.);
            fh_los_pos_ToT[iloscount]->GetXaxis()->SetTitle("X position / cm");
            fh_los_pos_ToT[iloscount]->GetYaxis()->SetTitle("Y position / cm");

            fh_los_dt_hits_ToT[iloscount] =
                new TH2F(Form("%s_dt_ToT", detName), Form("%s ToT dt ", detName), 100000, 0, 500, 500, 0, 200.);
            fh_los_dt_hits_ToT[iloscount]->GetXaxis()->SetTitle(
                "dt between two hits / micros");                             // dt between two hits / ns
            fh_los_dt_hits_ToT[iloscount]->GetYaxis()->SetTitle("ToT / ns"); // ToT / ns

            fh_los_tres_MCFD[iloscount] = new TH1F(
                Form("%s_time_res_MCFD", detName), Form("%s MCFD Time resolution - raw", detName), 4000, -4., 4.);
            fh_los_tres_MCFD[iloscount]->GetXaxis()->SetTitle("Time MCFD / ns");

            fh_los_tres_TAMEX[iloscount] = new TH1F(
                Form("%s_time_res_TAMEX", detName), Form("%s TAMEX Time resolution -raw ", detName), 4000, -4., 4.);
            fh_los_tres_TAMEX[iloscount]->GetXaxis()->SetTitle("Time TAMEX / ns");

            fh_los_tres_MTDC = new TH1F("los_time_res_MTDC", "LOS MTDC Time resolution - raw", 4000, -4, 4);
            fh_los_tres_MTDC->GetXaxis()->SetTitle("Time MTDC / ns");

            fh_los_tot[iloscount] =
                new TH2F(Form("%s_tot", detName), Form("%s ToT vs PMT", detName), 10, 0, 10, 3000, 0., 300.);
            fh_los_tot[iloscount]->GetXaxis()->SetTitle("PMT number");
            fh_los_tot[iloscount]->GetYaxis()->SetTitle("ToT / ns");

            fh_los_tot_mean[iloscount] =
                new TH1F(Form("%s_tot_mean", detName), Form("%s mean ToT", detName), 1500, 0., 300.);
            fh_los_tot_mean[iloscount]->GetYaxis()->SetTitle("Counts");
            fh_los_tot_mean[iloscount]->GetXaxis()->SetTitle("ToT / ns");

            fh_los_ihit_ToT[iloscount] =
                new TH2F(Form("%s_tot_ihit", detName), Form("%s ToT vs ihit", detName), 15, 0, 15, 9100, -10., 900.);
            fh_los_ihit_ToT[iloscount]->GetXaxis()->SetTitle("iHit");
            fh_los_ihit_ToT[iloscount]->GetYaxis()->SetTitle("ToT / ns");

            fh_los_dt_hits[iloscount] =
                new TH1F(Form("%s_los_dt_hits", detName), Form("%s los dt hits", detName), 100000, 0, 500);
            fh_los_dt_hits[iloscount]->GetXaxis()->SetTitle("dt between two hits / micros");
            fh_los_dt_hits[iloscount]->GetYaxis()->SetTitle("Counts");

            fh_dt_hits[iloscount] = new TH1F(Form("%s_dt_hits", detName), Form("%s dt hits", detName), 2000, 0, 500);
            fh_dt_hits[iloscount]->GetXaxis()->SetTitle("dt between two hits / micros");
            fh_dt_hits[iloscount]->GetYaxis()->SetTitle("Counts");

            fh_los_dt_first_ToT[iloscount] = new TH2F(
                Form("%s_dt_events_ToT", detName), Form("%s ToT vs deltat", detName), 5000, 0, 150000, 500, 0, 200.);
            fh_los_dt_first_ToT[iloscount]->GetXaxis()->SetTitle("dt between two events / ns");
            fh_los_dt_first_ToT[iloscount]->GetYaxis()->SetTitle("ToT / ns");

            cLos[iloscount]->Divide(4, 3);
            cLos[iloscount]->cd(1);
            fh_los_channels[iloscount]->Draw();
            cLos[iloscount]->cd(2);
            gPad->SetLogy();
            fh_los_multihit[iloscount]->Draw();
            cLos[iloscount]->cd(3);
            gPad->SetLogz();
            fh_los_tot[iloscount]->Draw("colz");
            cLos[iloscount]->cd(4);
            gPad->SetLogy();
            fh_los_tot_mean[iloscount]->Draw();
            cLos[iloscount]->cd(5);
            gPad->SetLogz();
            fh_los_ihit_ToT[iloscount]->Draw("colz");
            // cLos[iloscount]->cd(6);gPad->SetLogz();
            // fh_los_dt_first_ToT[iloscount]->Draw("colz");
            cLos[iloscount]->cd(6);
            gPad->SetLogz();
            fh_los_dt_hits_ToT[iloscount]->Draw("colz");
            cLos[iloscount]->cd(7);
            gPad->SetLogz();

            cLos[iloscount]->cd(8);
            gPad->SetLogy();
            fh_los_tres_MCFD[iloscount]->Draw();
            cLos[iloscount]->cd(9);
            gPad->SetLogz();
            fh_los_ToT_vs_ToT->Draw("colz");
            cLos[iloscount]->cd(10);
            gPad->SetLogz();
            fh_los_pos_ToT[iloscount]->Draw("colz");
            cLos[iloscount]->cd(11);
            gPad->SetLogz();
            fh_los_pos_MCFD[iloscount]->Draw("colz");
            cLos[iloscount]->cd(12);
            gPad->SetLogz();
            fh_los_pos_TAMEX[iloscount]->Draw("colz");
            cLos[iloscount]->cd(0);
            run->AddObject(cLos[iloscount]);
        }

        run->GetHttpServer()->RegisterCommand("Reset_LOS", Form("/Tasks/%s/->Reset_LOS_Histo()", GetName()));
    }

    //-----------------------------------------------------------------------

    if (fMappedItems.at(DET_BMON))
    {

        fh_TimePreviousEvent = new TH1F("TimePreviousEvent", "Time between 2 particles ", 300000, -3000, 3000);
        fh_TimePreviousEvent->GetXaxis()->SetTitle("time / µsec");
        fh_TimePreviousEvent->GetYaxis()->SetTitle("counts");

        TCanvas* cbmon = new TCanvas("Beam_Monitor", "Beam Monitors", 820, 10, 900, 900);
        Int_t Nbin_bmon = reset_time / read_time;

        fhTpat = new TH1F("Tpat", "Tpat", 20, 0, 20);
        fhTpat->GetXaxis()->SetTitle("Tpat value");

        fhTrigger = new TH1F("Trigger", "Trigger all", 20, 0, 20);
        fhTrigger->GetXaxis()->SetTitle("Trigger value");

        fh_spill_length = new TH1F("spill_length", "Spill ", Nbin_bmon, 0, reset_time);
        fh_spill_length->GetXaxis()->SetTitle("time / sec");

        fh_spill_times = new TH1F("CountsPerSpill", "counts per spill ", 120000, 0., 12.);
        fh_spill_times->GetXaxis()->SetTitle("time / sec");
        fh_spill_times->GetYaxis()->SetTitle("counts");

        fh_spill_times_zoom = new TH1F("CountsPerSpill_zoom", "counts per spill ", 50000, 0., 5.);
        fh_spill_times_zoom->GetXaxis()->SetTitle("time / sec");
        fh_spill_times_zoom->GetYaxis()->SetTitle("counts");

        fh_spill_times_pois = new TH1F("CountsPerSpill_pois", "counts per spill poisson", 50000, 0., 5.);
        fh_spill_times_pois->GetXaxis()->SetTitle("time / sec");
        fh_spill_times_pois->GetYaxis()->SetTitle("counts");
        // fh_spill_times_pois->SetLineColor(kRed);

        fh_spill_times_daq = new TH1F("CountsPerSpill_daq", "counts per spill daq", 50000, 0., 5.);
        fh_spill_times_daq->GetXaxis()->SetTitle("time / sec");
        fh_spill_times_daq->GetYaxis()->SetTitle("counts");
        // fh_spill_times_daq->SetLineColor(kRed);

        fh_spill_times_daq_first =
            new TH1F("CountsPerSpill_daq_first", "counts per spill daq first hit only", 50000, 0., 5.);
        fh_spill_times_daq_first->GetXaxis()->SetTitle("time / sec");
        fh_spill_times_daq_first->GetYaxis()->SetTitle("counts");
        // fh_spill_times_daq_first->SetLineColor(kRed);

        fh_frequencies = new TH1F("frequencies", "frequencies ", 10000., 0, 5000.);
        fh_frequencies->GetXaxis()->SetTitle("frequency / kHz");
        fh_frequencies->GetYaxis()->SetTitle("counts");

        fh_rate = new TH1F("rate", "accepted events relativ to poisson", 100, 0, 50);
        fh_rate->GetXaxis()->SetTitle("spill number");
        fh_rate->GetYaxis()->SetTitle("Percent");

        fh_TOFDOR = new TH1F("TOFDOR", "TOFDOR ", Nbin_bmon, 0, reset_time);
        fh_TOFDOR->GetXaxis()->SetTitle("time / sec");
        fh_TOFDOR->GetYaxis()->SetTitle("TOFDOR counts");

        fh_particle = new TH1F("particle_per_time", "Particle per time", 100, -50., 50.);
        fh_particle->GetXaxis()->SetTitle("particle per time");
        fh_particle->GetYaxis()->SetTitle("counts");

        fh_particle_pois = new TH1F("particle_per_time_pois", "Particle per time poisson", 100, -50., 50.);
        fh_particle_pois->GetXaxis()->SetTitle("particle per time");
        fh_particle_pois->GetYaxis()->SetTitle("counts");
        fh_particle_pois->SetLineColor(kRed);

        fh_TOFDOR_spill = new TH1F("TOFDOR_spill", "TOFDOR rate in kHz ", Nbin_bmon, 0, reset_time);
        fh_TOFDOR_spill->GetXaxis()->SetTitle("time / sec");
        fh_TOFDOR_spill->GetYaxis()->SetTitle("TOFDOR rate / kHz");

        cbmon->Divide(3, 3);
        //		cbmon->cd(1);gPad->SetLogy();
        //		fhTrigger->Draw();
        cbmon->cd(1);
        fh_spill_length->Draw();

        cbmon->cd(2);
        gPad->SetLogy();
        fh_dt_hits[0]->Draw();

        cbmon->cd(3);
        gPad->SetLogy();
        fh_los_dt_hits[0]->Draw();

        // cbmon->cd(4);
        // fh_spill_times_daq_first->Draw("hist");

        // cbmon->cd(5);
        // fh_spill_times_daq->Draw("hist");

        cbmon->cd(4);
        fh_spill_times->Draw("hist");
        // fh_spill_times_pois->Draw("same");

        cbmon->cd(5);
        fh_frequencies->Draw("hist");
        gPad->SetLogy();

        cbmon->cd(6);
        gPad->SetLogy();
        fh_particle->Draw("hist");
        fh_particle_pois->Draw("same");
        TLegend* leg_particle = new TLegend(0.7, 0.8, 0.9, 0.9);
        leg_particle->AddEntry(fh_particle, "Measured", "l");
        leg_particle->AddEntry(fh_particle_pois, "Poisson", "l");
        leg_particle->Draw("same");

        cbmon->cd(9);
        fh_rate->SetFillColor(3);
        fh_rate->GetYaxis()->SetRange(0, 110);
        fh_rate->Draw("B");

        //		cbmon->cd(9);
        //		fh_TOFDOR_spill->Draw("hist");

        cbmon->cd(0);

        run->AddObject(cbmon);

        run->GetHttpServer()->RegisterCommand("Reset_BMON", Form("/Tasks/%s/->Reset_BMON_Histo()", GetName()));
    }

    //-----------------------------------------------------------------------
    // Fiber Detectors 1-NOF_FIB_DET

    char canvName[255];
    TCanvas* FibCanvas[NOF_FIB_DET];

    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {

        std::stringstream fibhitobj;
        std::stringstream fibmapobj;
        std::stringstream FiName;
        std::string temp;

        if (fMappedItems.at(DET_FI_FIRST + ifibcount))
        {

            const char* chistName;
            const char* chistTitle;
            const char* detName;
            std::string tempTitle;
            std::string tempName;
            std::stringstream tempFibName;
            std::string tempFibNames;
            std::stringstream tempCanvName;

            detName = fDetectorNames[DET_FI_FIRST + ifibcount];

            cout << "I am creating canvas " << detName << endl;

            FibCanvas[ifibcount] = new TCanvas(detName, detName, 10, 10, 910, 910);

            // Channels:
            fh_channels_Fib[ifibcount] =
                new TH1F(Form("%s_channels", detName), Form("%s channels", detName), 520, 0., 520.);
            fh_channels_Fib[ifibcount]->GetXaxis()->SetTitle("Channel number");
            fh_channels_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Channels:
            fh_channels_single_Fib[ifibcount] =
                new TH1F(Form("%s_channels_single", detName), Form("%s channels of single PMTs", detName), 10, 0., 10.);
            fh_channels_single_Fib[ifibcount]->GetXaxis()->SetTitle("Channel number");
            fh_channels_single_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Fibers:
            fh_fibers_Fib[ifibcount] = new TH1F(
                Form("%s_fibers", detName), Form("%s fibers", detName), N_FIBER_PLOT_2019, 0., N_FIBER_PLOT_2019);
            fh_fibers_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_fibers_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Multiplicity (number of hit fibers):
            fh_mult_Fib[ifibcount] = new TH1F(Form("%s_mult", detName), Form("%s mult", detName), 100, 0., 100.);
            fh_mult_Fib[ifibcount]->GetXaxis()->SetTitle("Multiplicity");
            fh_mult_Fib[ifibcount]->GetYaxis()->SetTitle("Counts");

            // Multihit MAPMT:
            fh_multihit_m_Fib[ifibcount] = new TH2F(
                Form("%s_multihit_m", detName), Form("%s multihits MAPMT", detName), 520, 0., 520., 20, 0., 20.);
            fh_multihit_m_Fib[ifibcount]->GetXaxis()->SetTitle("MAPMT channel");
            fh_multihit_m_Fib[ifibcount]->GetYaxis()->SetTitle("Multihit");

            // Multihit SAPMT:
            fh_multihit_s_Fib[ifibcount] = new TH2F(
                Form("%s_multihit_s", detName), Form("%s multihits single PMT", detName), 16, 0., 16., 20, 0., 20.);
            fh_multihit_s_Fib[ifibcount]->GetXaxis()->SetTitle("Single PMT channel");
            fh_multihit_s_Fib[ifibcount]->GetYaxis()->SetTitle("Multihit");

            // ToT MAPMT:
            fh_ToT_m_Fib[ifibcount] = new TH2F(Form("%s_tot_m", detName),
                                               Form("%s ToT of MAPMT", detName),
                                               N_FIBER_PLOT_2019,
                                               0.,
                                               N_FIBER_PLOT_2019,
                                               2000,
                                               0.,
                                               200);
            fh_ToT_m_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_ToT_m_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // ToT SAPMT:
            fh_ToT_s_Fib[ifibcount] = new TH2F(Form("%s_tot_s", detName),
                                               Form("%s ToT of single PMT", detName),
                                               N_FIBER_PLOT_2019,
                                               0.,
                                               N_FIBER_PLOT_2019,
                                               2000,
                                               0.,
                                               200.);
            fh_ToT_s_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_ToT_s_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            // Time of fiber:
            fh_time_Fib[ifibcount] = new TH2F(Form("%s_TimevsFiber", detName),
                                              Form("%s Time vs Fiber", detName),
                                              N_FIBER_PLOT_2019,
                                              0.,
                                              N_FIBER_PLOT_2019,
                                              2048,
                                              -1024.,
                                              1024.);
            fh_time_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_time_Fib[ifibcount]->GetYaxis()->SetTitle("tMAPMT-tSPMT");

            // ToF LOS -> Fiber:
            fh_Fib_ToF[ifibcount] = new TH2F(Form("%s_tof", detName),
                                             Form("%s ToF LOS to Fiber", detName),
                                             N_FIBER_PLOT_2019,
                                             0.,
                                             N_FIBER_PLOT_2019,
                                             100000,
                                             0.,
                                             10000.);
            fh_Fib_ToF[ifibcount]->GetYaxis()->SetTitle("ToF / ns");
            fh_Fib_ToF[ifibcount]->GetXaxis()->SetTitle("Fiber number");

            // Not-calibrated position:
            fh_Fib_pos_xy[ifibcount] =
                new TH2F(Form("%s_pos", detName), Form("%s xy position", detName), 500, 0., 500, 1200, -600., 600.);
            fh_Fib_pos_xy[ifibcount]->GetXaxis()->SetTitle("x position (fiber)");
            fh_Fib_pos_xy[ifibcount]->GetYaxis()->SetTitle("y position (dt)");

            // hit fiber number vs. event number:
            fh_Fib_vs_Events[ifibcount] = new TH2F(Form("%s_fib_vs_event", detName),
                                                   Form("%s Fiber # vs. Event #", detName),
                                                   1000,
                                                   0,
                                                   5e6,
                                                   1100,
                                                   0.,
                                                   1100.);
            fh_Fib_vs_Events[ifibcount]->GetYaxis()->SetTitle("Fiber number");
            fh_Fib_vs_Events[ifibcount]->GetXaxis()->SetTitle("Event number");

            // ToT SAPMT:
            fh_ToT_single_Fib[ifibcount] = new TH2F(
                Form("%s_tot_single", detName), Form("%s ToT of single PMTs", detName), 10, 0., 10., 200, 0., 20.);
            fh_ToT_single_Fib[ifibcount]->GetXaxis()->SetTitle("Fiber number");
            fh_ToT_single_Fib[ifibcount]->GetYaxis()->SetTitle("ToT / ns");

            FibCanvas[ifibcount]->Divide(3, 4);
            FibCanvas[ifibcount]->cd(1);
            fh_channels_Fib[ifibcount]->Draw();
            FibCanvas[ifibcount]->cd(2);
            gPad->SetLogz();
            fh_multihit_m_Fib[ifibcount]->Draw("colz");
            FibCanvas[ifibcount]->cd(3);
            gPad->SetLogz();
            fh_multihit_s_Fib[ifibcount]->Draw("colz");
            FibCanvas[ifibcount]->cd(4);
            fh_fibers_Fib[ifibcount]->Draw();
            FibCanvas[ifibcount]->cd(5);
            fh_mult_Fib[ifibcount]->Draw();
            FibCanvas[ifibcount]->cd(6);
            gPad->SetLogz();
            fh_ToT_m_Fib[ifibcount]->Draw("colz");
            FibCanvas[ifibcount]->cd(7);
            gPad->SetLogz();
            fh_ToT_s_Fib[ifibcount]->Draw("colz");
            // FibCanvas[ifibcount]->cd(8);gPad->SetLogz();
            // fh_time_Fib[ifibcount]->Draw("colz");
            FibCanvas[ifibcount]->cd(8);
            gPad->SetLogz();
            fh_Fib_ToF[ifibcount]->Draw("colz");
            FibCanvas[ifibcount]->cd(9);
            fh_Fib_pos_xy[ifibcount]->Draw("colz");
            FibCanvas[ifibcount]->cd(10);
            fh_Fib_vs_Events[ifibcount]->Draw("colz");
            FibCanvas[ifibcount]->cd(11);
            fh_ToT_single_Fib[ifibcount]->Draw("colz");
            FibCanvas[ifibcount]->cd(12);
            fh_channels_single_Fib[ifibcount]->Draw();

            FibCanvas[ifibcount]->cd(0);
            run->AddObject(FibCanvas[ifibcount]);

        } // end if(Mapped)

    } // end for(ifibcount)
    if (NOF_FIB_DET > 0)
        run->GetHttpServer()->RegisterCommand("Reset_Fibers", Form("/Tasks/%s/->Reset_FIBERS_Histo()", GetName()));

    //---------------------------------------------------------------------------------------
    // Ptof detector

    if (fCalItems.at(DET_PTOF))
    {

        TCanvas* cPtof_plane = new TCanvas("Ptof_plane", "Ptof plane", 10, 10, 510, 510);
        cPtof_plane->Divide(1, 2);

        fh_ptof_channels = new TH1F("Ptof_channels", "Ptof channels", 65, 0., 65.);
        cPtof_plane->cd(1);
        fh_ptof_channels->Draw();

        fh_ptof_channels_cut = new TH1F("Ptof_channels_cut", "Ptof channels with cut", 65, 0., 65.);
        cPtof_plane->cd(2);
        fh_ptof_channels_cut->Draw();
        run->AddObject(cPtof_plane);

        TCanvas* cPtof_test = new TCanvas("Ptof_test", "Ptof test", 10, 10, 510, 510);
        cPtof_test->Divide(1, 1);

        fh_ptof_test1 = new TH1F("Ptof_test1", "Ptof test1", 1000, 0., 100.);
        fh_ptof_test2 = new TH1F("Ptof_test2", "Ptof test2", 1000, 0., 100.);
        cPtof_test->cd(1);
        fh_ptof_test1->Draw();
        //		cPtof_test->cd(2);
        //		fh_ptof_test2->Draw();

        run->AddObject(cPtof_test);

        TCanvas* cPtof1;
        cPtof1 = new TCanvas("Ptof1", "Ptof1", 10, 10, 500, 500);
        cPtof1->Divide(3, 3);
        for (Int_t j = 1; j < 10; j++)
        {
            char strName[255];
            sprintf(strName, "ToT_Bar_%d_PM_1", j);
            fh_ptof_TotPm1[j] = new TH1F(strName, "", 1000, 0., 100.);
            cPtof1->cd(j);
            fh_ptof_TotPm1[j]->Draw();
        }
        run->AddObject(cPtof1);

        TCanvas* cPtof2;
        cPtof2 = new TCanvas("Ptof2", "Ptof2", 10, 10, 500, 500);
        cPtof2->Divide(3, 3);
        for (Int_t j = 1; j < 10; j++)
        {
            char strName[255];
            sprintf(strName, "ToT_Bar_%d_PM_2", j);
            fh_ptof_TotPm2[j] = new TH1F(strName, "", 1000, 0., 100.);
            cPtof2->cd(j);
            fh_ptof_TotPm2[j]->Draw();
        }
        run->AddObject(cPtof2);
    }

    //---------------------------------------------------------------------------------------------------
    // TofD detector

    if (fMappedItems.at(DET_TOFD))
    {
        TCanvas* cTofd_planes = new TCanvas("TOFD_planes", "TOFD planes", 10, 10, 1100, 1000);
        cTofd_planes->Divide(5, 4);

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
            fh_tofd_TotPm[j] = new TH2F(strName3, strName4, 90, -45, 45, 3000, 0., 300.);
            fh_tofd_TotPm[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_TotPm[j]->GetYaxis()->SetTitle("ToT / ns");

            char strName7[255];
            sprintf(strName7, "tofd_ToF_plane_%d", j + 1);
            char strName8[255];
            sprintf(strName8, "TOFD-LOS ToF plane %d", j + 1);
            fh_tofd_ToF[j] = new TH2F(strName7, strName8, 15, 0, 15, 1000, 0, 100);
            fh_tofd_ToF[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_ToF[j]->GetYaxis()->SetTitle("ToF / ns");

            char strName9[255];
            sprintf(strName9, "tofd_multihit_plane_%d", j + 1);
            char strName10[255];
            sprintf(strName10, "Tofd multihit plane %d", j + 1);
            fh_tofd_multihit[j] = new TH2F(strName9, strName10, 45, 0., 45., 10, 0, 10);
            fh_tofd_multihit[j]->GetXaxis()->SetTitle("Bar number");
            fh_tofd_multihit[j]->GetYaxis()->SetTitle("Multihit");

            if (j < 3)
            {
                int jk = j + 1;
                char strName5[255];
                sprintf(strName5, "tofd_ToF_plane_%d_%d", jk, jk + 1);
                char strName6[255];
                sprintf(strName6, "Tofd ToF plane %d and %d ", jk, jk + 1);
                fh_tofd_dt[j] = new TH2F(strName5, strName6, 50, 0, 50, 800, -40., 40);
                fh_tofd_dt[j]->GetXaxis()->SetTitle("Bar number");
                fh_tofd_dt[j]->GetYaxis()->SetTitle("dt / ns");
            }
        }

        cTofd_planes->cd(1);
        fh_tofd_channels[0]->Draw();
        cTofd_planes->cd(2);
        gPad->SetLogz();
        fh_tofd_TotPm[0]->Draw("colz");
        cTofd_planes->cd(3);
        gPad->SetLogz();
        fh_tofd_ToF[0]->Draw("colz");
        cTofd_planes->cd(4);
        gPad->SetLogz();
        fh_tofd_multihit[0]->Draw("colz");

        cTofd_planes->cd(6);
        fh_tofd_channels[1]->Draw();
        cTofd_planes->cd(7);
        gPad->SetLogz();
        fh_tofd_TotPm[1]->Draw("colz");
        cTofd_planes->cd(8);
        gPad->SetLogz();
        fh_tofd_ToF[1]->Draw("colz");
        cTofd_planes->cd(9);
        gPad->SetLogz();
        fh_tofd_multihit[1]->Draw("colz");
        cTofd_planes->cd(10);
        gPad->SetLogz();
        fh_tofd_dt[0]->Draw("colz");

        cTofd_planes->cd(11);
        fh_tofd_channels[2]->Draw();
        cTofd_planes->cd(12);
        gPad->SetLogz();
        fh_tofd_TotPm[2]->Draw("colz");
        cTofd_planes->cd(13);
        gPad->SetLogz();
        fh_tofd_ToF[2]->Draw("colz");
        cTofd_planes->cd(14);
        gPad->SetLogz();
        fh_tofd_multihit[2]->Draw("colz");
        cTofd_planes->cd(15);
        gPad->SetLogz();
        fh_tofd_dt[1]->Draw("colz");

        cTofd_planes->cd(16);
        fh_tofd_channels[3]->Draw();
        cTofd_planes->cd(17);
        gPad->SetLogz();
        fh_tofd_TotPm[3]->Draw("colz");
        cTofd_planes->cd(18);
        gPad->SetLogz();
        fh_tofd_ToF[3]->Draw("colz");
        cTofd_planes->cd(19);
        gPad->SetLogz();
        fh_tofd_multihit[3]->Draw("colz");
        cTofd_planes->cd(20);
        gPad->SetLogz();
        fh_tofd_dt[2]->Draw("colz");

        cTofd_planes->cd(0);
        run->AddObject(cTofd_planes);

        run->GetHttpServer()->RegisterCommand("Reset_TOFD", Form("/Tasks/%s/->Reset_TOFD_Histo()", GetName()));
    }

    // -------------------------------------------------------------------------

    //------------------------------------------------------------------------
    // PSPX detector
    // get access to data
    if (fMappedItems.at(DET_PSPX) || fCalItems.at(DET_PSPX) || fHitItems.at(DET_PSPX))
    {
        TCanvas* cPspx_comp = new TCanvas("Pspx_comp", "Pspx Comparison", 10, 10, 1100, 1000);
        cPspx_comp->Divide(N_PSPX, 3);

        Int_t Emax = 1000000;

        if (fMappedItems.at(DET_PSPX))
        {
            // LOG(info) << "Init MappedPspx";

            for (UInt_t i = 0; i < N_PSPX; i++)
            {
                if (i % 2 == 0)
                { // even numbers = read out with energy filter
                    fh_pspx_channel_x[i] =
                        new TH1F(Form("pspx_%d_energyreadout_channel_x", i / 2),
                                 Form("Pspx %d Energy Readout: x Channel;x Position / Channel;Counts", i / 2 + 1),
                                 2 * N_STRIPS_PSPX,
                                 1,
                                 2 * N_STRIPS_PSPX + 1);
                    fh_pspx_channel_y[i] =
                        new TH1F(Form("pspx_%d_energyreadout_channel_y", i / 2),
                                 Form("Pspx %d Energy Readout: y Channel;y Position / Channel;Counts", i / 2 + 1),
                                 2 * N_STRIPS_PSPX,
                                 2 * N_STRIPS_PSPX + 1,
                                 4 * N_STRIPS_PSPX + 1);

                    fh_pspx_multiplicity_x[i] =
                        new TH1F(Form("pspx_%d_energyreadout_multiplicity_x", i / 2),
                                 Form("Pspx %d Energy Readout: x Multiplicity;Multiplicity;Counts", i / 2 + 1),
                                 10,
                                 0,
                                 10);
                    fh_pspx_multiplicity_y[i] =
                        new TH1F(Form("pspx_%d_energyreadout_multiplicity_y", i / 2),
                                 Form("Pspx %d Energy Readout: y Multiplicity;Multiplicity;Counts", i / 2 + 1),
                                 10,
                                 0,
                                 10);
                    fh_pspx_strips_position[i] = new TH2F(
                        Form("pspx_%d_energyreadout_strips", i),
                        Form("Pspx %d Energy Readout Position;x Position / Strips; y Position / Strips", i / 2 + 1),
                        N_STRIPS_PSPX,
                        1,
                        N_STRIPS_PSPX + 1,
                        N_STRIPS_PSPX,
                        N_STRIPS_PSPX + 1,
                        2 * N_STRIPS_PSPX + 1);
                }
                else
                { // odd numbers = read out with position filter
                    fh_pspx_channel_x[i] =
                        new TH1F(Form("pspx_%d_positionreadout_channel_x", i / 2),
                                 Form("Pspx %d Position Readout: x Channel;x Position / Channel;Counts", i / 2 + 1),
                                 2 * N_STRIPS_PSPX,
                                 1,
                                 2 * N_STRIPS_PSPX + 1);
                    fh_pspx_channel_y[i] =
                        new TH1F(Form("pspx_%d_positionreadout_channel_y", i / 2),
                                 Form("Pspx %d Position Readout: y Channel;y Position / Channel;Counts", i / 2 + 1),
                                 2 * N_STRIPS_PSPX,
                                 2 * N_STRIPS_PSPX + 1,
                                 4 * N_STRIPS_PSPX + 1);

                    fh_pspx_multiplicity_x[i] =
                        new TH1F(Form("pspx_%d_positionreadout_multiplicity_x", i / 2),
                                 Form("Pspx %d Position Readout: x Multiplicity;Multiplicity;Counts", i / 2 + 1),
                                 10,
                                 0,
                                 10);
                    fh_pspx_multiplicity_y[i] =
                        new TH1F(Form("pspx_%d_positionreadout_multiplicity_y", i / 2),
                                 Form("Pspx %d Position Readout: y Multiplicity;Multiplicity;Counts", i / 2 + 1),
                                 10,
                                 0,
                                 10);
                    fh_pspx_strips_position[i] = new TH2F(
                        Form("pspx_%d_positionreadout_strips", i),
                        Form("Pspx %d Position Readout Position;x Position / Strips; y Position / Strips", i / 2 + 1),
                        N_STRIPS_PSPX,
                        1,
                        N_STRIPS_PSPX + 1,
                        N_STRIPS_PSPX,
                        N_STRIPS_PSPX + 1,
                        2 * N_STRIPS_PSPX + 1);
                }
            }

            TCanvas* cPspx_strips = new TCanvas("Pspx_strips", "Pspx Strips", 10, 10, 1100, 1000);
            cPspx_strips->Divide(N_PSPX, 2);

            for (UInt_t i = 0; i < N_PSPX; i++)
            {
                cPspx_strips->cd(i + 1);
                fh_pspx_channel_x[i]->Draw();

                cPspx_strips->cd(i + 1 + N_PSPX);
                fh_pspx_channel_y[i]->Draw();
            }

            run->AddObject(cPspx_strips);

            TCanvas* cPspx_multiplicity = new TCanvas("Pspx_multiplicity", "Pspx Multiplicity", 10, 10, 1100, 1000);
            cPspx_multiplicity->Divide(N_PSPX, 2);

            for (UInt_t i = 0; i < N_PSPX; i++)
            {
                cPspx_multiplicity->cd(i + 1);
                fh_pspx_multiplicity_x[i]->Draw();

                cPspx_multiplicity->cd(i + 1 + N_PSPX);
                fh_pspx_multiplicity_y[i]->Draw();
            }

            run->AddObject(cPspx_multiplicity);

            // Fill cPspx_comp with Mapped level data
            for (UInt_t i = 0; i < N_PSPX; i++)
            {
                cPspx_comp->cd(i + 1);
                fh_pspx_strips_position[i]->Draw("colz");

                // Remove the current axis
                fh_pspx_strips_position[i]->GetYaxis()->SetLabelOffset(999);
                fh_pspx_strips_position[i]->GetYaxis()->SetTickLength(0);

                // Redraw the new inverted axis
                gPad->Update();
                TGaxis* newaxis = new TGaxis(gPad->GetUxmin(),
                                             gPad->GetUymax(),
                                             gPad->GetUxmin() - 0.001,
                                             gPad->GetUymin(),
                                             fh_pspx_strips_position[i]->GetYaxis()->GetXmin(),
                                             fh_pspx_strips_position[i]->GetYaxis()->GetXmax(),
                                             510,
                                             "+");
                newaxis->SetLabelOffset(0.003);
                newaxis->SetLabelSize(0.03);
                newaxis->SetTickLength(0.025);
                newaxis->Draw();
            }
        }

        if (fCalItems.at(DET_PSPX))
        {
            UInt_t nbins = 200;

            for (UInt_t i = 0; i < N_PSPX; i++)
            {
                if (i % 2 == 0)
                {
                    fh_pspx_cal_energy_frontback[i] = new TH2F(
                        Form("pspx_%d_energy_frontback", i / 2),
                        Form("Pspx %d Energy Back vs Front;Front Energy / arb.u.;Back Energy / arb.u.", i / 2 + 1),
                        nbins,
                        0,
                        Emax,
                        nbins,
                        -Emax,
                        0);
                }
                else
                {
                    fh_pspx_cal_energy_frontback[i] = new TH2F(
                        Form("pspx_%d_position_frontback", i / 2),
                        Form("Pspx %d Energy Back vs Front;Front Energy / arb.u.;Back Energy / arb.u.", i / 2 + 1),
                        nbins,
                        0,
                        Emax / 4,
                        nbins,
                        -Emax / 4,
                        0);
                }
            }

            /*fh_pspx_cal_energy_corr[i] =
                        new TH2F("pspx_%d_energy_corr",
                                 "Pspx %d Energy PSP2 vs. PSP3;PSP2 Energy / arb.u.;PSP3 Energy / arb.u.",
                                 nbins,
                                 0,
                                 Emax,
                                 nbins,
                                 0,
                                 Emax);
        */

            // Fill cPspx_comp with Cal level data
            for (UInt_t i = 0; i < N_PSPX; i++)
            {
                cPspx_comp->cd(i + 1 + N_PSPX); // i*2
                fh_pspx_cal_energy_frontback[i]->Draw("colz");
            }
        }

        if (fHitItems.at(DET_PSPX))
        {
            UInt_t nbins = 100;
            Double_t length = 10; // detector length, range of histogram

            for (UInt_t i = 0; i < (N_PSPX + 1) / 2; i++)
            {
                fh_pspx_hit_position[i] = new TH2F(Form("pspx_%d_position_cm", i),
                                                   Form("Pspx %d Position;x Position / cm; y Position / cm", i + 1),
                                                   nbins,
                                                   -length / 4.,
                                                   length / 4.,
                                                   nbins,
                                                   -length / 4.,
                                                   length / 4.);

                fh_pspx_hit_energy[i] = new TH1F(Form("pspx_%d_energy", i),
                                                 Form("Pspx %d Energy;Energy / arb.u.; Counts", i + 1),
                                                 nbins,
                                                 -Emax,
                                                 Emax);
            }

            TCanvas* cPspx_hit = new TCanvas("Pspx_hit", "Pspx Hit", 10, 10, 1100, 1000);
            cPspx_hit->Divide((N_PSPX + 1) / 2, 2);

            for (UInt_t i = 0; i < (N_PSPX + 1) / 2; i++)
            {
                cPspx_hit->cd(i + 1);
                fh_pspx_hit_position[i]->Draw("colz");

                cPspx_hit->cd(i + 1 + (N_PSPX + 1) / 2);
                fh_pspx_hit_energy[i]->Draw();
            }

            run->AddObject(cPspx_hit);

            // Fill cPspx_comp with Hit level data-((channel_y[i][0] + 1) / 2) + 3 * N_STRIPS_PSPX + 1)
            for (UInt_t i = 0; i < (N_PSPX + 1) / 2; i++)
            {
                cPspx_comp->cd(i * 2 + 2 * N_PSPX + 1); // supposed to be +2 if
                                                        // energy and position
                                                        // readout is used
                fh_pspx_hit_position[i]->Draw("colz");
            }
        }

        run->AddObject(cPspx_comp);

        run->GetHttpServer()->RegisterCommand("Reset_PSPX", Form("/Tasks/%s/->Reset_PSPX_Histo()", GetName()));
    }
    // -------------------------------------------------------------------------

    return kSUCCESS;
}

void R3BOnlineSpectraDec2019::Reset_LOS_Histo()
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
        }
        if (fNofLosDetectors == 2)
            fh_los_ToT_vs_ToT->Reset();
    }
}

void R3BOnlineSpectraDec2019::Reset_ROLU_Histo()
{
    fh_rolu_channels->Reset();
    fh_rolu_tot->Reset();
}

void R3BOnlineSpectraDec2019::Reset_BMON_Histo()
{
    fhTrigger->Reset();
    fhTpat->Reset();
    fh_spill_length->Reset();
    fh_spill_times->Reset();
    fh_spill_times_pois->Reset();
    fh_spill_times_daq->Reset();
    fh_spill_times_daq_first->Reset();
    fh_rate->Reset();
    fh_particle->Reset();
    fh_particle_pois->Reset();
    fh_frequencies->Reset();
    fh_TOFDOR->Reset();
    fh_TOFDOR_spill->Reset();
    for (Int_t iloscount = 0; iloscount < fNofLosDetectors; iloscount++)
    {
        fh_los_dt_hits[iloscount]->Reset();
        fh_dt_hits[iloscount]->Reset();
    }
}

void R3BOnlineSpectraDec2019::Reset_SCI8_Histo()
{
    fh_sci8_channels->Reset();
    fh_sci8_tres_MCFD->Reset();
    fh_sci8_tres_TAMEX->Reset();
    fh_sci8_tot->Reset();
    fh_sci8_tot_mean->Reset();
    fh_sci8_dt_hits->Reset();
    fh_sci8_dt_hits_l->Reset();
    fh_sci8_dt_hits_t->Reset();
    fh_sci8_multihit->Reset();
    fh_sci8_multihitVFTX->Reset();
    fh_sci8_multihitLEAD->Reset();
    fh_sci8_multihitTRAI->Reset();
    fh_tof_sci8->Reset();
}
void R3BOnlineSpectraDec2019::Reset_TOFD_Histo()
{
    for (int i = 0; i < N_PLANE_MAX_TOFD; i++)
    {
        fh_tofd_channels[i]->Reset();
        fh_tofd_multihit[i]->Reset();
        fh_tofd_ToF[i]->Reset();
        fh_tofd_TotPm[i]->Reset();
    }
    fh_tofd_dt[0]->Reset();
    fh_tofd_dt[1]->Reset();
    fh_tofd_dt[2]->Reset();
}
void R3BOnlineSpectraDec2019::Reset_FIBERS_Histo()
{
    if (fMappedItems.at(DET_PSPX))
    {
        for (UInt_t i = 0; i < N_PSPX; i++)
        {
            fh_pspx_channel_x[i]->Reset();
            fh_pspx_channel_y[i]->Reset();
            fh_pspx_multiplicity_x[i]->Reset();
            fh_pspx_multiplicity_y[i]->Reset();
            fh_pspx_strips_position[i]->Reset();
        }
    }
    if (fCalItems.at(DET_PSPX))
    {
        for (UInt_t i = 0; i < N_PSPX; i++)
        {
            fh_pspx_cal_energy_frontback[i]->Reset();
        }
    }
    if (fHitItems.at(DET_PSPX))
    {
        for (UInt_t i = 0; i < N_PSPX / 2; i++)
        {
            fh_pspx_hit_position[i]->Reset();
            fh_pspx_hit_energy[i]->Reset();
        }
    }
    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {
        if (fMappedItems.at(DET_FI_FIRST + ifibcount))
        {
            fh_channels_Fib[ifibcount]->Reset();
            fh_multihit_m_Fib[ifibcount]->Reset();
            fh_multihit_s_Fib[ifibcount]->Reset();
            fh_fibers_Fib[ifibcount]->Reset();
            fh_mult_Fib[ifibcount]->Reset();
            fh_ToT_m_Fib[ifibcount]->Reset();
            fh_ToT_s_Fib[ifibcount]->Reset();
            fh_ToT_single_Fib[ifibcount]->Reset();
            fh_Fib_ToF[ifibcount]->Reset();
            fh_Fib_pos_xy[ifibcount]->Reset();
            fh_Fib_vs_Events[ifibcount]->Reset();
            fh_channels_single_Fib[ifibcount]->Reset();
        }
    }
}

void R3BOnlineSpectraDec2019::Exec(Option_t* option)
{
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        // FairLogger::GetLogger()->Fatal(MESSAGE_ORIGIN, "FairRootManager not found");
        LOG(error) << "FairRootManager not found";
        return;
    }
    /*
        time = header->GetTimeStamp(); // /1.6; // divided by 1.6 for stand alone daq with Vulom generating time stamps
        if (time_begin == 0)
            time_begin = time;
        time -= time_begin;

        if (time_start == 0 && time > 0)
        {
            time_start = time;
            fNEvents_start = fNEvents;
        }

        if (header->GetTrigger() == 12)
        {
            spill_on = true;
            spill_on_sampler = true;
        }
    */
    if (fMappedItems.at(DET_SAMPLER))
    {
        // cout << "Samples:\n";
        auto det = fMappedItems.at(DET_SAMPLER);
        Int_t nHitsSamp = det->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHitsSamp; ihit++)
        {
            auto hit = dynamic_cast<R3BSamplerMappedData*>(det->At(ihit));
            // time is in steps of 10 ns
            // is is a 34 bit number, so max 1073741823
            samplerCurr = hit->GetTime();
            if (samplerBegin == 0)
                samplerBegin = samplerCurr;
            samplerCurr -= samplerBegin;

            // cout<<"Time - Sampler: "<<ihit<<"  "<< (double)time - (double)samplerCurr*10. <<endl;
            if (spill_on_sampler && ihit == 0)
            {
                samplerSpill = samplerCurr;
                spill_on_sampler = false;
            }

            long samp = (samplerCurr - samplerSpill); // time in 10 ns
            if (samp < 0)
                samp += 1073741823;

            if ((double)samp / 1e8 > fSpillLength && spillCounter > 0)
            {
                cout << "Missed spill end!" << endl;
                cout << samplerCurr << "  " << samplerPrev << "  " << samplerSpill << "  " << (double)samp / 1e8
                     << '\n';
            }

            fh_spill_times->Fill((double)samp / 1.e8);

            Double_t dt = ((double)(samplerCurr - samplerPrev)) / 100.;

            // cout << "sampler times: " << samplerCurr << "  " << samplerPrev << "  " << dt << endl;

            fh_dt_hits[0]->Fill(dt);
            fh_frequencies->Fill(1. / dt * 1.E3);
            // cout <<"freq " << 1./dt*1.E3 << endl;
            samplerPrev = samplerCurr;
        }
    }

    if (spill_on && time > 0)
    {
        time_spill_start = time; // spill start in nsec
        cout << "Spill start: " << time_spill_start << " " << endl;
        spill_on = false;
        spillCounter++;
    }

    if (header->GetTrigger() == 13)
        spill_off = true;

    if (spill_off)
    {
        time_spill_end = time; // spill end  in nsec
        cout << "Spill stop: " << time_spill_end << " " << endl;
        spill_off = false;

        Int_t pps = fh_spill_times->Integral(1, fh_spill_times->GetSize() - 2);
        Int_t pps_daq = fh_spill_times_daq->Integral(1, fh_spill_times_daq->GetSize() - 2);
        Int_t pps_daq_first = fh_spill_times_daq_first->Integral(1, fh_spill_times_daq_first->GetSize() - 2);
        if (time_spill_start > 0 && pps > 0 && spillCounter > 0)
        {

            cout << "pps: " << pps << "  " << pps_daq << "  " << pps_daq_first << endl;

            Double_t norm = (double)pps * 100.E-6 / fSpillLength; // expected particle per channel

            // calculate the accepted trigger
            Double_t a = pps * fDAQ_dead_time / (fSpillLength * 1e6);
            Double_t R_Poisson = (1 - exp(-a)) / a; // ratio of accepted trigger relativ to incoming particle

            Double_t accepted =
                ((double)pps / fSpillLength) / (1 + ((double)pps / fSpillLength) * fDAQ_dead_time * 1.e-6);

            cout << "Ratio pois: " << R_Poisson << "  " << accepted / (double)pps * fSpillLength << "  "
                 << (double)pps_daq_first / (double)pps << "  " << (double)pps_daq / (double)pps << endl;

            cout << "value " << (double)pps_daq_first / accepted / fSpillLength * 100. << endl;
            fh_rate->Fill(spillCounter, (double)pps_daq_first / accepted / fSpillLength * 100.);

            // now make histogram for random distribution
            for (Int_t i = 1; i < pps; i++)
            {
                Double_t num = (std::rand() / (float)RAND_MAX * fSpillLength);
                fh_spill_times_pois->Fill(num);
            }

            for (Int_t i = 1; i < fh_spill_times->GetSize() - 2; i++)
            {
                Int_t counts = fh_spill_times->GetBinContent(i);
                Int_t counts_pois = fh_spill_times_pois->GetBinContent(i);

                if (fh_spill_times->GetBinCenter(i) < fSpillLength && counts > 0)
                    fh_particle->Fill((double)counts - norm);
                if (fh_spill_times->GetBinCenter(i) < fSpillLength && counts_pois > 0)
                    fh_particle_pois->Fill((double)counts_pois - norm);

                // cout<<"Exp : " << i << " counts " << counts << " norm " << norm << " value " << (double)counts - norm
                // <<endl; cout<<"Pois: " << i << " counts " << counts_pois << " norm " << norm << " value " <<
                // (double)counts_pois - norm <<endl;
            }
        }
        fh_spill_times->Reset();
        fh_spill_times_pois->Reset();
        fh_spill_times_daq->Reset();
        fh_spill_times_daq_first->Reset();
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

    // cout<<"Trigger: "<<header->GetTrigger()<<endl;

    if (fMappedItems.at(DET_BMON))
    {
        fhTrigger->Fill(header->GetTrigger());
    }

    //   check for requested trigger (Todo: should be done globablly / somewhere else)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    if (fMappedItems.at(DET_BMON))
    {
        Int_t tpatbin;
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0)
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
        tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if (tpatvalue == 0)
            return;
    }

    if (fMappedItems.at(DET_BMON))
    {
        unsigned long IC;
        unsigned long SEETRAM, SEETRAM_raw;
        unsigned long TOFDOR;

        Bool_t bmon_read = false;

        auto det = fMappedItems.at(DET_BMON);
        Int_t nHitsbm = det->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nHitsbm; ihit++)
        {
            R3BBeamMonitorMappedData* hit = dynamic_cast<R3BBeamMonitorMappedData*>(det->At(ihit));
            if (!hit)
                continue;

            IC = hit->GetIC();                 // negative values if offset not high enough
            SEETRAM_raw = hit->GetSEETRAM();   // raw counts
            TOFDOR = hit->GetTOFDOR();         // only positive values possible
            SEETRAM = SEETRAM_raw * calib_SEE; // calibrated SEETRAM counts

            if (fNEvents == fNEvents_start)
            {
                see_mem = SEETRAM;
                ic_mem = IC;
                tofdor_mem = TOFDOR;
                time_mem = time_start;
                see_start = SEETRAM;
                ic_start = IC;
                tofdor_start = TOFDOR;
                time_prev_read = time_start;
            }

            if (time > 0)
            {

                // cout<<SEETRAM_raw<<"; "<<see_offset<<"; "<<calib_SEE<<", "<<SEETRAM<<", "<<see_mem<<endl;

                fh_spill_length->Fill((time - time_mem) / 1e9);

                // Spectra below are filled every read_time (secs)
                if (time_to_read == 0 && (time - time_prev_read) >= read_time * 1000000000)
                {
                    time_to_read = time;
                    bmon_read = true;
                }

                if (bmon_read)
                {

                    tdiff = double(time - time_mem) / 1.e9;
                    fNorm = 1.e-3 / (double(time - time_prev_read) / 1.e9); // kHz

                    // TOFDOR:
                    Int_t yTOFDOR = TOFDOR - tofdor_start;
                    fh_TOFDOR->Fill(tdiff, yTOFDOR);
                    fh_TOFDOR_spill->Fill(tdiff, (TOFDOR - tofdor_mem) * fNorm);
                    tofdor_mem = TOFDOR;

                    time_to_read = 0;
                    time_prev_read = time;
                    bmon_read = false;
                }

                if (spectra_clear)
                {
                    fh_spill_length->Reset();
                    fh_TOFDOR_spill->Reset("ICESM");
                    fh_frequencies->Reset("ICESM");
                    fh_TOFDOR->Reset("ICESM");
                    time_mem = time;
                    time_clear = -1.;
                    iclear_count = iclear_count + 1;
                    spectra_clear = false;
                    see_start = SEETRAM;
                    ic_start = IC;
                    tofdor_start = TOFDOR;
                }
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

    Int_t nParts;
    if (fCalItems.at(DET_ROLU))
    {

        auto det = fCalItems.at(DET_ROLU);
        nParts = det->GetEntriesFast();

        if (nParts < 1)
            return;

        Int_t iDet = 0;
        Double_t timeRolu_L[nParts][2][4];
        Double_t timeRolu_T[nParts][2][4];
        Double_t totRolu[nParts][2][4];
        for (int i = 0; i < nParts; i++)
        {
            for (int j = 0; j < 2; j++)
            {
                for (int k = 0; k < 4; k++)
                {
                    timeRolu_L[i][j][k] = 0. / 0.;
                    timeRolu_T[i][j][k] = 0. / 0.;
                    totRolu[i][j][k] = 0. / 0.;
                }
            }
        }

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

                //  cout<<"ROLU: "<<iPart<<", "<<iCha<<"; "<<timeRolu_L[iPart][iCha]<<",
                //  "<<timeRolu_T[iPart][iCha]<<endl;

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
    Double_t time_V[fNofLosDetectors][32][8]; // [det][multihit][pm]
    Double_t time_L[fNofLosDetectors][32][8];
    Double_t time_T[fNofLosDetectors][32][8];

    Double_t time_MTDC[32][8] = { { 0.0 / 0.0 } };
    Double_t LosTresMTDC[32] = { 0.0 / 0.0 };

    Double_t timeLosV[fNofLosDetectors][32];
    Double_t LosTresV[fNofLosDetectors][32];
    Double_t timeLosT[fNofLosDetectors][32];
    Double_t LosTresT[fNofLosDetectors][32];
    Double_t timeLos[fNofLosDetectors][32];
    Double_t totsum[fNofLosDetectors][32];
    Double_t tot[fNofLosDetectors][32][8];
    Double_t xT_cm[fNofLosDetectors][32];
    Double_t yT_cm[fNofLosDetectors][32];
    Double_t xToT_cm[fNofLosDetectors][32];
    Double_t yToT_cm[fNofLosDetectors][32];
    Double_t xV_cm[fNofLosDetectors][32];
    Double_t yV_cm[fNofLosDetectors][32];
    Double_t time_V_temp[32][8] = { { 0. / 0. } };

    for (int i = 0; i < fNofLosDetectors; i++)
    {
        for (int j = 0; j < 32; j++)
        {
            timeLosV[i][j] = 0.;
            LosTresV[i][j] = 0. / 0.;
            timeLosT[i][j] = 0.;
            LosTresT[i][j] = 0. / 0.;
            timeLos[i][j] = 0.;
            totsum[i][j] = 0.;
            xT_cm[i][j] = 0. / 0.;
            yT_cm[i][j] = 0. / 0.;
            xToT_cm[i][j] = 0. / 0.;
            yToT_cm[i][j] = 0. / 0.;
            xV_cm[i][j] = 0. / 0.;
            yV_cm[i][j] = 0. / 0.;
            for (int k = 0; k < 8; k++)
            {
                time_V[i][j][k] = 0. / 0.;
                time_L[i][j][k] = 0. / 0.;
                time_T[i][j][k] = 0. / 0.;
                tot[i][j][k] = 0. / 0.;
            }
        }
    }

    Int_t Multip;

    if (fMappedItems.at(DET_LOS))
    {
        auto det = fMappedItems.at(DET_LOS);
        Int_t nHits = det->GetEntriesFast();

        Multip = nHits;
        // cout<<"LOS: nHits"<<nHits<<endl;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BLosMappedData* hit = dynamic_cast<R3BLosMappedData*>(det->At(ihit));
            if (!hit)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hit->GetDetector(); // 1..
            Int_t iCha = hit->GetChannel();  // 1..
            Int_t iTyp = hit->GetType();     // 0,1,2,3

            // cout<<"LOS mapped: "<<ihit<<"; "<<iDet<<"; "<<iCha<<", "<<iTyp<<endl;

            if (iTyp < 3)
                fh_los_channels[iDet - 1]->Fill(iCha); // exclude MTDC data
        }
    }

    Int_t nPart;
    Int_t nPartLos[fNofLosDetectors];
    for (int i = 0; i < fNofLosDetectors; i++)
    {
        nPartLos[i] = 0;
    }

    if (fCalItems.at(DET_LOS))
    {
        Bool_t LOSID = false;

        auto det = fCalItems.at(DET_LOS);
        nPart = det->GetEntriesFast();

        Int_t iDet = 0;
        Double_t time_V_LOS1[32][8] = { { 0. } };
        Double_t time_V_LOS2[32][8] = { { 0. } };

        for (Int_t iPart = 0; iPart < nPart; iPart++)
        {
            /*
             * nPart is the number of particle passing through LOS detector in one event
             */
            R3BLosCalData* calData = dynamic_cast<R3BLosCalData*>(det->At(iPart));
            iDet = calData->GetDetector();
            Int_t nVftx = 0;
            Int_t nTamexL = 0;
            Int_t nTamexT = 0;

            for (Int_t iCha = 0; iCha < 8; iCha++)
            {
                if (iDet == 1)
                {
                    time_V_LOS1[nPartLos[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeV_ns(iCha))))
                    { // VFTX
                        time_V_LOS1[nPartLos[iDet - 1]][iCha] = calData->GetTimeV_ns(iCha);
                        nVftx++;
                    }
                    // cout<<"time_V_LOS1: "<<calData->GetTimeV_ns(iCha)<<endl;
                    time_L[iDet - 1][nPartLos[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeL_ns(iCha))))
                    { // TAMEX leading
                        time_L[iDet - 1][nPartLos[iDet - 1]][iCha] = calData->GetTimeL_ns(iCha);
                        nTamexL++;
                    }

                    time_T[iDet - 1][nPartLos[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeT_ns(iCha))))
                    { // TAMEX trailing
                        time_T[iDet - 1][nPartLos[iDet - 1]][iCha] = calData->GetTimeT_ns(iCha);
                        nTamexT++;
                    }
                    // cout<<"time Tamex trailing LOS1: "<<calData->GetTimeT_ns(iCha)<<endl;

                    time_MTDC[nPartLos[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeM_ns(iCha))))
                    { // MTDC
                        time_MTDC[nPartLos[iDet - 1]][iCha] = calData->GetTimeM_ns(iCha);
                    }
                }
                if (iDet == 2)
                {
                    time_V_LOS2[nPartLos[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeV_ns(iCha))))
                    { // VFTX
                        time_V_LOS2[nPartLos[iDet - 1]][iCha] = calData->GetTimeV_ns(iCha);
                        nVftx++;
                    }
                    // cout<<"time_V_LOS2: "<<calData->GetTimeV_ns(iCha)<<endl;
                    time_L[iDet - 1][nPartLos[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeL_ns(iCha))))
                    { // TAMEX leading
                        time_L[iDet - 1][nPartLos[iDet - 1]][iCha] = calData->GetTimeL_ns(iCha);
                        nTamexL++;
                    }

                    time_T[iDet - 1][nPartLos[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeT_ns(iCha))))
                    { // TAMEX trailing
                        time_T[iDet - 1][nPartLos[iDet - 1]][iCha] = calData->GetTimeT_ns(iCha);
                        nTamexT++;
                    }

                    time_MTDC[nPartLos[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeM_ns(iCha))))
                    { // MTDC
                        time_MTDC[nPartLos[iDet - 1]][iCha] = calData->GetTimeM_ns(iCha);
                    }
                }
            }
            // cout<<"n: "<<nVftx<<"  "<<nTamexL<<"  "<<nTamexT<<endl;
            //            if (nVftx == 8 && nTamexL == 8 && nTamexT == 8)
            if (nVftx > 6 && nTamexL > 6 && nTamexT > 6)
                nPartLos[iDet - 1]++;

            if (!calData)
            {
                cout << "LOS !calData" << endl;
                continue; // can this happen?
            }
        }

        // Sorting VFTX data:
        // cout<< "nhits: "<<nPart<<" los1: "<<nPartLos[0]<<" los2: "<<nPartLos[1]<<endl;

        // detector 1
        if (nPartLos[0] > 0)
        {
            std::qsort(time_V_LOS1,
                       nPartLos[0],
                       sizeof(*time_V_LOS1),
                       [](const void* arg1, const void* arg2) -> int
                       {
                           double const* lhs = static_cast<double const*>(arg1);
                           double const* rhs = static_cast<double const*>(arg2);

                           return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
                       });
            for (Int_t iPart = 0; iPart < nPartLos[0]; iPart++)
            {
                for (int ipm = 0; ipm < 8; ipm++)
                {
                    time_V[0][iPart][ipm] = time_V_LOS1[iPart][ipm];
                    // cout<<"Test "<<ipm <<"  "<<iPart<<"  "<< time_V_LOS1[iPart][ipm]<<endl;
                }
            }
        }

        // detector 2
        if (nPartLos[1] > 0)
        {
            std::qsort(time_V_LOS2,
                       nPartLos[1],
                       sizeof(*time_V_LOS2),
                       [](const void* arg1, const void* arg2) -> int
                       {
                           double const* lhs = static_cast<double const*>(arg1);
                           double const* rhs = static_cast<double const*>(arg2);

                           return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
                       });
            for (Int_t iPart = 0; iPart < nPartLos[1]; iPart++)
            {
                for (int ipm = 0; ipm < 8; ipm++)
                {
                    time_V[1][iPart][ipm] = time_V_LOS2[iPart][ipm];
                    // cout<<"Test "<<ipm <<"  "<<iPart<<"  "<< time_V_LOS2[iPart][ipm]<<endl;
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
            for (Int_t iPart = 0; iPart < nPartLos[iDet - 1]; iPart++)
            {
                Bool_t iLOSTypeMCFD = false;
                Bool_t iLOSTypeTAMEX = false;
                Bool_t iLOSType = false;

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
                    iLOSType = true;

                if (1 == 1)
                {

                    if (1 == 1) //(iLOSType)
                    {
                        int nPMT = 0;
                        int nPMV = 0;

                        for (int ipm = 0; ipm < 8; ipm++)
                        {

                            tot[iDet - 1][iPart][ipm] = 0. / 0.;

                            if (time_T[iDet - 1][iPart][ipm] > 0. && time_L[iDet - 1][iPart][ipm] > 0. &&
                                !(IS_NAN(time_T[iDet - 1][iPart][ipm])) && !(IS_NAN(time_L[iDet - 1][iPart][ipm])))
                            {
                                while (time_T[iDet - 1][iPart][ipm] - time_L[iDet - 1][iPart][ipm] <= 0.)
                                {
                                    time_T[iDet - 1][iPart][ipm] = time_T[iDet - 1][iPart][ipm] + 2048. * fClockFreq;
                                }

                                nPMT = nPMT + 1;
                                tot[iDet - 1][iPart][ipm] = time_T[iDet - 1][iPart][ipm] - time_L[iDet - 1][iPart][ipm];
                            }

                            if (tot[iDet - 1][iPart][ipm] != 0. && !(IS_NAN(tot[iDet - 1][iPart][ipm])))
                                totsum[iDet - 1][iPart] += tot[iDet - 1][iPart][ipm];

                            if (time_L[iDet - 1][iPart][ipm] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][ipm])))
                                timeLosT[iDet - 1][iPart] += time_L[iDet - 1][iPart][ipm];

                            /*
                                                    if(iDet == 2){
                                                        // Push VFTX times of LOS2
                                                        if( timeLos[iDet-2][iPart] > 0. &&
                               !(IS_NAN(timeLos[iDet-2][iPart])) > 0. && time_V[iDet-1][iPart][ipm] > 0. &&
                               !(IS_NAN(time_V[iDet-1][iPart][ipm]))){

                                                            while(time_V[iDet-1][iPart][ipm]-timeLos[iDet-2][iPart] <
                               8192.*5.){ time_V[iDet-1][iPart][ipm] =   time_V[iDet-1][iPart][ipm] + 8192.*5.;
                                                            }
                                                        }
                                                    }
                            */

                            // Calculate detector time
                            if (time_V[iDet - 1][iPart][ipm] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][ipm])))
                            {
                                timeLosV[iDet - 1][iPart] += time_V[iDet - 1][iPart][ipm];
                                nPMV = nPMV + 1;
                            }
                        }

                        totsum[iDet - 1][iPart] = totsum[iDet - 1][iPart] / nPMT;

                        timeLosV[iDet - 1][iPart] = timeLosV[iDet - 1][iPart] / nPMV;

                        timeLosT[iDet - 1][iPart] = timeLosT[iDet - 1][iPart] / nPMT;

                        timeLos[iDet - 1][iPart] = timeLosV[iDet - 1][iPart];

                        // cout<<"LOS: "<<iDet<<", "<<iPart<<"; "<<timeLos[iDet-1][iPart]<<endl;
                        // if(iDet == 1) cout<<"Time - Los1: "<<iPart<<"  "<<time - (int)timeLos[iDet-1][iPart]<<endl;

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

                        //   LosTresMTDC[iPart] =
                        //   ((time_MTDC[iPart][0]+time_MTDC[iPart][2]+time_MTDC[iPart][4]+time_MTDC[iPart][6]) -
                        //      (time_MTDC[iPart][1]+time_MTDC[iPart][3]+time_MTDC[iPart][5]+time_MTDC[iPart][7]))/4.;

                        LosTresMTDC[iPart] = ((time_MTDC[iPart][1] + time_MTDC[iPart][3]) -
                                              (time_MTDC[iPart][5] + time_MTDC[iPart][7])) /
                                             4.;
                        // LosTresMTDC[iPart] = (time_MTDC[iPart][1] + time_MTDC[iPart][7])/2.;

                        // Position from tamex:
                        xT_cm[iDet - 1][iPart] = ((time_L[iDet - 1][iPart][5] + time_L[iDet - 1][iPart][6]) / 2. -
                                                  (time_L[iDet - 1][iPart][1] + time_L[iDet - 1][iPart][2]) / 2.) *
                                                 (-1.);
                        yT_cm[iDet - 1][iPart] = ((time_L[iDet - 1][iPart][7] + time_L[iDet - 1][iPart][0]) / 2. -
                                                  (time_L[iDet - 1][iPart][3] + time_L[iDet - 1][iPart][4]) / 2.) *
                                                 (-1.);
                        xT_cm[iDet - 1][iPart] =
                            (xT_cm[iDet - 1][iPart] - flosOffsetXT[iDet - 1]) * flosVeffXT[iDet - 1];
                        yT_cm[iDet - 1][iPart] =
                            (yT_cm[iDet - 1][iPart] - flosOffsetYT[iDet - 1]) * flosVeffYT[iDet - 1];

                        // Position from VFTX:
                        xV_cm[iDet - 1][iPart] = ((time_V[iDet - 1][iPart][5] + time_V[iDet - 1][iPart][6]) / 2. -
                                                  (time_V[iDet - 1][iPart][1] + time_V[iDet - 1][iPart][2]) / 2.) *
                                                 (-1.);
                        yV_cm[iDet - 1][iPart] = ((time_V[iDet - 1][iPart][7] + time_V[iDet - 1][iPart][0]) / 2. -
                                                  (time_V[iDet - 1][iPart][3] + time_V[iDet - 1][iPart][4]) / 2.) *
                                                 (-1.);
                        xV_cm[iDet - 1][iPart] =
                            (xV_cm[iDet - 1][iPart] - flosOffsetXV[iDet - 1]) * flosVeffXV[iDet - 1];
                        yV_cm[iDet - 1][iPart] =
                            (yV_cm[iDet - 1][iPart] - flosOffsetYV[iDet - 1]) * flosVeffYV[iDet - 1];

                        // Position from ToT:
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

                        xToT_cm[iDet - 1][iPart] =
                            (xToT_cm[iDet - 1][iPart] - flosOffsetXQ[iDet - 1]) * flosVeffXQ[iDet - 1];
                        yToT_cm[iDet - 1][iPart] =
                            (yToT_cm[iDet - 1][iPart] - flosOffsetYQ[iDet - 1]) * flosVeffYQ[iDet - 1];

                        // cout<<"iDet "<<iDet<<" iPart "<<iPart<<"  "<<nPMV<<endl;
                        if (nPMV == 8)
                        {
                            // time difference between two hits
                            if (time0[iDet - 1] < 0)
                            {
                                time0[iDet - 1] = timeLosV[iDet - 1][iPart];
                                Double_t tt = (time0[iDet - 1] + time - time_spill_start) / 1.e9;
                                // if (iDet == 1)
                                // fh_spill_times_daq_first->Fill(tt);
                            }
                            time1[iDet - 1] = timeLosV[iDet - 1][iPart];
                            time_abs[iDet - 1] =
                                time1[iDet - 1] - time0[iDet - 1] + (double)(time - time_prev[iDet - 1]);
                            // cout<<"idet: "<<iDet<<"  "<<time_abs[iDet-1]<<endl;
                            // cout<<"Time: "<<time<<" time_prev: "<<time_prev<<" dt:
                            // "<<(double)(time-time_prev[iDet-1])<<endl; cout<<"time0: "<<time0[iDet-1]<<" time1:
                            // "<<time1[iDet-1]<<" dt: "<<time1[iDet-1] - time0[iDet-1]<<endl;
                            if (time_abs[iDet - 1] > 0 && time_abs[iDet - 1] < 5.E8)
                            {
                                fh_los_dt_hits_ToT[iDet - 1]->Fill(time_abs[iDet - 1] / 1000., totsum[iDet - 1][iPart]);
                                // cout<<"draw "<<time_abs[iDet-1]/1000.<<"  " <<totsum[iDet-1][iPart]<<endl;
                                fh_los_dt_hits[iDet - 1]->Fill(time_abs[iDet - 1] / 1000.);
                            }
                            time_prev[iDet - 1] = time;
                            time0[iDet - 1] = time1[iDet - 1];

                            // Plot Events over one spill
                            Double_t tt =
                                (time1[iDet - 1] - time0[iDet - 1] + (double)(time - time_spill_start)) / 1.e9;
                            //                            fh_spill_times_daq->Fill(tt);

                            // MH
                            /*
                                        // time difference between two events in one read-out window
                                        time2[iDet-1]=0.;
                                        if (time0[iDet-1] < 0 ) {
                                            time0[iDet-1] = timeLosV[iDet-1][iPart];
                                            Double_t tt = (time0[iDet-1] + time - time_spill_start)/1.e9;
                                            fh_spill_times_daq_first->Fill(tt);
                                        }
                                        if (time1[iDet-1] < 0 ) time1[iDet-1] = timeLosV[iDet-1][iPart];

                                        time2[iDet-1] = timeLosV[iDet-1][iPart] ;


                                        double_t time_abs[iDet-1] = time2[iDet-1] - time1[iDet-1] +
                            (double)(time-time_prev); if(time_abs[iDet-1] > 0 && time_abs[iDet-1] < 5.E8) {
                                            fh_los_dt_hits->Fill(time_abs[iDet-1] / 1000.);
                            //                fh_frequencies->Fill(1./time_abs*1.E6);
                                            //if(1./time_abs*1.E6<30.)
                                            //cout << "Time Test " << time1 << "  " << time2 << "  " << time2 - time1 <<
                            "  " << time << "  " << time_prev << "  " << time_abs<<"  "<<1./time_abs*1.E6<<endl;
                                        }

                                        time_prev[iDet-1] = time;
                                        time1[iDet-1]=time2[iDet-1];


                                        // Plot Events over one spill
                                        Double_t tt = (time2[iDet-1] - time0[iDet-1] + time - time_spill_start)/1.e9;
                                        fh_spill_times_daq->Fill(tt);
                            */
                            // end MH

                            if (foptcond == 0)
                            {
                                fh_los_tot_mean[iDet - 1]->Fill(totsum[iDet - 1][iPart]);
                                for (int ipm = 0; ipm < 8; ipm++)
                                {
                                    fh_los_tot[iDet - 1]->Fill(ipm + 1, tot[iDet - 1][iPart][ipm]);
                                }

                                fh_los_tres_MCFD[iDet - 1]->Fill(LosTresV[iDet - 1][iPart]);
                                fh_los_tres_TAMEX[iDet - 1]->Fill(LosTresT[iDet - 1][iPart]);
                                if (iDet == 1)
                                    fh_los_tres_MTDC->Fill(LosTresMTDC[iPart]);

                                fh_los_pos_MCFD[iDet - 1]->Fill(xV_cm[iDet - 1][iPart], yV_cm[iDet - 1][iPart]);
                                fh_los_pos_TAMEX[iDet - 1]->Fill(xT_cm[iDet - 1][iPart], yT_cm[iDet - 1][iPart]);
                                fh_los_pos_ToT[iDet - 1]->Fill(xToT_cm[iDet - 1][iPart], yToT_cm[iDet - 1][iPart]);
                            }

                            fh_los_ihit_ToT[iDet - 1]->Fill(iPart + 1, totsum[iDet - 1][iPart]);
                            fh_los_multihit[iDet - 1]->Fill(iPart + 1);
                        }
                    }
                }
                else
                {
                    cout << "Wrong detector ID for LOS!" << endl;
                }
            }
        }
    }

    if (fCalItems.at(DET_LOS) && fNofLosDetectors == 2)
    {
        for (Int_t iPart = 0; iPart < nPart; iPart++)
        {
            //	 fh_los_ToT_vs_ToT->Fill(totsum[0][iPart],totsum[1][iPart]);
        }
    }

    //----------------------------------------------------------------------
    // SCI8 detector
    //----------------------------------------------------------------------

    Double_t timeS8_V[10][2] = { { 0.0 / 0.0 } }; // [multihit][pm]
    Double_t timeS8_L[10][2] = { { 0.0 / 0.0 } };
    Double_t timeS8_T[10][2] = { { 0.0 / 0.0 } };
    Double_t timeSci8M[10] = { 0.0 };
    Double_t Sci8TresM[10] = { 0.0 / 0.0 };
    Double_t timeSci8T[10] = { 0.0 };
    Double_t Sci8TresT[10] = { 0.0 / 0.0 };
    Double_t timeSci8[10] = { 0.0 };
    Double_t totsumS8[10] = { 0.0 };
    Double_t totS8[10][8] = { { 0.0 / 0.0 } };

    Int_t MultipS8 = -1;

    if (fMappedItems.at(DET_SCI8))
    {
        auto det = fMappedItems.at(DET_SCI8);
        Int_t nHits = det->GetEntriesFast();

        MultipS8 = nHits;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSci8MappedData* hit = dynamic_cast<R3BSci8MappedData*>(det->At(ihit));
            if (!hit)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hit->GetDetector(); // 1..
            Int_t iCha = hit->GetChannel();  // 1..

            fh_sci8_channels->Fill(iCha);
        }
    }
    // assert(MultipS8 != -1);

    Int_t nPartS8;

    if (fCalItems.at(DET_SCI8))
    {
        auto det = fCalItems.at(DET_SCI8);
        nPartS8 = det->GetEntriesFast();

        fh_sci8_multihit->Fill(nPartS8);

        Int_t iDet = 0;
        Int_t nPartS8_VFTX[2] = { 0 };
        Int_t nPartS8_LEAD[2] = { 0 };
        Int_t nPartS8_TRAI[2] = { 0 };

        for (Int_t iPart = 0; iPart < nPartS8; iPart++)
        {

            /*
             * nPart is the number of particle passing through Sci8 detector in one event
             */

            R3BSci8CalData* calDataS8 = dynamic_cast<R3BSci8CalData*>(det->At(iPart));
            iDet = calDataS8->GetDetector();

            // VFTX Channels 1-2:
            if (!(IS_NAN(calDataS8->fTimeV_r_ns)))
            {
                timeS8_V[iPart][0] = calDataS8->fTimeV_r_ns;
                nPartS8_VFTX[0] += 1;
            }
            if (!(IS_NAN(calDataS8->fTimeV_l_ns)))
            {
                timeS8_V[iPart][1] = calDataS8->fTimeV_l_ns;
                nPartS8_VFTX[1] += 1;
            }

            // TAMEX Channels 1-2:
            if (!(IS_NAN(calDataS8->fTimeL_r_ns)))
            {
                timeS8_L[iPart][0] = calDataS8->fTimeL_r_ns;
                nPartS8_LEAD[0] += 1;
            }
            if (!(IS_NAN(calDataS8->fTimeT_r_ns)))
            {
                timeS8_T[iPart][0] = calDataS8->fTimeT_r_ns;
                nPartS8_TRAI[0] += 1;
            }
            if (!(IS_NAN(calDataS8->fTimeL_l_ns)))
            {
                timeS8_L[iPart][1] = calDataS8->fTimeL_l_ns;
                nPartS8_LEAD[1] += 1;
            }
            if (!(IS_NAN(calDataS8->fTimeT_l_ns)))
            {
                timeS8_T[iPart][1] = calDataS8->fTimeT_l_ns;
                nPartS8_TRAI[1] += 1;
            }

            if (iPart > 0 && MultipS8 % 2 == 0)
            {
                for (int k = 0; k < 2; k++)
                {
                    if (timeS8_V[iPart][k] > 0. && timeS8_V[iPart - 1][k] > 0. && !(IS_NAN(timeS8_V[iPart][k])) &&
                        !(IS_NAN(timeS8_V[iPart - 1][k])))
                    {
                        fh_sci8_dt_hits->Fill(timeS8_V[iPart][k] - timeS8_V[iPart - 1][k]);
                    }
                    if (timeS8_L[iPart][k] > 0. && timeS8_L[iPart - 1][k] > 0. && !(IS_NAN(timeS8_L[iPart][k])) &&
                        !(IS_NAN(timeS8_L[iPart - 1][k])))
                    {
                        fh_sci8_dt_hits_l->Fill(timeS8_L[iPart][k] - timeS8_L[iPart - 1][k]);
                    }
                    if (timeS8_T[iPart][k] > 0. && timeS8_T[iPart - 1][k] > 0. && !(IS_NAN(timeS8_T[iPart][k])) &&
                        !(IS_NAN(timeS8_T[iPart - 1][k])))
                    {
                        fh_sci8_dt_hits_t->Fill(timeS8_T[iPart][k] - timeS8_T[iPart - 1][k]);
                    }
                }
            }

            if (iDet == 1)
            {

                // if(iPart >= 0)
                if (1 == 1)
                {
                    int nPMT = 0;
                    int nPMV = 0;
                    int ilc = iPart;

                    for (int ipm = 0; ipm < 2; ipm++)
                    {

                        if (timeS8_T[iPart][ipm] > 0. && timeS8_L[iPart][ipm] > 0. && !(IS_NAN(timeS8_T[iPart][ipm])) &&
                            !(IS_NAN(timeS8_L[iPart][ipm])))
                        {
                            while (timeS8_T[iPart][ipm] - timeS8_L[iPart][ipm] < 0.)
                            {
                                timeS8_T[iPart][ipm] = timeS8_T[iPart][ipm] + 2048. * fClockFreq;
                            }

                            nPMT = nPMT + 1;
                            totS8[iPart][ipm] = timeS8_T[iPart][ipm] - timeS8_L[iPart][ipm];
                        }

                        totsumS8[iPart] += totS8[iPart][ipm];

                        if (totS8[iPart][ipm] != 0. && !(IS_NAN(totS8[iPart][ipm])))
                            fh_sci8_tot->Fill(ipm + 1, totS8[iPart][ipm]);

                        if (timeS8_L[iPart][ipm] > 0. && !(IS_NAN(timeS8_L[iPart][ipm])))
                            timeSci8T[iPart] += timeS8_L[iPart][ipm];

                        if (timeS8_V[iPart][ipm] > 0. && !(IS_NAN(timeS8_V[iPart][ipm])))
                        {
                            timeSci8M[iPart] += timeS8_V[iPart][ipm];
                            nPMV = nPMV + 1;
                        }
                    }

                    totsumS8[iPart] = totsumS8[iPart] / nPMT;
                    /*
                    if(totsum[iPart] < 88.) cout<<fNEvents<<"; "<<nPart<<"; "<<iPart<<",
                    "<<totsum[iPart]<<tot[iPart][0]<<
                    ", "<<tot[iPart][1]<<", " <<tot[iPart][2]<<", "<<tot[iPart][3]<<", "<<tot[iPart][4]<<", "
                    <<tot[iPart][5]<<", "<<tot[iPart][6]<<", "<<tot[iPart][7]<<endl;
                    */
                    timeSci8M[iPart] = timeSci8M[iPart] / nPMV;
                    timeSci8T[iPart] = timeSci8T[iPart] / nPMT;

                    timeSci8[iPart] = timeSci8M[iPart];

                    Double_t tof = (0.0 / 0.0);
                    // at first, take first LOS hit
                    if (!IS_NAN(timeLos[1][0]) && timeLos[1][0] > 0.)
                    {
                        if (timeSci8[iPart] - timeLos[1][0] < 8192. * 5. / 2.)
                            timeSci8[iPart] += 8192. * 5.;
                        if (timeSci8[iPart] - timeLos[1][0] > 8192. * 5. / 2.)
                            timeSci8[iPart] -= 8192. * 5.;
                        tof = timeSci8[iPart] - timeLos[1][0];
                    }

                    // cout<<"tof "<<tof<<"  "<<timeSci8[iPart]<<"  "<<timeLos[1][0]<<endl;
                    // loop over multihits of LOS
                    for (Int_t lhit = 1; lhit < nPartLos[1]; lhit++)
                    {
                        if (IS_NAN(timeLos[1][lhit]))
                            continue;
                        if (timeSci8[iPart] - timeLos[1][lhit] < 8192. * 5. / 2.)
                            timeSci8[iPart] += 8192. * 5.;
                        if (timeSci8[iPart] - timeLos[1][lhit] > 8192. * 5. / 2.)
                            timeSci8[iPart] -= 8192. * 5.;
                        Double_t tofi = timeSci8[iPart] - timeLos[1][lhit];
                        // if hit is closer to expected value then take it
                        if (abs(tofi - 7852.) < abs(tof - 7852.))
                            tof = tofi;
                        // cout<<"tofi "<<tof<<"  "<<timeSci8[iPart]<<"  "<<timeLos[1][0]<<endl;
                    }

                    fh_tof_sci8->Fill(tof);

                    if (nPMV == 2)
                        Sci8TresM[iPart] = (timeS8_V[iPart][1] - timeS8_V[iPart][0]);
                    if (nPMT == 2)
                        Sci8TresT[iPart] = (timeS8_L[iPart][1] - timeS8_L[iPart][0]);

                    if (nPMV == 2)
                        fh_sci8_tres_MCFD->Fill(Sci8TresM[iPart]);
                    if (nPMT == 2)
                        fh_sci8_tres_TAMEX->Fill(Sci8TresT[iPart]);
                    if (nPMT == 2)
                        fh_sci8_tot_mean->Fill(totsumS8[iPart]);
                }
            }
            else
            {
                cout << "Wrong detector ID for Sci8!" << endl;
            }
        }

        for (int ik = 0; ik < 2; ik++)
        {
            fh_sci8_multihitVFTX->Fill(ik + 1, nPartS8_VFTX[ik]);
            fh_sci8_multihitLEAD->Fill(ik + 1, nPartS8_LEAD[ik]);
            fh_sci8_multihitTRAI->Fill(ik + 1, nPartS8_TRAI[ik]);
        }
    }

    //----------------------------------------------------------------------
    // Fiber detectors
    //----------------------------------------------------------------------
    Double_t dtime = 0.0 / 0.0;
    Bool_t fib1a_cut = false;
    Bool_t fib1b_cut = false;
    Int_t iLosHit = 0;

    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {

        Int_t iFib = 0;

        auto detMapped = fMappedItems.at(DET_FI_FIRST + ifibcount);
        auto detHit = fHitItems.at(DET_FI_FIRST + ifibcount);

        if (detMapped)
        {
            Int_t nHits = detMapped->GetEntriesFast();
            std::vector<UInt_t> mapmt_num(512);
            std::vector<UInt_t> spmt_num(16);
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BBunchedFiberMappedData* hit = dynamic_cast<R3BBunchedFiberMappedData*>(detMapped->At(ihit));
                if (!hit)
                    continue;

                // channel numbers are stored 1-based (1..n)
                Int_t iCha = hit->GetChannel(); // 1..

                if (hit->IsMAPMT() && hit->IsLeading())
                {
                    fh_channels_Fib[ifibcount]->Fill(iCha); // Fill which clockTDC channel has events
                    ++mapmt_num.at(hit->GetChannel() - 1);  // multihit of a given clockTDC channel
                }

                if (!hit->IsMAPMT() && hit->IsLeading())
                {
                    fh_channels_single_Fib[ifibcount]->Fill(iCha); // Fill which single PMT channel has events
                    ++spmt_num.at(hit->GetChannel() - 1);          // multihit of a given PADI channel
                }
            }

            for (int i = 0; i < 512; ++i)
            {
                auto m = mapmt_num.at(i);
                if (m > 0)
                    fh_multihit_m_Fib[ifibcount]->Fill(i + 1, m); // multihit of a given clockTDC channel
            }

            for (int i = 0; i < 16; ++i)
            {
                auto s = spmt_num.at(i);
                if (s > 0)
                    fh_multihit_s_Fib[ifibcount]->Fill(i + 1, s); // multihit of a given PADI channel
            }
        }

        if (detHit)
        {
            Int_t nHits = detHit->GetEntriesFast();
            Double_t xpos = 0. / 0.;
            Double_t ypos = 0. / 0.;
            Double_t totMax = 0.;
            Double_t totMax_MA = 0.;
            Double_t tfib = 0. / 0., tof_fib = 0. / 0.;
            Double_t randx;
            Int_t iFibmax = -100;
            Int_t iFibMax_MA = -100;
            Double_t totMax_S = 0.;
            Double_t yposMax = -10000;
            Double_t tof_fibMax = -10000;
            // cout<<"new event: "<< nHits<<"  "<<ifibcount<<endl;

            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                Double_t tMAPMT = 0. / 0.;
                Double_t tSPMT = 0. / 0.;

                R3BBunchedFiberHitData* hit = dynamic_cast<R3BBunchedFiberHitData*>(detHit->At(ihit));
                if (!hit)
                    continue;

                iFib = hit->GetFiberId(); // 1..

                // times
                tMAPMT = hit->GetMAPMTTime_ns();
                tSPMT = hit->GetSPMTTime_ns();

                // if not both times present, go to next
                //                if (!(tMAPMT > 0.) || !(tSPMT > 0.))
                //                    continue;

                // loop over all hits of LOS1:
                for (Int_t lhit = 0; lhit < nPartLos[0]; lhit++)
                {
                    if (IS_NAN(timeLos[0][lhit]))
                        continue;

                    // if not multiplicity 8 or pile-up, go to next event
                    Bool_t pileup = false;
                    Int_t nPMT = 0;
                    for (int ipm = 0; ipm < 8; ipm++)
                    {
                        if (tot[0][iLosHit][ipm] > 0. && !(IS_NAN(tot[0][iLosHit][ipm])))
                            nPMT++;
                        if (tot[0][iLosHit][ipm] > 200.)
                            pileup = true;
                    }
                    //					if (pileup || nPMT != 8)
                    if (pileup)
                        continue;

                    dtime = tMAPMT - tSPMT;

                    // Not-calibrated ToF:
                    // tfib = (tMAPMT + tSPMT) / 2.;
                    // tfib = tSPMT;
                    tfib = tMAPMT;
                    if (tfib > 0. && !(IS_NAN(tfib)) && timeLos[0][lhit] > 0. && !(IS_NAN(timeLos[0][lhit])))
                        tof_fib = tfib - timeLos[0][lhit];

                    while (tof_fib < 0)
                    {
                        tof_fib += 2048 * 1000. / 150.;
                    }
                    while (tof_fib > 20000.)
                    {
                        tof_fib -= 2048 * 1000. / 150.;
                    }
                    cout << "Test " << tof_fib << "  " << tfib << "  " << hit->GetMAPMTToT_ns() << "  "
                         << hit->GetSPMTToT_ns() << "  " << lhit << "  " << timeLos[0][lhit] << endl;

                    // if ToF between LOS1 and fiber is not okay, go to next
                    //                    if (tof_fib < tof_fiber[ifibcount] - 10. || tof_fib > tof_fiber[ifibcount]
                    //                    + 10.)
                    //                        continue;

                    // Not-calibrated position:
                    randx = (std::rand() / (float)RAND_MAX);
                    //                if(iFib > 0) xpos = ((-n_fiber[ifibcount]/2.+iFib+(0.5-randx)))*0.21;
                    xpos = iFib;
                    ypos = dtime;

                    // if not resonable y-position, go to next
                    //                    if (ypos < 70. || ypos > 100.)
                    //                        continue;

                    // find fiber with most light which has resonable ToF and y-position
                    Double_t ToT_MA = hit->GetMAPMTToT_ns();
                    if (ToT_MA > totMax_MA)
                    {
                        totMax_MA = ToT_MA;
                        iFibMax_MA = iFib;
                        totMax_S = hit->GetSPMTToT_ns();
                        yposMax = ypos;
                        tof_fibMax = tof_fib;
                        //                        iLosHit = lhit;
                        // cout<<"new imax "<<iFibMax_MA<<" max "<<totMax_MA<<endl;
                    }

                    fh_fibers_Fib[ifibcount]->Fill(iFib);
                    fh_ToT_s_Fib[ifibcount]->Fill(iFib, hit->GetSPMTToT_ns());
                    fh_ToT_m_Fib[ifibcount]->Fill(iFib, hit->GetMAPMTToT_ns());
                    fh_time_Fib[ifibcount]->Fill(iFib, tMAPMT - tSPMT);
                    cout << "Result " << iFib << "  " << tof_fib << endl;
                    fh_Fib_ToF[ifibcount]->Fill(iFib, tof_fib);
                    fh_Fib_pos_xy[ifibcount]->Fill(xpos, ypos);
                    fh_Fib_vs_Events[ifibcount]->Fill(fNEvents, iFib);
                    fh_ToT_single_Fib[ifibcount]->Fill(1, hit->GetSPMTToT_ns());

                } // end for (lhit)
            }     // end for(ihit)
                  /*
                              fh_fibers_Fib[ifibcount]->Fill(iFibMax_MA);
                              fh_ToT_s_Fib[ifibcount]->Fill(iFibMax_MA, totMax_S);
                              fh_ToT_m_Fib[ifibcount]->Fill(iFibMax_MA, totMax_MA);
                              fh_time_Fib[ifibcount]->Fill(iFibMax_MA, yposMax);
                              fh_Fib_ToF[ifibcount]->Fill(iFibMax_MA, tof_fibMax);
                              fh_Fib_pos_xy[ifibcount]->Fill(xpos, ypos);
                              fh_Fib_vs_Events[ifibcount]->Fill(fNEvents, iFibMax_MA);
                              fh_ToT_single_Fib[ifibcount]->Fill(1, totMax_S);
                  */
            if (ifibcount == 0 && iFibMax_MA > 120 && iFibMax_MA < 136)
            {
                fib1a_cut = true;
            }
            if (ifibcount == 1 && iFibMax_MA > 120 && iFibMax_MA < 136)
            {
                fib1b_cut = true;
                // cout<<"Cut true "<< iFibMax_MA << endl;
            }
            if (nHits > 0)
                fh_mult_Fib[ifibcount]->Fill(nHits);
            // cout<<"Hits "<<nHits<<endl;

        } // end if(aHit[ifibcount])

    } // end for(ifibcount)

    // Plot LOS XY with condition on FIB position
    // xposfib = (-128.+iFibMax_MA+(0.5-randx))*0.021;

    if (foptcond && (fib1a_cut == true && fib1b_cut == true))
    {
        for (Int_t ilos = 0; ilos < fNofLosDetectors; ilos++)
        {
            if (iLosHit < nPartLos[ilos])
            {
                fh_los_pos_MCFD[ilos]->Fill(xV_cm[ilos][iLosHit], yV_cm[ilos][iLosHit]);
                fh_los_pos_TAMEX[ilos]->Fill(xT_cm[ilos][iLosHit], yT_cm[ilos][iLosHit]);
                fh_los_pos_ToT[ilos]->Fill(xToT_cm[ilos][iLosHit], yToT_cm[ilos][iLosHit]);
                fh_los_tres_MCFD[ilos]->Fill(LosTresV[ilos][iLosHit]);
                fh_los_tres_TAMEX[ilos]->Fill(LosTresT[ilos][iLosHit]);
                // fh_los_tres_MTDC->Fill(LosTresMTDC[iLosHit]);
                fh_los_tot_mean[ilos]->Fill(totsum[ilos][iLosHit]);
                for (int ipm = 0; ipm < 8; ipm++)
                {
                    fh_los_tot[ilos]->Fill(ipm + 1, tot[ilos][iLosHit][ipm]);
                }
                fh_los_ToT_vs_ToT->Fill(totsum[0][iLosHit], totsum[1][iLosHit]);
            }
        }
    }

    //----------------------------------------------------------------------
    // TOFD
    //----------------------------------------------------------------------

    Int_t NumPaddles[4] = { 0 };
    if (fMappedItems.at(DET_TOFD))
    {
        auto det = fMappedItems.at(DET_TOFD);
        Int_t nMapped = det->GetEntriesFast();
        Int_t iPlaneMem = 1, iBarMem = 0;
        for (Int_t imapped = 0; imapped < nMapped; imapped++)
        {
            auto mapped = dynamic_cast<R3BTofdMappedData const*>(det->At(imapped));
            if (!mapped)
                continue; // should not happen

            Int_t const iPlane = mapped->GetDetectorId(); // 1..n
            Int_t const iBar = mapped->GetBarId();        // 1..n
            Int_t const iSide = mapped->GetSideId();      // 1..n
            Int_t const iEdge = mapped->GetEdgeId();

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
                if (iSide == 1)
                    fh_tofd_channels[iPlane - 1]->Fill(-iBar - 1);
                if (iSide == 2)
                    fh_tofd_channels[iPlane - 1]->Fill(iBar);
            }
        }
    }

    if (fCalItems.at(DET_TOFD))
    {
        auto det = fCalItems.at(DET_TOFD);
        Int_t nCals = det->GetEntriesFast();

        Double_t tot1[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { 0. / 0. } } };
        Double_t tot2[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { 0. / 0. } } };
        Double_t t_paddle[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { 0. / 0. } } };
        Double_t t1l[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { 0. / 0. } } };
        Double_t t2l[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { 0. / 0. } } };
        Double_t t1t[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { 0. / 0. } } };
        Double_t t2t[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { 0. / 0. } } };
        Double_t ToF[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { 0. / 0. } } };

        Bool_t Bar_present[10][N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { { false } } };

        Int_t iBarMem = 0;
        Int_t jmult[N_PLANE_MAX_TOFD][N_PADDLE_MAX_TOFD] = { { 0 } };

        unsigned long long time0 = header->GetTimeStamp();
        double_t time1 = -1.;

        for (Int_t ical = 0; ical < nCals; ical++)
        {
            auto cal = dynamic_cast<R3BTofdCalData const*>(det->At(ical));
            if (!cal)
                continue; // should not happen

            Int_t const iPlane = cal->GetDetectorId(); // 1..n
            Int_t const iBar = cal->GetBarId();        // 1..n

            //      std::cout << iPlane << ' ' << iBar <<
            //          ',' << cal->GetTimeBL_ns() <<
            //          ' ' << cal->GetTimeBT_ns() <<
            //          ' ' << cal->GetTimeTL_ns() <<
            //          ' ' << cal->GetTimeTT_ns() << std::endl;

            // get all times of one bar
            /*
                         cout<<"TOFD Online: "<<fNEvents<<", "<<nCals<<"; "<<ical<<"; "<<iPlane<<", "<<iBar<<",
               "<<iBarMem<<", "<<", "<< cal->fTime1L_ns<<", "<<cal->fTime1T_ns<<", "<< cal->fTime2L_ns<<",
               "<<cal->fTime2T_ns<<endl;
            */
            /*
                  if(iBar != iBarMem )
                  {
                jmult[iPlane-1][iBar-1] = 0;
                if (!IS_NAN(cal->GetTimeBL_ns())) t1l[0][iPlane-1][iBar-1] = cal->GetTimeBL_ns();
                if (!IS_NAN(cal->GetTimeBT_ns())) t1t[0][iPlane-1][iBar-1] = cal->GetTimeBT_ns();
                if (!IS_NAN(cal->GetTimeTL_ns())) t2l[0][iPlane-1][iBar-1] = cal->GetTimeTL_ns();
                if (!IS_NAN(cal->GetTimeTT_ns())) t2t[0][iPlane-1][iBar-1] = cal->GetTimeTT_ns();
                Bar_present[0][iPlane-1][iBar-1] = true;
                  }
                  else
                  {
                jmult[iPlane-1][iBar-1] = jmult[iPlane-1][iBar-1] + 1;
                Int_t jm = jmult[iPlane-1][iBar-1];
                if (!IS_NAN(cal->GetTimeBL_ns())) t1l[jm][iPlane-1][iBar-1] = cal->GetTimeBL_ns();
                if (!IS_NAN(cal->GetTimeBT_ns())) t1t[jm][iPlane-1][iBar-1] = cal->GetTimeBT_ns();
                if (!IS_NAN(cal->GetTimeTL_ns())) t2l[jm][iPlane-1][iBar-1] = cal->GetTimeTL_ns();
                if (!IS_NAN(cal->GetTimeTT_ns())) t2t[jm][iPlane-1][iBar-1] = cal->GetTimeTT_ns();
                Bar_present[jm][iPlane-1][iBar-1] = true;
                  }
                  iBarMem = iBar;
            */

            Int_t jm = jmult[iPlane - 1][iBar - 1];
            /*		if (!IS_NAN(cal->GetTimeBL_ns())) t1l[jm][iPlane-1][iBar-1] = cal->GetTimeBL_ns();
                    if (!IS_NAN(cal->GetTimeBT_ns())) t1t[jm][iPlane-1][iBar-1] = cal->GetTimeBT_ns();
                    if (!IS_NAN(cal->GetTimeTL_ns())) t2l[jm][iPlane-1][iBar-1] = cal->GetTimeTL_ns();
                    if (!IS_NAN(cal->GetTimeTT_ns())) t2t[jm][iPlane-1][iBar-1] = cal->GetTimeTT_ns();*/
            Bar_present[jm][iPlane - 1][iBar - 1] = true;
            jmult[iPlane - 1][iBar - 1] = jmult[iPlane - 1][iBar - 1] + 1;

            double_t time2 = 0.;
            if (time1 < 0 && iPlane == 2)
                time1 = (t1l[jm][iPlane - 1][iBar - 1] + t2l[jm][iPlane - 1][iBar - 1]) / 2.;
            if (iPlane == 2)
                time2 = (t1l[jm][iPlane - 1][iBar - 1] + t2l[jm][iPlane - 1][iBar - 1]) / 2.;
            if (time1 > 0. && time2 > 0. && time2 > time1)
            {
                // cout<<"Time Test "<<time0<<"  "<<time1<< "   "<< time2 <<"  " <<time_previous_event <<endl;
                fh_TimePreviousEvent->Fill(time2 - time1);
                time2 = time1;
            }
        }

        /*
                for (Int_t ipl = 0; ipl < N_PLANE_MAX_TOFD; ipl++)
                for(Int_t ibr = 0; ibr < N_PADDLE_MAX_TOFD; ibr++)
                {
                fh_tofd_multihit[ipl]->Fill(ibr+1,jmult[ipl][ibr]+1);
                }
         */
        for (Int_t ipl = 0; ipl < N_PLANE_MAX_TOFD; ipl++)
            for (Int_t ibr = 0; ibr < N_PADDLE_MAX_TOFD; ibr++)
                for (Int_t jm = 0; jm < jmult[ipl][ibr]; jm++)
                {

                    if (Bar_present[jm][ipl][ibr])
                    {
                        Int_t iPlane = ipl + 1; // 1..n
                        Int_t iBar = ibr + 1;   // 1..n

                        fh_tofd_multihit[ipl]->Fill(ibr + 1, jmult[ipl][ibr]);

                        // calculate time over threshold and check if clock counter went out of range

                        while (t1t[jm][iPlane - 1][iBar - 1] < t1l[jm][iPlane - 1][iBar - 1])
                        {
                            t1t[jm][iPlane - 1][iBar - 1] += 2048. * fClockFreq;
                        }

                        while (t2t[jm][iPlane - 1][iBar - 1] < t2l[jm][iPlane - 1][iBar - 1])
                        {
                            t2t[jm][iPlane - 1][iBar - 1] += 2048. * fClockFreq;
                        }

                        int ilc = jm;

                        if (!(IS_NAN(timeLos[0][ilc])) && timeLos[0][ilc] > 0.)
                        {

                            while (t1l[jm][iPlane - 1][iBar - 1] < timeLos[0][ilc])
                            {
                                t1t[jm][iPlane - 1][iBar - 1] += 2048. * fClockFreq;
                                t1l[jm][iPlane - 1][iBar - 1] += 2048. * fClockFreq;
                            }
                            while (t2l[jm][iPlane - 1][iBar - 1] < timeLos[0][ilc])
                            {
                                t2t[jm][iPlane - 1][iBar - 1] += 2048. * fClockFreq;
                                t2l[jm][iPlane - 1][iBar - 1] += 2048. * fClockFreq;
                            }
                        }

                        // ToF
                        t_paddle[jm][iPlane - 1][iBar - 1] =
                            (t1l[jm][iPlane - 1][iBar - 1] + t2l[jm][iPlane - 1][iBar - 1]) / 2.;

                        //  between 2 bars in 2 planes
                        if (ipl > 0)
                        {
                            fh_tofd_dt[ipl - 1]->Fill(iBar,
                                                      t_paddle[jm][ipl][iBar - 1] - t_paddle[jm][ipl - 1][iBar - 1]);
                        }

                        if (!(IS_NAN(timeLos[0][ilc])) && timeLos[0][ilc] > 0.)
                        {
                            // between LOS and paddle
                            ToF[jm][iPlane - 1][iBar - 1] =
                                fmod(t_paddle[jm][iPlane - 1][iBar - 1] - timeLos[0][ilc] + 5 * 8192, 5 * 2048);
                            fh_tofd_ToF[iPlane - 1]->Fill(iBar, ToF[jm][iPlane - 1][iBar - 1]);
                        }
                        // ToT
                        tot1[jm][iPlane - 1][iBar - 1] = t1t[jm][iPlane - 1][iBar - 1] - t1l[jm][iPlane - 1][iBar - 1];
                        if (tot1[jm][iPlane - 1][iBar - 1] < 0)
                        {
                            cout << "Negative ToT " << tot1[jm][iPlane - 1] << ", for hit= " << jm << endl;
                            cout << "Times1: " << t1t[jm][iPlane - 1][iBar - 1] << " " << t1l[jm][iPlane - 1][iBar - 1]
                                 << endl;
                        }

                        tot2[jm][iPlane - 1][iBar - 1] = t2t[jm][iPlane - 1][iBar - 1] - t2l[jm][iPlane - 1][iBar - 1];
                        if (tot2[jm][iPlane - 1][iBar - 1] < 0)
                        {
                            cout << "Negative ToT " << tot2[jm][iPlane - 1][iBar - 1] << ", for hit= " << jm << endl;
                            cout << "Times2: " << t2t[jm][iPlane - 1][iBar - 1] << " " << t2l[jm][iPlane - 1][iBar - 1]
                                 << endl;
                        }

                        fh_tofd_TotPm[iPlane - 1]->Fill(iBar, tot2[jm][iPlane - 1][iBar - 1]);
                        fh_tofd_TotPm[iPlane - 1]->Fill(-iBar - 1, tot1[jm][iPlane - 1][iBar - 1]);
                    }
                }
    }

    //----------------------------------------------------------------------
    // PTOF
    //----------------------------------------------------------------------

    if (fCalItems.at(DET_PTOF))
    {
        auto det = fCalItems.at(DET_PTOF);

        Double_t tot1 = 0.;
        Double_t tot2 = 0.;
        Double_t t1l = 0.;
        Double_t t2l = 0.;
        Double_t t1t = 0.;
        Double_t t2t = 0.;
        Bool_t bar_quer1 = false;
        Bool_t bar_quer2 = false;

        Int_t nHits = det->GetEntriesFast();
        LOG(debug) << "nHits: " << nHits;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BPaddleCalData* hit = dynamic_cast<R3BPaddleCalData*>(det->At(ihit));

            if (!hit)
                continue; // should not happen

            Int_t iBar = hit->GetBar(); // 1..n
            // get all times of one bar
            t1l = hit->fTime1L_ns;
            t2l = hit->fTime2L_ns;
            t1t = hit->fTime1T_ns;
            t2t = hit->fTime2T_ns;

            if (!(t1l > 0 && t2l > 0 && t1t > 0 && t2t > 0))
                continue;

            fh_ptof_channels->Fill(iBar);
            LOG(debug) << "Bar: " << iBar;
            LOG(debug) << "times PM1: " << t1l << "  " << t1t << "  " << t1t - t1l;
            LOG(debug) << "times PM2: " << t2l << "  " << t2t << "  " << t2t - t2l;
            if (iBar == 7)
                bar_quer1 = true;
            if (iBar == 8)
                bar_quer2 = true;

            // calculate time over threshold and check if clock counter went out of range
            while (t1t - t1l < 0.)
            {
                t1t = t1t + 2048. * fClockFreq;
            }

            while (t2t - t2l < 0.)
            {
                t2t = t2t + 2048. * fClockFreq;
            }
            while (t1l - timeLos[0][0] < 0.)
            {
                t1t = t1t + 2048. * fClockFreq;
                t1l = t1l + 2048. * fClockFreq;
                t2t = t2t + 2048. * fClockFreq;
                t2l = t2l + 2048. * fClockFreq;
            }

            tot1 = t1t - t1l;
            if (tot1 < 0)
            {
                LOG(warn) << "Negative ToT " << tot1;
                LOG(warn) << "times1: " << t1t << " " << t1l;
            }

            tot2 = t2t - t2l;
            if (tot2 < 0)
            {
                LOG(warn) << "Negative ToT " << tot2;
                LOG(warn) << "times2: " << t2t << " " << t2l;
            }

            fh_ptof_TotPm1[iBar]->Fill(tot1);
            fh_ptof_TotPm2[iBar]->Fill(tot2);
            if (iBar == 2)
                fh_ptof_test1->Fill(sqrt(tot1 * tot1));
        }

        // once again

        nHits = det->GetEntriesFast();
        //		LOG(debug) << "nHits: " << nHits;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BPaddleCalData* hit = dynamic_cast<R3BPaddleCalData*>(det->At(ihit));

            if (!hit)
                continue; // should not happen

            Int_t iBar = hit->GetBar(); // 1..n
            // get all times of one bar
            t1l = hit->fTime1L_ns;
            t2l = hit->fTime2L_ns;
            t1t = hit->fTime1T_ns;
            t2t = hit->fTime2T_ns;

            if (!(t1l > 0 && t2l > 0 && t1t > 0 && t2t > 0))
                continue;

            // calculate time over threshold and check if clock counter went out of range
            while (t1t - t1l < 0.)
            {
                t1t = t1t + 2048. * fClockFreq;
            }

            while (t2t - t2l < 0.)
            {
                t2t = t2t + 2048. * fClockFreq;
            }
            while (t1l - timeLos[0][0] < 0.)
            {
                t1t = t1t + 2048. * fClockFreq;
                t1l = t1l + 2048. * fClockFreq;
                t2t = t2t + 2048. * fClockFreq;
                t2l = t2l + 2048. * fClockFreq;
            }

            tot1 = t1t - t1l;
            if (tot1 < 0)
            {
                LOG(warn) << "Negative ToT " << tot1;
                LOG(warn) << "times1: " << t1t << " " << t1l;
            }

            tot2 = t2t - t2l;
            if (tot2 < 0)
            {
                LOG(warn) << "Negative ToT " << tot2;
                LOG(warn) << "times2: " << t2t << " " << t2l;
            }

            if (bar_quer1 && bar_quer2)
            {
                //				fh_ptof_TotPm1[iBar]->Fill(tot1);
                //				fh_ptof_TotPm2[iBar]->Fill(tot2);
                //				if(iBar==2) fh_ptof_test2->Fill(sqrt(tot1*tot2));
                fh_ptof_channels_cut->Fill(iBar);
            }
        }
    }

    fNEvents += 1;
}

void R3BOnlineSpectraDec2019::FinishEvent()
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

void R3BOnlineSpectraDec2019::FinishTask()
{

    if (fMappedItems.at(DET_LOS))
    {
        for (Int_t iloscount = 0; iloscount < fNofLosDetectors; iloscount++)
        {
            fh_los_channels[iloscount]->Write();
            fh_los_tot[iloscount]->Write();
            fh_los_dt_hits_ToT[iloscount]->Write();
            fh_los_ihit_ToT[iloscount]->Write();
            fh_los_dt_hits[iloscount]->Write();
            fh_dt_hits[iloscount]->Write();
            fh_los_pos_MCFD[iloscount]->Write();
            fh_los_pos_TAMEX[iloscount]->Write();
            fh_los_pos_ToT[iloscount]->Write();
        }
    }

    if (fMappedItems.at(DET_BMON))
    {
        fhTpat->Write();
        fh_spill_length->Write();
        fhTrigger->Write();
        fh_spill_times->Write();
        fh_spill_times_pois->Write();
        fh_spill_times_daq->Write();
        fh_spill_times_daq_first->Write();
        fh_rate->Write();
        fh_frequencies->Write();
        fh_TOFDOR->Write();
        fh_particle->Write();
        fh_particle_pois->Write();
        fh_TOFDOR_spill->Write();
    }

    if (fCalItems.at(DET_TOFD))
    {
        fh_TimePreviousEvent->Write();

        for (Int_t i = 0; i < 4; i++)
        {
            fh_tofd_TotPm[i]->Write();
        }
    }
    for (Int_t ifibcount = 0; ifibcount < NOF_FIB_DET; ifibcount++)
    {
        if (fMappedItems.at(ifibcount + DET_FI_FIRST))
        {
            fh_channels_Fib[ifibcount]->Write();
            fh_fibers_Fib[ifibcount]->Write();
            fh_mult_Fib[ifibcount]->Write();
            fh_time_Fib[ifibcount]->Write();
            fh_multihit_m_Fib[ifibcount]->Write();
            fh_multihit_s_Fib[ifibcount]->Write();
            fh_ToT_m_Fib[ifibcount]->Write();
            fh_Fib_pos_xy[ifibcount]->Write();
            fh_ToT_s_Fib[ifibcount]->Write();
            fh_ToT_single_Fib[ifibcount]->Write();
            fh_Fib_vs_Events[ifibcount]->Write();
            fh_Fib_ToF[ifibcount]->Write();
            fh_channels_single_Fib[ifibcount]->Write();
        }
    }

    if (fMappedItems.at(DET_PSPX))
    {
        // LOG(info) << "Finish MappedPspx";

        for (UInt_t i = 0; i < N_PSPX; i++)
        {
            fh_pspx_channel_x[i]->Write();
            fh_pspx_channel_y[i]->Write();
            fh_pspx_multiplicity_x[i]->Write();
            fh_pspx_multiplicity_y[i]->Write();
            fh_pspx_strips_position[i]->Write();
        }
    }
    if (fCalItems.at(DET_PSPX))
    {
        for (UInt_t i = 0; i < N_PSPX; i++)
        {
            fh_pspx_cal_energy_frontback[i]->Write();
        }
    }
    if (fHitItems.at(DET_PSPX))
    {
        for (UInt_t i = 0; i < (N_PSPX + 1) / 2; i++)
        {
            fh_pspx_hit_energy[i]->Write();
            fh_pspx_hit_position[i]->Write();
        }
    }
}

ClassImp(R3BOnlineSpectraDec2019)
