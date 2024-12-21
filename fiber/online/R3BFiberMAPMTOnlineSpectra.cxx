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

// -------------------------------------------------------------
// -----      R3BFiberMAPMTOnlineSpectra source file       -----
// -----    Created 25/04/22 by J.L. Rodriguez-Sanchez     -----
// -------------------------------------------------------------

#include "R3BFiberMAPMTOnlineSpectra.h"
#include "R3BEventHeader.h"
#include "R3BFiberMAPMTCalData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFiberMappedData.h"
#include "R3BFiberMappingPar.h"
#include "R3BLogger.h"
#include "R3BShared.h"
#include "R3BTCalEngine.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TFolder.h"
#include "TGaxis.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TMath.h"
#include <TRandom3.h>
#include <algorithm>
#include <map>
#include <vector>

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;
namespace
{
    double c_fiber_coincidence_ns = 20; // nanoseconds.
    double c_tot_coincidence_ns = 100;  // nanoseconds
} // namespace

R3BFiberMAPMTOnlineSpectra::R3BFiberMAPMTOnlineSpectra()
    : R3BFiberMAPMTOnlineSpectra("FiberOnlineSpectra", 1)
{
}

R3BFiberMAPMTOnlineSpectra::R3BFiberMAPMTOnlineSpectra(const TString name, Int_t iVerbose)
    : FairTask(name + "OnlineSpectra", iVerbose)
    , fName(name)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fClockFreq(150.)
    , fClockPeriods(4096.)
    , fNEvents(0)
    , fChannelArray()
    , fMapPar(NULL)
    , fHitItems(NULL)
    , fCalItems(NULL)
    , fCalTriggerItems(NULL)
    , fMappedItems(NULL)
    , fNbfibersplot(520)
    , fNbfibers(512)
{
}

R3BFiberMAPMTOnlineSpectra::~R3BFiberMAPMTOnlineSpectra()
{
    if (fh_channels_Fib)
        delete fh_channels_Fib;
    if (fh_fibers_Fib)
        delete fh_fibers_Fib;
    if (fh_mult_Fib)
        delete fh_mult_Fib;
    if (fh_time_Fib)
        delete fh_time_Fib;
    if (fh_multihit_m_Fib)
        delete fh_multihit_m_Fib;
    if (fh_multihit_s_Fib)
        delete fh_multihit_s_Fib;
    if (fh_ToT_Fib)
        delete fh_ToT_Fib;
    if (fh_channels_single_Fib)
        delete fh_channels_single_Fib;
    if (fh_Fib_pos)
        delete fh_Fib_pos;
    if (fh_Fib_vs_Events)
        delete fh_Fib_vs_Events;
    if (fh_ToTup_vs_ToTdown)
        delete fh_ToTup_vs_ToTdown;
    if (fh_chan_corell)
        delete fh_chan_corell;
    if (fh_raw_tot_up)
        delete fh_raw_tot_up;
    if (fh_raw_tot_down)
        delete fh_raw_tot_down;
}

void R3BFiberMAPMTOnlineSpectra::SetParContainers()
{
    fMapPar = dynamic_cast<R3BFiberMappingPar*>(FairRuntimeDb::instance()->getContainer(fName + "MappingPar"));
    R3BLOG_IF(error, !fMapPar, "Couldn't get " << fName << "MappingPar");
    if (fMapPar)
    {
        fNbfibers = fMapPar->GetNbChannels();
        R3BLOG(info, fName << "MappingPar found with " << fNbfibers << " fibers");
    }
}

InitStatus R3BFiberMAPMTOnlineSpectra::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

InitStatus R3BFiberMAPMTOnlineSpectra::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    R3BLOG(info, "For firber " << fName);
    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(warn, header == nullptr, "EventHeader. not found");
    R3BLOG_IF(info, header, " EventHeader. found");

    // uncomment lines below when ucesb avaliable
    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // Get data levels
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject(fName + "Mapped"));
    R3BLOG_IF(fatal, fMappedItems == nullptr, fName + "Mapped not found");

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject(fName + "Cal"));
    R3BLOG_IF(fatal, fCalItems == nullptr, fName + "Cal not found");

    fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject(fName + "TriggerCal"));
    R3BLOG_IF(fatal, fCalTriggerItems == nullptr, fName + "TriggerCal not found");

    fHitItems = dynamic_cast<TClonesArray*>(mgr->GetObject(fName + "Hit"));
    R3BLOG_IF(warn, fHitItems == nullptr, fName + "Hit not found");

    //------------------------------------------------------------------------
    // create histograms
    //------------------------------------------------------------------------

    //-----------------------------------------------------------------------
    // Fiber Detector
    char canvName[255];
    auto mainfolder = new TFolder(fName, fName + " info");
    TCanvas* FibCanvas;

    if (fCalItems)
    {
        FibCanvas = new TCanvas(fName, fName, 10, 10, 910, 910);
        // Cal level
        // Channels:
        fh_channels_Fib = R3B::root_owned<TH1F>(
            fName + "Mapped_channels_up", fName + "Mapped channels up", fNbfibersplot, 1., fNbfibersplot + 1.);
        fh_channels_Fib->GetXaxis()->SetTitle("PMT channel number");
        fh_channels_Fib->GetYaxis()->SetTitle("Counts");
        fh_channels_Fib->SetFillColor(31);

        // Channels:
        fh_channels_single_Fib = R3B::root_owned<TH1F>(
            fName + "Mapped_channels_down", fName + "Mapped channels down", fNbfibersplot, 1., fNbfibersplot + 1.);
        fh_channels_single_Fib->GetXaxis()->SetTitle("PMT channel number");
        fh_channels_single_Fib->GetYaxis()->SetTitle("Counts");
        fh_channels_single_Fib->SetFillColor(31);

        // Multihit MAPMT:
        fh_multihit_m_Fib = R3B::root_owned<TH2F>(
            fName + "Cal_multihit_up", fName + "Cal multihits up", fNbfibersplot, 1., fNbfibersplot + 1., 20, 0., 20.);
        fh_multihit_m_Fib->GetXaxis()->SetTitle("Fiber number");
        fh_multihit_m_Fib->GetYaxis()->SetTitle("Multihit");

        // Multihit SAPMT:
        fh_multihit_s_Fib = R3B::root_owned<TH2F>(fName + "Cal_multihit_down",
                                                  fName + "Cal multihits down",
                                                  fNbfibersplot,
                                                  1.,
                                                  fNbfibersplot + 1.,
                                                  20,
                                                  0.,
                                                  20.);
        fh_multihit_s_Fib->GetXaxis()->SetTitle("Fiber number");
        fh_multihit_s_Fib->GetYaxis()->SetTitle("Multihit");

        // ToT raw up:
        fh_raw_tot_up = R3B::root_owned<TH2F>(
            fName + "Cal_totraw_up", fName + "Cal ToT raw up", fNbfibersplot, 1., fNbfibersplot + 1., 200, 0., 100.);
        fh_raw_tot_up->GetXaxis()->SetTitle("Fiber number");
        fh_raw_tot_up->GetYaxis()->SetTitle("ToT / ns");

        // ToT SAPMT:
        fh_raw_tot_down = R3B::root_owned<TH2F>(
            fName + "Cal_totraw_down", fName + "Cal ToT raw down", fNbfibersplot, 1., fNbfibersplot + 1, 200, 0., 100.);
        fh_raw_tot_down->GetXaxis()->SetTitle("Fiber number");
        fh_raw_tot_down->GetYaxis()->SetTitle("ToT / ns");

        fh_chan_corell = R3B::root_owned<TH2F>(
            fName + "Cal_chan_vs_chan", fName + "Cal ch_up vs ch_down", 512, 0., 512, 512, 0., 512.);
        fh_chan_corell->GetXaxis()->SetTitle("Fiber number down");
        fh_chan_corell->GetYaxis()->SetTitle("Fiber number up");

        fh_chan_dt_cal = R3B::root_owned<TH2F>(fName + "Cal_TimevsChannel",
                                               fName + "Cal Time bottom and top vs fiber number",
                                               2 * fNbfibersplot,
                                               -fNbfibersplot,
                                               fNbfibersplot,
                                               1700,
                                               -1500.,
                                               200.);
        fh_chan_dt_cal->GetXaxis()->SetTitle("Fiber number");
        fh_chan_dt_cal->GetYaxis()->SetTitle("time(up/down)-time_trigger / ns");

        // Hit level
        // Fibers:
        fh_fibers_Fib =
            R3B::root_owned<TH1F>(fName + "Hit_fibers", fName + "Hit fibers", fNbfibersplot, 1., fNbfibersplot + 1.);
        fh_fibers_Fib->GetXaxis()->SetTitle("Fiber number");
        fh_fibers_Fib->GetYaxis()->SetTitle("Counts");
        fh_fibers_Fib->SetFillColor(31);

        // Multiplicity (number of hit fibers):
        fh_mult_Fib = R3B::root_owned<TH1F>(fName + "Hit_mult", fName + "Hit mult", 150, 0., 150);
        fh_mult_Fib->GetXaxis()->SetTitle("Multiplicity");
        fh_mult_Fib->GetYaxis()->SetTitle("Counts");
        fh_mult_Fib->SetFillColor(31);

        // ToT :
        fh_ToT_Fib = R3B::root_owned<TH2F>(
            fName + "Hit_ToT_iFib", fName + "Hit ToT vs Fib", fNbfibersplot, 1., fNbfibersplot + 1, 200, 0., 100.);
        fh_ToT_Fib->GetXaxis()->SetTitle("Fiber number");
        fh_ToT_Fib->GetYaxis()->SetTitle("ToT / ns");

        // Time of fiber:
        fh_time_Fib = R3B::root_owned<TH2F>(fName + "Hit_TimevsFiber",
                                            fName + "Hit Time vs Fiber",
                                            fNbfibersplot,
                                            1.,
                                            fNbfibersplot + 1.,
                                            1000,
                                            -10000,
                                            10000);
        fh_time_Fib->GetXaxis()->SetTitle("Fiber number");
        fh_time_Fib->GetYaxis()->SetTitle("(tUp+tDown)/2 / ns");

        // calibrated position
        auto cxy = new TCanvas(fName + "_xy", fName + "_xy", 10, 10, 910, 910);
        fh_Fib_pos =
            R3B::root_owned<TH2F>(fName + "Hit_pos", fName + "Hit calibrated position", 300, -30, 30, 500, -500, 500);
        fh_Fib_pos->GetXaxis()->SetTitle("x position");
        fh_Fib_pos->GetYaxis()->SetTitle("y position");
        cxy->cd();
        fh_Fib_pos->Draw("colz");

        // hit fiber number vs. event number:
        fh_Fib_vs_Events = R3B::root_owned<TH2F>(fName + "Hit_fib_vs_event",
                                                 fName + "Hit Fiber # vs. Event #",
                                                 100,
                                                 0,
                                                 5e6,
                                                 fNbfibersplot,
                                                 1.,
                                                 fNbfibersplot + 1.);
        fh_Fib_vs_Events->GetYaxis()->SetTitle("Fiber number");
        fh_Fib_vs_Events->GetXaxis()->SetTitle("Event number");

        fh_ToTup_vs_ToTdown = R3B::root_owned<TH2F>(
            fName + "Hit_time_vs_event", fName + "Hit time # vs. Event #", 200, 0, 5e6, 2000, -1000., 1000.);
        fh_ToTup_vs_ToTdown->GetXaxis()->SetTitle("Event number");
        fh_ToTup_vs_ToTdown->GetYaxis()->SetTitle("Time  / ns");

        FibCanvas->Divide(4, 3);

        FibCanvas->cd(1);
        gPad->SetLogy();
        fh_channels_Fib->Draw();
        FibCanvas->cd(2);
        gPad->SetLogy();
        fh_channels_single_Fib->Draw();

        if (fCalItems)
        {
            FibCanvas->cd(3);
            gPad->SetLogz();
            fh_multihit_m_Fib->Draw("colz");
            FibCanvas->cd(4);
            gPad->SetLogz();
            fh_multihit_s_Fib->Draw("colz");
            FibCanvas->cd(5);
            gPad->SetLogz();
            fh_chan_corell->Draw("colz");
            FibCanvas->cd(6);
            gPad->SetLogz();
            fh_chan_dt_cal->Draw("colz");
            FibCanvas->cd(7);
            gPad->SetLogz();
            fh_raw_tot_up->Draw("colz");
            FibCanvas->cd(8);
            gPad->SetLogz();
            fh_raw_tot_down->Draw("colz");
        }
        if (fHitItems)
        {
            FibCanvas->cd(9);
            gPad->SetLogy();
            fh_fibers_Fib->Draw();
            FibCanvas->cd(10);
            gPad->SetLogy();
            fh_mult_Fib->Draw();
            FibCanvas->cd(12);
            gPad->SetLogz();
            fh_ToT_Fib->Draw("colz");

            FibCanvas->cd(11);
            gPad->SetLogz();
            fh_time_Fib->Draw("colz");
        }

        mainfolder->Add(FibCanvas);
        mainfolder->Add(cxy);
    } // end if(Mapped)

    run->AddObject(mainfolder);
    run->GetHttpServer()->RegisterCommand("Reset_" + fName, Form("/Objects/%s/->Reset_Histo()", GetName()));

    return kSUCCESS;
}
void R3BFiberMAPMTOnlineSpectra::Reset_Histo()
{
    R3BLOG(info, "For firber " << fName);
    if (fCalItems)
    {
        fh_channels_Fib->Reset();
        fh_mult_Fib->Reset();
        fh_multihit_m_Fib->Reset();
        fh_multihit_s_Fib->Reset();
        fh_chan_corell->Reset();
        fh_channels_single_Fib->Reset();
        fh_raw_tot_up->Reset();
        fh_raw_tot_down->Reset();
        fh_chan_dt_cal->Reset();
    }

    if (fHitItems)
    {
        fh_time_Fib->Reset();
        fh_fibers_Fib->Reset();
        fh_ToT_Fib->Reset();
        fh_Fib_vs_Events->Reset();
        fh_Fib_pos->Reset();
        fh_Fib_vs_Events->Reset();
        fh_ToTup_vs_ToTdown->Reset();
    }
}
void R3BFiberMAPMTOnlineSpectra::Exec(Option_t* option)
{
    fNEvents += 1;

    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    // fTpat = 1-16; fTpat_bit = 0-15
    if (fTpat1 > -1 && fTpat2 > -1)
    {
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
    }

    Int_t iFib = 0;
    Double_t dtime = 0.0 / 0.0, dtimeMax;
    Double_t yfib = -100000., xfib = -100000., tfibMax = 0. / 0.;
    Double_t rfib = 0.0 / 0.0;
    Int_t iFibMax;
    Double_t totMax;
    Double_t tofMax;
    Double_t tot_up = 0., tot_down = 0.;
    Int_t iCha = 0;
    Int_t iCha_up = 0, iCha_down = 0;
    Double_t tFib = 0. / 0.;
    Int_t iSide = 0;
    Double_t c_period = fClockPeriods * (1000. / fClockFreq);

    if (fMappedItems && fMappedItems->GetEntriesFast() > 0)
    {
        auto nMapp = fMappedItems->GetEntries();
        for (Int_t i = 0; i < nMapp; i++)
        {
            auto map_lead = dynamic_cast<R3BFiberMappedData const*>(fMappedItems->At(i));

            if (map_lead->IsLeading())
            {
                auto side_i = map_lead->GetSide() - 1;
                auto ch_i = map_lead->GetChannel() - 1;

                if (side_i == 1)
                {
                    fh_channels_Fib->Fill(ch_i); // Fill which channel has events
                }

                if (side_i == 0)
                {
                    fh_channels_single_Fib->Fill(ch_i); // Fill which channel has events
                }
            }
        }
    }

    if (fCalItems && fCalItems->GetEntriesFast() > 0)
    {
        UInt_t vmultihits_top[fNbfibers], vmultihits_bot[fNbfibers];
        for (Int_t i = 0; i < fNbfibers; i++)
        {
            vmultihits_top[i] = 0;
            vmultihits_bot[i] = 0;
        }
        // Resize per-channel info arrays.
        for (auto i = 0; i < 2; ++i)
        {
            fChannelArray[i].resize(fNbfibers);
        }

        double trig_time[8] = { 0 };

        //------ Collecting cal trigger hits --------
        size_t cal_num = fCalTriggerItems->GetEntriesFast();
        Double_t tl, tt; // lead and trail times of the trigger
        for (UInt_t j = 0; j < cal_num; ++j)
        {
            auto cur_cal = dynamic_cast<R3BFiberMAPMTCalData const*>(fCalTriggerItems->At(j));
            auto ch = cur_cal->GetChannel() - 1;
            tl = cur_cal->GetTime_ns();
            trig_time[ch] = tl;
        }

        Int_t nCals = fCalItems->GetEntries();

        for (auto side_i = 0; side_i < 2; ++side_i)
        {
            // Clear local helper containers.
            auto& array = fChannelArray[side_i];
            for (auto it = array.begin(); array.end() != it; ++it)
            {
                it->lead_list.clear();
            }
        }

        for (size_t j = 0; j < nCals; ++j)
        {
            auto cur_cal_lead = dynamic_cast<R3BFiberMAPMTCalData const*>(fCalItems->At(j));

            if (cur_cal_lead->IsLeading())
            {
                auto side_i = cur_cal_lead->GetSide() - 1;
                auto ch_i = cur_cal_lead->GetChannel() - 1;

                if (side_i == 1)
                {
                    //   fh_channels_Fib->Fill(ch_i); // Fill which channel has events
                    vmultihits_top[ch_i] += 1; // multihit of a given up killom channel
                }

                if (side_i == 0)
                {
                    //  fh_channels_single_Fib->Fill(ch_i); // Fill which channel has events
                    vmultihits_bot[ch_i] += 1; // multihit of a given down killom channel
                }

                auto time_trig = trig_time[fMapPar->GetTrigMap(side_i + 1, ch_i + 1)];
                auto time_ns =
                    fmod(cur_cal_lead->GetTime_ns() - time_trig + c_period + c_period / 2, c_period) - c_period / 2;

                if (side_i == 0)
                    fh_chan_dt_cal->Fill(-ch_i - 1, time_ns);
                if (side_i == 1)
                    fh_chan_dt_cal->Fill(ch_i + 1, time_ns);

                auto& channel = fChannelArray[side_i].at(ch_i);

                channel.lead_list.push_back(cur_cal_lead); // make list with leading times
            }
        }

        for (int i = 0; i < fNbfibers; ++i)
        {

            if (vmultihits_top[i] > 0)
                fh_multihit_m_Fib->Fill(i + 1,
                                        vmultihits_top[i]); // multihit of a given up killom channel

            if (vmultihits_bot[i] > 0)
                fh_multihit_s_Fib->Fill(i + 1,
                                        vmultihits_bot[i]); // multihit of a given down killom channel
        }

        for (size_t j = 0; j < nCals; ++j)
        {
            auto cur_cal_trail = dynamic_cast<R3BFiberMAPMTCalData const*>(fCalItems->At(j));

            if (cur_cal_trail->IsTrailing())
            {
                auto side_i = cur_cal_trail->GetSide() - 1;
                auto ch_i = cur_cal_trail->GetChannel() - 1;

                auto& channel = fChannelArray[side_i].at(ch_i);
                if (channel.lead_list.empty())
                {
                    continue;
                }
                auto lead = channel.lead_list.front();
                Double_t lead_raw = 0;
                Double_t trail_raw = 0;

                lead_raw = lead->GetTime_ns();
                trail_raw = cur_cal_trail->GetTime_ns();

                auto time_trig_lead = trig_time[fMapPar->GetTrigMap(lead->GetSide(), lead->GetChannel())];
                auto tlead_ns =
                    fmod(lead->GetTime_ns() - time_trig_lead + c_period + c_period / 2, c_period) - c_period / 2;

                auto time_trig_trail = trig_time[fMapPar->GetTrigMap(side_i + 1, ch_i + 1)];
                auto ttrail_ns =
                    fmod(cur_cal_trail->GetTime_ns() - time_trig_trail + c_period + c_period / 2, c_period) -
                    c_period / 2;

                auto tot_ns = fmod(ttrail_ns - tlead_ns + c_period + c_period / 2, c_period) - c_period / 2;

                if (tot_ns < c_tot_coincidence_ns && tot_ns > 0.)
                {
                    if (cur_cal_trail->GetSide() == 2)
                        fh_raw_tot_up->Fill(cur_cal_trail->GetChannel(), tot_ns);
                    if (cur_cal_trail->GetSide() == 1)
                        fh_raw_tot_down->Fill(cur_cal_trail->GetChannel(), tot_ns);
                }
            }
        }

        for (size_t j = 0; j < nCals; ++j)
        {
            auto cur_cal_top = dynamic_cast<R3BFiberMAPMTCalData const*>(fCalItems->At(j));
            if (cur_cal_top->IsLeading() && cur_cal_top->GetSide() == 2) // choose leading top times
            {
                auto cht_i = cur_cal_top->GetChannel() - 1;
                auto time_trig_top = trig_time[fMapPar->GetTrigMap(cur_cal_top->GetSide(), cht_i + 1)];
                auto time_top =
                    fmod(cur_cal_top->GetTime_ns() - time_trig_top + c_period + c_period / 2, c_period) - c_period / 2;
                for (size_t k = 0; k < nCals; ++k)
                {
                    auto cur_cal_bot = dynamic_cast<R3BFiberMAPMTCalData const*>(fCalItems->At(k));
                    if (cur_cal_bot->IsLeading() && cur_cal_bot->GetSide() == 1) // choose leading bottom times
                    {
                        auto chb_i = cur_cal_bot->GetChannel() - 1;
                        auto time_trig_bot = trig_time[fMapPar->GetTrigMap(cur_cal_bot->GetSide(), chb_i + 1)];
                        auto time_bot =
                            fmod(cur_cal_bot->GetTime_ns() - time_trig_bot + c_period + c_period / 2, c_period) -
                            c_period / 2;
                        auto dt = time_top - time_bot;
                        auto dt_mod = fmod(dt + c_period, c_period);

                        if (dt < 0)
                        {
                            // We're only interested in the short time-differences, so we
                            // want to move the upper part of the coarse counter range close
                            // to the lower range, i.e. we cut the middle of the range and
                            // glue zero and the largest values together.
                            dt_mod -= c_period;
                        }
                        if (std::abs(dt_mod) < c_fiber_coincidence_ns)
                            fh_chan_corell->Fill(cur_cal_bot->GetChannel(), cur_cal_top->GetChannel());
                    }
                }
            }
        }

        fChannelArray[0].clear();
        fChannelArray[1].clear();
    } // if Cal

    if (fHitItems && fHitItems->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitItems->GetEntriesFast();
        Double_t xposMax = 0. / 0.;
        Double_t yposMax = 0. / 0.;
        Double_t tfib = 0. / 0.;
        iFibMax = -1000;
        totMax = 0.;
        dtimeMax = 0. / 0.;

        // if (nHits > 1)
        //  return;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            Double_t tDown = 0. / 0.;
            Double_t tUp = 0. / 0.;

            auto hit = dynamic_cast<R3BFiberMAPMTHitData*>(fHitItems->At(ihit));
            if (!hit)
                continue;

            iFib = hit->GetFiberId(); // 1..

            // times
            tUp = hit->GetTopTime_ns();
            tDown = hit->GetBottomTime_ns();
            dtime = tUp - tDown;

            // (tUp+tDown)/2:
            tfib = hit->GetTime();

            // ToT from up and down MAPMT:
            Double_t ToT_up = hit->GetTopToT_ns();
            Double_t ToT_down = hit->GetBottomToT_ns();

            Double_t ToT = hit->GetEloss();
            Double_t xpos = hit->GetX();
            Double_t ypos = hit->GetY();

            fh_Fib_pos->Fill(xpos, ypos);

            if (ToT > totMax)
            {
                totMax = ToT;
                iFibMax = iFib;
                tfibMax = tfib;
                dtimeMax = dtime;
                xposMax = xpos;
                yposMax = ypos;
            }

            fh_fibers_Fib->Fill(iFib);
            fh_ToT_Fib->Fill(iFib, ToT);

            fh_ToTup_vs_ToTdown->Fill(fNEvents, tfib);

            fh_time_Fib->Fill(iFib, tfib);
            fh_Fib_vs_Events->Fill(fNEvents, iFib);

        } // end for(ihit)

        if (nHits > 0)
            fh_mult_Fib->Fill(nHits);

    } // end if(aHit)
}

void R3BFiberMAPMTOnlineSpectra::FinishEvent()
{
    if (fMappedItems)
    {
        fMappedItems->Clear();
    }
    if (fCalItems)
    {
        fCalItems->Clear();
    }
    if (fCalTriggerItems)
    {
        fCalTriggerItems->Clear();
    }
    if (fHitItems)
    {
        fHitItems->Clear();
    }
}

void R3BFiberMAPMTOnlineSpectra::FinishTask()
{
    if (fCalItems)
    {
        fh_channels_Fib->Write();
        fh_mult_Fib->Write();
        fh_multihit_m_Fib->Write();
        fh_multihit_s_Fib->Write();
        fh_channels_single_Fib->Write();
        fh_ToTup_vs_ToTdown->Write();
    }
    if (fHitItems)
    {
        fh_time_Fib->Write();
        fh_fibers_Fib->Write();
        fh_ToT_Fib->Write();
        fh_Fib_vs_Events->Write();
        fh_Fib_pos->Write();
        fh_Fib_vs_Events->Write();
        fh_chan_corell->Write();
        fh_raw_tot_up->Write();
        fh_raw_tot_down->Write();
    }
}

ClassImp(R3BFiberMAPMTOnlineSpectra)
