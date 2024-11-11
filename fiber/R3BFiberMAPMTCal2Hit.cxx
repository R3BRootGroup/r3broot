/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "mapping_fib23a_trig_hit.hh"
#include "mapping_fib23b_trig_hit.hh"
#include "mapping_fib30_trig.hh"
#include "mapping_fib31_trig.hh"
#include "mapping_fib32_trig.hh"
#include "mapping_fib33_trig.hh"

#include "R3BFiberMAPMTCal2Hit.h"
#include "R3BFiberMAPMTCalData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFiberMAPMTHitPar.h"
#include "R3BTCalEngine.h"

#include "TH1F.h"
#include "TH2F.h"
#include <TClonesArray.h>
#include <cassert>

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <iostream>

R3BFiberMAPMTCal2Hit::ToT::ToT(R3BFiberMAPMTCalData const* a_lead,
                               R3BFiberMAPMTCalData const* a_trail,
                               Double_t a_lead_ns,
                               Double_t a_tail_ns,
                               Double_t a_tot_ns)
    : lead(a_lead)
    , trail(a_trail)
    , lead_ns(a_lead_ns)
    , tail_ns(a_tail_ns)
    , tot_ns(a_tot_ns)
{
}

R3BFiberMAPMTCal2Hit::R3BFiberMAPMTCal2Hit(const char* a_name,
                                           Int_t a_verbose,
                                           Direction a_direction,
                                           UInt_t a_num_fibers,
                                           Bool_t a_is_calibrator)
    : FairTask(TString("R3B") + a_name + "Cal2Hit", a_verbose)
    , fName(a_name)
    , fClockFreq(150.)
    , fDirection(a_direction)
    , fNumFibers(a_num_fibers)
    , fIsCalibrator(a_is_calibrator)
    , fCalItems()
    , fCalTriggerItems()
    , fHitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fCalPar()
    , fNofHitPars()
    , fNofHitItems()
    , fChannelArray()
    , fnEvents(0)
    , ftofmin(TOF_MIN)
    , ftofmax(TOF_MAX)
    , fWrite(true)
    , fSimu(false)
    , fGate_ns(100.) // make it constructor parameter
{
    // Trigger mappin from included mapping files for each det.
    if (fName == "Fi30")
    {
        fib30_trig_map_setup();
        fTriggerMap[0] = g_fib30_trig_map[0];
        fTriggerMap[1] = g_fib30_trig_map[1];
    }
    if (fName == "Fi31")
    {
        fib31_trig_map_setup();
        fTriggerMap[0] = g_fib31_trig_map[0];
        fTriggerMap[1] = g_fib31_trig_map[1];
    }

    if (fName == "Fi32")
    {
        fib32_trig_map_setup();
        fTriggerMap[0] = g_fib32_trig_map[0];
        fTriggerMap[1] = g_fib32_trig_map[1];
    }
    if (fName == "Fi33")
    {
        fib33_trig_map_setup();
        fTriggerMap[0] = g_fib33_trig_map[0];
        fTriggerMap[1] = g_fib33_trig_map[1];
    }
    if (fName == "Fi23a")
    {
        fib23a_trig_map_hit_setup();
        fTriggerMap[0] = g_fib23a_trig_map_hit[0];
        fTriggerMap[1] = g_fib23a_trig_map_hit[1];
    }
    if (fName == "Fi23b")
    {
        fib23b_trig_map_hit_setup();
        fTriggerMap[0] = g_fib23b_trig_map_hit[0];
        fTriggerMap[1] = g_fib23b_trig_map_hit[1];
    }
}

R3BFiberMAPMTCal2Hit::~R3BFiberMAPMTCal2Hit()
{
    delete fHitItems;
    delete fCalPar;
}

InitStatus R3BFiberMAPMTCal2Hit::Init()
{
    std::cout << "Init:R3BFiberMAPMTCal2Hit" << std::endl;
    auto mgr = FairRootManager::Instance();
    if (!mgr)
    {
        LOG(error) << "FairRootManager not found.";
        return kERROR;
    }

    auto name = fName + "Cal";
    fCalItems = (TClonesArray*)mgr->GetObject(name);
    if (!fCalItems)
    {
        LOG(fatal) << "Branch " << name << " not found.";
        return kERROR;
    }
    auto name_mapmt_trig = fName + "TriggerCal";
    fCalTriggerItems = (TClonesArray*)mgr->GetObject(name_mapmt_trig);
    if (NULL == fCalTriggerItems)
        LOG(fatal) << "Branch " << name_mapmt_trig << " not found";

    maxevent = mgr->CheckMaxEventNo();

    // TClones branch with Hit items
    mgr->Register(fName + "Hit", "Land", fHitItems, kTRUE);

    // Resize per-channel info arrays.
    for (auto side_i = 0; side_i < 2; ++side_i)
    {
        fChannelArray[side_i].resize(fNumFibers);
    }

    if (!fIsCalibrator)
    {
        // Get calibration parameters if we're not a calibrator.
        auto container = fName + "HitPar";
        fHitPar = (R3BFiberMAPMTHitPar*)FairRuntimeDb::instance()->getContainer(container);
        if (!fHitPar)
        {
            LOG(error) << "Could not get " << container << " container.";
            fNofHitPars = 0;
        }
        else
        {
            fNofHitPars = fHitPar->GetNumModulePar();
            if (0 == fNofHitPars)
            {
                LOG(error) << "No Hit parameters in " << container << " container.";
                fHitPar = nullptr;
            }
        }
    }

    // create histograms
    TString chistName;
    TString chistTitle;
    // ToT bottom PMT raw:
    chistName = fName + "_ToT_bottom_raw";
    chistTitle = fName + " ToTbottom raw of fibers";
    fh_ToT_bottom_Fib_raw = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 200, 0., 100.);
    fh_ToT_bottom_Fib_raw->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_bottom_Fib_raw->GetYaxis()->SetTitle("ToT / ns");

    // ToT bottom gainmatched:
    chistName = fName + "_ToT_bottom";
    chistTitle = fName + " ToTbottom of fibers";
    fh_ToT_bottom_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 200, 0., 100.);
    fh_ToT_bottom_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_bottom_Fib->GetYaxis()->SetTitle("ToT / ns");

    // ToF Tofd -> Fiber:
    chistName = fName + "_time";
    chistTitle = fName + " Time";
    fh_Fib_ToF = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 2200, -1100., 1100.);
    fh_Fib_ToF->GetYaxis()->SetTitle("Time/ ns");
    fh_Fib_ToF->GetXaxis()->SetTitle("Fiber ID");

    // ToF Fiber for sync:
    chistName = fName + "_time_raw";
    chistTitle = fName + " Time raw";
    fh_Fib_ToF_raw = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 2400, -1000., 200.);
    fh_Fib_ToF_raw->GetYaxis()->SetTitle("Time/ ns");
    fh_Fib_ToF_raw->GetXaxis()->SetTitle("Fiber ID");

    // ToT top MAPMT raw:
    chistName = fName + "_ToT_top_raw";
    chistTitle = fName + " ToTtop raw of fibers";
    fh_ToT_top_Fib_raw = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 200, 0., 100.);
    fh_ToT_top_Fib_raw->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_top_Fib_raw->GetYaxis()->SetTitle("ToT / ns");

    // ToT top MAPMT gainmatched:
    chistName = fName + "_ToT_top";
    chistTitle = fName + " ToTtop of fibers";
    fh_ToT_top_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 200, 0., 100.);
    fh_ToT_top_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_top_Fib->GetYaxis()->SetTitle("ToT / ns");

    chistName = fName + "_time_top";
    chistTitle = fName + " time of top MAPMTs ";
    fh_time_top_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 2400, -1000., 200.);
    fh_time_top_Fib->GetXaxis()->SetTitle("single number");
    fh_time_top_Fib->GetYaxis()->SetTitle("time / ns");

    // time difference top-bottom:
    chistName = fName + "_dt";
    chistTitle = fName + " dt of fibers";
    fh_dt_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 4400, -1100., 1100.);
    fh_dt_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_dt_Fib->GetYaxis()->SetTitle("dt / ns");

    // time difference top-bottom for offsets:
    chistName = fName + "_dt_raw";
    chistTitle = fName + " dt_raw of fibers";
    fh_dt_Fib_raw = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 4400, -1100., 1100.);
    fh_dt_Fib_raw->GetXaxis()->SetTitle("Fiber number");
    fh_dt_Fib_raw->GetYaxis()->SetTitle("dt / ns");

    // time of MAPMT:
    chistName = fName + "_time_bottom";
    chistTitle = fName + " time of bottom MAPMTs";
    fh_time_bottom_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 2400, -1000., 200.);
    fh_time_bottom_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_time_bottom_Fib->GetYaxis()->SetTitle("time / ns");

    // histos to check hit parameters
    chistName = fName + "_check_tsync";
    chistTitle = fName + " check tsync for fiberNum 10";
    fh_time_check_tsync = new TH2F(chistName.Data(), chistTitle.Data(), 50000, 0, 1e6, 400, -100., 100.);
    fh_time_check_tsync->GetXaxis()->SetTitle("Event number");
    fh_time_check_tsync->GetYaxis()->SetTitle("dt / ns");

    return kSUCCESS;
}

InitStatus R3BFiberMAPMTCal2Hit::ReInit() { return kSUCCESS; }

void R3BFiberMAPMTCal2Hit::SetParContainers()
{
    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCalPar = (R3BFiberMAPMTHitPar*)FairRuntimeDb::instance()->getContainer(fName + "HitPar");
    if (fCalPar)
        cout << "R3BFiberCal2Hit::Init() container " << fName << "HitPar initialized. " << endl;
    if (!fCalPar)
    {
        LOG(error) << "R3BFiberCal2Hit::Init() Couldn't get " << fName << "HitPar. ";
    }
}

void R3BFiberMAPMTCal2Hit::Exec(Option_t* option)
{

    // cout<<"R3BFiberMAPMTCal2Hit::Exec entering"<<endl;

    if (fnEvents / 100000. == (int)fnEvents / 100000)
        std::cout << "\rEvents: " << fnEvents << " / " << maxevent << " (" << (int)(fnEvents * 100. / maxevent)
                  << " %) " << std::flush;

    for (auto side_i = 0; side_i < 2; ++side_i) // 0=top, 1=bottom
    {
        // Clear local helper containers.
        auto& array = fChannelArray[side_i];
        for (auto it = array.begin(); array.end() != it; ++it)
        {
            it->lead_list.clear();
            it->tot_list.clear();
        }
    }

    double trig_time[8];
    //------ Collecting cal trigger hits --------
    size_t cal_num_trig = fCalTriggerItems->GetEntriesFast();
    Double_t tl, tt; // lead and trile times of the trigger
    for (UInt_t j = 0; j < cal_num_trig; ++j)
    {
        auto cur_cal = (R3BFiberMAPMTCalData*)fCalTriggerItems->At(j);
        auto ch = cur_cal->GetChannel() - 1;
        tl = cur_cal->GetTime_ns();
        trig_time[ch] = tl;
    }

    size_t cal_num = fCalItems->GetEntriesFast();

    //   if(cal_num>0 ) cout<<"calNum: "<<fName<<", "<<cal_num<<", "<<endl;

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal_lead = (R3BFiberMAPMTCalData const*)fCalItems->At(j);

        if (cur_cal_lead->IsLeading())
        {
            auto side_i = cur_cal_lead->GetSide();
            auto ch_i = cur_cal_lead->GetChannel() - 1;
            auto& channel = fChannelArray[side_i].at(ch_i);
            channel.lead_list.push_back(cur_cal_lead);
        }
    }

    Double_t c_period = 4096. * (1000. / fClockFreq);

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal_trail = (R3BFiberMAPMTCalData const*)fCalItems->At(j);

        if (cur_cal_trail->IsTrailing())
        {
            auto side_i = cur_cal_trail->GetSide();

            auto ch_i = cur_cal_trail->GetChannel() - 1;

            auto& channel = fChannelArray[side_i].at(ch_i);
            if (channel.lead_list.empty())
            {
                continue;
            }
            auto lead = channel.lead_list.front();

            Double_t cur_cal_trig_ns = 0;
            Double_t lead_trig_ns = 0;
            Double_t lead_raw = 0;
            Double_t trail_raw = 0;
            Double_t trail_trig_ns = 0;
            Double_t lead_ns = 0. / 0., tot_ns = 0. / 0.;

            auto chlead_i = lead->GetChannel() - 1;
            cur_cal_trig_ns = trig_time[fTriggerMap[side_i][ch_i]];
            lead_trig_ns = trig_time[fTriggerMap[side_i][chlead_i]];

            auto cur_cal_ns =
                fmod(cur_cal_trail->GetTime_ns() - cur_cal_trig_ns + c_period + c_period / 2, c_period) - c_period / 2;
            lead_ns = fmod(lead->GetTime_ns() - lead_trig_ns + c_period + c_period / 2, c_period) - c_period / 2;
            tot_ns = fmod(cur_cal_ns - lead_ns + c_period + c_period / 2, c_period) - c_period / 2;

            lead_raw = lead->GetTime_ns();
            trail_raw = cur_cal_trail->GetTime_ns();
            trail_trig_ns = cur_cal_trig_ns;

            if (tot_ns < fGate_ns && tot_ns > 0.)
            {
                //	if(chlead_i == 0 && fName=="Fi31") cout<<"Side: "<<side_i<<", "<<lead_raw<<", "<<lead_trig_ns<<",
                //"<< 	lead_ns<<"; "<<fTriggerMap[side_i][chlead_i]<<endl;

                channel.tot_list.push_back(ToT(lead, cur_cal_trail, lead_ns, cur_cal_ns, tot_ns));
                channel.lead_list.pop_front();
            }
        }
    }

    auto const& down_array = fChannelArray[0];
    auto const& up_array = fChannelArray[1];

    for (auto it_down = down_array.begin(); down_array.end() != it_down; ++it_down) // over down channel 0...nmax-1;
    {
        auto const& down = *it_down;
        for (auto it_down_tot = down.tot_list.begin(); down.tot_list.end() != it_down_tot;
             ++it_down_tot) // over ihit(channel)
        {
            auto const& down_tot = *it_down_tot;
            auto down_sub_id = (down_tot.lead->GetChannel() - 1) / fNumFibers;
            auto fiber_down_ch = down_tot.lead->GetChannel();

            if ((fName == "Fi23a" || fName == "Fi23b") && (fiber_down_ch > 187 && fiber_down_ch < 197))
                continue;

            for (auto it_up = up_array.begin(); up_array.end() != it_up; ++it_up) // over up channel 0...nmax-1;
            {
                auto const& up = *it_up;
                for (auto it_up_tot = up.tot_list.begin(); up.tot_list.end() != it_up_tot;
                     ++it_up_tot) // over ihit(channel)
                {
                    auto const& up_tot = *it_up_tot;

                    // Check that the combo is inside one sub-det.

                    auto up_sub_id = (up_tot.lead->GetChannel() - 1) / fNumFibers;
                    if (down_sub_id != up_sub_id)
                    {
                        continue;
                    }

                    auto fiber_up_ch = up_tot.lead->GetChannel(); // 1...
                    if ((fName == "Fi23a" || fName == "Fi23b") && (fiber_up_ch > 187 && fiber_up_ch < 197))
                        continue;

                    Int_t fiber_id = -1000;
                    if (fiber_down_ch == fiber_up_ch)
                        fiber_id = fiber_down_ch;
                    else
                        continue;

                    // Calibrate hit fiber.
                    auto tot_down_raw = down_tot.tot_ns;
                    auto tot_up_raw = up_tot.tot_ns;
                    auto tot_down = down_tot.tot_ns;
                    auto tot_up = up_tot.tot_ns;
                    Double_t t_down = down_tot.lead_ns;
                    Double_t t_up = up_tot.lead_ns;
                    Double_t dtime = t_up - t_down;
                    Double_t tof = (t_up + t_down) / 2.;

                    // Fill histograms for gain match, offset and sync.
                    fh_ToT_bottom_Fib_raw->Fill(fiber_id, tot_down);
                    fh_ToT_top_Fib_raw->Fill(fiber_id, tot_up);
                    fh_dt_Fib_raw->Fill(fiber_id, dtime);
                    fh_Fib_ToF_raw->Fill(fiber_id, tof);

                    gainUp = 10.;
                    gainDown = 10.;
                    tsync = 0.;
                    offsetUp = 0.;
                    offsetDown = 0.;

                    if (!fIsCalibrator && fHitPar)
                    {
                        R3BFiberMAPMTHitModulePar* par = fHitPar->GetModuleParAt(fiber_id);
                        if (par)
                        {
                            gainUp = par->GetGainUp();
                            gainDown = par->GetGainDown();
                            tsync = par->GetSync();
                            offsetUp = par->GetOffsetUp();
                            offsetDown = par->GetOffsetDown();
                        }
                    }

                    tot_down *= 10. / gainDown;
                    tot_up *= 10. / gainUp;

                    //  tof -= tsync;
                    t_down -= offsetDown;
                    t_down -= tsync;
                    t_up -= offsetUp;
                    t_up -= tsync;
                    dtime = t_up - t_down;
                    tof = (t_up + t_down) / 2.;

                    // histogram after gain match, sync....
                    fh_dt_Fib->Fill(fiber_id, dtime);
                    fh_Fib_ToF->Fill(fiber_id, tof);
                    fh_ToT_bottom_Fib->Fill(fiber_id, tot_down);
                    fh_ToT_top_Fib->Fill(fiber_id, tot_up);
                    fh_time_bottom_Fib->Fill(fiber_id, t_down);
                    fh_time_top_Fib->Fill(fiber_id, t_up);

                    if (fiber_id == 10)
                        fh_time_check_tsync->Fill(fnEvents, tof);

                    Int_t numFibs = fNumFibers;
                    Double_t x = -10000.;
                    Double_t y = -10000.;
                    Double_t veff = 0.333; // cm/ns
                    Double_t randx = (std::rand() / (float)RAND_MAX);
                    if (fName == "Fi23a" || fName == "Fi23b")
                    {
                        Float_t fiber_thickness = 0.028000;
                        Int_t fiber_nbr = 384;
                        Float_t dead_layer = 0.9;
                        Float_t air_layer = 0.01 * 0.; // relative to fiber_thickness
                        Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);
                        if (fDirection == VERTICAL)
                        {
                            x = -detector_width / 2. + fiber_thickness * (1 + air_layer) / 2. +
                                double(fiber_id - 1.) * (1 + air_layer) * fiber_thickness;

                            // if(fSimu) x = x + 0.028;

                            // x = -detector_width / 2. +
                            //   (double(fiber_id - 1) + (double(fiber_id - 1.) * air_layer)) * fiber_thickness;
                            y = (t_down - t_up) * veff;
                        }
                        else
                        {
                            x = (t_down - t_up) * veff;

                            //   y = -detector_width / 2. +
                            //      (double(fiber_id - 1) + (double(fiber_id - 1) * air_layer)) * fiber_thickness;

                            y = -detector_width / 2. + fiber_thickness * (1 + air_layer) / 2. +
                                double(fiber_id - 1.) * (1 + air_layer) * fiber_thickness;

                            // if (fSimu) y = y + 0.028;
                        }
                    }
                    if (fName == "Fi30" || fName == "Fi31" || fName == "Fi32" || fName == "Fi33")
                    {
                        veff = 5. * 0.333;
                        Float_t fiber_thickness = 0.10000; // cm
                        if (fName == "Fi30")
                            fiber_thickness = 0.1034;
                        if (fName == "Fi31")
                            fiber_thickness = 0.1033;
                        if (fName == "Fi32")
                            fiber_thickness = 0.1024;
                        if (fName == "Fi33")
                            fiber_thickness = 0.1025;
                        Int_t fiber_nbr = 512;
                        Float_t dead_layer = 0.9;
                        Float_t air_layer = 0.01 * 0.; // relative to fiber_thickness
                        Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);

                        if (fDirection == VERTICAL)
                        {
                            //  x = -detector_width / 2. +
                            //    (double(fiber_id - 1) + (double(fiber_id - 1) * air_layer)) * fiber_thickness;

                            x = -detector_width / 2. + fiber_thickness * (1 + air_layer) / 2. +
                                double(fiber_id - 1.) * (1 + air_layer) * fiber_thickness;

                            y = (t_down - t_up) * veff;
                            // y = 0.;
                        }
                        else
                        {
                            x = (t_down - t_up) * veff;
                            // x = 0.;
                            //  y = -detector_width / 2. +
                            //    (double(fiber_id - 1) + (double(fiber_id - 1) * air_layer)) * fiber_thickness;

                            y = -detector_width / 2. + fiber_thickness * (1 + air_layer) / 2. +
                                double(fiber_id - 1.) * (1 + air_layer) * fiber_thickness;
                        }
                    }
                    // cout<<"Fiber y " << y << endl;
                    if (y < -60 || y > 60)
                    {
                        // continue;
                    }

                    Double_t eloss = sqrt(tot_down * tot_up);
                    //  eloss = 2.069*eloss-10.414;  // testing Q
                    // Z calib, run 773
                    // if (fName == "Fi30") eloss =
                    // -415.06629+174.67595*eloss-24.234506*eloss*eloss+1.1193663*eloss*eloss*eloss; if (fName ==
                    // "Fi31") eloss = -197.91146+90.109359*eloss-13.516177*eloss*eloss+0.68043217*eloss*eloss*eloss; if
                    // (fName == "Fi32") eloss =
                    // -156.00382+70.021718*eloss-10.350304*eloss*eloss+0.51455325*eloss*eloss*eloss; if (fName ==
                    // "Fi33") eloss = -298.60978+130.07812*eloss-18.699849*eloss*eloss+0.89852388*eloss*eloss*eloss;

                    Double_t t = tof;

                    if (fSimu)
                    { // for s494 exp
                        Double_t elossin = sqrt(eloss);
                        if (fName == "Fi30")
                        {
                            eloss = (6. - 2.) / (3.84711 - 1.21996) * elossin + 6. -
                                    (6. - 2.) / (3.84711 - 1.21996) * 3.84711;
                        }
                        if (fName == "Fi31")
                        {
                            eloss = (6. - 2.) / (3.84843 - 1.2226) * elossin + 6. -
                                    (6. - 2.) / (3.84843 - 1.2226) * 3.84843;
                        }
                        if (fName == "Fi32")
                        {
                            eloss = (6. - 2.) / (3.83080 - 1.21518) * elossin + 6. -
                                    (6. - 2.) / (3.83080 - 1.221518) * 3.83080;
                        }
                        if (fName == "Fi33")
                        {
                            eloss = (6. - 2.) / (3.82969 - 1.21722) * elossin + 6. -
                                    (6. - 2.) / (3.82969 - 1.21722) * 3.82969;
                        }
                        eloss = (int)(eloss + 0.5);
                    }

                    multi++;

                    if (!fIsCalibrator)
                    {
                        if (tof >= ftofmin && tof <= ftofmax)
                        {
                            new ((*fHitItems)[fNofHitItems++])
                                R3BFiberMAPMTHitData(0, x, y, eloss, t, fiber_id, t_down, t_up, tot_down, tot_up);
                        }
                    }
                }
            }
        }
    }

    fnEvents++;
}

void R3BFiberMAPMTCal2Hit::FinishEvent()
{
    fCalItems->Clear();
    fCalTriggerItems->Clear();
    fHitItems->Clear();
    fNofHitItems = 0;
}

void R3BFiberMAPMTCal2Hit::FinishTask()
{
    if (fWrite)
    {
        fh_ToT_bottom_Fib->Write();
        fh_ToT_top_Fib->Write();
        fh_ToT_bottom_Fib_raw->Write();
        fh_ToT_top_Fib_raw->Write();
        fh_dt_Fib->Write();
        fh_time_top_Fib->Write();
        fh_Fib_ToF->Write();
        fh_time_bottom_Fib->Write();
        fh_dt_Fib_raw->Write();
        fh_Fib_ToF_raw->Write();
        fh_time_check_tsync->Write();
    }

    if (fIsCalibrator)
    {

        R3BFiberMAPMTHitModulePar* mpar;

        UInt_t max = N_FIBER_MAX;
        if (fh_ToT_bottom_Fib->GetNbinsX() < N_FIBER_MAX)
            max = fh_ToT_bottom_Fib->GetNbinsX();

        for (UInt_t i = 1; i <= max; i++)
        {
            mpar = new R3BFiberMAPMTHitModulePar();
            mpar->SetFiber(i);
            fCalPar->AddModulePar(mpar);
        }

        for (UInt_t i = 1; i <= max; i++)
        {
            // gain bottom MAPMTs
            TH1D* proj = fh_ToT_bottom_Fib_raw->ProjectionY("", i + 1, i + 1, 0);
            for (UInt_t j = proj->GetNbinsX() - 2; j > 2; j--)
            {
                if (j == 2)
                {
                    // could not find maximum
                }
                if (proj->GetBinContent(j) > proj->GetMaximum() * 10. / 100.)
                {
                    R3BFiberMAPMTHitModulePar* par = fCalPar->GetModuleParAt(i);
                    par->SetGainDown(proj->GetBinCenter(j));

                    cout << fName << " fiberId: " << i << ",gainDown: " << proj->GetBinCenter(j) << endl;

                    break;
                }
            }
            // gain top MAPMT
            TH1D* proj1 = fh_ToT_top_Fib_raw->ProjectionY("", i + 1, i + 1, 0);
            for (UInt_t j = proj1->GetNbinsX() - 2; j > 2; j--)
            {
                if (j == 2)
                {
                    // could not find maximum
                }
                if (proj1->GetBinContent(j) > proj1->GetMaximum() * 10. / 100.)
                {
                    R3BFiberMAPMTHitModulePar* par1 = fCalPar->GetModuleParAt(i);
                    par1->SetGainUp(proj->GetBinCenter(j));

                    cout << fName << " fiberId: " << i << ",gainUp: " << proj1->GetBinCenter(j) << endl;

                    break;
                }
            }
            // time offset
            R3BFiberMAPMTHitModulePar* par2 = fCalPar->GetModuleParAt(i);
            TH1D* proj2 = fh_dt_Fib_raw->ProjectionY("", i + 1, i + 1, 0);
            par2->SetOffsetUp(0.5 * proj2->GetBinCenter(proj2->GetMaximumBin()));
            par2->SetOffsetDown(-0.5 * proj2->GetBinCenter(proj2->GetMaximumBin()));

            // tsync
            R3BFiberMAPMTHitModulePar* par3 = fCalPar->GetModuleParAt(i);
            TH1D* proj3 = fh_Fib_ToF_raw->ProjectionY("", i + 1, i + 1, 0);
            par3->SetSync(proj3->GetBinCenter(proj3->GetMaximumBin()));

            cout << fName << " fiberId: " << i << ", offset: " << 0.5 * proj2->GetBinCenter(proj2->GetMaximumBin())
                 << ", sync: " << proj3->GetBinCenter(proj3->GetMaximumBin()) << endl;
        }

        fCalPar->setChanged();
    }
}

ClassImp(R3BFiberMAPMTCal2Hit)
