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

#include "R3BFiberMAPMTCal2Hit.h"
#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BFiberMAPMTCalData.h"
#include "R3BFiberMAPMTHitData.h"
#include "R3BFiberMAPMTHitPar.h"
#include "R3BFiberMappingPar.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TClonesArray.h>

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
    , fDetId(1)
    , fClockFreq(150.)
    , fDirection(a_direction)
    , fNumFibers(a_num_fibers)
    , fIsCalibrator(a_is_calibrator)
    , fCalItems()
    , fCalTriggerItems()
    , fHitItems(new TClonesArray("R3BFiberMAPMTHitData"))
    , fMapPar(NULL)
    , fCalPar(NULL)
    , fHitPar(NULL)
    , fNofHitPars(0)
    , fChannelArray()
    , fnEvents(0)
    , ftofmin(-1000)
    , ftofmax(1000)
    , fWrite(false)
    , fGate_ns(100.)
    , fOnline(kFALSE)
    , fOrientation(STANDARD)
    , fTimeStitch(nullptr)
    , fHeader(nullptr)
{
    if (fName == "Fi23a")
        fDetId = 230;
    else if (fName == "Fi23b")
        fDetId = 231;
    else if (fName == "Fi30")
        fDetId = 30;
    else if (fName == "Fi31")
        fDetId = 31;
    else if (fName == "Fi32")
        fDetId = 32;
    else if (fName == "Fi33")
        fDetId = 33;
    else
    {
        R3BLOG(error, "Fiber " << fName << " not found");
    }
}

R3BFiberMAPMTCal2Hit::~R3BFiberMAPMTCal2Hit()
{
    if (fHitItems)
        delete fHitItems;
    if (fMapPar)
        delete fMapPar;
    if (fCalPar)
        delete fCalPar;
    if (fHitPar)
        delete fHitPar;
}

InitStatus R3BFiberMAPMTCal2Hit::Init()
{
    R3BLOG(info, "For firber " << fName);
    auto mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, !mgr, "FairRootManager not found.");

    fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(fatal, NULL == fHeader, "EventHeader. not found");

    auto name = fName + "Cal";
    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name));
    R3BLOG_IF(fatal, !fCalItems, "Branch " << name << " not found.");

    auto name_mapmt_trig = fName + "TriggerCal";
    fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name_mapmt_trig));
    R3BLOG_IF(fatal, NULL == fCalTriggerItems, "Branch " << name_mapmt_trig << " not found");

    // maxevent = mgr->CheckMaxEventNo();

    // TClones branch with Hit items
    mgr->Register(fName + "Hit", fName + " hit data", fHitItems, !fOnline);
    fHitItems->Clear();

    // Resize per-channel info arrays.
    for (auto side_i = 0; side_i < 2; ++side_i)
    {
        fChannelArray[side_i].resize(fNumFibers);
    }

    // create histograms
    TString chistName;
    TString chistTitle;
    // ToT bottom PMT raw:
    chistName = fName + "_ToT_bottom_raw";
    chistTitle = fName + " ToTbottom raw of fibers";
    fh_ToT_bottom_Fib_raw = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 200, 0., 100.);
    fh_ToT_bottom_Fib_raw->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_bottom_Fib_raw->GetYaxis()->SetTitle("ToT / ns");

    // ToT bottom gainmatched:
    chistName = fName + "_ToT_bottom";
    chistTitle = fName + " ToTbottom of fibers";
    fh_ToT_bottom_Fib = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 200, 0., 100.);
    fh_ToT_bottom_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_bottom_Fib->GetYaxis()->SetTitle("ToT / ns");

    // ToF Tofd -> Fiber:
    chistName = fName + "_time";
    chistTitle = fName + " Time";
    fh_Fib_ToF = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 100, -12000., 12000.);
    fh_Fib_ToF->GetYaxis()->SetTitle("Time/ ns");
    fh_Fib_ToF->GetXaxis()->SetTitle("Fiber ID");

    // ToF Fiber for sync:
    chistName = fName + "_time_raw";
    chistTitle = fName + " Time raw";
    fh_Fib_ToF_raw = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 100, -12000., 12000.);
    fh_Fib_ToF_raw->GetYaxis()->SetTitle("Time/ ns");
    fh_Fib_ToF_raw->GetXaxis()->SetTitle("Fiber ID");

    // ToT top MAPMT raw:
    chistName = fName + "_ToT_top_raw";
    chistTitle = fName + " ToTtop raw of fibers";
    fh_ToT_top_Fib_raw = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 200, 0., 100.);
    fh_ToT_top_Fib_raw->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_top_Fib_raw->GetYaxis()->SetTitle("ToT / ns");

    // ToT top MAPMT gainmatched:
    chistName = fName + "_ToT_top";
    chistTitle = fName + " ToTtop of fibers";
    fh_ToT_top_Fib = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 200, 0., 100.);
    fh_ToT_top_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_top_Fib->GetYaxis()->SetTitle("ToT / ns");

    chistName = fName + "_time_top";
    chistTitle = fName + " time of top MAPMTs ";
    fh_time_top_Fib =
        new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 2400, -15000., 12000.);
    fh_time_top_Fib->GetXaxis()->SetTitle("single number");
    fh_time_top_Fib->GetYaxis()->SetTitle("time / ns");

    // time difference top-bottom:
    chistName = fName + "_dt";
    chistTitle = fName + " dt of fibers";
    fh_dt_Fib = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 500, -100., 100.);
    fh_dt_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_dt_Fib->GetYaxis()->SetTitle("dt / ns");

    // time difference top-bottom for offsets:
    chistName = fName + "_dt_raw";
    chistTitle = fName + " dt_raw of fibers";
    fh_dt_Fib_raw = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 500, -100., 100.);
    fh_dt_Fib_raw->GetXaxis()->SetTitle("Fiber number");
    fh_dt_Fib_raw->GetYaxis()->SetTitle("dt / ns");

    // time of MAPMT:
    chistName = fName + "_time_bottom";
    chistTitle = fName + " time of bottom MAPMTs";
    fh_time_bottom_Fib =
        new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 5400, -5000., 11000.);
    fh_time_bottom_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_time_bottom_Fib->GetYaxis()->SetTitle("time / ns");

    // histos to check hit parameters
    chistName = fName + "_check_tsync";
    chistTitle = fName + " check tsync for fiberNum 10";
    fh_time_check_tsync = new TH2F(chistName.Data(), chistTitle.Data(), 50000, 0, 1e6, 400, -100., 100.);
    fh_time_check_tsync->GetXaxis()->SetTitle("Event number");
    fh_time_check_tsync->GetYaxis()->SetTitle("dt / ns");

    // Definition of a time stich object to correlate times coming from different systems
    fTimeStitch = new R3BCoarseTimeStitch();

    return kSUCCESS;
}

InitStatus R3BFiberMAPMTCal2Hit::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BFiberMAPMTCal2Hit::SetParContainers()
{
    fMapPar = dynamic_cast<R3BFiberMappingPar*>(FairRuntimeDb::instance()->getContainer(fName + "MappingPar"));
    if (!fMapPar)
    {
        R3BLOG(error, "Couldn't get " << fName << "MappingPar");
    }
    else
    {
        fNumFibers = fMapPar->GetNbChannels();
        R3BLOG(info, "Nb of fibers: " << fNumFibers);
    }
    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCalPar = dynamic_cast<R3BFiberMAPMTHitPar*>(FairRuntimeDb::instance()->getContainer(fName + "HitPar"));
    R3BLOG_IF(info, fCalPar, "Container " << fName << "HitPar initialized");
    R3BLOG_IF(error, !fCalPar, "Couldn't get " << fName << "HitPar");

    if (!fIsCalibrator)
    {
        // Get calibration parameters if we're not a calibrator.
        auto container = fName + "HitPar";
        fHitPar = dynamic_cast<R3BFiberMAPMTHitPar*>(FairRuntimeDb::instance()->getContainer(container));
        R3BLOG_IF(error, !fHitPar, "Could not get " << container << " container.");
    }
}

void R3BFiberMAPMTCal2Hit::Exec(Option_t* option)
{

    for (auto side_i = 0; side_i < 2; ++side_i)
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
        auto cur_cal = dynamic_cast<R3BFiberMAPMTCalData*>(fCalTriggerItems->At(j));
        auto ch = cur_cal->GetChannel() - 1;
        tl = cur_cal->GetTime_ns();
        trig_time[ch] = tl;
    }

    size_t cal_num = fCalItems->GetEntriesFast();

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal_lead = dynamic_cast<R3BFiberMAPMTCalData const*>(fCalItems->At(j));

        if (cur_cal_lead->IsLeading())
        {
            auto side_i = cur_cal_lead->GetSide() - 1;
            auto ch_i = cur_cal_lead->GetChannel() - 1;
            auto& channel = fChannelArray[side_i].at(ch_i);
            channel.lead_list.push_back(cur_cal_lead);
        }
    }

    for (size_t j = 0; j < cal_num; ++j)
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

            Double_t cur_cal_trig_ns = 0;
            Double_t lead_trig_ns = 0;
            Double_t lead_raw = 0;
            Double_t trail_raw = 0;
            Double_t trail_trig_ns = 0;
            Double_t lead_ns = 0. / 0., tot_ns = 0. / 0.;

            auto chlead_i = lead->GetChannel() - 1;

            cur_cal_trig_ns = trig_time[fMapPar->GetTrigMap(side_i + 1, ch_i + 1)];
            lead_trig_ns = trig_time[fMapPar->GetTrigMap(side_i + 1, chlead_i + 1)];

            auto cur_cal_ns =
                fTimeStitch->GetTime(cur_cal_trail->GetTime_ns() - cur_cal_trig_ns, "clocktdc", "clocktdc");
            lead_ns = fTimeStitch->GetTime(lead->GetTime_ns() - lead_trig_ns, "clocktdc", "clocktdc");
            tot_ns = fTimeStitch->GetTime(cur_cal_ns - lead_ns, "clocktdc", "clocktdc");

            lead_raw = lead->GetTime_ns();
            trail_raw = cur_cal_trail->GetTime_ns();
            trail_trig_ns = cur_cal_trig_ns;

            if (tot_ns < fGate_ns && tot_ns > 0.)
            {
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
                    // Double_t dtime = fTimeStitch->GetTime(t_up - t_down, "clocktdc", "clocktdc");
                    Double_t dtime = t_up - t_down;
                    Double_t tof =
                        fHeader ? fTimeStitch->GetTime((t_up + t_down) / 2. - fHeader->GetTStart(), "vftx", "clocktdc")
                                : (t_up + t_down) / 2.;

                    // Fill histograms for gain match, offset and sync.
                    if (fWrite)
                    {
                        fh_ToT_bottom_Fib_raw->Fill(fiber_id, tot_down);
                        fh_ToT_top_Fib_raw->Fill(fiber_id, tot_up);
                        fh_dt_Fib_raw->Fill(fiber_id, dtime);
                        fh_Fib_ToF_raw->Fill(fiber_id, tof);
                    }

                    gainUp = 10.;
                    gainDown = 10.;
                    tsync = 0.;
                    // offsetUp = 0.;
                    offsetDT = 0.;

                    if (!fIsCalibrator && fHitPar)
                    {
                        R3BFiberMAPMTHitModulePar* par = fHitPar->GetModuleParAt(fiber_id);
                        if (par)
                        {
                            gainUp = par->GetGainUp();
                            gainDown = par->GetGainDown();
                            // tsync = par->GetSync();
                            // offsetUp = par->GetOffsetUp();
                            offsetDT = par->GetOffsetDown();
                        }
                    }

                    tot_down *= 20. / gainDown;
                    tot_up *= 20. / gainUp;

                    //  tof -= tsync;
                    // t_down -= offsetDown;
                    // t_down -= tsync;
                    // t_up -= offsetUp;
                    // t_up -= tsync;
                    dtime -= offsetDT;
                    tof -= tsync;

                    // histogram after gain match, sync....
                    if (fWrite)
                    {
                        fh_dt_Fib->Fill(fiber_id, dtime);
                        fh_Fib_ToF->Fill(fiber_id, tof);
                        fh_ToT_bottom_Fib->Fill(fiber_id, tot_down);
                        fh_ToT_top_Fib->Fill(fiber_id, tot_up);
                        fh_time_bottom_Fib->Fill(fiber_id, t_down);
                        fh_time_top_Fib->Fill(fiber_id, t_up);
                    }

                    if (fiber_id == 10)
                        fh_time_check_tsync->Fill(fnEvents, tof);

                    Double_t x = -10000.;
                    Double_t y = -10000.;
                    // FIXME: veff should be taken from the hit-parameter container
                    Double_t veff = 12. / 2.; // cm/ns
                    Double_t randx = (std::rand() / (float)RAND_MAX);
                    if (fName == "Fi23a" || fName == "Fi23b")
                    {
                        Float_t fiber_thickness = 0.028;
                        Float_t air_layer = 0.01 * 0.; // relative to fiber_thickness
                        Float_t detector_width = fNumFibers * fiber_thickness * (1 + air_layer);
                        if (fDirection == VERTICAL)
                        {
                            x = -detector_width / 2. +
                                (double(fiber_id - 1) + (double(fiber_id - 1.) * air_layer)) * fiber_thickness;
                            y = (t_down - t_up) * veff;
                        }
                        else
                        {
                            x = (t_down - t_up) * veff;
                            y = -detector_width / 2. +
                                (double(fiber_id - 1) + (double(fiber_id - 1) * air_layer)) * fiber_thickness;
                        }
                    }
                    if (fName == "Fi30" || fName == "Fi31" || fName == "Fi32" || fName == "Fi33")
                    {
                        Float_t fiber_thickness = 0.10000; // cm
                        if (fName == "Fi30")
                            fiber_thickness = 0.1034;
                        if (fName == "Fi31")
                            fiber_thickness = 0.1033;
                        if (fName == "Fi32")
                            fiber_thickness = 0.1024;
                        if (fName == "Fi33")
                            fiber_thickness = 0.1025;

                        Float_t air_layer = 0.0; // 1 * 0.; // relative to fiber_thickness
                        Float_t detector_width = fNumFibers * fiber_thickness * (1.0 + air_layer);

                        // FIXME: This needs improvements because it depends on the fiber direction and orientation
                        if (fDirection == VERTICAL)
                        {
                            x = (fOrientation == STANDARD ? 1.0 : -1.0) *
                                (-1.0 * detector_width / 2.0 + (fiber_id - 1.0) * (1.0 + air_layer) * fiber_thickness);
                            y = (t_down - t_up) * veff;
                        }
                        else
                        {
                            x = (t_down - t_up) * veff;
                            y = (fOrientation == STANDARD ? 1.0 : -1.0) *
                                (-1.0 * detector_width / 2.0 + (fiber_id - 1.0) * (1.0 + air_layer) * fiber_thickness);
                        }
                    }

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

                    multi++;

                    if (!fIsCalibrator)
                    {
                        // if (tof >= ftofmin && tof <= ftofmax)
                        if (fabs(dtime) < fDTime_window)
                        {
                            new ((*fHitItems)[fHitItems->GetEntriesFast()]) R3BFiberMAPMTHitData(
                                fDetId, x, y, eloss, tof, fiber_id, t_down, t_up, tot_down, tot_up);
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
    if (fHitItems)
    {
        fHitItems->Clear();
    }
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
        Bool_t Redo = false; // Add a redo flag, to redo a calibration, if Maximum falls on noise.
        Double_t PercentOfMax = 0.5;

        R3BFiberMAPMTHitModulePar* mpar;

        for (UInt_t i = 1; i <= fNumFibers; i++)
        {
            mpar = new R3BFiberMAPMTHitModulePar();
            mpar->SetFiber(i);
            fCalPar->AddModulePar(mpar);
        }

        Int_t loopctr = 0;
        for (UInt_t i = 1; i <= fNumFibers; i++)
        {
            // gain bottom MAPMTs
            TH1D* proj = fh_ToT_bottom_Fib_raw->ProjectionY("", i, i, 0);
            Double_t PartMax = proj->GetMaximum() * PercentOfMax;
            Int_t LowerBound = proj->FindFirstBinAbove(PartMax, 1);
            loopctr = 0;

            for (UInt_t j = proj->GetNbinsX() - 2; j > 2; j--)
            {
                if (j == 2)
                {
                    // could not find maximum
                }
                if (proj->GetBinContent(j) > PartMax)
                {
                    if (j - LowerBound <= 90. /*Bins = 0.1 ns*/ && loopctr < 4)
                    {
                        // Redo, start from beginning, with lower threshold for "peak" acceptance.
                        j = proj->GetNbinsX() - 2;
                        PartMax = PartMax * 0.6; // lower threshold, do as many times until noise peak is no longer
                                                 // considered "Data-peak"
                        Redo = true;
                        loopctr++;
                    }
                    else
                    {
                        R3BFiberMAPMTHitModulePar* par = fCalPar->GetModuleParAt(i);
                        par->SetGainDown(proj->GetBinCenter(j));

                        R3BLOG(info, fName << " fiberId: " << i << ",gainDown: " << proj->GetBinCenter(j));
                        if (Redo == true)
                        {
                            PartMax = proj->GetMaximum() * PercentOfMax;
                            Redo = false;
                        }
                        break;
                    }
                }
            }
            // gain top MAPMT
            TH1D* proj1 = fh_ToT_top_Fib_raw->ProjectionY("", i, i, 0);
            Double_t PartMax1 = proj1->GetMaximum() * PercentOfMax;
            Int_t LowerBound1 = proj1->FindFirstBinAbove(PartMax1, 1);
            Redo = false;
            loopctr = 0;
            for (UInt_t j = proj1->GetNbinsX() - 2; j > 2; j--)
            {
                if (j == 2)
                {
                    // could not find maximum
                }
                if (proj1->GetBinContent(j) > PartMax1)
                {
                    if (j - LowerBound1 <= 90. /*ns*/ && loopctr < 4)
                    {
                        // Redo, start from beginning, with lower threshold for "peak" acceptance.
                        j = proj1->GetNbinsX() - 2;
                        PartMax1 = PartMax1 * 0.6; // lower threshold, do as many times until noise peak is no longer
                                                   // considered "Data-peak"
                        Redo = true;
                        loopctr++;
                    }
                    else
                    {
                        R3BFiberMAPMTHitModulePar* par1 = fCalPar->GetModuleParAt(i);
                        par1->SetGainUp(proj1->GetBinCenter(j));
                        R3BLOG(info, fName << " fiberId: " << i << ",gainUp: " << proj1->GetBinCenter(j));
                        break;
                    }
                }
            }

            // time offset
            R3BFiberMAPMTHitModulePar* par2 = fCalPar->GetModuleParAt(i);
            auto proj2 = fh_dt_Fib_raw->ProjectionY("", i, i, 0);
            // par2->SetOffsetUp(0.5 * proj2->GetBinCenter(proj2->GetMaximumBin()));
            auto Max = proj2->GetXaxis()->GetBinCenter(proj2->GetMaximumBin());
            TF1* fGauss = new TF1("fgaus", "gaus", Max - 8., Max + 8.);
            proj2->Fit("fgaus", "QR0");
            auto offsetdt = fGauss->GetParameter(1);
            par2->SetOffsetDown(offsetdt);
            // fGauss->Draw("LSAME");
            // proj2->Write();

            // tsync
            // R3BFiberMAPMTHitModulePar* par3 = fCalPar->GetModuleParAt(i);
            auto proj3 = fh_Fib_ToF_raw->ProjectionY("", i, i, 0);
            auto sync = proj3->GetXaxis()->GetBinCenter(proj3->GetMaximumBin());
            par2->SetSync(sync);

            R3BLOG(info, fName << " fiberId: " << i << ", offset_DT: " << offsetdt << ", sync: " << sync);
        }

        fCalPar->setChanged();
    }
}

ClassImp(R3BFiberMAPMTCal2Hit)
