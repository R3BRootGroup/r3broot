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

#include "R3BBunchedFiberCal2Hit.h"
#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberHitPar.h"
#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BFiberMappingPar.h"
#include "R3BLogger.h"
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

R3BBunchedFiberCal2Hit::ToT::ToT(R3BBunchedFiberCalData const* a_lead,
                                 R3BBunchedFiberCalData const* a_trail,
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

R3BBunchedFiberCal2Hit::R3BBunchedFiberCal2Hit(const char* a_name,
                                               Int_t a_verbose,
                                               enum R3BTCalEngine::CTDCVariant a_ctdc_variant,
                                               Direction a_direction,
                                               UInt_t a_sub_num,
                                               UInt_t a_mapmt_per_sub,
                                               UInt_t a_spmt_per_sub,
                                               Bool_t a_is_calibrator,
                                               Bool_t a_is_gain,
                                               Bool_t a_is_tsync)
    : FairTask(TString("R3B") + a_name + "Cal2Hit", a_verbose)
    , fName(a_name)
    , fClockFreq(R3BTCalEngine::CTDC_16_BWD_150 == a_ctdc_variant ? 150 : 250)
    , fDirection(a_direction)
    , fSubNum(a_sub_num)
    , fIsCalibrator(a_is_calibrator)
    , fIsGain(a_is_gain)
    , fIsTsync(a_is_tsync)
    , fCalItems()
    , fMAPMTCalTriggerItems()
    , fSPMTCalTriggerItems()
    , fHitItems(new TClonesArray("R3BBunchedFiberHitData"))
    , fMapPar(NULL)
    , fCalPar(NULL)
    , fHitPar(NULL)
    , fNofHitPars()
    , fNofHitItems()
    , fChannelArray()
    , fh_ToT_MA_Fib()
    , fh_ToT_Single_Fib()
    , fh_ToT_s_Fib()
    , fh_time_s_Fib()
    , fh_ToT_ToT()
    , fh_dt_Fib()
    , fh_time_Fib()
    , fh_Fib_ToF()
    , fnEvents(0)
    , fNumFibers(a_sub_num * a_mapmt_per_sub)
    , fOrientation(STANDARD)
    , fTimeStitch(nullptr)
    , fExpId(0)
    , fDetId(0)

{
    fChPerSub[0] = a_mapmt_per_sub;
    fChPerSub[1] = a_spmt_per_sub;

    if (fName == "Fi7")
    {
        fDetId = 7;
    }
    else if (fName == "Fi8")
    {
        fDetId = 8;
    }
    else if (fName == "Fi9")
    {
        fDetId = 9;
    }
    else if (fName == "Fi10")
    {
        fDetId = 10;
    }
    else if (fName == "Fi11")
    {
        fDetId = 11;
    }
    else if (fName == "Fi12")
    {
        fDetId = 12;
    }
    else if (fName == "Fi13")
    {
        fDetId = 13;
    }
}

R3BBunchedFiberCal2Hit::~R3BBunchedFiberCal2Hit()
{
    if (fHitItems)
        delete fHitItems;
    if (fCalPar)
        delete fCalPar;
    if (fMapPar)
        delete fMapPar;
}

InitStatus R3BBunchedFiberCal2Hit::Init()
{
    R3BLOG(info, "");
    auto mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, !mgr, "FairRootManager not found.");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (!header)
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    auto name = fName + "Cal";
    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name));
    R3BLOG_IF(fatal, !fCalItems, "Branch " << name << " not found.");

    auto name_mapmt_trig = fName + "TriggerCal";
    fMAPMTCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name_mapmt_trig));
    R3BLOG_IF(fatal, !fMAPMTCalTriggerItems, "Branch " << name_mapmt_trig << " not found.");

    auto name_spmt_trig = "BunchedFiberSPMTTrigCal";
    fSPMTCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name_spmt_trig));
    R3BLOG_IF(warning, !fSPMTCalTriggerItems, "Branch " << name_spmt_trig << " not found.");

    // maxevent = mgr->CheckMaxEventNo();
    // if (!fIsCalibrator)
    mgr->Register(fName + "Hit", "Fiber Hit Data", fHitItems, kTRUE);

    // Resize per-channel info arrays.
    for (auto side_i = 0; side_i < 2; ++side_i)
    {
        fChannelArray[side_i].resize(fSubNum * fChPerSub[side_i]);
    }

    for (int i = 0; i < fNumFibers; i++)
    {
        gain_temp[i] = 10.;
        tsync_temp[i] = 0.;
    }

    if (!fIsCalibrator)
    {
        // Get calibration parameters if we're not a calibrator.
        auto container = fName + "HitPar";
        fHitPar = dynamic_cast<R3BBunchedFiberHitPar*>(FairRuntimeDb::instance()->getContainer(container));

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

        if (fHitPar)
        {
            for (int i = 1; i <= fNumFibers; i++)
            {
                if (!fIsGain && fIsTsync) // gain already made, tsync will be done now
                {
                    R3BBunchedFiberHitModulePar* par = fHitPar->GetModuleParAt(i);
                    if (par && par->GetGainMA() > 0.)
                        gain_temp[i - 1] = par->GetGainMA();
                }
                if (fIsGain && !fIsTsync) // tsync already made, gain will be done now
                {
                    R3BBunchedFiberHitModulePar* par = fHitPar->GetModuleParAt(i);
                    if (par)
                        tsync_temp[i - 1] = par->GetSync();
                }
            }
        }
    }

    // create histograms
    TString chistName;
    TString chistTitle;
    // ToT MAPMT:
    chistName = fName + "_ToT_MAPMT";
    chistTitle = fName + " ToT of fibers";
    fh_ToT_MA_Fib = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 200, 0., 100.);
    fh_ToT_MA_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_MA_Fib->GetYaxis()->SetTitle("ToT / ns");

    // ToF Tofd -> Fiber:
    chistName = fName + "_tof";
    chistTitle = fName + " ToF";
    fh_Fib_ToF = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 44000, -1100., 1100.);
    fh_Fib_ToF->GetYaxis()->SetTitle("ToF / ns");
    fh_Fib_ToF->GetXaxis()->SetTitle("Fiber ID");

    // Test:
    chistName = fName + "_test";
    chistTitle = fName + " Tsync test";
    fh_Test = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 50, -1024., 1024.);
    fh_Test->GetYaxis()->SetTitle("Tsync / ns");
    fh_Test->GetXaxis()->SetTitle("Fiber ID");

    // ToT single PMT:
    chistName = fName + "_ToT_SAPMT";
    chistTitle = fName + " ToT of fibers";
    fh_ToT_Single_Fib = new TH2F(chistName.Data(),
                                 chistTitle.Data(),
                                 fNumFibers,
                                 1,
                                 fNumFibers + 1,
                                 500,
                                 0.,
                                 100.); // binning changed from 100 to 500
    fh_ToT_Single_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_Single_Fib->GetYaxis()->SetTitle("ToT / ns");

    // ToT SAPMT:
    for (Int_t i = 0; i < 4; i++)
    {
        char number[15];
        snprintf(number, sizeof(number), "%d", i);
        chistName = fName + "_ToT_SAPMT" + number;
        chistTitle = fName + " ToT of single PMTs " + number;
        fh_ToT_s_Fib[i] = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 500, 0., 100.);
        fh_ToT_s_Fib[i]->GetXaxis()->SetTitle("Fiber number");
        fh_ToT_s_Fib[i]->GetYaxis()->SetTitle("ToT / ns");
    }

    chistName = fName + "_time_SAPMT";
    chistTitle = fName + " time of single PMTs ";
    fh_time_s_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 6, 0, 6, 4000, -2000., 2000.);
    fh_time_s_Fib->GetXaxis()->SetTitle("single number");
    fh_time_s_Fib->GetYaxis()->SetTitle("time / ns");

    // ToT vs ToT SPMT:
    chistName = fName + "_ToT_ToT";
    chistTitle = fName + " ToT vs ToT of single PMTs";
    fh_ToT_ToT = new TH2F(chistName.Data(), chistTitle.Data(), 1000, 0., 100., 1000, 0., 100.);
    fh_ToT_ToT->GetXaxis()->SetTitle("SPMT 1");
    fh_ToT_ToT->GetYaxis()->SetTitle("SPMT 2");

    // time difference SPMT - MAPMT:
    chistName = fName + "_dt";
    chistTitle = fName + " dt of fibers";
    fh_dt_Fib = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 1000, -5000., 5000.);
    fh_dt_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_dt_Fib->GetYaxis()->SetTitle("dt / ns");

    // time of MAPMT:
    chistName = fName + "_time";
    chistTitle = fName + " time of fibers";
    fh_time_Fib = new TH2F(chistName.Data(), chistTitle.Data(), fNumFibers, 1, fNumFibers + 1, 4000, -2000., 2000.);
    fh_time_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_time_Fib->GetYaxis()->SetTitle("time / ns");

    // Definition of a time stich object to correlate times coming from different systems
    fTimeStitch = new R3BCoarseTimeStitch();

    return kSUCCESS;
}

InitStatus R3BBunchedFiberCal2Hit::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BBunchedFiberCal2Hit::SetParContainers()
{
    fMapPar = dynamic_cast<R3BFiberMappingPar*>(FairRuntimeDb::instance()->getContainer(fName + "MappingPar"));
    if (!fMapPar)
    {
        R3BLOG(error, "Couldn't get " << fName << "MappingPar");
    }
    else
    {
        fNumFibers = fMapPar->GetNbChannels();
        if (fName == "Fi10" || fName == "Fi11" || fName == "Fi12" || fName == "Fi13")
        {
            fNumFibers *= fSubNum;
        }
        R3BLOG(info, "Nb of fibers for " << fName << ": " << fNumFibers);
    }

    fCalPar = dynamic_cast<R3BBunchedFiberHitPar*>(FairRuntimeDb::instance()->getContainer(fName + "HitPar"));
    R3BLOG_IF(error, !fCalPar, "Couldn't get " << fName << "HitPar.");
}

void R3BBunchedFiberCal2Hit::Exec(Option_t* option)
{
    if (fExpId == 515 || (fExpId == 0 && header->GetExpId() == 515)) // If fExpId is not set, global ExpId will be used
    {
        S515();
    }
    else
    {
        Standard();
    }
}

void R3BBunchedFiberCal2Hit::S515()
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
    auto cal_num = fCalItems->GetEntriesFast();

    auto mapmt_trig_num = fMAPMTCalTriggerItems->GetEntriesFast();
    std::vector<R3BBunchedFiberCalData const*> mapmt_trig_table(fSubNum * fChPerSub[0]);
    // std::cout << mapmt_trig_table.size() << std::endl;
    for (size_t j = 0; j < mapmt_trig_num; ++j)
    {
        auto cal = dynamic_cast<R3BBunchedFiberCalData const*>(fMAPMTCalTriggerItems->At(j));
        auto idx = cal->GetChannel() - 1;
        if (idx >= mapmt_trig_table.size())
            mapmt_trig_table.resize(idx + 1);
        mapmt_trig_table.at(idx) = cal;
    }

    // Find multi-hit ToT for every channel.
    // The easiest safe way to survive ugly cases is to record all
    // leading edges per channel, and then pair up with whatever
    // trailing we have.
    // Not super efficient, but shouldn't crash if the data is not
    // perfect.
    unsigned n_lead = 0;
    unsigned n_trail = 0;
    int s_mult = 0;
    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = dynamic_cast<R3BBunchedFiberCalData const*>(fCalItems->At(j));
        // auto cur_cal = (R3BBunchedFiberCalData const*)fCalItems->At(j);
        if (cur_cal->IsLeading())
        {
            ++n_lead;
            // auto side_i = 0;
            auto side_i = cur_cal->IsMAPMT() ? 0 : 1;
            auto ch_i = cur_cal->GetChannel() - 1;
            auto& channel = fChannelArray[side_i].at(ch_i);
            channel.lead_list.push_back(cur_cal);
        }
        else
        {
            ++n_trail;
        }
    }

    if (n_lead == n_trail && n_lead == 0)
    {
        return;
    }

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = dynamic_cast<R3BBunchedFiberCalData const*>(fCalItems->At(j));
        // auto cur_cal = (R3BBunchedFiberCalData const*)fCalItems->At(j);
        if (!cur_cal->IsTrailing() || !cur_cal->IsMAPMT())
            continue;
        // auto side_i = 0;
        auto side_i = cur_cal->IsMAPMT() ? 0 : 1;
        auto ch_i = cur_cal->GetChannel() - 1;
        auto& channel = fChannelArray[side_i].at(ch_i);
        if (channel.lead_list.empty())
        {
            continue;
        }
        auto lead = channel.lead_list.front();

        Double_t cur_cal_trig_ns = 0;
        Double_t lead_trig_ns = 0;
        if (cur_cal->IsMAPMT() && fMapPar)
        {
            auto cur_cal_trig_i = fMapPar->GetTrigMap(side_i + 1, ch_i + 1); // fMAPMTTriggerMap[ch_i];
            auto lead_trig_i = fMapPar->GetTrigMap(side_i + 1, lead->GetChannel());
            // if(fName=="Fi3a") printf("3a trig curr %8u %8u lead %8u %8u  %8u\n", ch_i, cur_cal_trig_i,
            // lead->GetChannel() - 1, lead_trig_i, mapmt_trig_table.size()); if(fName=="Fi3b") printf("3b trig curr
            // %8u %8u lead %8u %8u\n", ch_i, cur_cal_trig_i, lead->GetChannel() - 1, lead_trig_i);
            if (cur_cal_trig_i < mapmt_trig_table.size() && lead_trig_i < mapmt_trig_table.size())
            {
                auto cur_cal_trig = mapmt_trig_table.at(cur_cal_trig_i);
                auto lead_trig = mapmt_trig_table.at(lead_trig_i);
                cur_cal_trig_ns = cur_cal_trig->GetTime_ns();
                lead_trig_ns = lead_trig->GetTime_ns();
            }
        }

        Double_t cur_cal_ns = 0.;
        Double_t lead_ns = 0.;
        Double_t tot_ns = 0.;
        if (side_i == 0)
        {
            cur_cal_ns = fTimeStitch->GetTime(cur_cal->GetTime_ns() - cur_cal_trig_ns, "clocktdc", "clocktdc");
            lead_ns = fTimeStitch->GetTime(lead->GetTime_ns() - lead_trig_ns, "clocktdc", "clocktdc");
            tot_ns = fTimeStitch->GetTime(cur_cal_ns - lead_ns, "clocktdc", "clocktdc");
        }
        else
        {
            cur_cal_ns = fTimeStitch->GetTime(cur_cal->GetTime_ns() - cur_cal_trig_ns, "tamex", "tamex");
            lead_ns = fTimeStitch->GetTime(lead->GetTime_ns() - lead_trig_ns, "tamex", "tamex");
            tot_ns = fTimeStitch->GetTime(cur_cal_ns - lead_ns, "tamex", "tamex");
        }

        auto fiber_id = (lead->GetChannel() - 1) * fChPerSub[1] + 1;

        // Apply calibration.
        Double_t gainMA = gain_temp[fiber_id - 1]; // 10.;
        Double_t gainS = gain_temp[fiber_id - 1];
        Double_t offset1 = 0.;
        Double_t offset2 = 0.;
        Double_t tsync = tsync_temp[fiber_id - 1]; // 0.;

        double tof = 0.;

        if (!fIsCalibrator && fHitPar)
        {
            R3BBunchedFiberHitModulePar* par = fHitPar->GetModuleParAt(fiber_id);
            if (par)
            {
                gainMA = par->GetGainMA();
                tsync = par->GetSync();
                gainS = par->GetGainS();
                offset1 = par->GetOffset1();
                offset2 = par->GetOffset2();
            }
        }

        fh_time_Fib->Fill(fiber_id, lead_ns);

        tof = lead_ns - tsync;
        auto caltot_ns = tot_ns * 10. / gainMA;

        // if (fIsCalibrator)
        {
            // histogram for offset determination

            // Fill histograms for gain match, and for debugging.
            fh_Fib_ToF->Fill(fiber_id, tof);
            if (caltot_ns > 0.)
                fh_ToT_MA_Fib->Fill(fiber_id, caltot_ns);
        }

        Double_t x = -10000.;
        Double_t y = -10000.;
        Double_t veff = 12. / 2.; // cm/ns

        if (fName == "Fi10" || fName == "Fi11" || fName == "Fi12" || fName == "Fi13")
        {
            Float_t fiber_thickness = 0.050000;
            Int_t fiber_nbr = 1024;
            Float_t dead_layer = 0.9;
            Float_t air_layer = 0.0; // relative to fiber_thickness
            Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);

            if (fDirection == VERTICAL)
            {
                x = (fOrientation == STANDARD ? 1.0 : -1.0) *
                    (-1.0 * detector_width / 2.0 + (fiber_id - 1.0) * (1.0 + air_layer) * fiber_thickness);
                y = 0;
            }
            else
            {
                x = 0.;
                y = (fOrientation == STANDARD ? 1.0 : -1.0) *
                    (-1.0 * detector_width / 2.0 + (fiber_id - 1.0) * (1.0 + air_layer) * fiber_thickness);
            }
        }

        new ((*fHitItems)[fHitItems->GetEntriesFast()])
            R3BBunchedFiberHitData(fDetId, x, y, caltot_ns, tof, fiber_id, lead_ns, 0., caltot_ns, 0.);
    }
}

void R3BBunchedFiberCal2Hit::Standard()
{
    multi = 0;
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
    for (Int_t i = 0; i < fNumFibers; i++)
    {
        energy[i] = 0.;
        counts[i] = 0;
    }

    size_t cal_num = fCalItems->GetEntriesFast();

    // Make direct mapping tables for trigger items.
    size_t mapmt_trig_num = fMAPMTCalTriggerItems->GetEntriesFast();
    std::vector<R3BBunchedFiberCalData const*> mapmt_trig_table(fSubNum * fChPerSub[0]);
    for (size_t j = 0; j < mapmt_trig_num; ++j)
    {
        auto cal = dynamic_cast<R3BBunchedFiberCalData const*>(fMAPMTCalTriggerItems->At(j));
        auto idx = cal->GetChannel() - 1;
        if (idx >= mapmt_trig_table.size())
            mapmt_trig_table.resize(idx + 1);
        mapmt_trig_table.at(idx) = cal;
    }

    // TODO: This will create a map for every fiber detector... Urg.
    // Also since it's shared between many detectors it must be dynamic, for now.
    size_t spmt_trig_num = fSPMTCalTriggerItems->GetEntriesFast();
    std::vector<R3BBunchedFiberCalData const*> spmt_trig_table;
    if (fSPMTCalTriggerItems)
    {
        for (size_t j = 0; j < spmt_trig_num; ++j)
        {
            auto cal = dynamic_cast<R3BBunchedFiberCalData const*>(fSPMTCalTriggerItems->At(j));
            auto idx = cal->GetChannel() - 1;
            if (idx >= spmt_trig_table.size())
                spmt_trig_table.resize(idx + 1);
            spmt_trig_table.at(idx) = cal;
        }
    }

    // Find multi-hit ToT for every channel.
    // The easiest safe way to survive ugly cases is to record all
    // leading edges per channel, and then pair up with whatever
    // trailing we have.
    // Not super efficient, but shouldn't crash if the data is not
    // perfect.
    unsigned n_lead = 0;
    unsigned n_trail = 0;
    int s_mult = 0;
    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = dynamic_cast<R3BBunchedFiberCalData const*>(fCalItems->At(j));
        // auto cur_cal = (R3BBunchedFiberCalData const*)fCalItems->At(j);
        if (cur_cal->IsLeading())
        {
            ++n_lead;
            auto side_i = cur_cal->IsMAPMT() ? 0 : 1;
            auto ch_i = cur_cal->GetChannel() - 1;
            auto& channel = fChannelArray[side_i].at(ch_i);
            channel.lead_list.push_back(cur_cal);
            if (side_i == 1)
                s_mult++;
        }
        else
        {
            ++n_trail;
        }
    }
    if (n_lead != n_trail)
    {
        //    return;
    }
    Int_t summmpt = 0;
    Int_t summsm1 = 0;
    Int_t summsm2 = 0;

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = dynamic_cast<R3BBunchedFiberCalData const*>(fCalItems->At(j));
        // auto cur_cal = (R3BBunchedFiberCalData const*)fCalItems->At(j);
        if (cur_cal->IsTrailing())
        {
            auto side_i = cur_cal->IsMAPMT() ? 0 : 1;
            auto ch_i = cur_cal->GetChannel() - 1;
            auto& channel = fChannelArray[side_i].at(ch_i);
            if (channel.lead_list.empty())
            {
                continue;
            }
            auto lead = channel.lead_list.front();

            Double_t cur_cal_trig_ns = 0;
            Double_t lead_trig_ns = 0;
            if (cur_cal->IsMAPMT() && fMapPar)
            {
                auto cur_cal_trig_i = fMapPar->GetTrigMap(side_i + 1, ch_i + 1); // fMAPMTTriggerMap[ch_i];
                auto lead_trig_i = fMapPar->GetTrigMap(side_i + 1, lead->GetChannel());
                // if(fName=="Fi3a") printf("3a trig curr %8u %8u lead %8u %8u  %8u\n", ch_i, cur_cal_trig_i,
                // lead->GetChannel() - 1, lead_trig_i, mapmt_trig_table.size()); if(fName=="Fi3b") printf("3b trig curr
                // %8u %8u lead %8u %8u\n", ch_i, cur_cal_trig_i, lead->GetChannel() - 1, lead_trig_i);
                if (cur_cal_trig_i < mapmt_trig_table.size() && lead_trig_i < mapmt_trig_table.size())
                {
                    auto cur_cal_trig = mapmt_trig_table.at(cur_cal_trig_i);
                    auto lead_trig = mapmt_trig_table.at(lead_trig_i);
                    cur_cal_trig_ns = cur_cal_trig->GetTime_ns();
                    lead_trig_ns = lead_trig->GetTime_ns();
                }
            }
            else if (cur_cal->IsSPMT() && fMapPar)
            {
                auto cur_cal_trig_i = fMapPar->GetTrigMap(side_i + 1, ch_i + 1); // fSPMTTriggerMap[ch_i];
                auto lead_trig_i = fMapPar->GetTrigMap(side_i + 1, lead->GetChannel());
                // if(fName=="Fi3a") printf("3a trig curr %8u %8u lead %8u %8u  %8u\n", ch_i, cur_cal_trig_i,
                // lead->GetChannel() - 1, lead_trig_i, spmt_trig_table.size()); if(fName=="Fi3b") printf("3b trig curr
                // %8u %8u lead %8u %8u\n", ch_i, cur_cal_trig_i, lead->GetChannel() - 1, lead_trig_i);
                if (cur_cal_trig_i < spmt_trig_table.size() && lead_trig_i < spmt_trig_table.size())
                {
                    auto cur_cal_trig = spmt_trig_table.at(cur_cal_trig_i);
                    auto lead_trig = spmt_trig_table.at(lead_trig_i);
                    cur_cal_trig_ns = cur_cal_trig->GetTime_ns();
                    lead_trig_ns = lead_trig->GetTime_ns();
                }
            }

            Double_t cur_cal_ns = 0.;
            Double_t lead_ns = 0.;
            Double_t tot_ns = 0.;
            if (side_i == 0)
            {
                cur_cal_ns = fTimeStitch->GetTime(cur_cal->GetTime_ns() - cur_cal_trig_ns, "clocktdc", "clocktdc");
                lead_ns = fTimeStitch->GetTime(lead->GetTime_ns() - lead_trig_ns, "clocktdc", "clocktdc");
                tot_ns = fTimeStitch->GetTime(cur_cal_ns - lead_ns, "clocktdc", "clocktdc");
            }
            else
            {
                cur_cal_ns = fTimeStitch->GetTime(cur_cal->GetTime_ns() - cur_cal_trig_ns, "tamex", "tamex");
                lead_ns = fTimeStitch->GetTime(lead->GetTime_ns() - lead_trig_ns, "tamex", "tamex");
                tot_ns = fTimeStitch->GetTime(cur_cal_ns - lead_ns, "tamex", "tamex");
            }

            if (side_i == 1)
                channel.tot_list.push_back(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));
            if (side_i == 0)
                channel.tot_list.push_front(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));

            channel.lead_list.pop_front();
        }
    }

    double s1 = 99.;
    double s2 = 99.;
    double s3 = 99.;
    double s4 = 99.;
    int single = 0;
    double tof = 0.;

    // Make every permutation to create fibers.
    auto const& mapmt_array = fChannelArray[0];
    auto const& spmt_array = fChannelArray[1];

    for (auto it_mapmt = mapmt_array.begin(); mapmt_array.end() != it_mapmt; ++it_mapmt) // over fiber_number 0...nmax-1
    {
        auto const& mapmt = *it_mapmt;
        for (auto it_mapmt_tot = mapmt.tot_list.begin(); mapmt.tot_list.end() != it_mapmt_tot;
             ++it_mapmt_tot) // over ihit(fiber)
        {
            auto const& mapmt_tot = *it_mapmt_tot;
            for (auto it_spmt = spmt_array.begin(); spmt_array.end() != it_spmt; ++it_spmt)
            {
                auto const& spmt = *it_spmt;
                for (auto it_spmt_tot = spmt.tot_list.begin(); spmt.tot_list.end() != it_spmt_tot; ++it_spmt_tot)
                {
                    auto const& spmt_tot = *it_spmt_tot;

                    // Check that the combo is inside one sub-det.
                    auto mapmt_sub_id = (mapmt_tot.lead->GetChannel() - 1) / fChPerSub[0];
                    auto spmt_sub_id = (spmt_tot.lead->GetChannel() - 1) / fChPerSub[1];
                    if (mapmt_sub_id != spmt_sub_id)
                    {
                        continue;
                    }

                    /*
                     * How to get a fiber ID for a fiber detector defined as:
                     *  SubNum = 2
                     *  MAPMT = 256
                     *  SPMT = 2
                     * This means we'll receive 512 MAPMT channels as 1..512, and 4 SPMT
                     * channels, but the first half (sub-detector) is completely
                     * decoupled from the second half. The sequence of all fibers in
                     * order is then, as (MAPMT,SPMT)-pairs:
                     *  (1,1), (1,2), (2,1), ... (256,1), (256,2),
                     *  (257,3), (257,4), (258,3), ... (512,3), (512,4)
                     */

                    // auto fiber_id = mapmt_tot.lead->GetChannel();
                    auto fiber_id = (mapmt_tot.lead->GetChannel() - 1) * fChPerSub[1] +
                                    ((spmt_tot.lead->GetChannel() - 1) % fChPerSub[1]) + 1;

                    single = spmt_tot.lead->GetChannel();

                    // Calibrate hit fiber.
                    auto tot_mapmt = mapmt_tot.tot_ns;
                    auto tot_spmt = spmt_tot.tot_ns;
                    Double_t t_mapmt = mapmt_tot.lead_ns;

                    Double_t t_mapmt1 = mapmt_tot.lead->GetTime_ns(); // MAPMT time without subtraction of trigger time

                    Double_t t_spmt = spmt_tot.lead_ns;
                    // only accept hits which are at the right time:
                    Bool_t simu = true;
                    if (!simu)
                    {
                        if (fName == "Fi3a" || fName == "Fi3b")
                        {
                            // s                            if (t_spmt < -400 || t_spmt > -270)
                            // s                                continue;
                            if (t_mapmt < -550 || t_mapmt > -350)
                                continue;
                        }
                        if (fName == "Fi10" || fName == "Fi11" || fName == "Fi12" || fName == "Fi13")
                        {
                            // s                            if (t_spmt < -300 || t_spmt > -160)
                            // s                                continue;
                            if (t_mapmt < -500 || t_mapmt > -300)
                                continue;
                        }
                    }

                    // if (fIsCalibrator)
                    {
                        fh_time_Fib->Fill(fiber_id, t_mapmt);
                        fh_time_s_Fib->Fill(single, t_spmt);
                    }

                    // Apply calibration.
                    Double_t gainMA = gain_temp[fiber_id - 1]; // 10.;
                    Double_t gainS = gain_temp[fiber_id - 1];
                    Double_t offset1 = 0.;
                    Double_t offset2 = 0.;
                    Double_t tsync = tsync_temp[fiber_id - 1]; // 0.;

                    if (!fIsCalibrator && fHitPar)
                    {
                        R3BBunchedFiberHitModulePar* par = fHitPar->GetModuleParAt(fiber_id);
                        if (par)
                        {
                            gainMA = par->GetGainMA();
                            tsync = par->GetSync();
                            gainS = par->GetGainS();
                            offset1 = par->GetOffset1();
                            offset2 = par->GetOffset2();
                        }
                    }

                    tot_mapmt *= 10. / gainMA;
                    tot_spmt *= 10. / gainS;
                    t_mapmt -= offset1;
                    t_spmt -= offset2;

                    if (tot_mapmt > 0. && tot_spmt > 0.)
                    {
                        tof = (t_mapmt + t_spmt) / 2.;
                    }
                    else
                    {
                        tof = t_mapmt;
                    }
                    tof -= tsync;

                    // if (fIsCalibrator)
                    {
                        // histogram for offset determination
                        fh_dt_Fib->Fill(fiber_id, t_spmt - t_mapmt);

                        // Fill histograms for gain match, and for debugging.
                        fh_Fib_ToF->Fill(fiber_id, tof);
                        fh_ToT_MA_Fib->Fill(fiber_id, tot_mapmt);
                        if (s_mult > 0)
                        {
                            fh_ToT_Single_Fib->Fill(fiber_id, tot_spmt);
                            fh_ToT_s_Fib[single - 1]->Fill(fiber_id, tot_spmt);
                        }
                    }

                    Double_t x = -10000.;
                    Double_t y = -10000.;
                    Double_t veff = 12. / 2.; // cm/ns

                    if (fName == "Fi10" || fName == "Fi11" || fName == "Fi12" || fName == "Fi13")
                    {
                        Float_t fiber_thickness = 0.050000;
                        Int_t fiber_nbr = 1024;
                        Float_t dead_layer = 0.9;
                        Float_t air_layer = 0.0; // relative to fiber_thickness
                        Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);

                        if (fDirection == VERTICAL)
                        {
                            x = (fOrientation == STANDARD ? 1.0 : -1.0) *
                                (-1.0 * detector_width / 2.0 + (fiber_id - 1.0) * (1.0 + air_layer) * fiber_thickness);
                            y = (t_spmt - t_mapmt) * veff;
                        }
                        else
                        {
                            x = (t_spmt - t_mapmt) * veff;
                            y = (fOrientation == STANDARD ? 1.0 : -1.0) *
                                (-1.0 * detector_width / 2.0 + (fiber_id - 1.0) * (1.0 + air_layer) * fiber_thickness);
                        }
                    }

                    if (fName == "Fi1a" || fName == "Fi1b" || fName == "Fi2a" || fName == "Fi2b" || fName == "Fi3a" ||
                        fName == "Fi3b")
                    {
                        Float_t fiber_thickness = 0.021000 * 2.; // s remove *2 when taking SPMT into analysis
                        Int_t fiber_nbr = 512 / 2;               // s remove /2 when taking SPMT into analysis
                        Float_t dead_layer = 0.9;
                        Float_t air_layer = 0.01; // relative to fiber_thickness
                        Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);
                        if (fDirection == VERTICAL)
                        {
                            x = -detector_width / 2. + fiber_thickness / 2. +
                                ((fiber_id - 1) + ((fiber_id - 1) * air_layer)) * fiber_thickness;
                            // s                            y = (t_spmt - t_mapmt) * 3.;
                            y = 0.;
                        }
                        else
                        {
                            // s                            x = (t_spmt - t_mapmt) * 3.;
                            x = 0.;
                            y = -detector_width / 2. + fiber_thickness / 2. +
                                ((fiber_id - 1) + ((fiber_id - 1) * air_layer)) * fiber_thickness;
                        }
                    }
                    // cout<<"Fiber y " << y << endl;
                    if (y < -60 || y > 60)
                    {
                        // continue;
                    }
                    if (tof < -20 || tof > 20)
                    {
                        // continue;
                    }

                    Double_t eloss = 0.;
                    if (tot_mapmt > 0. && tot_spmt > 0.)
                    {
                        eloss = sqrt(tot_mapmt * tot_spmt);
                    }
                    else
                    {
                        eloss = tot_mapmt;
                    }
                    // if (fIsCalibrator)
                    {
                        fh_ToT_ToT->Fill(tot_mapmt, tot_spmt);
                    }

                    energy[fiber_id - 1] = eloss;
                    counts[fiber_id - 1] = counts[fiber_id - 1] + 1;
                    multi++;

                    // if (!fIsCalibrator)
                    new ((*fHitItems)[fHitItems->GetEntriesFast()])
                        R3BBunchedFiberHitData(0, x, y, eloss, tof, fiber_id, t_mapmt, t_spmt, tot_mapmt, tot_spmt);
                }
            }
        }
    }

    fnEvents++;
    return;
}

void R3BBunchedFiberCal2Hit::FinishEvent()
{
    if (fHitItems)
        fHitItems->Clear();
}

void R3BBunchedFiberCal2Hit::FinishTask()
{

    fh_ToT_MA_Fib->Write();
    fh_ToT_Single_Fib->Write();
    fh_dt_Fib->Write();
    fh_time_Fib->Write();
    fh_Fib_ToF->Write();
    fh_Test->Write();
    fh_time_s_Fib->Write();

    /*

        fh_multi->Write();
    */

    for (Int_t i = 0; i < 4; i++)
    {
        fh_ToT_s_Fib[i]->Write();
    }
    fh_ToT_ToT->Write();

    if (fIsCalibrator)
    {

        R3BBunchedFiberHitModulePar* mpar;

        UInt_t max = fNumFibers;
        if (fh_ToT_MA_Fib->GetNbinsX() < fNumFibers)
            max = fh_ToT_MA_Fib->GetNbinsX();

        for (UInt_t i = 1; i <= max; i++)
        {
            mpar = new R3BBunchedFiberHitModulePar();
            mpar->SetFiber(i);
            fCalPar->AddModulePar(mpar);
        }

        // tsync
        if (fIsTsync)
        {
            for (UInt_t i = 1; i <= max; i++)
            {
                // R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);

                if (!fIsGain)
                {
                    R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);
                    par->SetGainMA(gain_temp[i - 1]);
                }
                R3BBunchedFiberHitModulePar* par1 = fCalPar->GetModuleParAt(i);
                TH1D* proj = fh_Fib_ToF->ProjectionY("", i, i, 0);
                par1->SetSync(proj->GetBinCenter(proj->GetMaximumBin()));

                cout << "MA fiber: " << fName << ", " << i << " tsync: " << proj->GetBinCenter(proj->GetMaximumBin())
                     << ", " << par1->GetSync() << ", " << par1->GetGainMA() << endl;
            }
        }
        // time offset

        for (UInt_t i = 1; i <= max; i++)
        {
            TH1D* proj = fh_dt_Fib->ProjectionY("", i, i, 0);
            R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);
            par->SetOffset1(0.5 * proj->GetBinCenter(proj->GetMaximumBin()));
            par->SetOffset2(-0.5 * proj->GetBinCenter(proj->GetMaximumBin()));

            cout << "MA fiber: " << i << " offset: " << 0.5 * proj->GetBinCenter(proj->GetMaximumBin()) << endl;
        }

        if (fIsGain)
        {
            // MAPMT
            for (UInt_t i = 1; i <= max; i++)
            {
                TH1D* proj = fh_ToT_MA_Fib->ProjectionY("", i, i, 0);

                Int_t binmax = proj->GetMaximumBin();
                auto totsync = proj->GetXaxis()->GetBinCenter(binmax);

                R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);
                par->SetGainMA(totsync);
                cout << "MAPMT fiber: " << i << " par: " << totsync << endl;

                /* for (UInt_t j = proj->GetNbinsX() ; j > 0; j--)
                 {
                     if (j == 2)
                     {
                         // could not find maximum
                     }
                     if (proj->GetBinContent(j) > proj->GetMaximum() / 100.)
                     {
                         R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);
                         par->SetGainMA(proj->GetBinCenter(j));
                         cout << "MAPMT fiber: " << i << " par: " << proj->GetBinCenter(j) << endl;
                         // par->SetGainMA(j - 1);
                         //	cout << "MA fiber: " << i << " par: " << proj->GetBinCenter(j) << endl;
                         break;
                     }
                 }*/
            }

            // SPMT
            /* for (UInt_t i = 1; i <= max; i++)
             {
                 TH1D* proj = fh_ToT_Single_Fib->ProjectionY("", i, i, 0);
                 for (UInt_t j = proj->GetNbinsX(); j > 0; j--)
                 {
                     if (j == 2)
                     {
                         // could not find maximum
                     }

                     if (proj->GetBinContent(j) > proj->GetMaximum() / 10.)
                     {
                         R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);
                         par->SetGainS(proj->GetBinCenter(j));
                         // cout<<"S fiber: "<< i<<" par: "<<proj->GetBinCenter(j)<<endl;
                         // par->SetGainS(j - 1);
                         cout << "SPMT fiber: " << i << " par: " << proj->GetBinCenter(j) << endl;
                         break;
                     }
                 }
             }*/

            // Spmt-test
            for (UInt_t i = 1; i <= max; i++)
            {
                TH1D* proj = fh_ToT_Single_Fib->ProjectionY("", i, i, 0);
                Int_t binmax = proj->GetMaximumBin();
                auto s_tot = proj->GetXaxis()->GetBinCenter(binmax);

                R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);
                par->SetGainS(s_tot);

                cout << "SPMT fiber: " << i << " par: " << s_tot << endl;
            }
        }

        fCalPar->setChanged();
    }
}

ClassImp(R3BBunchedFiberCal2Hit);
