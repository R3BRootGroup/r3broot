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

#include "R3BBunchedFiberCal2Hit_s494.h"
#include "R3BBunchedFiberCalData.h"
#include "R3BBunchedFiberHitData.h"
#include "R3BBunchedFiberHitPar.h"
#include "R3BTCalEngine.h"
#include "TH1F.h"
#include "TH2F.h"
#include <TClonesArray.h>
#include <cassert>

#include "R3BTofdHitData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <iostream>

R3BBunchedFiberCal2Hit_s494::ToT::ToT(R3BBunchedFiberCalData const* a_lead,
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

R3BBunchedFiberCal2Hit_s494::R3BBunchedFiberCal2Hit_s494(const char* a_name,
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
    , fHitItems(new TClonesArray("R3BBunchedFiberHitData"))
    , fCalPar()
    , fHitPar()
    , fNofHitPars()
    , fNofHitItems()
    , fChannelArray()
    , fh_ToT_MA_Fib()
    , fh_ToT_SA_Fib()
    , fh_time_SA_Fib()
    , fh_dt_Fib()
    , fh_time_MA_Fib()
    , fh_Fib_ToF()
    , fnEvents(0)

{
    fChPerSub[0] = a_mapmt_per_sub;
    fChPerSub[1] = a_spmt_per_sub;
}

R3BBunchedFiberCal2Hit_s494::~R3BBunchedFiberCal2Hit_s494()
{
    if (fHitItems)
        delete fHitItems;
    if (fCalPar)
        delete fCalPar;
}

InitStatus R3BBunchedFiberCal2Hit_s494::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    cout << "Init:R3BBunchedFiberCal2Hit_s494" << endl;
    auto mgr = FairRootManager::Instance();
    if (!mgr)
    {
        LOG(error) << "FairRootManager not found.";
        return kERROR;
    }

    auto name = fName + "Cal";
    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name));
    if (!fCalItems)
    {
        LOG(fatal) << "Branch " << name << " not found.";
        return kERROR;
    }
    auto name_mapmt_trig = fName + "TriggerCal";
    fMAPMTCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name_mapmt_trig));
    if (NULL == fMAPMTCalTriggerItems)
        LOG(fatal) << "Branch " << name_mapmt_trig << " not found";

    maxevent = mgr->CheckMaxEventNo();

    mgr->Register(fName + "Hit", fName + " hit data", fHitItems, kTRUE);
    // Resize per-channel info arrays.
    for (auto side_i = 0; side_i < 2; ++side_i)
    {
        fChannelArray[side_i].resize(fSubNum * fChPerSub[side_i]);
    }

    ///////

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
    }

    //    }
    // create histograms
    TString chistName;
    TString chistTitle;
    // ToT MAPMT:
    chistName = fName + "_ToT_MA";
    chistTitle = fName + " ToT of fibers";
    //    fh_ToT_MA_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 2100, 0., 2100., 100, 0., 41.666667);
    fh_ToT_MA_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 250, 0., 100.);
    fh_ToT_MA_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_MA_Fib->GetYaxis()->SetTitle("ToT / ns");

    // ToF Tofd -> Fiber:
    chistName = fName + "_time";
    chistTitle = fName + " Time";
    fh_Fib_ToF = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 44000, -1100., 1100.);
    fh_Fib_ToF->GetYaxis()->SetTitle("Time/ ns");
    fh_Fib_ToF->GetXaxis()->SetTitle("Fiber ID");

    // Test:
    chistName = fName + "_test";
    chistTitle = fName + " Tsync test";
    fh_Test = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 50, -1024., 1024.);
    fh_Test->GetYaxis()->SetTitle("Tsync / ns");
    fh_Test->GetXaxis()->SetTitle("Fiber ID");

    // ToT single PMT:
    chistName = fName + "_ToT_SA";
    chistTitle = fName + " ToT of fibers";
    fh_ToT_SA_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 250, 0., 100.);
    fh_ToT_SA_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_ToT_SA_Fib->GetYaxis()->SetTitle("ToT / ns");

    chistName = fName + "_time_SA";
    chistTitle = fName + " time of single PMTs ";
    fh_time_SA_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 4000, -2000., 2000.);
    fh_time_SA_Fib->GetXaxis()->SetTitle("single number");
    fh_time_SA_Fib->GetYaxis()->SetTitle("time / ns");

    // time difference SPMT - MAPMT:
    chistName = fName + "_dt";
    chistTitle = fName + " dt of fibers";
    fh_dt_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 1000, -5000., 5000.);
    fh_dt_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_dt_Fib->GetYaxis()->SetTitle("dt / ns");

    // time of MAPMT:
    chistName = fName + "_time_MA";
    chistTitle = fName + " time of fibers";
    fh_time_MA_Fib = new TH2F(chistName.Data(), chistTitle.Data(), 512, 0, 512, 4000, -2000., 2000.);
    fh_time_MA_Fib->GetXaxis()->SetTitle("Fiber number");
    fh_time_MA_Fib->GetYaxis()->SetTitle("time / ns");

    cout << "R3BBunchedFiberCal2Hit: Spectra created!" << endl;

    return kSUCCESS;
}

InitStatus R3BBunchedFiberCal2Hit_s494::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BBunchedFiberCal2Hit_s494::SetParContainers()
{
    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCalPar = dynamic_cast<R3BBunchedFiberHitPar*>(FairRuntimeDb::instance()->getContainer(fName + "HitPar"));
    if (!fCalPar)
    {
        LOG(error) << "R3BFiberCal2Hit_s494::Init() Couldn't get " << fName << "HitPar. ";
    }
}

void R3BBunchedFiberCal2Hit_s494::Exec(Option_t* option)
{
    // cout<<"R3BBunchedFiberCal2Hit: in exec"<<endl;
    //	if(fnEvents/10000.==(int)fnEvents/10000) cout<<"Events: "<<fnEvents<<"         \r"<<std::flush;
    multi = 0;
    if (fnEvents / 100000. == (int)fnEvents / 100000)
        std::cout << "\rEvents: " << fnEvents << " / " << maxevent << " (" << (int)(fnEvents * 100. / maxevent)
                  << " %) " << std::flush;

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
    for (Int_t i = 0; i < 1025; i++)
    {
        energy[i] = 0.;
        counts[i] = 0;
    }
    //   cout<<"multihit"<<endl;
    // fSide: 0-bottom, 1-top, 2-trigger
    size_t cal_num = fCalItems->GetEntriesFast();

    // Make direct mapping tables for trigger items.
    size_t mapmt_trig_num = fMAPMTCalTriggerItems->GetEntries();
    std::vector<R3BBunchedFiberCalData const*> mapmt_trig_table(fSubNum * fChPerSub[0] / 128);
    for (size_t j = 0; j < mapmt_trig_num; ++j)
    {
        auto cal = dynamic_cast<R3BBunchedFiberCalData const*>(fMAPMTCalTriggerItems->At(j));
        mapmt_trig_table.at(cal->GetChannel() - 1) = cal;
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
        if (cur_cal->IsLeading())
        {
            ++n_lead;
            auto side_i = cur_cal->GetSide(); // cur_cal->IsMAPMT() ? 0 : 1;
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

    // puts("Event");
    // bool do_print = false;

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = dynamic_cast<R3BBunchedFiberCalData const*>(fCalItems->At(j));
        if (cur_cal->IsTrailing())
        {
            auto side_i = cur_cal->GetSide(); // cur_cal->IsMAPMT() ? 0 : 1;
            // Double_t c_period = 0 == side_i ? 4096. * (1000. / fClockFreq) : 2048. * (1000. / 200.);
            Double_t c_period = 4096. * (1000. / fClockFreq);

            auto ch_i = cur_cal->GetChannel() - 1;
            auto& channel = fChannelArray[side_i].at(ch_i);
            if (channel.lead_list.empty())
            {
                continue;
            }
            auto lead = channel.lead_list.front();

            Double_t cur_cal_trig_ns = 0;
            Double_t lead_trig_ns = 0;
            /*
                        if (cur_cal->IsMAPMT() && fMAPMTTriggerMap)
                        {

                            auto cur_cal_trig_i = fMAPMTTriggerMap[ch_i];
                            auto lead_trig_i = fMAPMTTriggerMap[lead->GetChannel() - 1];
                            // if(fName=="Fi30") printf("30 trig curr %8u %8u lead %8u %8u  %8u\n", ch_i,
               cur_cal_trig_i,
                            //  lead->GetChannel() - 1, lead_trig_i, mapmt_trig_table.size());
                            // if(fName=="Fi30") printf("30 trig curr %8u %8u lead %8u %8u\n", ch_i, cur_cal_trig_i,
                            // lead->GetChannel() - 1, lead_trig_i);
                            if (cur_cal_trig_i < mapmt_trig_table.size() && lead_trig_i < mapmt_trig_table.size())
                            {

                                auto cur_cal_trig = mapmt_trig_table.at(cur_cal_trig_i);
                                auto lead_trig = mapmt_trig_table.at(lead_trig_i);
                                cur_cal_trig_ns = cur_cal_trig->GetTime_ns();
                                lead_trig_ns = lead_trig->GetTime_ns();
                            }

                        }
                        else if (cur_cal->IsSPMT() && fMAPMTTriggerMap)
                        {
                            auto cur_cal_trig_i = fMAPMTTriggerMap[ch_i];
                            auto lead_trig_i = fMAPMTTriggerMap[lead->GetChannel() - 1];
                            // if(fName=="Fi3a") printf("3a trig curr %8u %8u lead %8u %8u  %8u\n", ch_i,
               cur_cal_trig_i,
                            // lead->GetChannel() - 1, lead_trig_i, spmt_trig_table.size()); if(fName=="Fi3b")
               printf("3b trig curr
                            // %8u %8u lead %8u %8u\n", ch_i, cur_cal_trig_i, lead->GetChannel() - 1, lead_trig_i);
                            if (cur_cal_trig_i < mapmt_trig_table.size() && lead_trig_i < mapmt_trig_table.size())
                            {
                                auto cur_cal_trig = mapmt_trig_table.at(cur_cal_trig_i);
                                auto lead_trig = mapmt_trig_table.at(lead_trig_i);
                                cur_cal_trig_ns = cur_cal_trig->GetTime_ns();
                                lead_trig_ns = lead_trig->GetTime_ns();
                            }
                        }


                        auto cur_cal_ns =
                            fmod(cur_cal->GetTime_ns() - cur_cal_trig_ns + c_period + c_period / 2, c_period) - c_period
               / 2; auto lead_ns = fmod(lead->GetTime_ns() - lead_trig_ns + c_period + c_period / 2, c_period) -
               c_period / 2; auto tot_ns = fmod(cur_cal_ns - lead_ns + c_period + c_period / 2, c_period) - c_period /
               2;
            */
            auto cur_cal_ns = fmod(cur_cal->GetTime_ns() + c_period + c_period / 2, c_period) - c_period / 2;
            auto lead_ns = fmod(lead->GetTime_ns() + c_period + c_period / 2, c_period) - c_period / 2;
            auto tot_ns = fmod(cur_cal_ns - lead_ns + c_period + c_period / 2, c_period) - c_period / 2;

            //   cout<<fName<<","<<lead_ns<<", "<<tot_ns<<endl;

            if (tot_ns > 0 && tot_ns < 1000)
            {
                // if (side_i==0) tot_ns -= 9.;//taken out for the moment MH
                /*
                  if (side_i == 1)
                      channel.tot_list.push_back(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));
                  if (side_i == 0)
                      channel.tot_list.push_front(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));
                */
                // In s494 all FibDet are readd out with MAPMT->CTDC and no TAMEX, thus nothing has to be inversed
                channel.tot_list.push_back(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));
                channel.lead_list.pop_front();
            }
        }
    }
    // if (do_print) for (size_t j = 0; j < cal_num; ++j) {
    // auto cur_cal = (R3BBunchedFiberCalData const*)fCalItems->At(j);
    // if (cur_cal->IsMAPMT())
    // printf("%3u: %3u %u %f\n", j, cur_cal->GetChannel(), cur_cal->IsLeading(), cur_cal->GetTime_ns());
    //}

    double tot_mapmt_max = -1;
    double tot_spmt_max = -1;
    int tot_mapmt_max_fiber_id = 0;
    int tot_spmt_max_fiber_id = 0;
    int tofd_id = 0;
    double tof;
    double t_tofd;
    // cout<<"permutations"<<endl;

    //   cout << "new Event ********************* " << fName << endl;
    // Make every permutation to create fibers.
    auto const& mapmt_array = fChannelArray[0];
    auto const& spmt_array = fChannelArray[1];

    for (auto it_mapmt = mapmt_array.begin(); mapmt_array.end() != it_mapmt;
         ++it_mapmt) // over MA channel 0...nmax-1; nmax=fSubNum*fChPerSub[0]
    {
        auto const& mapmt = *it_mapmt;
        for (auto it_mapmt_tot = mapmt.tot_list.begin(); mapmt.tot_list.end() != it_mapmt_tot;
             ++it_mapmt_tot) // over ihit(channel)
        {
            auto const& mapmt_tot = *it_mapmt_tot;
            auto mapmt_sub_id = (mapmt_tot.lead->GetChannel() - 1) / fChPerSub[0];
            auto fiber_MA_ch = mapmt_tot.lead->GetChannel();

            for (auto it_spmt = spmt_array.begin(); spmt_array.end() != it_spmt;
                 ++it_spmt) // over SA channel 0...nmax-1; nmax=fSubNum*fChPerSub[1]
            {
                auto const& spmt = *it_spmt;
                for (auto it_spmt_tot = spmt.tot_list.begin(); spmt.tot_list.end() != it_spmt_tot;
                     ++it_spmt_tot) // over ihit(channel)
                {
                    auto const& spmt_tot = *it_spmt_tot;

                    // Check that the combo is inside one sub-det.

                    auto spmt_sub_id = (spmt_tot.lead->GetChannel() - 1) / fChPerSub[1];
                    if (mapmt_sub_id != spmt_sub_id)
                    {
                        continue;
                    }

                    auto fiber_SA_ch = spmt_tot.lead->GetChannel();
                    Int_t fiber_id = -1000;
                    if (fiber_MA_ch == fiber_SA_ch)
                        fiber_id = fiber_MA_ch;
                    else
                        continue;

                    /*
                        auto fiber_id = (mapmt_tot.lead->GetChannel() - 1) * fChPerSub[1] +
                                        ((spmt_tot.lead->GetChannel() - 1) % fChPerSub[1]) + 1;
                     */

                    // cout << "mapmt: " << mapmt_tot.lead->GetChannel() << " fiber_SA_ch: " << fiber_SA_ch << " fiber:
                    // " << fiber_id
                    //     << endl;
                    // TODO: Use it_sub->direction to find real life coordinates.

                    // Fix fiber installation mistakes.
                    fiber_id = FixMistake(fiber_id);

                    // Calibrate hit fiber.
                    Double_t tot_mapmt_raw = mapmt_tot.tot_ns;
                    Double_t tot_spmt_raw = spmt_tot.tot_ns;
                    auto tot_mapmt = mapmt_tot.tot_ns;
                    auto tot_spmt = spmt_tot.tot_ns;
                    Double_t t_mapmt = mapmt_tot.lead_ns;
                    Double_t t_spmt = spmt_tot.lead_ns;

                    Double_t t_mapmt1 = mapmt_tot.lead->GetTime_ns(); // MAPMT time without subtraction of trigger time
                    Double_t t_spmt1 = spmt_tot.lead->GetTime_ns();   // SAPMT time without subtraction of trigger time

                    // cout << "ToT fiber_SA_ch: " << tot_spmt << "  ToT multi: " << tot_mapmt << endl;
                    // cout << "Time fiber_SA_ch: " << t_spmt << "  Time multi: " << t_mapmt
                    //	<< "  ts: " << t_spmt1 << " tm: " << t_mapmt1 <<  endl;
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

                    fh_time_MA_Fib->Fill(fiber_id, t_mapmt);
                    fh_time_SA_Fib->Fill(fiber_SA_ch, t_spmt);

                    // cout << fName << " tS " << t_spmt << " tM: " << t_mapmt << " totS: " <<
                    //	tot_spmt << " totM: " << tot_mapmt << endl;
                    gainMA = 10.;
                    gainSA = 10.;
                    tsync = 0.;
                    offset1 = 0.;
                    offset2 = 0.;

                    if (!fIsCalibrator && fHitPar)
                    {
                        R3BBunchedFiberHitModulePar* par = fHitPar->GetModuleParAt(fiber_id);
                        if (par)
                        {
                            gainMA = par->GetGainMA();
                            gainSA = par->GetGainS();
                            tsync = par->GetSync();
                            offset1 = par->GetOffset1();
                            offset2 = par->GetOffset2();
                        }
                    }

                    tot_mapmt *= 10. / gainMA;
                    tot_spmt *= 10. / gainSA;
                    t_mapmt += offset1;
                    t_spmt += offset2;
                    tof = (t_mapmt + t_spmt) / 2.;
                    // tof = t_mapmt;
                    tof -= tsync;

                    // histogram for offset determination
                    fh_dt_Fib->Fill(fiber_id, t_spmt - t_mapmt);

                    // Fill histograms for gain match, and for debugging.
                    fh_Fib_ToF->Fill(fiber_id, tof);
                    fh_ToT_MA_Fib->Fill(fiber_id, tot_mapmt);
                    fh_ToT_SA_Fib->Fill(fiber_id, tot_spmt);

                    // Int_t numFibs = fSubNum * fChPerSub[0];
                    Int_t numFibs = fSubNum * fChPerSub[0] * fChPerSub[1];
                    Double_t x = -10000.;
                    Double_t y = -10000.;
                    Double_t veff = 12. / 2.; // cm/ns
                    Double_t randx = (std::rand() / (float)RAND_MAX);
                    if (fName == "Fi23a" || fName == "Fi23b")
                    {
                        Float_t fiber_thickness = 0.028000;
                        Int_t fiber_nbr = 384;
                        Float_t dead_layer = 0.9;
                        Float_t air_layer = 0.01; // relative to fiber_thickness
                        Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);
                        if (fDirection == VERTICAL)
                        {
                            x = -detector_width / 2. + fiber_thickness * (0.5 - randx) +
                                ((fiber_id - 1) + ((fiber_id - 1) * air_layer)) * fiber_thickness;
                            y = (t_spmt - t_mapmt) * veff;
                        }
                        else
                        {
                            x = (t_spmt - t_mapmt) * veff;

                            y = -detector_width / 2. + fiber_thickness * (0.5 - randx) +
                                ((fiber_id - 1) + ((fiber_id - 1) * air_layer)) * fiber_thickness;
                        }
                    }
                    if (fName == "Fi30" || fName == "Fi31" || fName == "Fi32" || fName == "Fi33")
                    {
                        Float_t fiber_thickness = 0.10000; // cm
                        Int_t fiber_nbr = 512;
                        Float_t dead_layer = 0.9;
                        Float_t air_layer = 0.01; // relative to fiber_thickness
                        Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);

                        if (fDirection == VERTICAL)
                        {
                            x = -detector_width / 2. + fiber_thickness * (0.5 - randx) +
                                ((fiber_id - 1) + ((fiber_id - 1) * air_layer)) * fiber_thickness;
                            y = (t_spmt - t_mapmt) * veff;
                            // y = 0.;
                        }
                        else
                        {
                            x = (t_spmt - t_mapmt) * veff;
                            // x = 0.;
                            y = -detector_width / 2. + fiber_thickness * (0.5 - randx) +
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
                    Double_t eloss = sqrt(tot_mapmt * tot_spmt);
                    Double_t t = tof;

                    energy[fiber_id] = eloss;
                    counts[fiber_id] = counts[fiber_id] + 1;
                    multi++;

                    //  cout << "save fiber " << fName<< "  "  << fiber_id << " pos " << x << endl;
                    if (!fIsCalibrator)
                    {
                        new ((*fHitItems)[fHitItems->GetEntriesFast()])
                            R3BBunchedFiberHitData(0, x, y, eloss, t, fiber_id, t_mapmt, t_spmt, tot_mapmt, tot_spmt);

                        // cout << fName << " x: " << x << " y: " << y << endl;
                    }
                }
            }
        }
    }

    fnEvents++;

    // cout<<"end exec"<<endl;
}

void R3BBunchedFiberCal2Hit_s494::FinishEvent()
{
    fHitItems->Clear();
    fNofHitItems = 0;
}

void R3BBunchedFiberCal2Hit_s494::FinishTask()
{
    fh_ToT_MA_Fib->Write();
    fh_ToT_SA_Fib->Write();
    fh_dt_Fib->Write();
    fh_time_SA_Fib->Write();
    fh_Fib_ToF->Write();
    fh_Test->Write();
    fh_time_MA_Fib->Write();

    /*

        fh_multi->Write();
    */

    if (fIsCalibrator)
    {

        R3BBunchedFiberHitModulePar* mpar;

        UInt_t max = N_FIBER_MAX;
        if (fh_ToT_MA_Fib->GetNbinsX() < N_FIBER_MAX)
            max = fh_ToT_MA_Fib->GetNbinsX();

        for (UInt_t i = 1; i <= max; i++)
        {
            mpar = new R3BBunchedFiberHitModulePar();
            mpar->SetFiber(i);
            fCalPar->AddModulePar(mpar);
        }

        for (UInt_t i = 1; i <= max; i++)
        {
            //   R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);

            // gain MA
            TH1D* proj = fh_ToT_MA_Fib->ProjectionY("", i + 1, i + 1, 0);
            for (UInt_t j = proj->GetNbinsX() - 2; j > 2; j--)
            {
                if (j == 2)
                {
                    // could not find maximum
                }
                if (proj->GetBinContent(j) > proj->GetMaximum() * 10. / 100.)
                {
                    R3BBunchedFiberHitModulePar* par = fCalPar->GetModuleParAt(i);
                    par->SetGainMA(proj->GetBinCenter(j));

                    if (i < 512)
                        cout << fName << " Fiber: " << i << " gain MA: " << proj->GetBinCenter(j) << endl;
                    break;
                }
            }
            // gain SA
            TH1D* proj1 = fh_ToT_SA_Fib->ProjectionY("", i + 1, i + 1, 0);
            for (UInt_t j = proj1->GetNbinsX() - 2; j > 2; j--)
            {
                if (j == 2)
                {
                    // could not find maximum
                }
                if (proj1->GetBinContent(j) > proj1->GetMaximum() * 10. / 100.)
                {
                    R3BBunchedFiberHitModulePar* par1 = fCalPar->GetModuleParAt(i);
                    par1->SetGainS(proj->GetBinCenter(j));

                    if (i < 512)
                        cout << fName << " Fiber: " << i << " gain SA: " << proj1->GetBinCenter(j) << endl;
                    break;
                }
            }
            // time offset
            R3BBunchedFiberHitModulePar* par2 = fCalPar->GetModuleParAt(i);
            TH1D* proj2 = fh_dt_Fib->ProjectionY("", i + 1, i + 1, 0);
            par2->SetOffset1(0.5 * proj2->GetBinCenter(proj2->GetMaximumBin()));
            par2->SetOffset2(-0.5 * proj2->GetBinCenter(proj2->GetMaximumBin()));

            // tsync
            R3BBunchedFiberHitModulePar* par3 = fCalPar->GetModuleParAt(i);
            TH1D* proj3 = fh_Fib_ToF->ProjectionY("", i + 1, i + 1, 0);
            par3->SetSync(proj3->GetBinCenter(proj3->GetMaximumBin()));

            if (i < 512)
                cout << fName << " Fiber: " << i << ", tsync" << proj3->GetBinCenter(proj3->GetMaximumBin())
                     << ", offset: " << 0.5 * proj2->GetBinCenter(proj2->GetMaximumBin()) << endl;
        }

        fCalPar->setChanged();
    }
}

void R3BBunchedFiberCal2Hit_s494::MAPMTTriggerMapSet(unsigned const* a_map, size_t a_map_bytes)
{
    size_t const c_length = a_map_bytes / sizeof *a_map;
    if (c_length != fSubNum * fChPerSub[0])
    {
        LOG(fatal) << "MAPMT trigger map length=" << c_length << " not compatible with sub#=" << fSubNum
                   << " mapmt_per_sub=" << fChPerSub[0] << "!";
    }
    fMAPMTTriggerMap = a_map;
}

ClassImp(R3BBunchedFiberCal2Hit_s494);
