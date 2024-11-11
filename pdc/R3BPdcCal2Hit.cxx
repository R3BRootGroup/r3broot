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
// ------------------------------------------------------------
// -----                 R3BPdcCal2Hit                   -----
// -----            Created June 2020 by M.Heil            -----
// ------------------------------------------------------------

#include "R3BPdcCal2Hit.h"
#include "R3BPdcCalData.h"
#include "R3BPdcHitData.h"
#include "R3BPdcHitModulePar.h"
#include "R3BPdcHitPar.h"

#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"

#include "TClonesArray.h"
#include "TMath.h"
#include "TSystem.h"

#include <iostream>
#include <sstream>
#include <string>

#include "mapping_pdc_trig.hh"

extern "C"
{
//#include "ext_data_client.h"
#include "ext_h101_pdc.h"
}
#define LENGTH(x) (sizeof x / sizeof *x)
#define M_PI 3.14159265358979323846

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

R3BPdcCal2Hit::ToT::ToT(R3BPdcCalData const* a_lead,
                        R3BPdcCalData const* a_trail,
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

R3BPdcCal2Hit::R3BPdcCal2Hit(const char* name, Bool_t a_is_calibrator, Bool_t a_is_sync, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems()
    , fCalTriggerItems()
    , fHitItems(new TClonesArray("R3BPdcHitData"))
    , fIsCalibrator(a_is_calibrator)
    , fIsSync(a_is_sync)
    , fNofHitItems()
    , fNofHitPars()
    , fHitPar()
    , fnEvents(0)
    , fClockFreq(1. / CTDC_16_CLOCK_MHZ * 1000.)
{

    pdc_trig_map_setup();
}

R3BPdcCal2Hit::~R3BPdcCal2Hit()
{

    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
    // delete fHitPar;
}

InitStatus R3BPdcCal2Hit::Init()
{

    LOG(info) << "PDC at Init() start" << endl;

    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = (R3BEventHeader*)mgr->GetObject("R3BEventHeader");
    fCalItems = (TClonesArray*)mgr->GetObject("PdcCal");
    if (NULL == fCalItems)
        LOG(fatal) << "Branch PdcCal not found";
    fCalTriggerItems = (TClonesArray*)mgr->GetObject("PdcTriggerCal");
    if (NULL == fCalTriggerItems)
        LOG(fatal) << "Branch PdcTriggerCal not found";
    maxevent = mgr->CheckMaxEventNo();

    // request storage of Hit data in output tree
    mgr->Register("PdcHit", "Land", fHitItems, kTRUE);

    if ((fIsCalibrator && !fIsSync) || !fIsCalibrator)
    {
        fHitPar = (R3BPdcHitPar*)FairRuntimeDb::instance()->getContainer("PdcHitPar");
        if (!fHitPar)
        {
            LOG(error) << "Could not get access to PdcHitPar-Container.";
            fNofHitPars = 0;
            // return kFATAL;
        }
        else
        {
            fNofHitPars = fHitPar->GetNumModulePar();
            if (fNofHitPars == 0)
            {
                LOG(error) << "There are no Hit parameters in container PdcHitPar";
                fHitPar = nullptr;
                // return kFATAL;
            }
        }
    }

    if (fIsCalibrator && !fIsSync && fHitPar)
    {
        for (Int_t iplane = 0; iplane < N_PLANE_MAX_PDC; iplane++)
        {
            for (Int_t iwire = 1; iwire <= N_WIRE_MAX; iwire++)
            {
                Int_t ii = iplane * N_WIRE_MAX + iwire;
                R3BPdcHitModulePar* par = fHitPar->GetModuleParAt(iplane * N_WIRE_MAX + iwire);
                if (par)
                {
                    tsync_mem[ii - 1] = par->GetSync();
                    // cout<<"Tsync_mem: "<<iplane<<", "<<iwire<<"; "<<ii<<", "<<tsync_mem[ii-1] <<endl;
                }
            }
        }
    }

    auto plane_num = LENGTH(EXT_STR_h101_PDC_onion::PDC_P);
    fPlaneArray.resize(plane_num);
    auto ch_num = LENGTH(EXT_STR_h101_PDC_onion::PDC_P[0].TLCMI);
    for (auto it = fPlaneArray.begin(); fPlaneArray.end() != it; ++it)
    {
        it->resize(ch_num);
    }

    char chistName[255];
    char chistTitle[255];
    for (Int_t j = 0; j < N_PLANE_MAX_PDC; j++)
    {
        sprintf(chistName, "pdc_time_raw_plane_%d", j + 1);
        sprintf(chistTitle, "PDC raw times plane %d", j + 1);
        fh_time_raw[j] = new TH2F(chistName, chistTitle, 150, 0, 150., 2400, -1400., 600.);
        fh_time_raw[j]->GetXaxis()->SetTitle("wire");
        fh_time_raw[j]->GetYaxis()->SetTitle("time / ns");

        sprintf(chistName, "pdc_tsync_plane_%d", j + 1);
        sprintf(chistTitle, "PDC tsync plane %d", j + 1);
        fh_tsync[j] = new TH2F(chistName, chistTitle, 150, 0, 150., 500, -500., 500.);
        fh_tsync[j]->GetXaxis()->SetTitle("wire");
        fh_tsync[j]->GetYaxis()->SetTitle("tsync / ns");

        sprintf(chistName, "pdc_radius_%d", j + 1);
        sprintf(chistTitle, "PDC r1+r2 %d", j + 1);
        fh_radius[j] = new TH2F(chistName, chistTitle, 150, 0., 150., 200, 0., 20.);
        fh_radius[j]->GetXaxis()->SetTitle("wire");
        fh_radius[j]->GetYaxis()->SetTitle("r1+r2 in mm");

        sprintf(chistName, "pdc_radius_dt_%d", j + 1);
        sprintf(chistTitle, "PDC distance between two wires %d", j + 1);
        fh_radius_dt[j] = new TH2F(chistName, chistTitle, 500, -1000, 1000., 200, 0., 20.);
        fh_radius_dt[j]->GetYaxis()->SetTitle("ri in mm");
        fh_radius_dt[j]->GetXaxis()->SetTitle("ti in ns");

        sprintf(chistName, "pdc_mult_%d", j + 1);
        sprintf(chistTitle, "PDC multiplicity %d", j + 1);
        fh_mult[j] = new TH2F(chistName, chistTitle, 150, 0., 150., 10, 0., 10.);
        fh_mult[j]->GetXaxis()->SetTitle("wire");
        fh_mult[j]->GetYaxis()->SetTitle("multiplicity");

        sprintf(chistName, "fh_time_drift_%d", j + 1);
        sprintf(chistTitle, "PDC drift time difference %d", j + 1);
        fh_time_drift[j] = new TH2F(chistName, chistTitle, 150, 0., 150., 1200, -500., 500.);
        fh_time_drift[j]->GetXaxis()->SetTitle("wire");
        fh_time_drift[j]->GetYaxis()->SetTitle("dt in ns");

        sprintf(chistName, "pdc_fish_%d", j + 1);
        sprintf(chistTitle, "PDC fish %d", j + 1);
        fh_fish[j] = new TH2F(chistName, chistTitle, 600, -1000, 1000., 600, -1000., 1000.);
        fh_fish[j]->GetXaxis()->SetTitle("t1+t2");
        fh_fish[j]->GetYaxis()->SetTitle("t1-t2");

        if (method2)
        {
            sprintf(chistName, "pdc_check_%d", j + 1);
            sprintf(chistTitle, "PDC check %d", j + 1);
            fh_check[j] = new TH2F(chistName, chistTitle, 200, 0, 20., 600, -1000., 1000.);
            fh_check[j]->GetXaxis()->SetTitle("r1+r2");
            fh_check[j]->GetYaxis()->SetTitle("t1-t2");
        }
        else
        {
            sprintf(chistName, "pdc_check_%d", j + 1);
            sprintf(chistTitle, "PDC check %d", j + 1);
            fh_check[j] = new TH2F(chistName, chistTitle, 600, -1000., 1000., 200, -10, 10.);
            fh_check[j]->GetXaxis()->SetTitle("t1+t2");
            fh_check[j]->GetYaxis()->SetTitle("r2-r1");
        }

        sprintf(chistName, "pdc_r1_vs_r2_%d", j + 1);
        sprintf(chistTitle, "PDC r1 vs r2 %d", j + 1);
        fh_rvspos[j] = new TH2F(chistName, chistTitle, 100, 0, 10., 100, 0., 10.); // 200, 0, 20., 3000, 200., 800.);
        fh_rvspos[j]->GetXaxis()->SetTitle("r1 / mm");
        fh_rvspos[j]->GetYaxis()->SetTitle("r2 / mm");

        sprintf(chistName, "pdc_angle_%d", j + 1);
        sprintf(chistTitle, "PDC angle %d", j + 1);
        fh_angle[j] = new TH1F(chistName, chistTitle, 2000, -1, 1);
        fh_angle[j]->GetXaxis()->SetTitle("Reconstructed angle / rad");
        fh_angle[j]->GetYaxis()->SetTitle("Counts");

        sprintf(chistName, "pdc_xtc_error_%d", j + 1);
        sprintf(chistTitle, "PDC xtc error %d", j + 1);
        fh_error[j] = new TH1F(chistName, chistTitle, N_WIRE_MAX, 0, N_WIRE_MAX);
        fh_error[j]->GetXaxis()->SetTitle("WIRE NUMBER");
        fh_error[j]->GetYaxis()->SetTitle("error");
    }
    if (!fIsSync)
    {
        TCanvas* cTrigg = new TCanvas("PDCcal2hit", "PDCcal2hit", 50, 50, 1000, 1000);
        cTrigg->Divide(5, 4);
        for (Int_t j = 0; j < N_PLANE_MAX_PDC; j++)
        {
            cTrigg->cd(1 + j * (N_PLANE_MAX_PDC + 1));
            gPad->SetLogz();
            fh_fish[j]->Draw("colz");
            cTrigg->cd(2 + j * (N_PLANE_MAX_PDC + 1));
            gPad->SetLogz();
            fh_check[j]->Draw("colz");
            cTrigg->cd(3 + j * (N_PLANE_MAX_PDC + 1));
            gPad->SetLogz();
            fh_radius_dt[j]->Draw("colz");
            cTrigg->cd(4 + j * (N_PLANE_MAX_PDC + 1));
            gPad->SetLogz();
            fh_radius[j]->Draw("colz");
            cTrigg->cd(5 + j * (N_PLANE_MAX_PDC + 1));
            gPad->SetLogz();
            fh_rvspos[j]->Draw("colz");
        }
        cTrigg->cd(0);
    }

    sprintf(chistName, "pdc_running_sum");
    sprintf(chistTitle, "PDC running sum");
    fh_running_sum = new TH1F(chistName, chistTitle, 3000, -500, 2000);
    fh_running_sum->GetXaxis()->SetTitle("time");
    fh_running_sum->GetYaxis()->SetTitle("radius");

    sprintf(chistName, "pdc_xtc");
    sprintf(chistTitle, "PDC xtc");
    fh_xtc = new TH1F(chistName, chistTitle, 3000, -500, 2000);
    fh_xtc->GetXaxis()->SetTitle("time");
    fh_xtc->GetYaxis()->SetTitle("radius");

    sprintf(chistName, "pdc_wire_cor");
    sprintf(chistTitle, "PDC wire correlation");
    fh_wire_cor = new TH2F(chistName, chistTitle, 150, 0, 150., 150, 0., 150.);
    fh_wire_cor->GetXaxis()->SetTitle("wire");
    fh_wire_cor->GetYaxis()->SetTitle("wire");

    return kSUCCESS;
}

InitStatus R3BPdcCal2Hit::ReInit() { return kSUCCESS; }

// Note that the container may still be empty at this point.
void R3BPdcCal2Hit::SetParContainers()
{
    fHitPar = (R3BPdcHitPar*)FairRuntimeDb::instance()->getContainer("PdcHitPar");
    if (!fHitPar)
    {
        LOG(error) << "R3BPdcCal2Hit::Init() Could not get access to PdcHitPar-Container.";
        fNofHitPars = 0;
        return;
    }
    else
    {
        LOG(info) << "R3BPdcCal2Hit::Init() container PdcHitPar initialized" << endl;
    }
    LOG(info) << "PDC SetParContainers() end" << endl;
}

void R3BPdcCal2Hit::Exec(Option_t* option)
{
    if (0 == fnEvents % 10000)
    {
        std::cout << "\rEvents: " << fnEvents << " / " << maxevent << " (" << (int)(fnEvents * 100. / maxevent)
                  << " %) " << std::flush;
    }
    for (auto it = fPlaneArray.begin(); fPlaneArray.end() != it; ++it)
    {
        auto& ch_arr = *it;
        for (auto it2 = ch_arr.begin(); ch_arr.end() != it2; ++it2)
        {
            it2->lead_list.clear();
            it2->tot_list.clear();
        }
    }

    // Make direct mapping tables for trigger items.
    size_t trig_num = fCalTriggerItems->GetEntries();
    // cout << "trig num: " << trig_num << endl;
    if (trig_num < 4)
    {
        // LOG(error) << "No trigger times found!!!" << endl;
        // return;
    }

    std::vector<R3BPdcCalData const*> trig_table(LENGTH(EXT_STR_h101_PDC_onion::PDC_P));
    for (size_t j = 0; j < trig_num; ++j)
    {
        auto cal = (R3BPdcCalData const*)fCalTriggerItems->At(j);
        // cout << "cal: " << cal << " wire: " << cal->GetWireId() << endl;
        trig_table.at(cal->GetWireId() - 1) = cal;
    }

    // Read in data
    // Find multi-hit ToT for every channel.
    // The easiest safe way to survive ugly cases is to record all
    // leading edges per channel, and then pair up with whatever
    // trailing we have.
    // Not super efficient, but shouldn't crash if the data is not
    // perfect.
    size_t cal_num = fCalItems->GetEntriesFast();
    unsigned n_lead = 0;
    unsigned n_trail = 0;
    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = (R3BPdcCalData const*)fCalItems->At(j);
        if (cur_cal->GetEdgeId() == 1)
        {
            ++n_lead;
            auto plane_i = cur_cal->GetPlaneId() - 1;
            auto wire_i = cur_cal->GetWireId() - 1;
            auto& channel = fPlaneArray.at(plane_i).at(wire_i);
            channel.lead_list.push_back(cur_cal);
            // cout << "leading wire:" << wire_i <<  " time: " << cur_cal->GetTime_ns() << endl;
        }
        else
        {
            ++n_trail;
        }
    }
    // if (n_lead != n_trail+1 && n_lead != 0) // one more leading for trigger
    if (n_lead != n_trail)
    {
        // cout << "leading: " << n_lead << " trailing: " <<  n_trail << endl;
        // cout << "Number of leading edges not equal to number of trailing edges!" << endl;
        return;
    }

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = (R3BPdcCalData const*)fCalItems->At(j);
        if (cur_cal->GetEdgeId() == 2)
        {
            Double_t c_period = 4096. * (1000. / fClockFreq);

            auto plane_i = cur_cal->GetPlaneId() - 1;
            auto wire_i = cur_cal->GetWireId() - 1;
            auto& channel = fPlaneArray.at(plane_i).at(wire_i);
            // cout << "trailing wire:" << wire_i <<  " time: " << cur_cal->GetTime_ns() << endl;
            if (channel.lead_list.empty())
            {
                // cout << "lead list empty for plane, wire " << plane_i + 1 << "  " << wire_i + 1 << endl;

                continue;
            }
            auto lead = channel.lead_list.front();

            Double_t cur_cal_trig_ns = 0;
            Double_t lead_trig_ns = 0;

            assert(plane_i == lead->GetPlaneId() - 1);
            assert(wire_i == lead->GetWireId() - 1);

            // cout << "Plane: " << plane_i << endl;
            // cout << "wire: " << wire_i << endl;

            // Trigger time is the same for leading and trailing,
            // but we still need to subtract it for when we want the
            // edge times.

            auto cur_cal_trig_i = g_pdc_trig_map[plane_i][wire_i];
            // cout << "cur_cal_trig_i: " << cur_cal_trig_i << endl;
            // if(cur_cal_trig_i == 1) cur_cal_trig_i=0; // please remove again MH

            auto cur_cal_trig = trig_table.at(cur_cal_trig_i);
            // cout << "cur_cal_trig: " << cur_cal_trig << endl;

            cur_cal_trig_ns = cur_cal_trig->GetTime_ns();
            // cout << "trigger time: " << cur_cal_trig_ns << endl;
            auto cur_cal_ns =
                fmod(cur_cal->GetTime_ns() - cur_cal_trig_ns + c_period + c_period / 2, c_period) - c_period / 2;
            auto lead_ns =
                fmod(lead->GetTime_ns() - cur_cal_trig_ns + c_period + c_period / 2, c_period) - c_period / 2;
            auto tot_ns = fmod(cur_cal_ns - lead_ns + c_period + c_period / 2, c_period) - c_period / 2;

            /*
            cout << plane_i << ' ' << wire_i << '\n';
            cout << "Lead: " << lead_ns << " = " << lead->GetTime_ns() << " - " << cur_cal_trig_ns << endl;
            cout << "Trail: " << cur_cal_ns << " = " << cur_cal->GetTime_ns() << " - " << cur_cal_trig_ns << endl;
            cout << "ToT: " << tot_ns << " = " << cur_cal_ns << " - " << lead_ns << endl;
            */
            if (tot_ns > 0 && tot_ns < 1000)
            {
                channel.tot_list.push_back(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));
                channel.lead_list.pop_front();
            }
        }
    }

    // cout << "***************** new event ****************" << endl;
    Int_t ID_old = 0;
    Int_t ID_new = 0;
    Double_t t_old = 0.;
    Double_t t_new = 0.;
    Double_t tsync_old = 0. / 0.;
    Double_t tsync_new = 0. / 0.;
    Double_t tsync_pdc = 0. / 0.;
    Double_t r_old = 0.;
    Double_t r_new = 0.;
    Double_t tot_old = 0.;
    Double_t tot_new = 0.;
    Int_t wire_old = 0;
    Int_t wire_new = 0;

    for (auto it = fPlaneArray.begin(); fPlaneArray.end() != it; ++it)
    {
        auto& ch_arr = *it;
        for (auto it2 = ch_arr.begin(); ch_arr.end() != it2; ++it2)
        {
            auto const& channel = *it2;
            Int_t mult = 0;
            Bool_t first = true;
            Double_t firstLeading = 0.;
            Double_t lastTrailing = 0.;
            Double_t TotMax = 0.;
            Double_t sumTot = 0.;

            for (auto it_pdc_tot = channel.tot_list.begin(); channel.tot_list.end() != it_pdc_tot;
                 ++it_pdc_tot) // over ihit(fiber)
            {
                auto const& pdc_tot = *it_pdc_tot;
                auto plane = pdc_tot.lead->GetPlaneId();
                auto wire = pdc_tot.lead->GetWireId();
                // cout << "Cal2Hit plane: " << plane << " wire: " << wire << endl;
                /*
                                cout << "wire before: " << wire << endl;
                                if(wire%8 == 1) wire=wire;
                                else if(wire%8 == 2) wire=wire+1;
                                else if(wire%8 == 3) wire=wire+2;
                                else if(wire%8 == 4) wire=wire+3;
                                else if(wire%8 == 5) wire=wire-3;
                                else if(wire%8 == 6) wire=wire-2;
                                else if(wire%8 == 7) wire=wire-1;
                                else if(wire%8 == 0) wire=wire;
                                cout << "wire after: " << wire << endl;
                */
                auto tot_pdc = pdc_tot.tot_ns;
                Double_t t_pdc = pdc_tot.lead_ns;

                mult++;

                // make wire correlation plot
                /*
                                auto& ch_arr1 = *it;
                                for (auto it3 = ch_arr1.begin(); ch_arr1.end() != it3; ++it3)
                                {
                                    auto const& channel1 = *it3;

                                    for (auto it_pdc_tot1 = channel1.tot_list.begin(); channel1.tot_list.end() !=
                   it_pdc_tot1;
                                         ++it_pdc_tot1) // over ihit(fiber)
                                    {
                                        auto const& pdc_tot1 = *it_pdc_tot1;
                                        auto plane1 = pdc_tot1.lead->GetPlaneId();
                                        auto wire1 = pdc_tot1.lead->GetWireId();

                                        //cout << "wire1 before: " << wire1 << endl;
                                        //if(wire1%8 == 1) wire=wire;
                                        //else if(wire1%8 == 2) wire1=wire1+1;
                                        //else if(wire1%8 == 3) wire1=wire1+2;
                                        //else if(wire1%8 == 4) wire1=wire1+3;
                                        //else if(wire1%8 == 5) wire1=wire1-3;
                                        //else if(wire1%8 == 6) wire1=wire1-2;
                                        //else if(wire1%8 == 7) wire1=wire1-1;
                                        //else if(wire1%8 == 0) wire1=wire1;

                                        auto tot_pdc1 = pdc_tot1.tot_ns;
                                        Double_t t_pdc1 = pdc_tot1.lead_ns;

                                        if(wire != wire1) fh_wire_cor->Fill(wire, wire1);


                                    }
                                }
                */

                Double_t r = 0.;
                Int_t ID = 0;
                ID = plane;
                Double_t tsync = 0.;

                if (!fIsSync && fHitPar)
                {
                    // cout << "reading sync parameter" << endl;
                    R3BPdcHitModulePar* par = fHitPar->GetModuleParAt((ID - 1) * N_WIRE_MAX + wire);
                    if (par)
                    {
                        tsync = par->GetSync();
                    }
                }
                tsync_pdc = t_pdc - tsync;
                // cout<<"Tsync: "<<ID<<", "<<wire<<", "<<tsync<<", "<<t_pdc<<", "<<tsync_pdc<<endl;

                // cout  << "plane: " << plane <<  " time: " << t_pdc << " wire: " << wire << " ToT: " << tot_pdc <<
                // endl;

                // cuts
                if (tot_pdc < 35.)
                {
                    // continue;
                }
                if (!fIsSync && (tsync_pdc < -50. || tsync_pdc > 250.))
                {
                    // continue;
                }

                // store first hit in wire
                // if(1>0) // v48
                // if(first)  //v49
                if (tot_pdc > TotMax) // v50
                {
                    TotMax = tot_pdc;
                    ID_new = ID;
                    t_new = t_pdc;
                    tot_new = tot_pdc;
                    wire_new = wire;
                    tsync_new = tsync_pdc;
                    firstLeading = t_pdc;
                    first = false;
                }
                lastTrailing = tsync_pdc + tot_pdc;
                sumTot += tot_pdc;

                Double_t x = 0.;
                Double_t y = 0.;

                if (ID == 1 || ID == 3)
                {
                    x = wire * wire_distance;
                    y = 0.;
                }
                if (ID == 2 || ID == 4)
                {
                    y = wire * wire_distance;
                    x = 0;
                }

                // MH
                new ((*fHitItems)[fNofHitItems++]) R3BPdcHitData(t_pdc, x, y, tot_pdc, ID, wire);

                // cout << "cal2hit plane: " << ID << " wire: " << wire << " ToT: " << tot_pdc << " t: " << t_pdc <<
                // endl; cout << firstLeading << "  " << lastTrailing << endl;
            }
            Double_t totalToT = lastTrailing - firstLeading;
            if (mult > 0)
            {
                // tot_new = sumTot;
                fh_mult[ID_new - 1]->Fill(wire_new, mult);
                // cout << "Test1: " << t_new << "  " << tot_new << "  " << sumTot << endl;
            }

            if (ID_new == ID_old && wire_new - 1 == wire_old)
            {
                // cout  << "found Pair: " << ID_new << "  " << wire_new << "  " << wire_old << endl;
                Double_t pos = 0.;
                Double_t x = 0.;
                Double_t y = 0.;
                Double_t dtmin = 0.;
                Double_t dtmax = 0.;
                Double_t tmin = 0.;
                Double_t tmax = 0.;
                Double_t r = -1.;
                Double_t r1 = -1.;
                Double_t r2 = -1.;

                fh_time_raw[ID_new - 1]->Fill(wire_new, t_new);
                fh_time_raw[ID_new - 1]->Fill(wire_old, t_old);

                fh_tsync[ID_new - 1]->Fill(wire_new, tsync_new);
                fh_tsync[ID_new - 1]->Fill(wire_old, tsync_old);

                fh_time_drift[ID_new - 1]->Fill(wire_new, tsync_new - tsync_old);
                Double_t tdrift = tsync_new - tsync_old;

                if (!fIsCalibrator && fHitPar)
                {
                    Int_t xtc_npoints = 0;
                    Double_t x_xtc[300] = { 0 };
                    Double_t t_xtc[300] = { 0 };
                    Double_t apar, bpar;

                    if (!method2)
                    {
                        // Version MH

                        R3BPdcHitModulePar* par = fHitPar->GetModuleParAt((ID_new - 1) * N_WIRE_MAX + wire_new);
                        if (par)
                        {
                            dtmin = par->GetdTmin();
                            dtmax = par->GetdTmax();
                            // cout  << "Parameter: " << dtmin << "  " << dtmax << endl;
                        }
                        // calculate radius
                        Double_t m = wire_distance / (dtmax - dtmin);
                        Double_t b = wire_distance / 2. - m * dtmax;
                        r = m * (tsync_new - tsync_old) + b;
                        // cout  << "dt: " << (t_new - t_old) << " Radius: " << r << endl;

                        fh_radius_dt[ID_new - 1]->Fill(wire_new, r);

                        // end version MH

                        // Version AK
                        /*
                        R3BPdcHitModulePar* par = fHitPar->GetModuleParAt((ID_new - 1) * N_WIRE_MAX + wire_new);
                        if (par)
                        {
                            dtmin = par->GetTmin();
                            dtmax = par->GetTmax();
                            xtc_npoints = par->GetNPoints();
                            xtc_x = par->GetXT_xArray();
                            xtc_t = par->GetXT_tArray();
                            for (Int_t ic = 0; ic < xtc_npoints; ic++)
                            {
                                x_xtc[ic] = *(xtc_x + ic);
                                t_xtc[ic] = *(xtc_t + ic);
                            }
                            LOG(debug) << "Parameter: " << dtmin << "  " << dtmax << endl;

                            Double_t ttemp1 = 0, ttemp2 = 0;
                            if (xtc_npoints < 1)
                                r = 0. / 0.;
                            if (tdrift < dtmin && dtmin < dtmax)
                            {
                                r = -wire_distance / 2.; //MH sollte das nicht durch 2 dividiert werden?
                            }
                            else if (tdrift > dtmax && dtmin < dtmax)
                            {
                                r = wire_distance;
                            }
                            else
                            {
                                for (Int_t ic = 0; ic < xtc_npoints - 1; ic++)
                                {

                                    if (tdrift == t_xtc[ic])
                                        r = x_xtc[ic];
                                    else if (tdrift > t_xtc[ic] && tdrift < t_xtc[ic + 1])
                                    {
                                        Double_t dx = x_xtc[ic + 1] - x_xtc[ic];
                                        Double_t dt = t_xtc[ic + 1] - t_xtc[ic];

                                        ttemp1 = t_xtc[ic];
                                        ttemp2 = t_xtc[ic + 1];

                                        apar = dx / dt;
                                        bpar = x_xtc[ic] - apar * t_xtc[ic];

                                        r = apar * tdrift + bpar;
                                    }
                                }
                            }

                            fh_radius_dt[ID_new - 1]->Fill(tdrift, r);
                            cout << "Test tdrift = : " << tdrift << "  r = " << r << endl;
                        }
                        */
                        // end version AK
                    }

                    if (method2)
                    {
                        R3BPdcHitModulePar* par1 = fHitPar->GetModuleParAt((ID_old - 1) * N_WIRE_MAX + wire_old);
                        if (par1)
                        {
                            dtmin = par1->GetTmin();
                            dtmax = par1->GetTmax();
                            xtc_npoints = par1->GetNPoints();
                            if (xtc_npoints > 0)
                                cout << "npoints: " << xtc_npoints << endl;
                            xtc_x = par1->GetXT_xArray();
                            xtc_t = par1->GetXT_tArray();
                            for (Int_t ic = 0; ic < xtc_npoints; ic++)
                            {
                                x_xtc[ic] = *(xtc_x + ic);
                                t_xtc[ic] = *(xtc_t + ic);
                            }

                            LOG(debug) << "Parameter: " << dtmin << "  " << dtmax << endl;
                        }

                        if (xtc_npoints < 1)
                            r1 = -1;
                        if (tsync_old < dtmin && dtmin < dtmax)
                        {
                            r1 = 0;
                        }
                        else if (tsync_old > dtmax && dtmin < dtmax)
                        {
                            r1 = wire_distance;
                        }
                        else
                        {
                            for (Int_t ic = 0; ic < xtc_npoints - 1; ic++)
                            {
                                if (tsync_old == t_xtc[ic])
                                    r1 = x_xtc[ic];
                                else if (tsync_old > t_xtc[ic] && tsync_old < t_xtc[ic + 1])
                                {
                                    Double_t dx = x_xtc[ic + 1] - x_xtc[ic];
                                    Double_t dt = t_xtc[ic + 1] - t_xtc[ic];

                                    apar = dx / dt;
                                    bpar = x_xtc[ic] - apar * t_xtc[ic];

                                    r1 = apar * tsync_old + bpar;
                                    fh_radius_dt[ID_new - 1]->Fill(tsync_old, r1);
                                }
                            }

                            LOG(debug) << "Radius1: " << r1 << endl;
                        }
                        // cout<<"Radius 1: "<<r1<<", "<<tsync_old<<", "<<xtc_npoints<<", "<<dtmin<<"; "<<dtmax<<endl;

                        R3BPdcHitModulePar* par2 = fHitPar->GetModuleParAt((ID_new - 1) * N_WIRE_MAX + wire_new);
                        if (par2)
                        {
                            dtmin = par2->GetTmin();
                            dtmax = par2->GetTmax();
                            xtc_npoints = par2->GetNPoints();
                            xtc_x = par2->GetXT_xArray();
                            xtc_t = par2->GetXT_tArray();
                            for (Int_t ic = 0; ic < xtc_npoints; ic++)
                            {
                                x_xtc[ic] = *(xtc_x + ic);
                                t_xtc[ic] = *(xtc_t + ic);
                            }
                            LOG(debug) << "Parameter: " << dtmin << "  " << dtmax << endl;
                        }
                        Double_t ttemp1 = 0, ttemp2 = 0;
                        if (xtc_npoints < 1)
                            r2 = -1;
                        if (tsync_new < dtmin && dtmin < dtmax)
                        {
                            r2 = 0.;
                        }
                        else if (tsync_new > dtmax && dtmin < dtmax)
                        {
                            r2 = wire_distance;
                        }
                        else
                        {
                            for (Int_t ic = 0; ic < xtc_npoints - 1; ic++)
                            {

                                if (tsync_new == t_xtc[ic])
                                    r2 = x_xtc[ic];
                                else if (tsync_new > t_xtc[ic] && tsync_new < t_xtc[ic + 1])
                                {
                                    Double_t dx = x_xtc[ic + 1] - x_xtc[ic];
                                    Double_t dt = t_xtc[ic + 1] - t_xtc[ic];

                                    ttemp1 = t_xtc[ic];
                                    ttemp2 = t_xtc[ic + 1];

                                    apar = dx / dt;
                                    bpar = x_xtc[ic] - apar * t_xtc[ic];

                                    r2 = apar * tsync_new + bpar;
                                    fh_radius_dt[ID_new - 1]->Fill(tsync_new, r2);
                                }
                            }

                            LOG(debug) << "Radius2: " << r2 << endl;
                        }

                        //  cout<<"Radius 2: "<<r2<<", "<<tsync_new<<", "<<xtc_npoints<<", "<<dtmin<<"; "<<dtmax<<",
                        //  "<<ttemp1<<", "<<ttemp2<<endl;

                        r = r1 + r2;
                    }
                }

                if (!fIsCalibrator && fHitPar)
                {
                    if (!method2)
                        pos = (wire_old * wire_distance + wire_distance / 2. - r);
                    if (method2 && r1 > 0. && r2 > 0.)
                    {
                        pos = ((wire_old * wire_distance + r1) + (wire_new * wire_distance - r2)) / 2.;
                        double arg = (r1 + r2) / (2. * wire_distance);
                        double angle = (arg >= 1) ? -M_PI / 3 : (M_PI / 6 - asin(arg));
                        fh_angle[ID_new - 1]->Fill(angle);
                    }
                }
                else
                {
                    pos = ((wire_new * wire_distance) + (wire_old * wire_distance)) / 2.;
                }

                if (ID_new == 1 || ID_new == 3)
                {
                    x = pos;
                    y = 0;
                }
                else
                {
                    x = 0;
                    y = pos;
                }

                /*
                 if((tsync_new - tsync_old) >=0 ) cout<<"dt positiv => t_old= "<<tsync_old<<"; t_new= "<<tsync_new<<",
                 r_old= "<<r1<<
                 ", r_new= "<<r2<<", w_old= "<<wire_old<<", w_new= "<<wire_new<<endl;

                 if((tsync_new - tsync_old) < 0 ) cout<<"dt negativ => t_old= "<<tsync_old<<"; t_new= "<<tsync_new<<",
                 r_old= "<<r1<<
                 ", r_new= "<<r2<<", w_old= "<<wire_old<<", w_new= "<<wire_new<<endl;
                */

                fh_fish[ID_new - 1]->Fill((tsync_new + tsync_old), tsync_new - tsync_old);
                fh_radius[ID_new - 1]->Fill(wire_new, r);
                if (method2)
                {
                    fh_check[ID_new - 1]->Fill(r, tsync_new - tsync_old);
                    fh_rvspos[ID_new - 1]->Fill(r1, r2);
                }
                else
                {
                    fh_check[ID_new - 1]->Fill(tsync_new + tsync_old, r);
                }

                if (!fIsCalibrator)
                {
                    // if (tsync_new > -100 && tsync_new < 300 && tsync_old > -100 && tsync_old < 300)
                    {
                        // MH
                        //                        new ((*fHitItems)[fNofHitItems++]) R3BPdcHitData(
                        //                            (tsync_new + tsync_old) / 2., x, y, sqrt(tot_new * tot_old),
                        //                            ID_new, wire_new);

                        // cout << "Hit level ID: " << ID_new << " x: " << x << " y: " << y << " ToT: " <<
                        //	sqrt(tot_new * tot_old) << " t: " << (tsync_new + tsync_old) / 2. << " wire: " << wire_new
                        //<< endl;
                    }
                }
            }
            ID_old = ID_new;
            t_old = t_new;
            r_old = r_new;
            tot_old = tot_new;
            wire_old = wire_new;
            tsync_old = tsync_new;
        }
    }

    // Sort hits im time.

    fnEvents++;
}

void R3BPdcCal2Hit::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
        fNofHitItems = 0;
    }
}

void R3BPdcCal2Hit::FinishTask()
{

    for (Int_t j = 0; j < N_PLANE_MAX_PDC; j++)
    {
        fh_time_raw[j]->Write();
        fh_radius[j]->Write();
        fh_radius_dt[j]->Write();
        fh_mult[j]->Write();
        fh_time_drift[j]->Write();
        fh_fish[j]->Write();
        fh_check[j]->Write();
        fh_rvspos[j]->Write();
        fh_tsync[j]->Write();
        fh_angle[j]->Write();
    }
    fh_wire_cor->Write();

    if (fIsCalibrator)
    {
        cout << "IN CALIBRATOR!" << endl;
        cout << "Test" << endl;

        if (fIsSync)
        {

            R3BPdcHitModulePar* mpar;

            for (Int_t plane = 0; plane < N_PLANE_MAX_PDC; plane++)
            {
                for (UInt_t i = 1; i <= N_WIRE_MAX; i++)
                {
                    mpar = new R3BPdcHitModulePar();
                    mpar->SetWire(plane * N_WIRE_MAX + i);
                    // cout << "New module: " << plane*N_WIRE_MAX + i << endl;
                    fHitPar->AddModulePar(mpar);
                }
            }
            // tsync
            for (Int_t plane = 0; plane < N_PLANE_MAX_PDC; plane++)
            {
                for (UInt_t i = 1; i <= N_WIRE_MAX; i++)
                {
                    Double_t tsync = -10000.;
                    TH1D* proj0 = fh_time_raw[plane]->ProjectionY("", i + 1, i + 1, 0);
                    // one can look for the highest bin an use this for sync
                    // Double_t tsync = proj0->GetBinCenter(proj0->GetMaximumBin());

                    // but one can also look at the sharp increasing edge and find the 50% value
                    Int_t numBins = proj0->GetNbinsX() - 2;
                    for (Int_t j = 1; j < numBins; j++)
                    {
                        if ((proj0->GetBinContent(j) / proj0->GetMaximum()) > 0.5 && tsync < -9999)
                        {
                            tsync = proj0->GetBinCenter(j);
                        }
                    }

                    R3BPdcHitModulePar* par0 = fHitPar->GetModuleParAt(plane * N_WIRE_MAX + i);
                    par0->SetSync(tsync);

                    // cout << "Plane: " << plane << "  " << plane * N_WIRE_MAX + i << endl;
                    // cout << "Set tsync  parameter of wire: " << i << " tsync: " << tsync << endl;
                }
            }
        }
        if (!fIsSync)
        {
            // Bool_t method2 = false; -> in header file included
            if (!method2)
            {
                // short version MH

                for (Int_t plane = 0; plane < N_PLANE_MAX_PDC; plane++)
                {
                    for (UInt_t i = 1; i <= N_WIRE_MAX; i++)
                    {
                        R3BPdcHitModulePar* par0 = fHitPar->GetModuleParAt(plane * N_WIRE_MAX + i);
                        Double_t tsync = par0->GetSync();

                        // find maximum value of drift times
                        TH1D* proj = fh_time_drift[plane]->ProjectionY("", i + 1, i + 1, 0);
                        // Double_t wire_distance = 7.8; // in mm
                        Double_t tmin = -10000.;
                        Double_t tmax = -10000.;
                        Int_t numBins = proj->GetNbinsX() - 2;
                        if (proj->GetSum() < 100)
                            continue;
                        for (UInt_t j = 1; j < numBins; j++)
                        {
                            // cout << "sum: " << sum << "  " << proj->GetBinContent(j)/proj->GetSum() * 8. << endl;
                            // cout << "Bin: " << proj->GetBinCenter(j) <<  "  " << proj->GetBinContent(j) << endl;
                            // cout << "Maximum: " << proj->GetMaximum() << endl;
                            // cout << "Test1: " << proj->GetBinContent(j)/proj->GetMaximum() << endl;
                            // cout << "Test2: " << proj->GetBinContent(numBins-j)/proj->GetMaximum() << endl;
                            if (proj->GetBinContent(j) / proj->GetMaximum() * 100. > 5. && tmin < -9999)
                            {
                                tmin = proj->GetBinCenter(j);
                            }
                            if (proj->GetBinContent(numBins - j) / proj->GetMaximum() * 100. > 5. && tmax < -9999)
                            {
                                tmax = proj->GetBinCenter(numBins - j);
                            }
                        }

                        // tmin -= 70.;
                        // tmax += 70.;

                        R3BPdcHitModulePar* par2 = fHitPar->GetModuleParAt(plane * N_WIRE_MAX + i);
                        par2->SetdTmin(tmin);
                        par2->SetdTmax(tmax);
                        par2->SetSync(tsync);

                        //  cout << "Plane: " << plane << "  " << plane * N_WIRE_MAX + i << endl;
                        //  cout << "Set parameter of wire: " << i << " tmin: " << tmin << " tmax: " << tmax
                        //       << " tsync: " << tsync << endl;
                    }
                }
                // end of short version

                // long version from AK
                /*
                if (!fIsSync)
                {
                    fHitPar->Clear();
                }

                cout << "Entering method1" << endl;

                TCanvas* cfit = new TCanvas("cfit", "fit", 50, 50, 500, 500);
                cfit->Clear();
                cfit->Divide(2, 2);

                Double_t sum;
                Double_t tmin = 0;
                Double_t tmax = 0;
                Double_t total;
                Int_t rebind_sum = 0;
                Int_t index_low = 0;
                Int_t index_high = 0;
                Int_t iwire = 0;
                Int_t rebinned_points = 0;
                int index_start;
                int distance = 15;
                for (Int_t plane = 0; plane < N_PLANE_MAX_PDC; plane++)
                {
                    cout << "Starting plane " << plane << endl;
                    long sumcheck = fh_time_drift[plane]->GetSum();

                    for (UInt_t i = 1; i <= N_WIRE_MAX; i++)
                    {

                        rebinned_points = 0;
                        tmin = 0;
                        tmax = 0;
                        rebind_sum = 0;
                        index_start = 0;
                        for (Int_t icount = 0; icount < xtc_points; icount++)
                        {
                            rebinned_x_mm[icount] = -1;
                            rebinned_t_ns[icount] = -10000;
                        }
                        //  cout<<"Starting wire "<<i<<endl;
                        iwire = i;

                        // find maximum value of drift times
                        TH1D* proj = fh_time_drift[plane]->ProjectionY("", i + 1, i + 1, 0);

                        //cfit->cd(1);
                        //  proj->Draw();

                        Int_t maxContent = proj->GetMaximum();
                        Int_t sumproj = proj->GetSum();
                        // make running sum of drift times
                        if (maxContent > 15.)
                        {
                            //   cout<<"For plane "<<plane<<", wire "<<i<<", maxContent "<<maxContent<<endl;

                            sum = 0;
                            for (Int_t icount = 0; icount < 4000; icount++)
                            {
                                t_ns[icount] = -10000;
                                x_mm[icount] = -1;
                            }

                            total = 0.;
                            sum = 0.;

                            Double_t yloghigh, yloglow, xhigh, xlow, slope, curveStart_ns = -1000;
                            Double_t tlow = 1000, thigh = -1000;

                            for (int j = 1; j < proj->GetNbinsX() - distance; j++)
                            {

                                yloghigh = log(proj->GetBinContent(j + distance));
                                yloglow = log(proj->GetBinContent(j));
                                xhigh = proj->GetBinCenter(j + distance);
                                xlow = proj->GetBinCenter(j);
                                slope = (yloghigh - yloglow) / (xhigh - xlow);
                                if (slope > 0. && proj->GetBinContent(j) > 0 && proj->GetBinContent(j + distance) > 0)
                                {
                                    index_start = j;
                                    tlow = xlow;
                                    curveStart_ns = tlow;
                                    break;
                                }
                            }
                            Double_t nbins = double(proj->GetNbinsX());
                            Double_t binwidth = proj->GetBinWidth(0);
                            Double_t firstbin = proj->GetBinCenter(0);
                            Double_t lastbin = firstbin + nbins * binwidth;
                            Int_t centrebin = int(nbins / 2.);
                            for (int j = proj->GetNbinsX(); j > centrebin + distance; j--)
                            {

                                yloghigh = log(proj->GetBinContent(j));
                                yloglow = log(proj->GetBinContent(j - distance));
                                xhigh = proj->GetBinCenter(j);
                                xlow = proj->GetBinCenter(j - distance);
                                slope = (yloghigh - yloglow) / (xhigh - xlow);
                                if (slope < 0. && proj->GetBinContent(j) > 0 && proj->GetBinContent(j - distance) > 0)
                                {
                                    index_high = j;
                                    thigh = xhigh;
                                    break;
                                }
                            }

                            for (UInt_t j = index_start; j < index_high + 1; j++) // proj->GetNbinsX(); j++)
                            {
                                total += proj->GetBinContent(j);
                            }

                            Double_t float_index = 0;
                            Int_t npoints = 0;
                            Double_t time_ns = -1000;
                            // for (UInt_t j = 1; j < proj->GetNbinsX(); j++)
                            for (UInt_t j = index_start; j < index_high + 1; j++)
                            {
                                time_ns = proj->GetBinCenter(j);
                                if (time_ns >= curveStart_ns)
                                {
                                    if (npoints == 0)
                                    {
                                        x_mm[npoints] = 0;
                                        t_ns[npoints] = time_ns; // curveStart_ns;
                                        npoints = 1;
                                    }
                                    else if (npoints > 0)
                                    {
                                        sum = sum + proj->GetBinContent(j);
                                        // add current position to xt curve
                                        x_mm[npoints] = sum / total * wire_distance;
                                        t_ns[npoints] = time_ns;
                                        npoints++;
                                    }
                                }
                            }

                            for (UInt_t j = 0; j < npoints; j++)
                            {
                                fh_running_sum->Fill(t_ns[j], x_mm[j]);
                            }

                            // +++++++++++++++++++++++++++++++++
                            //    Rebin (remove points)
                            // +++++++++++++++++++++++++++++++++
                            Double_t max_err = float(sumproj) / float(sumcheck); // 1./maxContent;// 0.1;
                            if (max_err < 1. / maxContent)
                                max_err = 1. / maxContent;

                            Double_t max_abs_err;
                            rebinned_points = 0;
                            // rebinned_x_mm[rebinned_points] = x_mm[0];
                            // rebinned_t_ns[rebinned_points] = t_ns[0];
                            UInt_t rebinned_counter = 0;

                            Double_t time_interval = thigh - tlow + 1.;
                            Int_t idelta = int(time_interval / double(xtc_points) + 1.);

                            cout << "*** " << plane << ", " << i << "; " << tlow << ", " << thigh << "; "
                                 << time_interval << ", " << idelta << ", " << npoints << endl;

                            if (double(npoints) / double(idelta) >= double(xtc_points))
                            {
                                cout << "WARNING more than 300 points! " << npoints << ", " << idelta << "; "
                                     << double(npoints) / double(idelta) << ", " << npoints / xtc_points << endl;
                                idelta = npoints / xtc_points;
                                cout << "new idelta = " << idelta << endl;
                            }

                            for (Int_t index = 1; index < npoints; index++)
                            {
                                //	cout<<"Index: "<<index_start<<", "<<index_high<<"; "<<index<<"; "<<index %
                                // idelta<<endl;
                                if (0 == index % idelta)
                                {
                                    rebinned_x_mm[rebinned_points] = x_mm[index];
                                    rebinned_t_ns[rebinned_points] = t_ns[index];
                                    rebinned_points++;
                                }
                            }

                            if (rebinned_points > 300)
                            {
                                rebinned_points = 300;
                            }
                            rebind_sum = rebinned_points;
                            for (Int_t irp = 0; irp < rebind_sum; irp++)
                            {
                                //  cout<<"Rebined: "<<rebind_sum<<", "<<irp<<", t_ns: "<<rebinned_t_ns[irp]<<"; x_mm:
                                //  "<<rebinned_x_mm[irp]<<endl;
                                fh_xtc->Fill(rebinned_t_ns[irp], rebinned_x_mm[irp]);
                            }
                            tmin = rebinned_t_ns[0];
                            //    tmax = rebinned_t_ns[rebind_sum - 1];
                            tmax = thigh;

                            Double_t err = 0;
                            Double_t err_check = 0;
                            Int_t isum = 0;

                            for (Int_t j = 1; j < npoints; j++)
                            {
                                Double_t temp_x = x_mm[j];
                                Double_t temp_t = t_ns[j];

                                for (UInt_t k = 0; k < rebind_sum - 1; k++)
                                {
                                    if (temp_t >= rebinned_t_ns[k] && temp_t < rebinned_t_ns[k + 1])
                                    {

                                        Double_t x1 = rebinned_x_mm[k];
                                        Double_t t1 = rebinned_t_ns[k];
                                        Double_t x2 = rebinned_x_mm[k + 1];
                                        Double_t t2 = rebinned_t_ns[k + 1];

                                        Double_t dx = x2 - x1;
                                        Double_t dt = t2 - t1;

                                        Double_t apar = dx / dt;
                                        Double_t bpar = x1 - apar * t1;

                                        Double_t xc = apar * temp_t + bpar;

                                        err = err + (xc - temp_x) * (xc - temp_x);
                                        isum++;
                                    }
                                }
                            }
                            err_check = sqrt(err / float(isum));
                            fh_error[plane]->Fill(i, err_check);
                            // cout<<"Err: "<<err_check<<", "<<err<<", "<<isum<<endl;

                            //cfit->cd(2);
                            //fh_running_sum->SetMarkerStyle(21);
                            //fh_running_sum->Draw("hist p ");
                            //cfit->cd(3);
                            //fh_xtc->SetMarkerStyle(23);
                            //fh_xtc->Draw("hist p");

                            //fh_running_sum->Clear();
                            //fh_xtc->Clear();
                        }
                        R3BPdcHitModulePar* par3 = fHitPar->GetModuleParAt(plane * N_WIRE_MAX + iwire);
                        par3->SetTmin(tmin);
                        par3->SetTmax(tmax);
                        par3->SetSync(tsync_mem[plane * N_WIRE_MAX + iwire - 1]);
                        par3->SetNPoints(rebind_sum);

                        for (Int_t icount = 0; icount < xtc_points; icount++)
                        {
                            xt_tarray[icount] = -10000;
                            xt_xarray[icount] = -1;
                        }

                        for (Int_t icount = 0; icount < rebind_sum; icount++)
                        {
                            xt_tarray[icount] = rebinned_t_ns[icount];
                            xt_xarray[icount] = rebinned_x_mm[icount];
                        }

                        par3->SetXT_xArray(xt_xarray, rebind_sum);
                        par3->SetXT_tArray(xt_tarray, rebind_sum);

                        cout << "Plane: " << plane + 1 << ", wire: " << i << ", channel: " << plane * N_WIRE_MAX + iwire
                             << ", tmin: " << tmin << ", tmax: " << tmax
                             << ", tsync: " << tsync_mem[plane * N_WIRE_MAX + iwire - 1]
                             << ", xtc points: " << rebind_sum << endl;
                    }
                    cfit->cd(plane + 1);
                    fh_error[plane]->SetMarkerStyle(22);
                    fh_error[plane]->Draw("hist p pmc");
                }
                // end of long version
                */
            }

            if (method2)
            {
                if (!fIsSync)
                {
                    // fHitPar->Clear();
                }

                TCanvas* cfit = new TCanvas("cfit", "fit", 50, 50, 500, 500);
                cfit->Clear();
                cfit->Divide(2, 2);

                Double_t sum;
                Double_t tmin = 0;
                Double_t tmax = 0;
                Int_t npoints = 0;
                Double_t total;
                Int_t index_low;
                Int_t index_high;
                Int_t iwire = 0;
                Int_t rebinned_points = 0;
                for (Int_t plane = 0; plane < N_PLANE_MAX_PDC; plane++)
                {

                    for (UInt_t i = 1; i <= N_WIRE_MAX; i++)
                    {

                        rebinned_points = 0;
                        tmin = 0;
                        tmax = 0;
                        npoints = 0;
                        for (Int_t icount = 0; icount < xtc_points; icount++)
                        {
                            t_ns[icount] = -10000;
                            x_mm[icount] = -1;
                        }

                        for (Int_t icount = 0; icount < xtc_points; icount++)
                        {
                            rebinned_x_mm[icount] = -1;
                            rebinned_t_ns[icount] = -10000;
                        }
                        //  cout<<"Starting wire "<<i<<endl;
                        iwire = i;

                        // find maximum value of drift times
                        TH1D* proj = fh_tsync[plane]->ProjectionY("", i + 1, i + 1, 0);

                        /*  cfit->cd(1);
                          proj->Draw();*/

                        uint maxContent = proj->GetMaximum();

                        // make running sum of drift times
                        if (maxContent > 2.)
                        {
                            //   cout<<"For plane "<<plane<<", wire "<<i<<", maxContent "<<maxContent<<endl;

                            total = 0.;
                            sum = 0.;
                            index_low = 0;
                            index_high = 0;

                            for (int j = 1; j < proj->GetNbinsX(); j++)
                            {
                                if (proj->GetBinContent(j) <= 0.1 * maxContent)
                                {
                                    index_low = j;
                                    // cout << "index low: " << j << ", " << index_low << ", " << proj->GetBinContent(j)
                                    //     << ", " << 0.1 * maxContent << ", " << maxContent << endl;
                                }
                                if (proj->GetBinContent(j) >= 0.9 * maxContent)
                                {
                                    index_high = j;
                                    if (plane == 1)
                                        cout << "index low: " << index_low << ", " << proj->GetBinContent(index_low)
                                             << ", " << 0.1 * maxContent << ", " << maxContent << endl;
                                    if (plane == 1)
                                        cout << "index high: " << index_high << ", " << proj->GetBinContent(j) << ". "
                                             << 0.9 * maxContent << "; " << maxContent << endl;

                                    break;
                                }
                            }


                            int content_low = proj->GetBinContent(index_low);
                            int content_high = proj->GetBinContent(index_high);
                            Double_t tlow = proj->GetBinCenter(index_low);
                            Double_t thigh = proj->GetBinCenter(index_high);
                            // imagine a line through the two points where rising edge reaches 10% and 90%
                            // from its maximum. We use the intersection of this line with the x-axis
                            // as the "start" of the curve.
                            Double_t curveStart_ns = (tlow == thigh) ? tlow
                                                                     : tlow - float(content_low) * (thigh - tlow) /
                                                                                  (content_high - content_low);

                            Double_t findex_start = (index_low == index_high)
                                                        ? index_low
                                                        : float(index_low) - float(content_low) *
                                                                                 (index_high - index_low) /
                                                                                 (content_high - content_low);

                            Int_t index_start = findex_start;
                            Int_t index_end = findex_start + 240. / proj->GetBinWidth(0);

                            for (UInt_t j = index_start; j < index_end; j++)
                            {
                                total += proj->GetBinContent(j);
                            }

                            Double_t float_index = 0;
                            npoints = 0;
                            Double_t time_ns = -1000;
                            for (UInt_t j = 1; j < index_end; j++)
                            {
                                time_ns = proj->GetBinCenter(j);
                                if (time_ns > curveStart_ns)
                                {
                                    if (npoints == 0)
                                    {
                                        x_mm[npoints] = 0;
                                        t_ns[npoints] = time_ns; // curveStart_ns;
                                        npoints++;
                                    }
                                    else
                                    {
                                        sum = sum + proj->GetBinContent(j);
                                        // add current position to xt curve
                                        x_mm[npoints] = sum / total * wire_distance;
                                        t_ns[npoints] = time_ns;
                                        // cout << "Test point: " << npoints << " sum: " << sum << " x: " <<
                                        // x_mm[npoints] << " t: " << t_ns[npoints] << endl;
                                        npoints++;
                                    }
                                }
                            }

                            for (UInt_t j = 0; j < npoints; j++)
                            {
                                fh_running_sum->Fill(t_ns[j], x_mm[j]);
                                fh_xtc->Fill(t_ns[j], x_mm[j]);
                            }

                            /*
                                                        // +++++++++++++++++++++++++++++++++
                                                        //    Rebin (remove points)
                                                        // +++++++++++++++++++++++++++++++++
                                                        Double_t max_err = 0.015;

                                                        Double_t max_abs_err;
                                                        rebinned_points = 0;
                                                        rebinned_x_mm[rebinned_points] = x_mm[0];
                                                        rebinned_t_ns[rebinned_points] = t_ns[0];
                                                        UInt_t rebinned_counter = 0;

                                                        for (Int_t index = 1; index < npoints; index++)
                                                        {
                                                            if (rebinned_counter + 2 > npoints)
                                                            {
                                                                break;
                                                            }
                                                            for (Int_t ij = rebinned_counter + 2; ij < npoints; ij++)
                                                            {

                                                                if (rebinned_counter > 600)
                                                                    max_err = 0.002;
                                                                Double_t temp_x = x_mm[ij];
                                                                Double_t temp_t = t_ns[ij];

                                                                Double_t start_x = rebinned_x_mm[rebinned_points - 1];
                                                                Double_t start_t = rebinned_t_ns[rebinned_points - 1];

                                                                Double_t dx = temp_x - start_x;
                                                                Double_t dt = temp_t - start_t;

                                                                Double_t apar = dx / dt;
                                                                Double_t bpar = start_x - apar * start_t;

                                                                Double_t tc;
                                                                Double_t xc;
                                                                Int_t isum = 0;
                                                                Double_t err = 0;
                                                                for (Int_t iz = rebinned_counter; iz < ij; iz++)
                                                                {
                                                                    tc = t_ns[iz];
                                                                    xc = apar * tc + bpar;

                                                                    err = err + (xc - x_mm[iz]) * (xc - x_mm[iz]);
                                                                    isum++;
                                                                }

                                                                Double_t err_check = sqrt(err / float(isum));

                                                                if (err_check > max_err)
                                                                {
                                                                    rebinned_points++;
                                                                    rebinned_x_mm[rebinned_points] = x_mm[ij];
                                                                    rebinned_t_ns[rebinned_points] = t_ns[ij];
                                                                    rebinned_counter = ij;

                                                                    break;
                                                                }
                                                            }

                                                            fh_xtc->Fill(rebinned_t_ns[rebinned_points],
                               rebinned_x_mm[rebinned_points]);
                                                        }

                                                        tmin = rebinned_t_ns[0];
                                                        tmax = rebinned_t_ns[rebinned_points];

                                                        Double_t err = 0;
                                                        Double_t err_check = 0;
                                                        Int_t isum = 0;

                                                        for (UInt_t j = 0; j < npoints; j++)
                                                        {
                                                            Double_t temp_x = x_mm[j];
                                                            Double_t temp_t = t_ns[j];

                                                            for (UInt_t k = 0; k < rebinned_points - 1; k++)
                                                            {
                                                                if (temp_t >= rebinned_t_ns[k] && temp_t <
                               rebinned_t_ns[k + 1])
                                                                {

                                                                    Double_t x1 = rebinned_x_mm[k];
                                                                    Double_t t1 = rebinned_t_ns[k];
                                                                    Double_t x2 = rebinned_x_mm[k + 1];
                                                                    Double_t t2 = rebinned_t_ns[k + 1];

                                                                    Double_t dx = x2 - x1;
                                                                    Double_t dt = t2 - t1;

                                                                    Double_t apar = dx / dt;
                                                                    Double_t bpar = x1 - apar * t1;

                                                                    Double_t xc = apar * temp_t + bpar;

                                                                    err = err + (xc - temp_x) * (xc - temp_x);
                                                                    isum++;
                                                                }
                                                            }
                                                        }
                                                        err_check = sqrt(err / float(isum));
                                                        fh_error[plane]->Fill(i, err_check);
                                                        // cout<<"Err: "<<err_check<<", "<<err<<", "<<isum<<endl;
                            */
                        }

                        tmin = -1000;
                        tmax = -1000;
                        if (npoints > 0)
                        {
                            tmin = t_ns[0];
                            tmax = t_ns[npoints - 1];
                        }
                        if (plane == 1 && npoints > 0)
                            cout << "reading par: " << plane << "  " << iwire << endl;
                        R3BPdcHitModulePar* par3 = fHitPar->GetModuleParAt(plane * N_WIRE_MAX + iwire);
                        par3->SetTmin(tmin);
                        par3->SetTmax(tmax);
                        if (plane == 1 && npoints > 0)
                            cout << "Tmin, Tmax: " << tmin << "  " << tmax << endl;
                        par3->SetSync(tsync_mem[plane * N_WIRE_MAX + iwire - 1]);
                        par3->SetNPoints(npoints);
                        par3->SetXT_xArray(x_mm, npoints);
                        par3->SetXT_tArray(t_ns, npoints);
                        /*
                                                for (Int_t icount = 0; icount < xtc_points; icount++)
                                                {
                                                    xt_tarray[icount] = -10000;
                                                    xt_xarray[icount] = -1;
                                                }

                                                for (Int_t icount = 0; icount < rebinned_points; icount++)
                                                {
                                                    xt_tarray[icount] = rebinned_t_ns[icount];
                                                    xt_xarray[icount] = rebinned_x_mm[icount];
                                                }

                                                par3->SetXT_xArray(xt_xarray, rebinned_points);
                                                par3->SetXT_tArray(xt_tarray, rebinned_points);

                        */

                        if (plane == 1 && npoints > 0)
                            cout << "Plane: " << plane + 1 << ", wire: " << i
                                 << ", channel: " << plane * N_WIRE_MAX + iwire << ", tmin: " << tmin
                                 << ", tmax: " << tmax << ", tsync: " << tsync_mem[plane * N_WIRE_MAX + iwire - 1]
                                 << ", xtc points: " << npoints << endl;
                    }
                    // cfit->cd(plane + 1);
                    // fh_error[plane]->SetMarkerStyle(22);
                    // fh_error[plane]->Draw("hist p pmc");
                }
            }
        }

        fHitPar->setChanged();
        for (Int_t j = 1; j < 2; j++)
        {
            for (UInt_t i = 1; i <= N_WIRE_MAX; i++)
            {
                R3BPdcHitModulePar* par4 = fHitPar->GetModuleParAt(j * N_WIRE_MAX + i);
                // cout<<"Testing par for wire: "<<i<<", tsync: "<< par4->GetSync()<<", tmin: "<<par4->GetTmin()<<",
                // tmax: "<<par4->GetTmax()<<endl;
                if (i == 25)
                {
                    Int_t np = par4->GetNPoints();
                    xtc_x = par4->GetXT_xArray();
                    xtc_t = par4->GetXT_tArray();
                    for (Int_t ic = 0; ic < par4->GetNPoints(); ic++)
                    {
                        Double_t xx = *(xtc_x + ic);
                        Double_t tt = *(xtc_t + ic);
                        cout << "x: " << xx << ", t: " << tt << endl;
                    }
                }
            }
        }

        fh_running_sum->Write();
        fh_xtc->Write();
        //        for (Int_t j = 0; j < N_PLANE_MAX_PDC; j++)
        //        {
        //            fh_error[j]->Write();
        //        }
    }
}

ClassImp(R3BPdcCal2Hit)
