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
// -----                 R3BTofdCal2Hit                   -----
// -----            Created June 2020 by M.Heil            -----
// ------------------------------------------------------------

#include "R3BPdcCal2Hit.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BPdcCalData.h"
#include "R3BPdcHitData.h"
#include "R3BTCalEngine.h"
#include <FairRootManager.h>
// #include "R3BPdcHitModulePar.h"
// #include "R3BPdcHitPar.h"

#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"

#include "TClonesArray.h"
#include "TMath.h"

#include <iostream>
#include <sstream>
#include <string>

#include "mapping_pdc_trig.hh"

extern "C"
{
// #include "ext_data_client.h"
#include "ext_h101_pdc.h"
}
#define LENGTH(x) (sizeof x / sizeof *x)

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

R3BPdcCal2Hit::R3BPdcCal2Hit()
    : FairTask("PdcCal2Hit", 1)
    , fCalItems(NULL)
    , fCalTriggerItems(NULL)
    , fHitItems(new TClonesArray("R3BPdcHitData"))
    , fNofHitItems(0)
    , fNofHitPars(0)
    //    , fHitPar(NULL)
    , fnEvents(0)
    , fClockFreq(1. / CTDC_16_CLOCK_MHZ * 1000.)
{
}

R3BPdcCal2Hit::R3BPdcCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems(NULL)
    , fCalTriggerItems(NULL)
    , fHitItems(new TClonesArray("R3BPdcHitData"))
    , fNofHitItems(0)
    , fNofHitPars(0)
    //    , fHitPar(NULL)
    , fnEvents(0)
    , fClockFreq(1. / CTDC_16_CLOCK_MHZ * 1000.)
{
}

R3BPdcCal2Hit::~R3BPdcCal2Hit()
{
    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

InitStatus R3BPdcCal2Hit::Init()
{
    pdc_trig_map_setup();

    /*
        fHitPar = (R3BPdcHitPar*)FairRuntimeDb::instance()->getContainer("PdcHitPar");
        if (!fHitPar)
        {
            LOG(error) << "Could not get access to PdcHitPar-Container.";
            fNofHitPars = 0;
            return kFATAL;
        }

        fNofHitPars = fHitPar->GetNumModulePar();
        if (fNofHitPars == 0)
        {
            LOG(error) << "There are no Hit parameters in container PdcHitPar";
            return kFATAL;
        }
    */

    auto plane_num = LENGTH(EXT_STR_h101_PDC_onion::PDC_P);
    fPlaneArray.resize(plane_num);
    auto ch_num = LENGTH(EXT_STR_h101_PDC_onion::PDC_P[0].TLCMI);
    for (auto it = fPlaneArray.begin(); fPlaneArray.end() != it; ++it)
    {
        it->resize(ch_num);
    }

    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (header == nullptr)
    {
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
        R3BLOG(warn, "R3BEventHeader was found instead of EventHeader.");
    }

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("PdcCal"));
    if (NULL == fCalItems)
        LOG(fatal) << "Branch PdcCal not found";
    fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("PdcTriggerCal"));
    if (NULL == fCalTriggerItems)
        LOG(fatal) << "Branch PdcTriggerCal not found";
    maxevent = mgr->CheckMaxEventNo();

    // request storage of Hit data in output tree
    mgr->Register("PdcHit", "Land", fHitItems, kTRUE);

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BPdcCal2Hit::SetParContainers()
{
    /*
        fHitPar = (R3BPdcHitPar*)FairRuntimeDb::instance()->getContainer("PdcHitPar");
        if (!fHitPar)
        {
            LOG(error) << "Could not get access to PdcHitPar-Container.";
            fNofHitPars = 0;
            return;
        }
    */
}

InitStatus R3BPdcCal2Hit::ReInit()
{
    SetParContainers();
    return kSUCCESS;
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
    std::vector<R3BPdcCalData const*> trig_table(LENGTH(EXT_STR_h101_PDC_onion::PDC_P));
    for (size_t j = 0; j < trig_num; ++j)
    {
        auto cal = dynamic_cast<R3BPdcCalData const*>(fCalTriggerItems->At(j));
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
        auto cur_cal = dynamic_cast<R3BPdcCalData const*>(fCalItems->At(j));
        if (cur_cal->GetEdgeId() == 1)
        {
            ++n_lead;
            auto plane_i = cur_cal->GetPlaneId() - 1;
            auto wire_i = cur_cal->GetWireId() - 1;
            auto& channel = fPlaneArray.at(plane_i).at(wire_i);
            channel.lead_list.push_back(cur_cal);
            // cout << "leading wire:" << ch_i <<  " time: " << cur_cal->GetTime_ns() << endl;
        }
        else
        {
            ++n_trail;
        }
    }
    if (n_lead != n_trail)
    {
        cout << "Number of leading edges not equal to number of trailing edges!" << endl;
        //    return;
    }

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = dynamic_cast<R3BPdcCalData const*>(fCalItems->At(j));
        if (cur_cal->GetEdgeId() == 2)
        {
            Double_t c_period = 4096. * (1000. / fClockFreq);

            auto plane_i = cur_cal->GetPlaneId() - 1;
            auto wire_i = cur_cal->GetWireId() - 1;
            auto& channel = fPlaneArray.at(plane_i).at(wire_i);
            // cout << "trailing wire:" << ch_i <<  " time: " << cur_cal->GetTime_ns() << endl;
            if (channel.lead_list.empty())
            {
                cout << "lead list empty for plane, wire " << plane_i + 1 << "  " << wire_i + 1 << endl;
                continue;
            }
            auto lead = channel.lead_list.front();

            Double_t cur_cal_trig_ns = 0;
            Double_t lead_trig_ns = 0;

            assert(plane_i == lead->GetPlaneId() - 1);
            assert(wire_i == lead->GetWireId() - 1);

            // Trigger time is the same for leading and trailing,
            // but we still need to subtract it for when we want the
            // edge times.

            auto cur_cal_trig_i = g_pdc_trig_map[plane_i][wire_i];
            auto cur_cal_trig = trig_table.at(cur_cal_trig_i);
            cur_cal_trig_ns = cur_cal_trig->GetTime_ns();

            auto cur_cal_ns =
                fmod(cur_cal->GetTime_ns() - cur_cal_trig_ns + c_period + c_period / 2, c_period) - c_period / 2;
            auto lead_ns =
                fmod(lead->GetTime_ns() - cur_cal_trig_ns + c_period + c_period / 2, c_period) - c_period / 2;
            auto tot_ns = fmod(cur_cal_ns - lead_ns + c_period + c_period / 2, c_period) - c_period / 2;

            /*
                        cout << plane_i << ' ' << wire_i << '\n';
                        cout << "Lead: " << lead_ns << " = " << lead->GetTime_ns() << " - " << cur_cal_trig_ns << endl;
                        cout << "Trail: " << cur_cal_ns << " = " << cur_cal->GetTime_ns() << " - " << cur_cal_trig_ns <<
               endl; cout << "ToT: " << tot_ns << " = " << cur_cal_ns << " - " << lead_ns << endl;
            */
            if (tot_ns > 0 && tot_ns < 1000)
            {
                channel.tot_list.push_back(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));
                channel.lead_list.pop_front();
            }
        }
    }
    // cout << "***************** new event ****************" << endl;
    for (auto it = fPlaneArray.begin(); fPlaneArray.end() != it; ++it)
    {
        auto& ch_arr = *it;
        for (auto it2 = ch_arr.begin(); ch_arr.end() != it2; ++it2)
        {
            auto const& channel = *it2;
            for (auto it_pdc_tot = channel.tot_list.begin(); channel.tot_list.end() != it_pdc_tot;
                 ++it_pdc_tot) // over ihit(fiber)
            {
                auto const& pdc_tot = *it_pdc_tot;
                auto plane = pdc_tot.lead->GetPlaneId();
                auto wire = pdc_tot.lead->GetWireId();
                auto tot_pdc = pdc_tot.tot_ns;
                Double_t t_pdc = pdc_tot.lead_ns;

                // cout << "plane: " << plane <<  " time: " << t_pdc << " wire: " << wire << " ToT: " << tot_pdc <<
                // endl;

                Double_t x = 0.;
                Double_t y = 0.;
                Int_t ID = 0;

                if (plane & 1)
                {
                    x = wire;
                    y = 0;
                }
                else
                {
                    x = 0;
                    y = wire;
                }
                ID = plane;

                // cout << "Hit level ID: " << ID << " x: " << x << " y: " << y << " ToT: " << tot_pdc << " t: " <<
                // t_pdc << endl;
                new ((*fHitItems)[fNofHitItems++]) R3BPdcHitData(t_pdc, x, y, tot_pdc, ID, wire);
            }
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

void R3BPdcCal2Hit::FinishTask() {}

ClassImp(R3BPdcCal2Hit);
