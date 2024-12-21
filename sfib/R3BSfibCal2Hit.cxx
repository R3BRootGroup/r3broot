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

#include "R3BSfibCal2Hit.h"
#include "FairLogger.h"
#include "R3BSfibCalData.h"
#include "R3BSfibHitData.h"
#include <FairRootManager.h>
#include <TClonesArray.h>
#include <iostream>

R3BSfibCal2Hit::ToT::ToT(R3BSfibCalData const* a_lead,
                         R3BSfibCalData const* a_trail,
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

R3BSfibCal2Hit::R3BSfibCal2Hit(Int_t a_verbose, enum R3BTCalEngine::CTDCVariant a_ctdc_variant)
    : FairTask("R3BSfibCal2Hit", a_verbose)
    , fnEvents()
    , fClockFreq(R3BTCalEngine::CTDC_16_BWD_150 == a_ctdc_variant ? 150 : 250)
    , fCalItems()
    , fTopCalTriggerItems()
    , fBotCalTriggerItems()
    , fHitItems(new TClonesArray("R3BSfibHitData"))
    , fTopTriggerMap()
    , fBotTriggerMap()
    , fNofHitItems()
    , fChannelArray()
{
}

R3BSfibCal2Hit::~R3BSfibCal2Hit() { delete fHitItems; }

InitStatus R3BSfibCal2Hit::Init()
{
    cout << "Init::R3BSfibCal2Hit" << endl;
    auto mgr = FairRootManager::Instance();
    if (!mgr)
    {
        LOG(error) << "FairRootManager not found.";
        return kERROR;
    }

    auto name = "SfibCal";
    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject(name));
    if (!fCalItems)
    {
        LOG(fatal) << "Branch " << name << " not found.";
        return kERROR;
    }

    mgr->Register("SfibHit", "Land", fHitItems, kTRUE);
    for (auto side_i = 0; side_i < 2; ++side_i)
    {
        fChannelArray[side_i].resize(256);
    }

    return kSUCCESS;
}

InitStatus R3BSfibCal2Hit::ReInit() { return kSUCCESS; }

void R3BSfibCal2Hit::SetParContainers() {}

void R3BSfibCal2Hit::Exec(Option_t* option)
{
    if (fnEvents / 100000. == (int)fnEvents / 100000)
        std::cout << "\rEvents: " << fnEvents << std::flush;

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

    size_t cal_num = fCalItems->GetEntriesFast();

    std::vector<R3BSfibCalData const*> top_trig_table(256 / 128);
    std::vector<R3BSfibCalData const*> bot_trig_table(256 / 128);
#if 0
    // Make direct mappingo tables for trigger items.
    size_t top_trig_num = fTopCalTriggerItems->GetEntries();
    for (size_t j = 0; j < top_trig_num; ++j)
    {
        auto cal = (R3BSfibCalData const*)fTopCalTriggerItems->At(j);
        top_trig_table.at(cal->GetChannel() - 1) = cal;
    }

    // TODO: This will create a map for every fiber detector... Urg.
    // Also since it's shared between many detectors it must be dynamic, for now.
    size_t bot_trig_num = fBotCalTriggerItems->GetEntries();
    for (size_t j = 0; j < bot_trig_num; ++j)
    {
        auto cal = (R3BSfibCalData const*)fBotCalTriggerItems->At(j);
        bot_trig_table.at(cal->GetChannel() - 1) = cal;
    }
#endif

    // Find multi-hit ToT for every channel.
    // The easiest safe way to survive ugly cases is to record all
    // leading edges per channel, and then pair up with whatever
    // trailing we have.
    // Not super efficient, but shouldn't crash if the data is not
    // perfect.
    unsigned n_lead = 0;
    unsigned n_trail = 0;
    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = dynamic_cast<R3BSfibCalData const*>(fCalItems->At(j));
        if (cur_cal->IsLeading())
        {
            ++n_lead;
            auto side_i = cur_cal->IsTop() ? 1 : 0;
            auto ch_i = cur_cal->GetChannel() - 1;
            auto& channel = fChannelArray[side_i].at(ch_i);
            channel.lead_list.push_back(cur_cal);
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

    for (size_t j = 0; j < cal_num; ++j)
    {
        auto cur_cal = dynamic_cast<R3BSfibCalData const*>(fCalItems->At(j));
        if (cur_cal->IsTrailing())
        {
            auto side_i = cur_cal->IsTop() ? 1 : 0;
            Double_t c_period = 0 == side_i ? 4096. * (1000. / fClockFreq) : 2048. * (1000. / 200.);

            auto ch_i = cur_cal->GetChannel() - 1;
            auto& channel = fChannelArray[side_i].at(ch_i);
            if (channel.lead_list.empty())
            {
                continue;
            }
            auto lead = channel.lead_list.front();

            Double_t cur_cal_trig_ns = 0;
            Double_t lead_trig_ns = 0;
            if (cur_cal->IsTop() && fTopTriggerMap)
            {
                auto cur_cal_trig_i = fTopTriggerMap[ch_i];
                auto lead_trig_i = fTopTriggerMap[lead->GetChannel() - 1];
                if (cur_cal_trig_i < top_trig_table.size() && lead_trig_i < top_trig_table.size())
                {
                    auto cur_cal_trig = top_trig_table.at(cur_cal_trig_i);
                    auto lead_trig = top_trig_table.at(lead_trig_i);
                    cur_cal_trig_ns = cur_cal_trig->GetTime_ns();
                    lead_trig_ns = lead_trig->GetTime_ns();
                }
            }
            else if (!cur_cal->IsTop() && fBotTriggerMap)
            {
                auto cur_cal_trig_i = fBotTriggerMap[ch_i];
                auto lead_trig_i = fBotTriggerMap[lead->GetChannel() - 1];
                if (cur_cal_trig_i < bot_trig_table.size() && lead_trig_i < bot_trig_table.size())
                {
                    auto cur_cal_trig = bot_trig_table.at(cur_cal_trig_i);
                    auto lead_trig = bot_trig_table.at(lead_trig_i);
                    cur_cal_trig_ns = cur_cal_trig->GetTime_ns();
                    lead_trig_ns = lead_trig->GetTime_ns();
                }
            }

            auto cur_cal_ns =
                fmod(cur_cal->GetTime_ns() - cur_cal_trig_ns + c_period + c_period / 2, c_period) - c_period / 2;
            auto lead_ns = fmod(lead->GetTime_ns() - lead_trig_ns + c_period + c_period / 2, c_period) - c_period / 2;
            auto tot_ns = fmod(cur_cal_ns - lead_ns + c_period + c_period / 2, c_period) - c_period / 2;

            if (tot_ns > 0 && tot_ns < 1000)
            {
                if (side_i == 1)
                    channel.tot_list.push_back(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));
                if (side_i == 0)
                    channel.tot_list.push_front(ToT(lead, cur_cal, lead_ns, cur_cal_ns, tot_ns));

                channel.lead_list.pop_front();
            }
        }
    }

    // Make every permutation to create fibers.
    auto const& bot_array = fChannelArray[0];
    auto const& top_array = fChannelArray[1];

    for (auto it_top = top_array.begin(); top_array.end() != it_top; ++it_top)
    {
        auto const& top = *it_top;
        for (auto it_top_tot = top.tot_list.begin(); top.tot_list.end() != it_top_tot; ++it_top_tot)
        {
            auto const& top_tot = *it_top_tot;
            auto top_i = top_tot.lead->GetChannel() - 1;
            for (auto it_bot = bot_array.begin(); bot_array.end() != it_bot; ++it_bot)
            {
                auto const& bot = *it_bot;
                for (auto it_bot_tot = bot.tot_list.begin(); bot.tot_list.end() != it_bot_tot; ++it_bot_tot)
                {
                    auto const& bot_tot = *it_bot_tot;
                    auto bot_i = bot_tot.lead->GetChannel() - 1;

                    // Check that the combo is inside one block of 8x8 sorting.
                    auto top_sub_i = top_i / 64;
                    auto bot_sub_i = bot_i / 64;
                    if (top_sub_i != bot_sub_i)
                        continue;

                    auto fiber_id = (bot_i & 7) + (top_i * 8);

                    auto tot_top = top_tot.tot_ns;
                    Double_t t_top = top_tot.lead_ns;

                    Double_t t_top1 = top_tot.lead->GetTime_ns();

                    new ((*fHitItems)[fNofHitItems++])
                        R3BSfibHitData(fiber_id, top_tot.lead_ns, bot_tot.lead_ns, top_tot.tot_ns, bot_tot.tot_ns);
                }
            }
        }
    }
    fnEvents++;
}

void R3BSfibCal2Hit::FinishEvent()
{
    fHitItems->Clear();
    fNofHitItems = 0;
}

void R3BSfibCal2Hit::FinishTask() {}

ClassImp(R3BSfibCal2Hit)
