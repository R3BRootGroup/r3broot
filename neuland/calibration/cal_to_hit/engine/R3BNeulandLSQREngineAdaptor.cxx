/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2023 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BNeulandLSQREngineAdaptor.h"
#include "R3BDataMonitor.h"
#include "R3BNeulandHitPar.h"
#include <FairRuntimeDb.h>
#include <R3BLogger.h>
#include <R3BNeulandParDirCreator.h>
#include <R3BShared.h>
#include <algorithm>

namespace
{
    void fill_module_parameters(const std::vector<R3BNeulandHitModulePar>& old_module_pars,
                                R3B::Neuland::Cal2HitPar& new_module_pars)
    {
        for (const auto& par : old_module_pars)
        {
            auto new_par = R3B::Neuland::HitModulePar{};
            // GetModuleId return 1 based module Num
            new_par.module_num = par.GetModuleId();
            new_par.tDiff.value = par.GetTDiff();
            new_par.tSync.value = par.GetTSync();
            new_par.effectiveSpeed.value = par.GetEffectiveSpeed();
            new_par.lightAttenuationLength.value = par.GetLightAttenuationLength();
            new_par.lightAttenuationFactor.value =
                std::exp(R3B::Neuland::TotalBarLength / new_par.lightAttenuationLength.value / 2);

            // Getters accept 1 as the left side and 2 as the right side

            const auto left_index = 1 + toIndex(R3B::Side::left);
            const auto right_index = 1 + toIndex(R3B::Side::right);

            new_par.pedestal.left() = par.GetPedestal(left_index);
            new_par.pedestal.right() = par.GetPedestal(right_index);

            new_par.energyGain.left().value = par.GetEnergyGain(left_index);
            new_par.energyGain.right().value = par.GetEnergyGain(right_index);

            new_par.PMTSaturation.left().value = par.GetPMTSaturation(left_index);
            new_par.PMTSaturation.right().value = par.GetPMTSaturation(right_index);

            new_par.PMTThreshold.left().value = par.GetPMTThreshold(left_index);
            new_par.PMTThreshold.right().value = par.GetPMTThreshold(right_index);
            new_module_pars.AddModulePar(std::move(new_par));
        }
    }
} // namespace

namespace R3B::Neuland::Calibration
{
    void LSQREngineAdaptor::Init()
    {
        const auto module_size = GetModuleSize();
        if (module_size == 0)
        {
            throw R3B::runtime_error("module size is 0 during the initialization!");
        }

        // To set the number of planes in legacy hitCalibrationEngine
        auto hit_par_temp = R3BNeulandHitPar{};
        auto number_of_planes = static_cast<int>(module_size) / BarsPerPlane;
        hit_par_temp.SetNumberOfPlanes(number_of_planes);
        cosmic_tracker_ = std::make_unique<CosmicTracker>(number_of_planes);
        hit_cal_engine_.Init(&hit_par_temp);
    }

    void LSQREngineAdaptor::AddSignals(const CalData& signals)
    {
        for (const auto& [plane_num, plane_cal_data] : signals)
        {
            for (const auto& [bar_num, bar_signal] : plane_cal_data.bar_cal_data)
            {
                add_bar_signal(bar_signal, Side::left);
                add_bar_signal(bar_signal, Side::right);
            }
        }
    }

    void LSQREngineAdaptor::add_bar_signal(const BarCalData& barSignal, Side side)
    {
        // change to 0 based indexing
        const auto module_id = static_cast<int>(barSignal.module_num - 1);

        const auto& calData = (side == Side::left) ? barSignal.left : barSignal.right;

        for (const auto& calSignal : calData)
        {
            // 1 is the right side and 0 is the left side
            hit_cal_engine_.Set(module_id,
                                static_cast<int>(toIndex(side)),
                                calSignal.leading_time.value - calSignal.trigger_time.value,
                                static_cast<int>(calSignal.time_over_threshold.value));
            if (hit_cal_engine_.IsValid(module_id))
            {
                cosmic_tracker_->AddPoint(module_id, hit_cal_engine_.GetPosition(module_id));
            }
        }
    }

    void LSQREngineAdaptor::EndOfEvent(unsigned int event_num)
    {
        const auto added_points_num = cosmic_tracker_->GetBarIDs().size();
        const auto& track = cosmic_tracker_->GetTrack();
        if (track.Interactions.empty() or added_points_num < 6)
        {
            R3BLOG(debug, "Failed to get the cosmic track!");
        }
        else
        {
            hit_cal_engine_.Add(track, event_num);
        }
    }

    void LSQREngineAdaptor::Calibrate(Cal2HitPar& hit_par)
    {
        hit_par.Reset();
        // auto* new_dir = ParDirCreator{}.mkdir(DEFAULT_HIST_MONITOR_DIR);
        auto* new_dir = ParDirCreator{}.mkdir("Histograms");
        fill_module_parameters(hit_cal_engine_.Calibrate(new_dir), hit_par);
    }

} // namespace R3B::Neuland::Calibration
