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

#include "R3BNeulandMapToCalTask.h"
#include "R3BNeulandCalData2.h"
#include "R3BNeulandCommon.h"
#include <FairRootManager.h>
#include <FairRuntimeDb.h>
#include <R3BEventHeader.h>
#include <R3BException.h>
#include <R3BLogger.h>
#include <TH1D.h>
#include <TH2D.h>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include <range/v3/all.hpp>

namespace
{
    // NOLINTBEGIN
    auto GetWalkCorrection(double value) -> double
    {
        value = value * 17.2278 / 162.464;
        const auto veronika3 = -4.29359 + 17.3841 / sqrt(value) + 0.073181; // veronika 3
        const auto veronika2 = 1.59667 + 78.5274 * pow(value, -1.97051) - 4.00192 / value - 0.125473 * value +
                               0.00130958 * value * value; // veronika 2
        const auto result = (value > 25.0000) ? veronika3 : (veronika3 + veronika2) / 2.;
        return -1. * result;
    }
    // NOLINTEND

} // namespace

namespace R3B::Neuland
{

    Map2CalTask::Map2CalTask()
        : Map2CalTask{ "NeulandMapToCalTask", 1 }
    {
    }

    Map2CalTask::Map2CalTask(std::string_view name, int iVerbose)
        : CalibrationTask(name.data(), iVerbose)
    {
    }

    void Map2CalTask::HistogramInit(DataMonitor& histograms)
    {
        const auto total_bar_num = plane_num_ * BarsPerPlane;
        constexpr auto BINSIZE_TOT = 1000;
        constexpr auto BINSIZE_LT = 1000;
        constexpr auto BINSIZE_TIME = 8000;
        constexpr auto HIST_MAX_TIME = 6000.;
        constexpr auto HIST_MIN_TIME = -5000.;
        constexpr auto HIST_MAX_TOT = 600.;
        constexpr auto HIST_MAX_LT = 2000.;

        histograms.add_hist<TH1D>("module_num", "Counts with module ids", total_bar_num, -0.5, total_bar_num + 0.5);
        histograms.add_hist<TH1D>("ToT", "Time over threshold", BINSIZE_TOT * 2, 0., BINSIZE_TOT + 0.);
        histograms.add_hist<TH1D>("LeadingTime", "Leading time", BINSIZE_LT, 0., HIST_MAX_LT);
        histograms.add_hist<TH1D>("TriggerTime", "Trigger time", BINSIZE_LT, 0., HIST_MAX_LT);
        histograms.add_hist<TH1D>("Leading_minus_trigger", "Time", BINSIZE_TIME, HIST_MIN_TIME, HIST_MAX_TIME);
        histograms.add_hist<TH1I>(
            "Bar_hitNum_l", "Number of hits on the left PMT per bar", total_bar_num, 0.5, 0.5 + total_bar_num);
        histograms.add_hist<TH1I>(
            "Bar_hitNum_r", "Number of hits on the right PMT per bar", total_bar_num, 0.5, 0.5 + total_bar_num);

        histograms.add_hist<TH2D>("TimeLVsBar",
                                  "Time_vs_bar_left",
                                  total_bar_num,
                                  0.5,
                                  0.5 + total_bar_num,
                                  BINSIZE_TIME,
                                  HIST_MIN_TIME,
                                  HIST_MAX_TIME);
        histograms.add_hist<TH2D>("TimeRVsBar",
                                  "Time_vs_bar_right",
                                  total_bar_num,
                                  0.5,
                                  0.5 + total_bar_num,
                                  BINSIZE_TIME,
                                  HIST_MIN_TIME,
                                  HIST_MAX_TIME);

        histograms.add_hist<TH2D>(
            "ToTLVsBar", "ToT_vs_bar_left", total_bar_num, 0.5, 0.5 + total_bar_num, BINSIZE_TOT, 0., HIST_MAX_TOT);
        histograms.add_hist<TH2D>(
            "ToTRVsBar", "ToT_vs_bar_right", total_bar_num, 0.5, 0.5 + total_bar_num, BINSIZE_TOT, 0., HIST_MAX_TOT);
    }

    void Map2CalTask::ExtraInit(FairRootManager* /*rootMan*/)
    {
        mappedData_.init();
        if (not base_par_->GetTrigIDMap().empty())
        {
            R3BLOG(debug, "NeuLAND data has trig time!");
            trigMappedData_.init();
        }
        else
        {
            R3BLOG(debug, "NeuLAND data has no trig time!");
        }
        calData_.init(not IsHistDisabled());
        set_pmt_num();
        set_ct_freq();
    }

    void Map2CalTask::set_pmt_num()
    {
        if (plane_num_ = base_par_->GetNumOfPlanes(); plane_num_ == 0)
        {
            R3BLOG(warn, "plane number obtained from the calibration parameters is 0!");
        }
        else
        {
            total_pmt_nums_ = 2 * plane_num_ * BarsPerPlane;
            R3BLOG(info, fmt::format("total number of PMTs set to be {}", total_pmt_nums_).c_str());
        }

        if (total_pmt_nums_ == 0)
        {
            throw R3B::runtime_error("The number of PMTs cannot be zero!");
        }
    }

    void Map2CalTask::set_ct_freq()
    {
        if (auto ct_freq = calibrationPar_->GetSlowClockFrequency(); ct_freq == 0)
        {
            R3BLOG(warn, "Coarse time frequency obtained from parameteers is 0! Use default value.");
            coarse_time_frequency_ = COARSE_TIME_CLOCK_FREQUENCY_MHZ;
        }
        else
        {
            coarse_time_frequency_ = ct_freq;
            max_coarse_time_ = 1000. * coarse_time_max_num_ / coarse_time_frequency_;
            R3BLOG(info, fmt::format("Coarse time frequency set to be {} MHz", ct_freq).c_str());
        }
    }

    auto Map2CalTask::CheckConditions() const -> bool
    {
        auto signal_size = mappedData_.size();
        R3BLOG(debug2,
               fmt::format("Minimal signal size: {}. Current signal size: {}. Number of PMTs: {}",
                           signal_min_size_,
                           signal_size,
                           total_pmt_nums_));
        if (signal_size < signal_min_size_)
        {
            R3BLOG(debug2,
                   fmt::format("condition of the minimal size is not met with current paddle signal size. Skip the "
                               "current event."));
            return false;
        }
        return is_pulse_mode_ ? signal_size > total_pmt_nums_ : signal_size < total_pmt_nums_ / 2;
    }

    void Map2CalTask::TriggeredExec() { calibrate(); }

    auto Map2CalTask::convert_to_real_time(R3BTCalPar2* calPar,
                                           SingleEdgeSignal signal,
                                           FTType ftType,
                                           unsigned int module_num) const -> ValueError<double>
    {
        const auto& modulePar = calPar->GetParamAt(module_num);
        const auto fineTime = modulePar.GetFineTime(ftType, signal.fine);

        // signal time in ns
        const auto coarseTime = 1000. / coarse_time_frequency_ * static_cast<float>(signal.coarse + 1);
        return { coarseTime - fineTime.value, fineTime.error };
    }

    auto Map2CalTask::get_trigger_time(unsigned int module_num, Side side) const -> ValueError<double>
    {

        const auto& triggerMap = base_par_->GetTrigIDMap();
        if (triggerMap.empty())
        {
            return { 0., 0. };
        }

        auto triggerIDPair = triggerMap.find(module_num);
        if (triggerIDPair == triggerMap.end())
        {
            const auto eventNum = GetEventHeader()->GetEventno();
            R3BLOG(error,
                   fmt::format("Can't find the trigger ID for the module ID {} in the event {}", module_num, eventNum));
            return { -1., 0 };
        }
        const auto triggerID = (side == Side::left) ? triggerIDPair->second.first : triggerIDPair->second.second;
        const auto trigData = std::find_if(trigMappedData_.begin(),
                                           trigMappedData_.end(),
                                           [&triggerID](const auto& ele) { return ele.first == triggerID; });
        if (trigData == trigMappedData_.end())
        {
            const auto eventNum = GetEventHeader()->GetEventno();
            R3BLOG(error,
                   fmt::format("No such trigger ID {} in mappedTrigData for moduleNum {} in the event {}!",
                               triggerID,
                               module_num,
                               eventNum));
            R3BLOG(error,
                   fmt::format("Available trigIDs: {}", fmt::join(trigMappedData_.get() | ranges::views::keys, ", ")));
            return { -1., 0 };
        }
        return convert_to_real_time(calibrationTrigPar_, trigData->second.signal, FTType::trigger, trigData->first);
    }

    auto Map2CalTask::get_tot(DoubleEdgeSignal pmtSignal,
                              unsigned int module_num,
                              Side side) const -> ValueError<double>
    {
        const auto leadFType = (side == Side::left) ? FTType::leftleading : FTType::rightleading;
        const auto trailFType = (side == Side::left) ? FTType::lefttrailing : FTType::righttrailing;
        const auto leadingT = convert_to_real_time(calibrationPar_, pmtSignal.leading, leadFType, module_num);
        const auto trailingT = convert_to_real_time(calibrationPar_, pmtSignal.trailing, trailFType, module_num);
        const auto time_over_thres = trailingT - leadingT;
        R3BLOG(debug3, fmt::format("leading :{} trailing: {}", leadingT.value, trailingT.value));
        return (time_over_thres.value > 0) ? time_over_thres : time_over_thres + max_coarse_time_;
    }

    void Map2CalTask::overflow_correct(R3B::Neuland::CalDataSignal& calSignal) const
    {
        // trigger time should come after signal time
        if (calSignal.leading_time.value > calSignal.trigger_time.value)
        {
            calSignal.leading_time -= max_coarse_time_;
        }
    }

    auto Map2CalTask::doubleEdgeSignal_to_calSignal(const DoubleEdgeSignal& double_edge_signal,
                                                    R3B::Side side,
                                                    unsigned int module_num) const -> CalDataSignal
    {
        auto calDataSignal = CalDataSignal{};
        const auto ftType = (side == R3B::Side::left) ? FTType::leftleading : FTType::rightleading;
        calDataSignal.time_over_threshold = get_tot(double_edge_signal, module_num, side);
        const auto walk_correction =
            (is_walk_enabled_) ? GetWalkCorrection(calDataSignal.time_over_threshold.value) : 0.;
        calDataSignal.leading_time =
            convert_to_real_time(calibrationPar_, double_edge_signal.leading, ftType, module_num) + walk_correction;
        calDataSignal.trigger_time = get_trigger_time(module_num, side);
        overflow_correct(calDataSignal);
        R3BLOG(debug, fmt::format("Adding a new cal signal: {}", calDataSignal));
        return calDataSignal;
    }

    auto Map2CalTask::mapBarSignal_to_calSignals(const MapBarSignal& map_bar_signals,
                                                 unsigned int module_num,
                                                 R3B::Side side) const -> std::vector<CalDataSignal>
    {
        const auto& signals = (side == Side::left) ? map_bar_signals.left : map_bar_signals.right;
        auto calSignals = std::vector<CalDataSignal>{};
        calSignals.reserve(signals.size());
        std::transform(signals.begin(),
                       signals.end(),
                       std::back_inserter(calSignals),
                       [module_num, side, this](const auto& dESignal)
                       { return doubleEdgeSignal_to_calSignal(dESignal, side, module_num); });
        return calSignals;
    }

    void Map2CalTask::fill_cal_data(BarCalData& cal, const MapBarSignal& signals)
    {
        cal.left = mapBarSignal_to_calSignals(signals, cal.module_num, Side::left);
        cal.right = mapBarSignal_to_calSignals(signals, cal.module_num, Side::right);

        if (not IsHistDisabled())
        {
            histogram_monitor(cal, Side::left);
            histogram_monitor(cal, Side::right);
        }
    }

    void Map2CalTask::calibrate()
    {
        R3BLOG(debug2, fmt::format("mapped Data size: {}", mappedData_.size()));
        auto& cal_data = calData_.get();
        for (const auto& planeSignals : mappedData_)
        {
            const auto planeNum = static_cast<int>(planeSignals.plane_num);
            auto [plane_cal, is_plane_success] = cal_data.try_emplace(planeNum, PlaneCalData{ planeNum });
            for (const auto& [barNum, barSignals] : planeSignals.bars)
            {
                if (not IsHistDisabled())
                {
                    GetHistMonitor().get("module_num")->Fill(barNum);
                }
                const auto module_num = static_cast<int>(Neuland_PlaneBar2ModuleNum(planeNum, barNum));
                auto [bar_cal, is_bar_success] =
                    plane_cal->second.bar_cal_data.try_emplace(module_num, BarCalData{ module_num });
                fill_cal_data(bar_cal->second, barSignals);
            }
        }
    }

    void Map2CalTask::histogram_monitor(const BarCalData& cal, Side side)
    {
        auto& histograms = GetHistMonitor();
        const auto module_num = cal.module_num;
        const auto& calSignals = (side == Side::left) ? cal.left : cal.right;

        for (const auto& signal : calSignals)
        {
            const auto lTime = signal.leading_time.value;
            const auto time_over_thres = signal.time_over_threshold.value;
            const auto triggerTime = signal.trigger_time.value;
            histograms.get("ToT")->Fill(time_over_thres);
            histograms.get("LeadingTime")->Fill(lTime);
            histograms.get("TriggerTime")->Fill(triggerTime);
            histograms.get("Leading_minus_trigger")->Fill(lTime - triggerTime);
            if (side == Side::left)
            {
                histograms.get("TimeLVsBar")->Fill(module_num, lTime - triggerTime);
                histograms.get("ToTLVsBar")->Fill(module_num, time_over_thres);
                histograms.get("Bar_hitNum_l")->Fill(module_num);
            }
            else
            {
                histograms.get("TimeRVsBar")->Fill(module_num, lTime - triggerTime);
                histograms.get("ToTRVsBar")->Fill(module_num, time_over_thres);
                histograms.get("Bar_hitNum_r")->Fill(module_num);
            }
        }
    }

    void Map2CalTask::FinishEvent() {}
} // namespace R3B::Neuland
