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

#include "R3BDigitizingTamex.h"

#include <FairRuntimeDb.h>
#include <cmath>

#include "R3BException.h"
#include "R3BNeulandHitModulePar.h"
#include "R3BNeulandHitPar.h"
#include <FairRunAna.h>
#include <R3BLogger.h>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <utility>

namespace R3B::Digitizing::Neuland::Tamex
{
    // some declarations for static functions:
    namespace
    {
        template <class T>
        auto CheckOverlapping(const T& peak, std::vector<T>& peaks) -> decltype(peaks.begin());
        template <class T>
        void ReOverlapping(typename std::vector<T>::iterator v_iter, std::vector<T>& peaks);
        template <class T>
        void RemovePeakAt(typename std::vector<T>::iterator v_iter, std::vector<T>& peaks);
    } // namespace

    // global variables for default options:
    const size_t TmxPeaksInitialCapacity = 10;
    R3BNeulandHitPar* Channel::neuland_hit_par_ = nullptr; // NOLINT

    Params::Params(TRandom3& rnd)
        : fRnd{ &rnd }
    {
    }

    Params::Params(const Params& other)
    {
        if (other.fRnd == nullptr)
        {
            throw std::runtime_error(
                "R3BDigitizingTamex: copy constructor of Params cannot take nullptr of random generator!");
        }
        *this = other;
    }

    PMTPeak::PMTPeak(Digitizing::Channel::Hit pmtHit, const Channel& channel)
        : time_(pmtHit.time)
    {
        auto par = channel.GetParConstRef();
        // apply saturation coefficent
        qdc_ = pmtHit.light / (1. + par.fSaturationCoefficient * pmtHit.light);
    };

    auto PMTPeak::operator+=(const PMTPeak& other) -> PMTPeak&
    {
        qdc_ += other.qdc_;
        time_ = (time_ < other.time_) ? time_ : other.time_;
        return *this;
    }

    FQTPeak::FQTPeak(const PMTPeak& pmtPeak, Channel* channel)
        : qdc_(pmtPeak.GetQDC())
        , leading_edge_time_(pmtPeak.GetLETime())
        , channel_ptr_(channel)
    {
        if (channel_ptr_ == nullptr)
        {
            LOG(fatal) << "channel is not bound to FQTPeak object!";
        }
        const auto& par = channel->GetParConstRef();

        // calculate the time and the width of the signal
        width_ = QdcToWidth(qdc_, par);
        trailing_edge_time_ = leading_edge_time_ + width_;
    }

    auto FQTPeak::operator==(const FQTPeak& other) const -> bool
    {
        if (other.leading_edge_time_ == 0 && leading_edge_time_ == 0)
        {
            LOG(warn) << "the times of both PMT signals are 0!";
        }
        return (leading_edge_time_ <= (other.leading_edge_time_ + other.width_)) &&
               (other.leading_edge_time_ <= (leading_edge_time_ + width_));
    }

    void FQTPeak::operator+=(const FQTPeak& other)
    {
        if (channel_ptr_ == nullptr)
        {
            throw R3B::logic_error("channel is not bound to FQTPeak object!");
        }
        leading_edge_time_ =
            (leading_edge_time_ < other.leading_edge_time_) ? leading_edge_time_ : other.leading_edge_time_;
        trailing_edge_time_ =
            (trailing_edge_time_ > other.trailing_edge_time_) ? trailing_edge_time_ : other.trailing_edge_time_;
        width_ = trailing_edge_time_ - leading_edge_time_;
        qdc_ = WidthToQdc(width_, channel_ptr_->GetParConstRef());
    }

    Channel::Channel(ChannelSide side, PeakPileUpStrategy strategy, const Params& par)
        : Digitizing::Channel{ side }
        , pileup_strategy_{ strategy }
        , par_{ par }
    {
        pmt_peaks_.reserve(TmxPeaksInitialCapacity);
    }

    Channel::Channel(ChannelSide side, PeakPileUpStrategy strategy, TRandom3& rnd)
        : Channel{ side, strategy, Params{ rnd } }
    {
    }

    void Channel::GetHitPar(const std::string& hitParName)
    {
        if (hitParName.empty())
        {
            LOG(info) << "DigitizingTamex: Using default parameter for Tamex Channels.";
            return;
        }
        auto* run = FairRunAna::Instance();
        auto* rtdb = run->GetRuntimeDb();
        neuland_hit_par_ = dynamic_cast<R3BNeulandHitPar*>(rtdb->findContainer(hitParName.c_str()));
        if (neuland_hit_par_ != nullptr)
        {
            LOG(info) << "DigitizingTamex: HitPar " << hitParName
                      << " has been found in the root file. Using calibration values in rootfile.";
        }
        else
        {
            LOG(info) << "DigitizingTamex: HitPar " << hitParName << " cannot be found. Using default values.";
            neuland_hit_par_ = nullptr;
        }
    }

    void Channel::AttachToPaddle(Digitizing::Paddle* /*v_paddle*/)
    {
        auto* paddle = GetPaddle();
        if (paddle == nullptr)
        {
            return;
        }

        if (CheckPaddleIDInHitPar())
        {
            neuland_hit_module_par_ = neuland_hit_par_->GetModuleParAt(GetPaddle()->GetPaddleID() - 1);
            if (CheckPaddleIDInHitModulePar())
            {
                auto side = GetSide();
                par_.fSaturationCoefficient = neuland_hit_module_par_->GetPMTSaturation(static_cast<int>(side));
                par_.fEnergyGain = neuland_hit_module_par_->GetEnergyGain(static_cast<int>(side));
                par_.fPedestal = neuland_hit_module_par_->GetPedestal(static_cast<int>(side));
                par_.fPMTThresh = neuland_hit_module_par_->GetPMTThreshold(static_cast<int>(side));
                par_.fQdcMin = 1 / par_.fEnergyGain;
            }
        }
    }

    auto Channel::CheckPaddleIDInHitModulePar() const -> bool
    {
        auto is_valid = false;
        if (neuland_hit_module_par_ == nullptr || GetPaddle() == nullptr)
        {
            return false;
        }

        if (GetPaddle()->GetPaddleID() != neuland_hit_module_par_->GetModuleId())
        {
            LOG(warn) << "Channel::SetHitModulePar:Wrong paddleID for the parameters!";
            is_valid = false;
        }
        else
        {
            is_valid = true;
        }
        return is_valid;
    }

    auto Channel::CheckPaddleIDInHitPar() const -> bool
    {
        auto is_valid = false;
        if (neuland_hit_par_ == nullptr)
        {
            return false;
        }
        if (not neuland_hit_par_->hasChanged())
        {
            R3BLOG(warn, "Can't setup parameter in the root file correctly!.");
            return false;
        }

        auto PaddleId_max = neuland_hit_par_->GetNumModulePar();
        if (GetPaddle()->GetPaddleID() > PaddleId_max)
        {
            LOG(warn) << "Paddle id " << GetPaddle()->GetPaddleID() << " exceeds the id " << PaddleId_max
                      << " in the parameter file!";
            is_valid = false;
        }
        else
        {
            is_valid = true;
        }

        return is_valid;
    }

    void Channel::AddHit(Hit newHit)
    {
        if (newHit.time < par_.fTimeMin || newHit.time > par_.fTimeMax)
        {
            return;
        }
        InvalidateSignals();
        InvalidateTrigTime();
        pmt_peaks_.emplace_back(newHit, *this);
    }

    auto Channel::CreateSignal(const FQTPeak& peak) const -> Signal
    {
        auto peakQdc = peak.GetQDC();
        auto peakTime = peak.GetLETime();
        auto qdc = ToQdc(peakQdc);

        auto signal = Signal{};
        signal.qdcUnSat = ToUnSatQdc(qdc);
        signal.qdc = qdc;
        signal.tdc = ToTdc(peakTime);
        signal.side = this->GetSide();
        LOG(debug) << "R3BDigitizingTamex: Create a signal with qdc " << signal.qdc << " and tdc " << signal.tdc
                   << std::endl;
        return signal;
    }

    template <typename Peak>
    void Channel::PeakPileUp(/* inout */ std::vector<Peak>& peaks)
    {
        if (peaks.size() <= 1)
        {
            return;
        }

        std::sort(peaks.begin(), peaks.end(), std::less{});
        for (auto front_peak = peaks.begin(); front_peak != peaks.end(); ++front_peak)
        {
            auto end_peak = std::remove_if(front_peak + 1,
                                           peaks.end(),
                                           [&front_peak](auto& peak)
                                           {
                                               if (*front_peak == peak)
                                               {
                                                   (*front_peak) += peak;
                                                   return true;
                                               }
                                               return false;
                                           });
            peaks.erase(end_peak, peaks.end());
        }
    }

    void Channel::PeakPileUpWithDistance(/* inout */ std::vector<FQTPeak>& peaks, double distance)
    {
        if (peaks.empty())
        {
            return;
        }
        std::sort(peaks.begin(), peaks.end(), std::less{});

        for (auto front_peak = peaks.begin(); front_peak != peaks.end(); ++front_peak)
        {
            auto last_leading_time = front_peak->GetLETime();
            auto end_peak = std::remove_if(front_peak + 1,
                                           peaks.end(),
                                           [&distance, &front_peak, &last_leading_time](FQTPeak& peak)
                                           {
                                               if ((peak - last_leading_time) < distance)
                                               {
                                                   front_peak->AddQDC(peak.GetQDC());
                                                   last_leading_time = peak.GetLETime();
                                                   return true;
                                               }
                                               return false;
                                           });
            peaks.erase(end_peak, peaks.end());
        }
    }

    void Channel::PeakPileUpInTimeWindow(/* inout */ std::vector<FQTPeak>& peaks, double time_window)
    {
        if (peaks.empty())
        {
            return;
        }
        std::sort(peaks.begin(), peaks.end(), std::less{});

        auto& front_peak = peaks.front();
        std::for_each(peaks.begin() + 1,
                      peaks.end(),
                      [&front_peak, time_window](FQTPeak& peak)
                      {
                          if ((peak - front_peak) < time_window)
                          {
                              front_peak.AddQDC(peak.GetQDC());
                          }
                      });
        peaks.erase(peaks.begin() + 1, peaks.end());
    }

    void Channel::FQTPeakPileUp(/* inout */ std::vector<FQTPeak>& peaks)
    {
        switch (pileup_strategy_)
        {
            case PeakPileUpStrategy::width:
                PeakPileUp(peaks);
                break;
            case PeakPileUpStrategy::distance:
                PeakPileUpWithDistance(peaks, par_.fPileUpDistance);
                break;
            case PeakPileUpStrategy::time_window:
                PeakPileUpInTimeWindow(peaks, par_.fPileUpTimeWindow);
                break;
            default:
                break;
        }
    }

    template <typename Peak>
    void Channel::ApplyThreshold(std::vector<Peak>& peaks)
    {
        // apply threshold on energy using c++ erase-remove idiom:
        auto it_end =
            std::remove_if(peaks.begin(),
                           peaks.end(),
                           [this](const auto& peak) { return peak.GetQDC() < this->GetParConstRef().fPMTThresh; });
        peaks.erase(it_end, peaks.end());
    }

    auto Channel::ConstructFQTPeaks(std::vector<PMTPeak>& pmtPeaks) -> std::vector<FQTPeak>
    {
        auto FQTPeaks = std::vector<FQTPeak>{};
        FQTPeaks.reserve(pmtPeaks.size());

        // sorting pmt peaks according to time:
        std::sort(pmtPeaks.begin(), pmtPeaks.end());

        PeakPileUp(pmtPeaks);
        ApplyThreshold(pmtPeaks);
        for (auto const& peak : pmtPeaks)
        {
            FQTPeaks.emplace_back(peak, this);
        }
        return FQTPeaks;
    }

    auto Channel::ConstructSignals() -> Signals
    {
        fqt_peaks_ = ConstructFQTPeaks(pmt_peaks_);
        // signal pileup:
        FQTPeakPileUp(fqt_peaks_);

        // construct Channel signals:
        auto signals = std::vector<Signal>{};
        signals.reserve(fqt_peaks_.size());

        for (const auto& peak : fqt_peaks_)
        {
            signals.emplace_back(CreateSignal(peak));
        }
        return signals;
    }

    auto Channel::GetFQTPeaks() -> const std::vector<FQTPeak>&
    {

        if (!Is_ValidSignals())
        {
            ConstructSignals();
        }
        return fqt_peaks_;
    }

    auto Channel::GetPMTPeaks() -> const std::vector<PMTPeak>&
    {
        if (!Is_ValidSignals())
        {
            ConstructSignals();
        }
        return pmt_peaks_;
    }

    auto Channel::ToQdc(double qdc) const -> double
    {
        // apply energy smearing
        qdc = par_.fRnd->Gaus(qdc, par_.fEResRel * qdc);
        return qdc;
    }

    auto Channel::ToTdc(double time) const -> double { return time + par_.fRnd->Gaus(0., par_.fTimeRes); }

    auto Channel::ToUnSatQdc(double qdc) const -> double
    {
        // Apply reverse saturation
        if (par_.fExperimentalDataIsCorrectedForSaturation)
        {
            qdc = qdc / (1 - par_.fSaturationCoefficient * qdc);
        }
        // Apply reverse attenuation
        return qdc;
    }
} // namespace R3B::Digitizing::Neuland::Tamex
