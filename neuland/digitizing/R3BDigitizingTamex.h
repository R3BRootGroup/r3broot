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

#ifndef R3BROOT_DIGITIZINGTAMEX_H
#define R3BROOT_DIGITIZINGTAMEX_H
/**
 * NeuLAND Tamex digitizing module with multiple hits per event
 * @author Yanzhao Wang
 * @since July 15, 2022
 *
 *
 */

#include "R3BDigitizingChannel.h"
#include "R3BDigitizingPaddle.h"
#include "TRandom3.h"
#include "Validated.h"
#include <optional>

class R3BNeulandHitPar;
class R3BNeulandHitModulePar;

namespace R3B::Digitizing::Neuland::Tamex
{
    class Channel;
    struct Params
    {
        // NOLINTBEGIN
        double fPMTThresh = 1.;                // [MeV]
        double fSaturationCoefficient = 0.012; // Saturation coefficient of PMTs
        bool fExperimentalDataIsCorrectedForSaturation = true;
        double fTimeRes = 0.15; // time + Gaus(0., fTimeRes) [ns]
        double fEResRel = 0.05; // Gaus(e, fEResRel * e) []
        double fEnergyGain = 15.0;
        double fPedestal = 14.0;
        double fTimeMax = 1000.;          // ns
        double fTimeMin = 1.;             // ns
        double fPileUpTimeWindow = 1000.; // ns
        double fPileUpDistance = 100.;    // ns
        double fQdcMin = 0.067;
        TRandom3* fRnd = nullptr;
        // NOLINTEND

        explicit Params(TRandom3&);

        // rule of 5
        Params(Params&&) = delete;
        auto operator=(const Params&) -> Params& = default;
        auto operator=(Params&&) -> Params& = delete;
        ~Params() = default;
        Params(const Params& other);
    };

    enum class PeakPileUpStrategy
    {
        width,
        distance,
        time_window
    };

    class PMTPeak
    {
      public:
        PMTPeak() = default;
        PMTPeak(Digitizing::Channel::Hit pmtHit, const Channel&);
        auto operator<(const PMTPeak& rhs) const -> bool { return (time_ < rhs.time_); }
        auto operator==(const PMTPeak& rhs) const -> bool { return std::abs(time_ - rhs.time_) < peakWidth; }
        auto operator+=(const PMTPeak& other) -> PMTPeak&;
        [[nodiscard]] auto GetQDC() const -> double { return qdc_; }
        [[nodiscard]] auto GetLETime() const -> double { return time_; }
        static constexpr double peakWidth = 15.0; // ns

      private:
        double qdc_ = 0.0;
        double time_ = 0.0;
    };

    class FQTPeak
    {
      public:
        FQTPeak(const PMTPeak& pmtPeak, Channel* channel);
        FQTPeak() = default;

        // Getters:
        [[nodiscard]] auto GetWidth() const -> double { return width_; }
        [[nodiscard]] auto GetQDC() const -> double { return qdc_; }
        [[nodiscard]] auto GetLETime() const -> double { return leading_edge_time_; }
        [[nodiscard]] auto GetTETime() const -> double { return trailing_edge_time_; }

        auto operator==(const FQTPeak& other) const -> bool;
        void operator+=(const FQTPeak& other);
        auto operator-(const FQTPeak& other) const -> double { return leading_edge_time_ - other.leading_edge_time_; }
        auto operator-(double time) const -> double { return leading_edge_time_ - time; };
        auto operator>(const FQTPeak& other) const -> bool { return leading_edge_time_ - other.leading_edge_time_ > 0; }
        auto operator<(const FQTPeak& other) const -> bool { return leading_edge_time_ - other.leading_edge_time_ < 0; }

        void AddQDC(double qdc) { qdc_ += qdc; }

        template <typename Par>
        static auto WidthToQdc(double width, const Par& par) -> double
        {
            return std::max(1., width - par.fPedestal) / par.fEnergyGain;
        }

        template <typename Par>
        static auto QdcToWidth(double qdc, const Par& par) -> double
        {
            auto width = 0.0;
            if (qdc > par.fQdcMin)
            {
                width = qdc * par.fEnergyGain + par.fPedestal;
            }
            else
            {
                width = qdc * par.fEnergyGain * (par.fPedestal + 1);
            }
            return width;
        }
        explicit operator Digitizing::Channel::Signal() const;

      private:
        double width_ = 0.0;              // the temperal width of the TmxPeak in [ns]
        double qdc_ = 0.0;                // the qdc value in [MeV] (without threshold)
        double leading_edge_time_ = 0.0;  // leading edge of the TmxPeak in [ns]
        double trailing_edge_time_ = 0.0; // tailing edge of the TmxPeak
        Channel* channel_ptr_ = nullptr;
    };

    class Channel : public Digitizing::Channel
    {
      public:
        Channel(ChannelSide, PeakPileUpStrategy strategy, TRandom3&);
        Channel(ChannelSide, PeakPileUpStrategy strategy, const Params&);
        explicit Channel(ChannelSide side, PeakPileUpStrategy strategy = PeakPileUpStrategy::width)
            : Channel(side, strategy, GetDefaultRandomGen())
        {
        }
        // Setters:
        void SetPileUpStrategy(PeakPileUpStrategy strategy) { pileup_strategy_ = strategy; }

        // Getters:
        auto GetPar() -> Tamex::Params& { return par_; }
        auto GetParConstRef() const -> const Tamex::Params& { return par_; }
        auto GetFQTPeaks() -> const std::vector<FQTPeak>&;
        auto GetPMTPeaks() -> const std::vector<PMTPeak>&;

        void AddHit(Hit /*hit*/) override;
        auto CreateSignal(const FQTPeak& peak) const -> Signal;
        static void GetHitPar(const std::string& hitParName);

      private:
        PeakPileUpStrategy pileup_strategy_ = PeakPileUpStrategy::width;
        std::vector<PMTPeak> pmt_peaks_;
        std::vector<FQTPeak> fqt_peaks_;
        static R3BNeulandHitPar* neuland_hit_par_; // NOLINT
        R3BNeulandHitModulePar* neuland_hit_module_par_ = nullptr;
        Tamex::Params par_;

        // private virtual functions
        auto ConstructSignals() -> Signals override;
        void AttachToPaddle(Digitizing::Paddle* paddle) override;

        // private non-virtual functions
        auto CheckPaddleIDInHitPar() const -> bool;
        auto CheckPaddleIDInHitModulePar() const -> bool;
        void SetHitModulePar(int PaddleId);
        auto ToQdc(double) const -> double;
        auto ToTdc(double) const -> double;
        auto ToUnSatQdc(double) const -> double;
        template <typename Peak>
        void ApplyThreshold(/* inout */ std::vector<Peak>&);
        auto ConstructFQTPeaks(std::vector<PMTPeak>& pmtPeaks) -> std::vector<FQTPeak>;
        template <typename Peak>
        static void PeakPileUp(/* inout */ std::vector<Peak>& peaks);

        static void PeakPileUpWithDistance(/* inout */ std::vector<FQTPeak>& peaks, double distance);
        static void PeakPileUpInTimeWindow(/* inout */ std::vector<FQTPeak>& peaks, double time_window);
        void FQTPeakPileUp(/* inout */ std::vector<FQTPeak>& peaks);
    };

} // namespace R3B::Digitizing::Neuland::Tamex
#endif // R3BROOT_DIGITIZINGTAMEX_H
