#pragma once

#include "R3BTrigIDMappingFinder.h"
#include <FairParSet.h>
#include <R3BNeulandParSet.h>
#include <R3BShared.h>
#include <cmath>
#include <unordered_map>

namespace R3B
{
    struct FTChannel2TimeRelation
    {
        double hist_overflow{};
        std::vector<ValueError<double>> value_error;
        ClassDefNV(FTChannel2TimeRelation, 1);
    };

    enum class FTType
    {
        leftleading,
        rightleading,
        lefttrailing,
        righttrailing,
        trigger,
    };

    class TCalVFTXModulePar
    {
      public:
        using ValueErrors = std::vector<ValueError<double>>;
        TCalVFTXModulePar() = default;
        explicit TCalVFTXModulePar(unsigned int module_num)
            : module_num_{ module_num }
        {
        }

        // Setters:
        void Set_correlation(std::map<FTType, FTChannel2TimeRelation> correlation)
        {
            finetime_correlation_ = std::move(correlation);
        }

        // Getters:
        [[nodiscard]] auto GetFineTime(FTType ftType, uint channel_Num) const -> ValueError<double>;
        [[nodiscard]] auto GetModuleNum() const { return module_num_; }

      private:
        unsigned int module_num_ = 0;
        std::map<FTType, FTChannel2TimeRelation> finetime_correlation_;

      public:
        ClassDefNV(TCalVFTXModulePar, 1);
    };

    class Map2CalPar : public Neuland::ParSet
    {
      public:
        using TrigIDMap = Neuland::TrigIDMappingFinder::TrigIDMap;
        explicit Map2CalPar(std::string_view name = "NeulandMap2CalPar",
                            std::string_view title = "Neuland Map2Cal calibration",
                            std::string_view context = "TestDefaultContext",
                            bool own = true);
        // ~R3BTCalPar2() { LOG(debug) << " R3BTCalPar2 deleted!"; }
        void AddModuleParam(unsigned int module_num, TCalVFTXModulePar modulePar);

        // Setters:
        void SetTrigEnabled(bool is_enabled) { is_trig_enabled_ = is_enabled; }
        void SetSlowClockFrequency(float frq) { slow_clock_frequency_ = frq; }
        void SetDetParName(std::string_view detParName) { fairDetParName_ = detParName; }
        void SetModuleParam(unsigned int module_num, TCalVFTXModulePar par)
        {
            moduleParams_.insert_or_assign(module_num, std::move(par));
        }

        void SetParam(std::unordered_map<unsigned int, TCalVFTXModulePar> moduleParams)
        {
            moduleParams_ = std::move(moduleParams);
        }

        // Getters:
        [[nodiscard]] auto GetSlowClockFrequency() const -> const auto& { return slow_clock_frequency_; }
        [[nodiscard]] auto GetParams() const -> const auto& { return moduleParams_; }
        [[nodiscard]] auto GetParamAt(unsigned int module_num) const -> const TCalVFTXModulePar&;
        auto HasTrigEnabled() const -> bool { return is_trig_enabled_; }

      private:
        float slow_clock_frequency_ = 0; // GHz
        bool is_trig_enabled_ = true;
        std::unordered_map<unsigned int, TCalVFTXModulePar> moduleParams_;
        std::string fairDetParName_ = "FairGenericParIo";

        void clear() override { moduleParams_.clear(); }

      public:
        ClassDefOverride(Map2CalPar, 1);
    };
} // namespace R3B

using R3BTCalPar2 = R3B::Map2CalPar;
