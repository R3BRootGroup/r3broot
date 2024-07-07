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

#pragma once

#include "R3BShared.h"
#include <R3BNeulandParSet.h>
#include <unordered_map>

namespace R3B::Neuland
{
    struct HitModulePar
    {
        unsigned int module_num = 0;               // 1-based
        ValueError<double> tDiff;                  // ns
        ValueError<double> tSync;                  // ns
        ValueError<double> effectiveSpeed;         // cm/ns
        ValueError<double> offset_effective_c;         // cm
        ValueError<double> lightAttenuationLength; // cm, 1/alpha
        ValueError<double> lightAttenuationFactor; // exp(alpha*L/2)
        LRPair<int> pedestal;                      //
        LRPair<ValueError<double>> energyGain;     //
        LRPair<ValueError<double>> PMTSaturation;  //
        LRPair<ValueError<double>> PMTThreshold;   //
        ClassDefNV(HitModulePar, 1);
    };

    class Cal2HitPar : public ParSet
    {
      public:
        explicit Cal2HitPar(std::string_view name = "NeulandHitPar",
                            std::string_view title = "Neuland Cal2Hit calibration",
                            std::string_view context = "TestDefaultContext",
                            Bool_t own = kTRUE);

        Cal2HitPar(const Cal2HitPar&) = default;
        Cal2HitPar(Cal2HitPar&&) = default;
        auto operator=(const Cal2HitPar&) -> Cal2HitPar& = default;
        auto operator=(Cal2HitPar&&) -> Cal2HitPar& = default;
        ~Cal2HitPar() override = default;

        void CalculateEnergyCutoff();
        void Reset() { clear(); }

        // setters:
        void SetDistanceToTarget(double distance) { distance_to_target_ = distance; }
        void SetEnergyCutoff(double cutoff) { energy_cut_ = cutoff; }
        void SetGlobalTimeOffset(double offset) { global_time_offset_ = offset; }
        void AddModulePar(HitModulePar module_par)
        {
            const auto mNum = module_par.module_num;
            module_pars_.insert_or_assign(mNum, std::move(module_par));
        }

        // getter:
        auto GetDistanceToTarget() const { return distance_to_target_; }
        auto GetEnergyCutoff() const { return energy_cut_; }
        auto GetGlobalTimeOffset() const { return global_time_offset_; }
        [[deprecated("Use GetDistancesToFirstPlane instead")]] auto GetDistanceToFirstPlane(
            unsigned int plane_num) const
        {
            return distances_to_first_plane_.at(plane_num);
        }
        auto GetDistancesToFirstPlane() const -> const auto& { return distances_to_first_plane_; }
        auto GetNumModulePar() const { return module_pars_.size(); }
        auto GetModuleParAt(unsigned int module_num) const -> const HitModulePar&;
        auto GetModulePars() const -> const auto& { return module_pars_; }
        auto GetListOfModulePar() const -> const std::unordered_map<unsigned int, R3B::Neuland::HitModulePar>&
        {
            return module_pars_;
        }
        auto GetListOfModuleParRef() -> auto& { return module_pars_; }

      private:
        double global_time_offset_ = 0.; // in ns
        double distance_to_target_ = 0.; // in cm
        double energy_cut_ = 0.;         // in MeV
        std::vector<double> distances_to_first_plane_;
        std::unordered_map<unsigned int, R3B::Neuland::HitModulePar> module_pars_;
        void clear() override
        {
            global_time_offset_ = 0.;
            distance_to_target_ = 0.;
            energy_cut_ = 0.;
            distances_to_first_plane_.clear();
            module_pars_.clear();
        }

      public:
        ClassDefOverride(Cal2HitPar, 1);
    };

} // namespace R3B::Neuland

using R3BNeulandHitPar2 = R3B::Neuland::Cal2HitPar;
