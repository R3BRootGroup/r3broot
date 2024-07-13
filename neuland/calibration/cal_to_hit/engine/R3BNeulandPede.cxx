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

#include "R3BNeulandPede.h"
#include "Utilities.h"
#include <R3BLogger.h>

constexpr auto DEFAULT_RES_FILENAME = "millepede.res";
constexpr auto PEDE_STEER_FILENAME = "neuland_steer.txt";
constexpr auto PARAMETER_FILENAME = "neuland_pars.txt";

namespace R3B::Neuland::Calibration
{
    namespace
    {
        void calculate_time_offset(R3B::Neuland::Cal2HitPar& cal_to_hit_par)
        {
            auto& module_pars = cal_to_hit_par.GetListOfModuleParRef();
            for (auto& [module_num, module_par] : module_pars)
            {
                if (module_par.effectiveSpeed.value != 0)
                {
                    module_par.tDiff = module_par.offset_effective_c / module_par.effectiveSpeed;
                }
            }
        }
    } // namespace

    void PedeHandler::init()
    {
        par_result_.set_filename(DEFAULT_RES_FILENAME);
        pede_launcher_.set_steer_filename(PEDE_STEER_FILENAME);
        pede_launcher_.set_parameter_filename(PARAMETER_FILENAME);
    }

    void PedeHandler::calibrate()
    {
        R3BLOG(info, "Launching pede algorithm...");
        pede_launcher_.sync_launch();
        pede_launcher_.end();

        par_result_.read();
    }

    void PedeHandler::fill_module_parameters(Neuland::Cal2HitPar& cal_to_hit_par)
    {
        const auto& pars = par_result_.get_pars();
        for (const auto& [par_id, par] : pars)
        {
            const auto [module_num, global_label] = to_module_num_label(par_id, max_number_of_modules_);
            auto& module_pars = cal_to_hit_par.GetListOfModuleParRef();

            auto& par_ref = module_pars.emplace(module_num, HitModulePar{}).first->second;
            switch (global_label)
            {
                case GlobalLabel::tsync:
                    par_ref.tSync += ValueErrorD{ par.value, par.error } * scale_factor_;
                    break;
                case GlobalLabel::offset_effective_c:
                    // The value here is the product of tDiff and effectiveSped. Real tDiff will be calculated later
                    par_ref.offset_effective_c += (ValueErrorD{ par.value, par.error } * scale_factor_);
                    break;
                case GlobalLabel::effective_c:
                    par_ref.effectiveSpeed += ValueErrorD{ par.value, par.error };
                    break;
                default:
                    throw std::runtime_error("An error occured with unrecognized global tag");
            }
        }
        calculate_time_offset(cal_to_hit_par);
    }

    void PedeHandler::init_steer_writer(const Cal2HitPar& /*cal_to_hit_par*/)
    {
        auto steer_writer = SteerWriter{};
        steer_writer.set_filepath(PEDE_STEER_FILENAME);
        steer_writer.set_parameter_file(PARAMETER_FILENAME);
        steer_writer.set_data_filepath(data_filename_);
        steer_writer.add_method(SteerWriter::Method::inversion,
                                std::make_pair(static_cast<float>(pede_interartion_number_), pede_error_threshold_));
        steer_writer.add_other_options(std::vector<std::string>{ "outlierdownweighting", "4" });
        // steer_writer.add_other_options(
        //     std::vector<std::string>{ "scaleerrors", fmt::format("{}", error_scale_factor_) });

        // const auto& module_pars = cal_to_hit_par.GetModulePars();
        // for (const auto& [module_num, module_par] : module_pars)
        // {
        //     steer_writer.add_parameter_default(
        //         to_global_label_id(static_cast<int>(module_num), GlobalLabel::effective_c),
        //         std::make_pair(static_cast<float>(module_par.effectiveSpeed.value), 0.F));
        //     steer_writer.add_parameter_default(
        //         to_global_label_id(static_cast<int>(module_num), GlobalLabel::offset_effective_c),
        //         std::make_pair(
        //             static_cast<float>(module_par.effectiveSpeed.value * module_par.tDiff.value / SCALE_FACTOR),
        //             0.F));
        // }
        steer_writer.write();
    }
} // namespace R3B::Neuland::Calibration
