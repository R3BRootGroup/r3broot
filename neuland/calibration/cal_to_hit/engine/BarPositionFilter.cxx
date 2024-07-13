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

#include "BarPositionFilter.h"
#include <Math/WrappedMultiTF1.h>
#include <R3BException.h>
#include <R3BLogger.h>

namespace R3B::Neuland::Calibration
{
    BarPositionFitter::BarPositionFitter()
    {
        static constexpr auto init_values = std::array{ 1., 0. };
        fitter_.Config().SetParamsSettings(static_cast<int>(init_values.size()), init_values.data());
        fitter_.Config().SetUpdateAfterFit(false);
    }
    auto BarPositionFitter::fit_positions() -> bool
    {
        if (not fit_status_.has_value())
        {
            fit_status_ = fit_with(horizontal_bars_, y_pos_fitting_function_) and
                          fit_with(vertical_bars_, x_pos_fitting_function_);
        }
        return fit_status_.value();
    }

    void BarPositionFitter::clear()
    {
        fit_status_.reset();
        horizontal_bars_.clear();
        vertical_bars_.clear();
        y_pos_fitting_function_.SetParameter(0, 0.);
        y_pos_fitting_function_.SetParameter(1, 0.);
        x_pos_fitting_function_.SetParameter(0, 0.);
        x_pos_fitting_function_.SetParameter(1, 0.);
    }

    void BarPositionFitter::add_one_bar_signal(const BarPosition& bar_position)
    {
        if (fit_status_.has_value())
        {
            fit_status_.reset();
        }
        auto& bar_positions = bar_position.is_horizontal ? horizontal_bars_ : vertical_bars_;
        bar_positions.displacements.push_back(bar_position.displacement);
        bar_positions.displacement_errors.push_back(bar_position.displacement_error);
        bar_positions.pos_z.push_back(bar_position.pos_z);
    }

    auto BarPositionFitter::get_prediction(bool is_horizontal, double pos_z) const -> double
    {
        if (not fit_status_.has_value())
        {
            throw R3B::logic_error("Cannot make prediction without fitting the data first!");
        }
        return is_horizontal ? x_pos_fitting_function_(pos_z) : y_pos_fitting_function_(pos_z);
    }

    auto BarPositionFitter::fit_with(const BarPositions& positions, TF1& func) -> bool
    {
        if (positions.size() < 2)
        {
            return false;
        }
        auto fit_data = ROOT::Fit::BinData{ positions.size(),
                                            positions.pos_z.data(),
                                            positions.displacements.data(),
                                            nullptr,
                                            positions.displacement_errors.data() };
        fitter_.SetFunction(ROOT::Math::WrappedMultiTF1{ func, static_cast<unsigned int>(func.GetNdim()) }, false);
        const auto is_success = fitter_.Fit(fit_data);
        if (not is_success)
        {
            R3BLOG(error, "Fitting failed!");
        }
        return true;
    }
} // namespace R3B::Neuland::Calibration
