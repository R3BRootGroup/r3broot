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
        return fit_with(horizontal_bars_, y_pos_fitting_function_) and
               fit_with(vertical_bars_, x_pos_fitting_function_);
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
