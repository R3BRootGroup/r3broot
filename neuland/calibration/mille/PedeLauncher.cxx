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

#define BOOST_PROCESS_USE_STD_FS 1

#include "PedeLauncher.h"
#include <boost/process.hpp>
#include <filesystem>
#include <fmt/format.h>
#include <fmt/ranges.h>

namespace bp = boost::process;
namespace fs = std::filesystem;

namespace R3B::Millepede
{
    void Launcher::launch()
    {
#ifdef PEDEDIR
        const auto exe_string = fmt::format("{}/{}", PEDEDIR, executable_);
        auto launch_args = std::vector<std::string>{ "-i", steer_filename_ };
        auto exe_path = std::filesystem::path(exe_string);
        if (not std::filesystem::exists(exe_path))
        {
            throw std::runtime_error(fmt::format("Command {:?} cannot be found or executed!", exe_string));
        }

        try
        {
            auto pede_program =
                bp::child{ exe_path,
                           bp::args(launch_args),
                           bp::std_out > stdout,
                           bp::std_err > stderr,
                           ios_,
                           bp::on_exit = [](int pid, const std::error_code& err)
                           { fmt::print("Child process returns {} with error code: {}\n", pid, err.message()); } };
            ios_.run();
            pede_program.wait();
            fmt::print("Child process with PID {} ended.\n", pede_program.id());
            ios_.reset();
        }
        catch (std::exception& ex)
        {
            throw std::runtime_error(fmt::format("Error occured when launching pede with \"{} {}\": \n {}",
                                                 bp::search_path(executable_).string(),
                                                 fmt::join(launch_args, " "),
                                                 ex.what()));
        }
#endif
    }

    void Launcher::end()
    {
        // auto old_result_path = fs::current_path() / fs::path{ DEFAULT_RES_FILENAME };
        // auto new_result_path = fs::current_path() / fs::path{ parameter_file_ };

        // if (fs::exists(old_result_path))
        // {
        //     fmt::print(
        //         "Info: copy the file {:?} to the file {:?}!\n", old_result_path.string(), new_result_path.string());
        //     fs::copy(old_result_path, new_result_path, fs::copy_options::overwrite_existing);
        // }
        // else
        // {
        //     fmt::print("Error: millepede.res doesn't exist!\n");
        // }
    }

} // namespace R3B::Millepede
