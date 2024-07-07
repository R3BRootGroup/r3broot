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

#include "R3BDataMonitor.h"
#include <FairRootFileSink.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>
#include <fmt/chrono.h>
#include <fmt/format.h>

namespace
{
    auto create_datatime_rootfile(std::string_view base_filename) -> std::unique_ptr<TFile>
    {
        auto time = std::chrono::system_clock::now();
        const auto file_name = fmt::format("{}-{:%Y-%m-%d-%Hh-%Mm-%Ss}.root", base_filename, time);
        auto* root_file = TFile::Open(file_name.data(), "RECREATE");
        return std::unique_ptr<TFile>{ root_file };
    }

} // namespace

namespace R3B
{
    auto DataMonitor::get(const std::string& histName) -> TH1*
    {
        if (auto hist = histograms_.find(histName); hist != histograms_.end())
        {
            return hist->second.get();
        }
        throw R3B::logic_error(fmt::format("Histogram with the name {} doesn't exist!", histName));
    }

    auto DataMonitor::get_canvas(const std::string& histName) -> DataMonitorCanvas&
    {
        if (auto canvas = canvases_.find(histName); canvas != canvases_.end())
        {
            return canvas->second;
        }
        throw R3B::logic_error(fmt::format("Canvas with the name {} doesn't exist!", histName));
    }

    void DataMonitor::save_to_sink(std::string_view folderName, FairSink* sinkFile)
    {
        auto* hist_dir = get_hist_dir(sinkFile);
        auto* new_dir = folderName.empty() ? hist_dir : hist_dir->mkdir(folderName.data(), "", true);
        if (new_dir == nullptr)
        {
            throw R3B::runtime_error(
                fmt::format("Failed to create a sub directory {} for the histrogams!", folderName));
        }
        R3BLOG(info,
               fmt::format("Saving figures to the directory {:?} in the root file {:?}",
                           folderName,
                           new_dir->GetFile()->GetName()));

        write_all(new_dir);
        // old_dir->cd();
    }

    auto DataMonitor::get_hist_dir(FairSink* sinkFile) -> TDirectory*
    {
        auto* rootSinkFile = dynamic_cast<FairRootFileSink*>(sinkFile);
        if (rootSinkFile == nullptr)
        {
            throw R3B::logic_error("Cannot save the histograms as the output file is not a root file!");
        }
        auto* rootFile = rootSinkFile->GetRootFile();
        auto* hist_dir = rootFile->mkdir(DEFAULT_HIST_MONITOR_DIR, "", true);
        if (hist_dir == nullptr)
        {
            throw R3B::runtime_error("Cannot create a directory for the histrogams!");
        }
        return hist_dir;
    }

    void DataMonitor::draw_canvases()
    {
        for (auto& [canvas_name, canvas] : canvases_)
        {
            canvas.draw();
        }
    }

    void DataMonitor::register_canvases(FairRun* run)
    {
        if (auto* run_online = dynamic_cast<FairRunOnline*>(run); run_online != nullptr)
        {
            for (auto& [name, canvas] : canvases_)
            {
                run_online->AddObject(canvas.get_canvas());
            }
        }
    }

    void DataMonitor::reset_all_hists()
    {
        for (auto& [name, hist] : histograms_)
        {
            hist->Reset();
        }
    }

    void DataMonitor::save_to_file(std::string_view filename)
    {
        if (filename.empty())
        {
            filename = save_filename_;
        }
        auto rootfile = create_datatime_rootfile(filename);
        R3BLOG(info, fmt::format("Saving histograms to {}", rootfile->GetName()));
        write_all(rootfile.get());
    }

    void DataMonitor::write_all(TDirectory* dir)
    {
        for (auto& [name, hist] : histograms_)
        {
            if (hist->GetEntries() == 0)
            {
                R3BLOG(warn, fmt::format("Histogram {} is empty while written to the file!", hist->GetName()));
            }
            dir->WriteObject(hist.get(), hist->GetName());
        }
        for (auto& [name, graph] : graphs_)
        {
            dir->WriteObject(graph.get(), graph->GetName(), "overwrite");
        }
    }
} // namespace R3B
