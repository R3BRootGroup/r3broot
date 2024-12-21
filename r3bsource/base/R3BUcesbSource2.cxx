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

#include "R3BUcesbSource2.h"
#include <FairRootManager.h>
#include <FairRun.h>
#include <R3BEventHeader.h>
#include <R3BException.h>
#include <R3BLogger.h>
#include <R3BUcesbDecl.h>
#include <boost/core/span.hpp>
#include <ext_data_client.h>
#include <fmt/format.h>

namespace R3B
{
    UcesbSource::UcesbSource(std::string_view lmdfile_name,
                             std::string_view ntuple_options,
                             std::string_view ucesb_path,
                             EventStructType* event_struct,
                             size_t event_struct_size)
        : event_struct_size_{ event_struct_size }
        , event_struct_{ event_struct }
        , lmdfile_name_{ lmdfile_name }
        , ntuple_options_{ ntuple_options }
        , ucesb_path_{ ucesb_path }
    {
    }

    bool UcesbSource::Init()
    {
        init_runID();
        init_ucesb();
        return true;
    }

    UcesbSource::~UcesbSource() { ucesb_server_launcher_.Close(); }

    void UcesbSource::init_ucesb()
    {
        auto command_string = fmt::format("{0} {1} --ntuple={2},STRUCT,-", ucesb_path_, lmdfile_name_, ntuple_options_);
        if (max_event_num_ > 0)
        {
            command_string = fmt::format("{} --max-events={}", command_string, max_event_num_);
        }
        R3BLOG(info, fmt::format("Calling ucesb with command: {}", command_string));

        ucesb_server_launcher_.Launch(std::move(command_string));
    }

    bool UcesbSource::InitUnpackers()
    {
        if (auto* frm = FairRootManager::Instance(); frm != nullptr)
        {
            R3BLOG(debug, "Checking the register of R3BEventHeader");
            if (event_header_ = dynamic_cast<R3BEventHeader*>(frm->GetObject("EventHeader.")); event_header_ == nullptr)
            {
                throw R3B::runtime_error("EventHeader. was not defined properly!");
            }
        }
        R3BLOG(debug, "EventHeader. was defined properly");

        init_readers();
        setup_ucesb();

        return true;
    }

    void UcesbSource::setup_ucesb()
    {
        // TODO: convert to std::bitset
        // could be initialzed in type UcesbMap. But C++ doesn't allow static cast of enum class pointer to its
        // underlying type
        auto is_struct_map_success = uint32_t{};
        if (ucesb_client_.setup(
                nullptr, 0, ucesb_client_struct_info_.Get(), &is_struct_map_success, event_struct_size_) == 0)
        {
            ucesb_client_struct_info_.CheckStructMapping(this);
        }
        else
        {
            R3BLOG(error, "ext_data_clnt::setup() failed");
            const auto* msg = (ucesb_client_.last_error() == nullptr) ? UCESB_NULL_STR_MSG : ucesb_client_.last_error();
            throw R3B::runtime_error(fmt::format("UCESB error: {}", msg));
        }
    }

    int UcesbSource::ReadEvent(unsigned int /*eventID*/)
    {
        auto ret_val = ucesb_client_.fetch_event(event_struct_, event_struct_size_);
        if (ret_val > 0)
        {
            ForEachReader([](auto& reader) { reader->R3BRead(); });
        }
        else if (ret_val == 0)
        {
            R3BLOG(info, "Reached the maximal event num on the ucesb server.");
            // ending event loop here
            return 1;
        }
        else
        {
            R3BLOG(error, "ext_data_clnt::fetch_event() failed");
            const auto* msg = (ucesb_client_.last_error() == nullptr) ? UCESB_NULL_STR_MSG : ucesb_client_.last_error();
            throw R3B::runtime_error(fmt::format("UCESB error: {}", msg));
        }

        return 0;
    }

    void print_uint32_with_size(const uint32_t* data, ssize_t size)
    {
        // TODO: use ranges library instead of reinterpret_cast
        constexpr auto column_size = 8;
        using SubDataType = const std::array<uint32_t, column_size>;

        auto data_span = std::span<SubDataType>(reinterpret_cast<SubDataType*>(data), size / column_size);
        R3BLOG(info, "Raw data:");
        auto index = uint32_t{};
        for (const auto& row_data : data_span)
        {
            fmt::print("RAW{0:04x}: {1:08x}\n", index, fmt::join(row_data, " "));
            index += column_size;
        }
    }

    void UcesbSource::print_raw_data()
    {
        // TODO: what's the best way to deal with this void** monstrosity
        const void* raw_data = nullptr;
        ssize_t raw_data_size = 0;
        auto ret_val = ucesb_client_.get_raw_data(&raw_data, &raw_data_size);

        if (ret_val == 0)
        {
            if (raw_data != nullptr)
            {
                const auto* data = reinterpret_cast<const uint32_t*>(raw_data);
                print_uint32_with_size(data, raw_data_size);
            }
        }
        else
        {
            R3BLOG(error, "ext_data_clnt::get_raw_data()");
            throw R3B::runtime_error("Failed to get raw data.");
        }
    }

    void UcesbSource::init_runID()
    {
        auto* run = FairRun::Instance();
        if (run == nullptr)
        {
            throw R3B::runtime_error("FairRun is not available!");
        }

        if (run_id_ != 0)
        {
            R3BLOG(info, fmt::format("Setting the run ID of the FairRun to be {} from FairSource", run_id_));
            run->SetRunId(run_id_);
        }
        else if (auto run_id = run->GetRunId(); run_id != 0)
        {
            R3BLOG(info, fmt::format("Setting the run ID of the FairSource to be {} from FairRun", run_id));
            run_id_ = run_id;
        }
        else
        {
            R3BLOG(warn, "Run ID of neither FairRun nor FairSource is set!");
        }
    }

    void UcesbSource::FillEventHeader(FairEventHeader* feh) { feh->SetRunId(run_id_); }

    int UcesbSource::CheckMaxEventNo(int EvtEnd)
    {
        max_event_num_ = (EvtEnd == 0) ? max_event_num_ : EvtEnd;
        return static_cast<int>(max_event_num_);
    }

    // readers looping methods:
    void UcesbSource::init_readers()
    {
        ForEachReader(
            [this](auto& reader)
            {
                if (not reader->Init(ucesb_client_struct_info_.Get()))
                {
                    const auto* msg =
                        (ucesb_client_.last_error() == nullptr) ? UCESB_NULL_STR_MSG : ucesb_client_.last_error();
                    throw R3B::runtime_error(fmt::format("UCESB error: {}", msg));
                }
            });
    }

    bool UcesbSource::ReInitUnpackers()
    {
        ForEachReader(
            [](auto& reader)
            {
                if (!reader->ReInit())
                {
                    throw R3B::runtime_error("ReInit of a reader failed.");
                }
            });
        return true;
    }
} // namespace R3B
