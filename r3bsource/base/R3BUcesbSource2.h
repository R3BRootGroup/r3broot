#pragma once

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

#include "R3BReader.h"
#include <FairSource.h>
#include <R3BUcesbLauncher.h>
#include <R3BUcesbMappingFlag.h>
#include <R3BUcesbStructInfo.h>
#include <ext_data_clnt.hh>

struct EXT_STR_h101_t;
using EventStructType = EXT_STR_h101_t;

namespace R3B
{
    class UcesbSource : public FairSource
    {
      public:
        UcesbSource() = default;
        UcesbSource(std::string_view lmdfile_name,
                    std::string_view ntuple_options,
                    std::string_view ucesb_path,
                    EventStructType* event_struct,
                    size_t event_struct_size);
        // rule of five:
        ~UcesbSource() override;
        UcesbSource(const UcesbSource&) = delete;
        UcesbSource(UcesbSource&&) = delete;
        UcesbSource& operator=(const UcesbSource&) = delete;
        UcesbSource& operator=(UcesbSource&&) = delete;

        // setters:
        void SetMaxEvents(unsigned int max_event_num) { max_event_num_ = max_event_num; }
        void SetLMDFileName(std::string_view lmdfile_name) { lmdfile_name_ = lmdfile_name; }
        void SetNTupleOptions(std::string_view ntuple_options) { ntuple_options_ = ntuple_options; }
        void SetUcesbPath(std::string_view ucesb_path) { ucesb_path_ = ucesb_path; }
        void SetEventStructSize(size_t event_size) { event_struct_size_ = event_size; }
        // non-owning
        void SetEventStruct(EventStructType* event_struct) { event_struct_ = event_struct; }
        void SetRawDataPrint(bool print_raw_data) { has_raw_data_printing_ = print_raw_data; }
        void SetRunID(unsigned int run_id) { run_id_ = run_id; }
        void AllowExtraMap(UcesbMap flag) { ucesb_client_struct_info_.SetExtraMapFlags(flag); }

        template <typename ReaderType>
        auto AddReader(std::unique_ptr<ReaderType> reader) -> ReaderType*;
        template <typename ReaderType, typename... Args>
        auto AddReader(Args&&... args) -> ReaderType*;
        // TODO: C++20 concepts
        template <typename UnaryOp>
        void ForEachReader(UnaryOp&& opt);

        template <typename Predicate>
        auto FindReaderIf(Predicate&& pred) -> R3BReader*;

        // deprecate the old API because of bad memory managerment
        [[deprecated("Please use smart pointer method to add a reader")]] auto* AddReader(R3BReader* a_reader)
        {
            return AddReader(std::unique_ptr<R3BReader>(a_reader));
        }

      private:
        bool has_raw_data_printing_ = false;
        unsigned int run_id_ = 0;
        unsigned int max_event_num_ = 0;
        size_t event_struct_size_ = 0;
        EventStructType* event_struct_ = nullptr; // non-owning
        R3BEventHeader* event_header_ = nullptr;  // non-owning
        std::vector<std::unique_ptr<R3BReader>> readers_;
        std::string lmdfile_name_;
        std::string ntuple_options_;
        std::string ucesb_path_;

        ext_data_clnt ucesb_client_;
        UcesbStructInfo ucesb_client_struct_info_;
        UcesbServerLauncher ucesb_server_launcher_ = UcesbServerLauncher{ &ucesb_client_ };

        // private non-virtual methods:

        void init_runID();
        void init_ucesb();
        void init_ucesb_bp();
        void init_readers();
        void setup_ucesb();
        void print_raw_data();

        // private virtual methods:
        bool Init() override;
        bool InitUnpackers() override;
        bool ReInitUnpackers() override;
        void Close() override {}
        void SetParUnpackers() override
        {
            ForEachReader([](auto& reader) { reader->SetParContainers(); });
        }
        void Reset() override
        {
            ForEachReader([](auto& reader) { reader->Reset(); });
        }
        void FillEventHeader(FairEventHeader* feh) override;
        int ReadEvent(unsigned int eventID = 0) override;
        int CheckMaxEventNo(int EvtEnd = 0) override;
        bool SpecifyRunId() override { return true; }
        Source_Type GetSourceType() override { return kONLINE; }

      public:
        ClassDefInlineOverride(R3B::UcesbSource, 1);
    };

    template <typename ReaderType>
    auto UcesbSource::AddReader(std::unique_ptr<ReaderType> reader) -> ReaderType*
    {
        auto& reader_ref = readers_.emplace_back(std::move(reader));
        return static_cast<ReaderType*>(reader_ref.get());
    }

    template <typename ReaderType, typename... Args>
    auto UcesbSource::AddReader(Args&&... args) -> ReaderType*
    {
        static_assert(std::is_base_of_v<R3BReader, ReaderType>, "The reader type must be derived from R3BReader!");
        return AddReader(std::make_unique<ReaderType>(std::forward<Args>(args)...));
    }

    template <typename UnaryOp>
    void UcesbSource::ForEachReader(UnaryOp&& opt)
    {
        for (auto& reader : readers_)
        {
            opt(reader);
        }
    }

    template <typename Predicate>
    auto UcesbSource::FindReaderIf(Predicate&& pred) -> R3BReader*
    {
        auto res = std::find_if(readers_.begin(), readers_.end(), [&pred](auto& reader) { return pred(reader.get()); });
        return (res == readers_.end()) ? nullptr : res->get();
    }
} // namespace R3B

using R3BUcesbSource2 = R3B::UcesbSource;
