/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
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
#include <FairFileSourceBase.h>
#include <TObjString.h>
#include <chrono>
#include <optional>

class FairRootManager;
class TChain;
class TFolder;

class R3BEventProgressPrinter
{
  public:
    R3BEventProgressPrinter() = default;
    void SetRunID(unsigned int runID) { run_id_ = runID; }
    void SetMaxEventNum(unsigned int max_event_num) { max_event_num_ = max_event_num; }
    void SetRefreshRate_Hz(float rate);
    void ShowProgress(uint64_t event_num);

  private:
    uint64_t max_event_num_ = 0;
    float refresh_rate_ = 2.; // Hz
    std::chrono::milliseconds refresh_period_{ static_cast<int>(1000. / refresh_rate_) };
    unsigned int run_id_ = 0;
    std::chrono::time_point<std::chrono::steady_clock> previous_t_ = std::chrono::steady_clock::now();
    uint64_t previous_event_num_ = 0;

    void Print(uint64_t event_num, double speed_per_ms);
};

class R3BInputRootFiles
{
  public:
    using Strings = std::vector<std::string>;
    R3BInputRootFiles() = default;
    auto AddFileName(std::string name, bool is_tree_file = false) -> std::optional<std::string>;
    void SetInputFileChain(TChain* chain);
    void RegisterTo(FairRootManager*);
    [[nodiscard]] auto is_empty() const -> bool { return fileNames_.empty(); }

    [[nodiscard]] auto is_friend() const -> bool { return is_friend_; }
    void Make_as_friend() { is_friend_ = true; }
    void SetFriend(R3BInputRootFiles& friendFiles);
    // Getters:
    [[nodiscard]] auto GetBranchListRef() const -> const auto& { return branchList_; }
    [[nodiscard]] auto GetBaseFileName() const -> const auto& { return fileNames_.front(); }
    [[nodiscard]] auto GetTreeName() const -> const auto& { return treeName_; }
    [[nodiscard]] auto GetFolderName() const -> const auto& { return folderName_; }
    [[nodiscard]] auto GetTitle() const -> const auto& { return title_; }
    [[nodiscard]] auto GetEntries() const -> int64_t;
    [[nodiscard]] auto GetChain() const -> TChain* { return rootChain_; }
    [[nodiscard]] auto GetInitialRunID() const { return initial_RunID_; }

    // Setters:
    void SetTreeName(std::string_view treeName) { treeName_ = treeName; }
    void SetTitle(std::string_view title) { title_ = title; }
    void SetFileHeaderName(std::string_view fileHeader) { fileHeader_ = fileHeader; }
    void SetRunID(uint run_id) { initial_RunID_ = run_id; }

    // rule of five:
    ~R3BInputRootFiles() = default;
    R3BInputRootFiles(const R3BInputRootFiles&) = delete;
    R3BInputRootFiles(R3BInputRootFiles&&) = default;
    R3BInputRootFiles& operator=(const R3BInputRootFiles&) = delete;
    R3BInputRootFiles& operator=(R3BInputRootFiles&&) = default;

  private:
    bool is_friend_ = false;
    uint initial_RunID_ = 0;
    // TODO: title of each file group seems not necessary. Consider to remove it in the future.
    std::string title_;
    std::string treeName_ = "evt";
    std::string folderName_;
    std::string fileHeader_;
    Strings fileNames_;
    Strings branchList_;
    std::vector<TObjString> timeBasedBranchList_;
    std::vector<R3B::unique_rootfile> validRootFiles_;
    std::vector<TFolder*> validMainFolders_;
    TChain* rootChain_ = nullptr;

    void Intitialize(std::string_view filename, bool is_tree_file = false);
    auto ValidateFile(const std::string& filename, bool is_tree_file = false) -> bool;
    static auto ExtractMainFolder(TFile*) -> std::optional<TKey*>;
    auto ExtractRunId(TFile* rootFile) -> std::optional<uint>;
    void register_branch_name();
};

class R3BFileSource2 : public FairFileSourceBase
{
  public:
    R3BFileSource2();
    explicit R3BFileSource2(std::string file, std::string_view title = "InputRootFile");
    R3BFileSource2(std::vector<std::string> fileNames, std::string_view title);
    explicit R3BFileSource2(std::vector<std::string> fileNames);

    void AddFile(std::string file_name, bool is_tree_file = false);
    void AddFile(std::vector<std::string> file_names, bool is_tree_file = false);
    void AddFriend(std::string, bool is_tree_file = false);
    void AddFriend(std::vector<std::string> file_names, bool is_tree_file = false);

    [[nodiscard]] auto GetEventEnd() const { return event_end_; }

    // setters:
    void SetFileHeaderName(std::string_view fileHeaderName) { inputDataFiles_.SetFileHeaderName(fileHeaderName); }
    // Set event print refresh rate in Hz
    void SetEventPrintRefreshRate(float rate) { event_progress_.SetRefreshRate_Hz(rate); }
    void SetInitRunID(int run_id)
    {
        inputDataFiles_.SetRunID(run_id);
        SetRunId(run_id);
    }

  private:
    int event_end_ = 0;
    R3BInputRootFiles inputDataFiles_;
    R3BEventProgressPrinter event_progress_;
    FairEventHeader* evtHeader_ = nullptr;
    std::vector<R3BInputRootFiles> inputFriendFiles_;
    std::vector<std::string> dataFileNames_;
    std::vector<std::string> friendFileNames_;

    Bool_t Init() override;
    Int_t ReadEvent(UInt_t eventID = 0) override;
    void Close() override {}
    void Reset() override {}
    Bool_t InitUnpackers() override { return kTRUE; }
    Bool_t ReInitUnpackers() override { return kTRUE; }
    Source_Type GetSourceType() override { return kFILE; }
    void SetParUnpackers() override {}
    Int_t CheckMaxEventNo(Int_t EvtEnd = 0) override;
    void ReadBranchEvent(const char* BrName) override;
    void ReadBranchEvent(const char* BrName, Int_t Entry) override;
    void FillEventHeader(FairEventHeader* evtHeader) override;
    Bool_t ActivateObject(TObject** obj, const char* BrName) override;
    Bool_t ActivateObjectAny(void** obj, const std::type_info& info, const char* BrName) override;
    // WTF is this?
    Bool_t SpecifyRunId() override { return true; }

  public:
    ClassDefOverride(R3BFileSource2, 0) // NOLINT
};
