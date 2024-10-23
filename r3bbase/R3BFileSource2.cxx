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

#include "R3BFileSource2.h"

#include "R3BException.h"
#include "R3BLogger.h"
#include <FairEventHeader.h>
#include <FairFileHeader.h>
#include <FairRootManager.h>
#include <FairRun.h>
#include <TBranchElement.h>
#include <TClonesArray.h>
#include <TFolder.h>
#include <TKey.h>
#include <fmt/chrono.h>
#include <fmt/color.h>
#include <fmt/core.h>
#include <vector>

namespace
{
    constexpr auto DEFAULT_TITLE = "InputRootFile";

    template <typename ContainerType, typename DataType>
    auto Vector2TContainer(std::vector<DataType>& vec) -> std::unique_ptr<ContainerType>
    {
        using RawType = std::remove_reference_t<std::remove_pointer_t<DataType>>;
        static_assert(std::is_base_of_v<TObject, RawType>);
        auto list = std::make_unique<ContainerType>();
        for (auto& iter : vec)
        {
            if constexpr (std::is_pointer_v<DataType>)
            {
                list->Add(iter);
            }
            else
            {
                list->Add(&iter);
            }
        }
        return list;
    }

    template <typename StringType = std::string>
    auto GetBranchList(TFile* rootFile, std::string_view listName) -> std::vector<StringType>
    {
        auto branchList = std::vector<StringType>{};
        if (auto* list = dynamic_cast<TList*>(rootFile->Get(listName.data())); list != nullptr)
        {
            for (const auto& str : TRangeDynCast<TObjString>(list))
            {
                branchList.emplace_back(str->GetString().Data());
            }
        }
        else
        {
            throw R3B::logic_error(
                fmt::format("No branch list named {0} in input file {1}", listName, rootFile->GetName()));
        }
        return branchList;
    }

    template <typename UnaryFunc>
    void loop_through_branches(TFile* root_file, std::string_view tree_name, UnaryFunc&& action)
    {
        auto* tree = root_file->Get<TTree>(tree_name.data());
        auto* branches = tree->GetListOfBranches();
        for (auto* branch : TRangeDynCast<TBranchElement>(branches))
        {
            action(branch);
        }
    }

    template <typename StringType = std::string>
    auto GetBranchListFromTree(TFile* root_file, std::string_view tree_name) -> std::vector<StringType>
    {
        auto branch_name_list = std::vector<StringType>{};
        loop_through_branches(root_file,
                              tree_name,
                              [&branch_name_list](auto* branch) { branch_name_list.emplace_back(branch->GetName()); });
        return branch_name_list;
    }

    auto get_tca_data_class(TBranchElement* branch) -> std::string
    {
        TClonesArray* buffer = nullptr;
        branch->SetAddress(&buffer);
        branch->GetEntry(0);
        branch->SetAddress(nullptr);
        if (buffer != nullptr)
        {
            auto class_name = std::string{ buffer->GetClass()->GetName() };
            R3BLOG(debug,
                   fmt::format("Determine the class name {:?} of the branch {:?}", class_name, branch->GetName()));
            return class_name;
        }
        R3BLOG(warn, fmt::format("Cannot determine the class name of the branch {:?}", branch->GetName()));
        return std::string{ "TObject" };
    }

    void add_branches_to_folder(TFolder* folder, TFile* root_file, std::string_view tree_name)
    {
        loop_through_branches(root_file,
                              tree_name,
                              [folder](auto* branch)
                              {
                                  auto class_name = std::string_view{ branch->GetClassName() };
                                  if (class_name == "TClonesArray")
                                  {
                                      const auto data_class = get_tca_data_class(branch);
                                      auto tca_obj = std::make_unique<TClonesArray>(data_class.data());
                                      tca_obj->SetName(branch->GetName());
                                      folder->Add(tca_obj.release());
                                  }
                              });
    }

    auto HasBranchList(TFile* rootFile, const std::vector<std::string>& branchList) -> bool
    {
        auto const newBranchList = GetBranchList(rootFile, "BranchList");
        auto view1 = std::vector<std::string_view>(branchList.begin(), branchList.end());
        auto view2 = std::vector<std::string_view>(newBranchList.begin(), newBranchList.end());

        std::sort(view1.begin(), view1.end());
        std::sort(view2.begin(), view2.end());
        return view1 == view2;
    }

    template <typename ContainerType>
    auto GetDataFromAnyFolder(TFile* rootFile, const ContainerType& folderNames) -> std::optional<TKey*>
    {
        for (auto const& name : folderNames)
        {
            R3BLOG(debug, "looking for " + name);
            auto* dataFolder = dynamic_cast<TKey*>(rootFile->FindKey(name.c_str()));
            if (dataFolder != nullptr)
            {
                R3BLOG(debug, name + " has been found!");
                return dataFolder;
            }
        }
        return {};
    }

    auto Get_TChain_FromFairRM(FairRootManager* rootMan) -> TChain*
    {
        auto const chainTitle = "/" + std::string{ FairRootManager::GetFolderName() };
        auto inChain = std::make_unique<TChain>(FairRootManager::GetTreeName(), chainTitle.c_str());
        R3BLOG(debug, "Chain created");
        LOG(info) << "chain name: " << FairRootManager::GetTreeName();
        rootMan->SetInChain(inChain.release());
        return FairRootManager::Instance()->GetInChain();
    }
} // namespace

void R3BEventProgressPrinter::SetRefreshRate_Hz(float rate)
{
    if (rate <= 0.)
    {
        throw R3B::logic_error(fmt::format("Refresh rate {} must be a positive floating point value", rate));
    }

    refresh_rate_ = rate;
    refresh_period_ = std::chrono::milliseconds(static_cast<int>(1000. / rate));
}

void R3BEventProgressPrinter::ShowProgress(uint64_t event_num)
{
    if (event_num < 1)
    {
        return;
    }
    if (max_event_num_ == 0)
    {
        throw R3B::logic_error("Maximal event number has not been set up!");
    }
    const auto now_t = std::chrono::steady_clock::now();
    const auto time_spent = std::chrono::ceil<std::chrono::milliseconds>(now_t - previous_t_);
    if (time_spent > refresh_period_)
    {
        const auto processed_events = event_num - previous_event_num_;
        const auto events_per_millisecond =
            static_cast<double>(processed_events) / static_cast<double>(time_spent.count());
        Print(event_num, events_per_millisecond);

        previous_t_ = now_t;
        previous_event_num_ = event_num;
    }
}

void R3BEventProgressPrinter::Print(uint64_t event_num, double speed_per_ms)
{
    if (speed_per_ms <= 0.)
    {
        return;
    }
    const auto event_num_str =
        fmt::format(fg(fmt::terminal_color::bright_green) | fmt::emphasis::bold, "{:^5d}k", event_num / 1000);
    const auto speed_str = fmt::format(fg(fmt::color::white), "{:^6.1F}k/s", speed_per_ms);
    const auto progress_str = fmt::format(fg(fmt::terminal_color::bright_yellow) | fmt::emphasis::bold,
                                          "{:^6.2F}",
                                          100. * static_cast<double>(event_num) / static_cast<double>(max_event_num_));
    const auto time_left_ms =
        std::chrono::milliseconds{ (max_event_num_ - event_num) / static_cast<int>(std::ceil(speed_per_ms)) };
    fmt::print("Events processed: {0} ({1})  Progress: {2}% (time left: {3:%H h %M m %S s})   Run ID: {4}\r",
               event_num_str,
               speed_str,
               progress_str,
               std::chrono::ceil<std::chrono::seconds>(time_left_ms),
               run_id_);
    std::cout << std::flush;
}

auto R3BInputRootFiles::AddFileName(std::string fileName, bool is_tree_file) -> std::optional<std::string>
{
    auto const msg = fmt::format("Adding {} to file source\n", fileName);
    R3BLOG(info, msg);
    if (fileNames_.empty())
    {
        Intitialize(fileName, is_tree_file);
        register_branch_name();
    }
    if (!ValidateFile(fileName, is_tree_file))
    {
        return fileName;
    }
    fileNames_.emplace_back(std::move(fileName));
    return {};
}

void R3BInputRootFiles::register_branch_name()
{

    for (auto const& branchName : branchList_)
    {
        FairRootManager::Instance()->AddBranchToList(branchName.c_str());
    }
}

void R3BInputRootFiles::SetInputFileChain(TChain* chain)
{
    if (rootChain_ != nullptr)
    {
        throw R3B::logic_error("TChain has already been created!");
    }
    rootChain_ = chain;
    for (auto const& filename : fileNames_)
    {
        rootChain_->AddFile(filename.c_str(), TTree::kMaxEntries, treeName_.c_str());
    }
}
void R3BInputRootFiles::RegisterTo(FairRootManager* rootMan)
{
    if (is_friend_)
    {
        return;
    }

    if (validMainFolders_.empty())
    {
        throw R3B::runtime_error("There is no main folder to be registered!");
    }

    if (!is_friend_)
    {
        auto listOfFolders = Vector2TContainer<TObjArray>(validMainFolders_);
        R3BLOG(debug, fmt::format("Set {} main folder(s) to FairRootManager.", listOfFolders->GetEntries()));
        rootMan->SetListOfFolders(listOfFolders.release());
        rootMan->SetTimeBasedBranchNameList(Vector2TContainer<TList>(timeBasedBranchList_).release());
        SetInputFileChain(Get_TChain_FromFairRM(rootMan));
    }
}

auto R3BInputRootFiles::ExtractMainFolder(TFile* rootFile) -> std::optional<TKey*>
{
    auto const folderNames =
        std::array<std::string, 4>{ FairRootManager::GetFolderName(), "r3broot", "cbmout", "cbmroot" };

    return GetDataFromAnyFolder(rootFile, folderNames);
}

auto R3BInputRootFiles::ValidateFile(const std::string& filename, bool is_tree_file) -> bool
{
    auto rootFile = R3B::make_rootfile(filename.c_str());

    if (is_tree_file)
    {
        if (!is_friend_)
        {
            auto folder = std::make_unique<TFolder>("r3broot", "r3broot");
            add_branches_to_folder(folder.get(), rootFile.get(), treeName_);
            validRootFiles_.push_back(std::move(rootFile));
            validMainFolders_.push_back(folder.release());
        }
        return true;
    }

    auto folderKey = ExtractMainFolder(rootFile.get());
    auto res = folderKey.has_value() && HasBranchList(rootFile.get(), branchList_);
    if (res)
    {
        if (!folderName_.empty() && (folderKey.value()->GetName() != folderName_))
        {
            R3BLOG(warn, "Different folder name!");
        }
        if (!is_friend_)
        {
            validRootFiles_.push_back(std::move(rootFile));
            validMainFolders_.push_back((folderKey.value())->ReadObject<TFolder>());
        }
    }
    return res;
}

auto R3BInputRootFiles::ExtractRunId(TFile* rootFile) -> std::optional<uint>
{
    //
    auto* header = rootFile->Get<FairFileHeader>(fileHeader_.c_str());
    if (header == nullptr)
    {
        return {};
    }
    auto runID = header->GetRunId();
    return runID;
}

void R3BInputRootFiles::Intitialize(std::string_view filename, bool is_tree_file)
{
    auto file = R3B::make_rootfile(filename.data());

    if (is_tree_file)
    {
        branchList_ = GetBranchListFromTree(file.get(), treeName_);
        return;
    }

    if (const auto runID = ExtractRunId(file.get()); runID.has_value() && runID.value() != 0)
    {
        auto const msg = fmt::format(R"(Successfully extract RunID "{}" from root file "{}")", runID.value(), filename);
        R3BLOG(debug, msg);
        initial_RunID_ = runID.value();
    }
    else
    {
        auto const msg = fmt::format("Failed to extract RunID from root file \"{}\"", filename);
        R3BLOG(error, msg);
    }

    if (auto folderKey = ExtractMainFolder(file.get()); folderKey.has_value())
    {
        folderName_ = folderKey.value()->GetName();
    }
    else
    {
        throw R3B::logic_error(fmt::format("Cannot find main folder from the root file {}!", filename));
    }

    branchList_ = GetBranchList(file.get(), "BranchList");

    if (timeBasedBranchList_ = GetBranchList<TObjString>(file.get(), "TimeBasedBranchList");
        timeBasedBranchList_.empty())
    {
        LOG(warn) << "No time based branch list in input file";
    }
}

void R3BInputRootFiles::SetFriend(R3BInputRootFiles& friendFiles)
{
    if (is_friend_)
    {
        throw R3B::logic_error("Can not set friendFiles with another friendFile!");
    }
    auto chain = std::make_unique<TChain>(friendFiles.GetTitle().c_str(), friendFiles.GetFolderName().c_str());
    friendFiles.SetInputFileChain(chain.get());
    rootChain_->AddFriend(chain.release());
}

[[nodiscard]] auto R3BInputRootFiles::GetEntries() const -> int64_t
{
    if (rootChain_ == nullptr)
    {
        throw R3B::logic_error("Can't get entries before being initialized!");
    }
    return rootChain_->GetEntries();
}

R3BFileSource2::R3BFileSource2(std::vector<std::string> fileNames, std::string_view title)
{
    LOG(debug) << "Creating a new R3BFileSource!";
    inputDataFiles_.SetTitle(title);
    inputDataFiles_.SetFileHeaderName("FileHeader");
    for (auto& name : fileNames)
    {
        if (name.empty())
        {
            continue;
        }
        AddFile(std::move(name));
    }
}

R3BFileSource2::R3BFileSource2(std::string file, std::string_view title)
    : R3BFileSource2(std::vector<std::string>{ std::move(file) }, title)
{
}

R3BFileSource2::R3BFileSource2(std::vector<std::string> fileNames)
    : R3BFileSource2(std::move(fileNames), DEFAULT_TITLE)
{
}

R3BFileSource2::R3BFileSource2()
    : R3BFileSource2(std::string{})
{
}

void R3BFileSource2::AddFile(std::string file_name, bool is_tree_file)
{
    if (auto const res = inputDataFiles_.AddFileName(std::move(file_name), is_tree_file); res.has_value())
    {
        if (not dataFileNames_.empty())
        {

            R3BLOG(
                error,
                fmt::format(
                    "Root file {0} is incompatible with the first root file {1}", res.value(), dataFileNames_.front()));
        }
        else
        {
            R3BLOG(error, fmt::format("Failed to add the first root file {:?}", file_name));
        }
    }
    dataFileNames_.emplace_back(file_name);
}

void R3BFileSource2::AddFile(std::vector<std::string> file_names, bool is_tree_file)
{
    for (auto& file_name : file_names)
    {
        AddFile(std::move(file_name), is_tree_file);
    }
}

void R3BFileSource2::AddFriend(std::vector<std::string> file_names, bool is_tree_file)
{
    for (auto& file_name : file_names)
    {
        AddFriend(std::move(file_name), is_tree_file);
    }
}

void R3BFileSource2::AddFriend(std::string file_name, bool is_tree_file)
{
    //
    auto rootfile = R3B::make_rootfile(file_name.c_str());
    auto friendGroup = std::find_if(inputFriendFiles_.begin(),
                                    inputFriendFiles_.end(),
                                    [&rootfile](const auto& friends)
                                    { return HasBranchList(rootfile.get(), friends.GetBranchListRef()); });
    if (friendGroup == inputFriendFiles_.end())
    {
        auto newFriendGroup = R3BInputRootFiles{};
        newFriendGroup.Make_as_friend();
        inputFriendFiles_.push_back(std::move(newFriendGroup));
        friendGroup = --inputFriendFiles_.end();
        friendGroup->SetTitle(fmt::format("FriendTree_{}", inputFriendFiles_.size()));
    }
    auto res = friendGroup->AddFileName(file_name, is_tree_file);
    if (res.has_value())
    {
        R3BLOG(error,
               fmt::format("Friend file {0} is incompatible with the first friend file {1}",
                           res.value(),
                           friendGroup->GetBaseFileName()));
    }
    else
    {
        // TODO: really need it?
        friendFileNames_.emplace_back(std::move(file_name));
    }
}

Bool_t R3BFileSource2::Init()
{
    if (inputDataFiles_.is_empty())
    {
        throw R3B::logic_error{ "No input file available!" };
    }

    inputDataFiles_.RegisterTo(FairRootManager::Instance());

    for (auto& friendGroup : inputFriendFiles_)
    {
        inputDataFiles_.SetFriend(friendGroup);
    }

    event_progress_.SetMaxEventNum(inputDataFiles_.GetEntries());
    event_progress_.SetRunID(inputDataFiles_.GetInitialRunID());

    return true;
}

void R3BFileSource2::FillEventHeader(FairEventHeader* evtHeader)
{
    if (evtHeader == nullptr)
    {
        throw R3B::logic_error("Filled event header is empty!");
    }
    evtHeader_ = evtHeader;

    // Set runID for event header:
    auto const init_runID = inputDataFiles_.GetInitialRunID();

    if (init_runID == 0)
    {
        throw R3B::logic_error("RunId is not being set!");
    }

    if (init_runID != GetRunId())
    {
        R3BLOG(
            warn,
            fmt::format("runID {} being set is different from the runID {} in the data file!", GetRunId(), init_runID));
    }
    SetRunId(init_runID); // NOLINT
    evtHeader->SetRunId(init_runID);
}

Int_t R3BFileSource2::CheckMaxEventNo(Int_t EvtEnd)
{
    event_end_ = (EvtEnd <= 0) ? inputDataFiles_.GetEntries() : EvtEnd; // NOLINT
    R3BLOG(info, fmt::format("Setting printing event max to {}", event_end_));
    event_progress_.SetMaxEventNum(event_end_);
    return event_end_;
}

void R3BFileSource2::ReadBranchEvent(const char* BrName)
{
    auto const currentEventID = evtHeader_->GetMCEntryNumber();
    ReadBranchEvent(BrName, currentEventID);
}

void R3BFileSource2::ReadBranchEvent(const char* BrName, Int_t entryID)
{
    auto const read_bytes = inputDataFiles_.GetChain()->FindBranch(BrName)->GetEntry(entryID);
    if (read_bytes == 0)
    {
        LOG(warn) << fmt::format("Failed to read the data of the event {0} from the branch {1}", entryID, BrName);
    }
}

Int_t R3BFileSource2::ReadEvent(UInt_t eventID)
{
    auto* chain = inputDataFiles_.GetChain();
    if (fair::Logger::GetConsoleSeverity() == fair::Severity::info)
    {
        event_progress_.ShowProgress(eventID);
    }

    auto read_bytes = chain->GetEntry(eventID);
    if (read_bytes == 0)
    {
        LOG(warn) << fmt::format("Failed to read the data of the event {0} from the source", eventID);
        return 1;
    }
    return 0;
}

Bool_t R3BFileSource2::ActivateObject(TObject** obj, const char* BrName)
{
    auto* chain = inputDataFiles_.GetChain();
    chain->SetBranchStatus(BrName, true);
    chain->SetBranchAddress(BrName, obj);
    return kTRUE;
}

Bool_t R3BFileSource2::ActivateObjectAny(void** obj, const std::type_info& info, const char* BrName)
{
    auto* chain = inputDataFiles_.GetChain();
    if (chain != nullptr)
    {
        return ActivateObjectAnyImpl(chain, obj, info, BrName);
    }
    return kFALSE;
}

ClassImp(R3BFileSource2);
