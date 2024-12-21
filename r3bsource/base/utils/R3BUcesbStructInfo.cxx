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

#include "R3BUcesbStructInfo.h"
#include "R3BUcesbDecl.h"
#include <R3BException.h>
#include <R3BLogger.h>
#include <R3BUcesbSource2.h>

namespace R3B
{
    // return pointer of the reader if it requires the item. If not required, return nullptr
    inline auto check_struct_item_requried(std::string_view item_name, UcesbSource* source) -> R3BReader*
    {
        return source->FindReaderIf([item_name](R3BReader* reader)
                                    { return reader->MismappedItemRequired(item_name); });
    }

    bool UcesbStructInfo::check_struct_item(ext_data_structure_item* item, UcesbSource* source)
    {
        auto map_flag = static_cast<UcesbMap>(item->_map_success);
        auto is_match_ok = ((map_flag & ~(map_success_conditions_)) == UcesbMap::zero);
        if (is_match_ok)
        {
            return true;
        }
        auto* required_reader = check_struct_item_requried(item->_var_name, source);
        if (required_reader == nullptr || required_reader->AllowExtraCondition(map_flag, map_success_conditions_))
        {
            // no reader requries this item or reader accepts extra flag
            return true;
        }

        R3BLOG(error,
               fmt::format("Failed to map the item {} required from {} due to the map flag: {}",
                           item->_var_name,
                           required_reader->GetName(),
                           map_flag));
        return false;
    }

    void UcesbStructInfo::CheckStructMapping(UcesbSource* source)
    {
        auto* struct_info = static_cast<ext_data_structure_info*>(struct_info_);
        auto is_checking_ok = true;
        for (auto* item = struct_info->_items; item != nullptr; item = item->_next_off_item)
        {
            is_checking_ok &= check_struct_item(item, source);
        }

        if (not is_checking_ok)
        {
            terminate();
        }
    }

    void UcesbStructInfo::terminate()
    {
        R3BLOG(error, "ext_data_clnt::setup() failed to map all items:");

        if (fair::Logger::GetConsoleSeverity() < fair::Severity::info)
        {
            ext_data_struct_info_print_map_success(static_cast<ext_data_structure_info*>(struct_info_),
                                                   stderr,
                                                   static_cast<UcesbMapUType>(map_success_conditions_));
        }
        throw R3B::runtime_error("ext_data_clnt::setup() mapping failure may cause unexpected analysis results "
                                 "due to missing data members. Unpacker needs fixing.\n\n\n");
    }

    // void UcesbStructInfo::check_struct_info_mapping_old(const UcesbMap& is_map_success)
    // {
    //     if ((is_map_success & ~(map_success_conditions_)) != UcesbMap::zero)
    //     {
    //         R3BLOG(error, "ext_data_clnt::setup() failed to map all items:");
    //         ext_data_struct_info_print_map_success(static_cast<ext_data_structure_info*>(ucesb_client_struct_info_),
    //                                                stderr,
    //                                                static_cast<UcesbMapUType>(map_success_conditions_));
    //         throw R3B::runtime_error("ext_data_clnt::setup() mapping failure may cause unexpected analysis results "
    //                                  "due to missing data members. Unpacker needs fixing.");
    //     }
    // }
} // namespace R3B
