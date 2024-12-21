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

#include "R3BUcesbMappingFlag.h"
#include <R3BReader.h>
#include <ext_data_struct_info.hh>
#include <vector>

struct ext_data_structure_item;
namespace R3B
{
    class UcesbSource;
    class UcesbStructInfo
    {
      public:
        UcesbStructInfo() = default;
        auto Get() -> ext_data_struct_info* { return &struct_info_; }
        void CheckStructMapping(UcesbSource* source);
        void SetExtraMapFlags(UcesbMap extr_flag) { map_success_conditions_ |= extr_flag; }
        void SetMapSuccessCondition(UcesbMap condition) { map_success_conditions_ = condition; }

      private:
        void terminate();
        bool check_struct_item(ext_data_structure_item* item, UcesbSource* source);
        ext_data_struct_info struct_info_;
        UcesbMap map_success_conditions_ = UcesbMap::ok;

        // void check_struct_info_mapping_old(const UcesbMap& is_map_success);
    };
} // namespace R3B
