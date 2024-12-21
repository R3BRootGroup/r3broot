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

#include <bitset>
#include <cstdint>
#include <ext_data_client.h>
#include <fmt/format.h>
#include <ostream>

constexpr auto UCESB_MAP_BITSIZE = 32;

namespace R3B
{
    enum class UcesbMap : uint32_t
    {
        zero = 0x0000U,
        match = EXT_DATA_ITEM_MAP_MATCH,
        no_dest = EXT_DATA_ITEM_MAP_NO_DEST,
        not_found = EXT_DATA_ITEM_MAP_NOT_FOUND,
        type_mismatch = EXT_DATA_ITEM_MAP_TYPE_MISMATCH,
        ctrl_mismatch = EXT_DATA_ITEM_MAP_CTRL_MISMATCH,
        array_fewer = EXT_DATA_ITEM_MAP_ARRAY_FEWER,
        array_more = EXT_DATA_ITEM_MAP_ARRAY_MORE,
        not_done = EXT_DATA_ITEM_MAP_NOT_DONE,
        ok = EXT_DATA_ITEM_MAP_OK,
        ok_no_dest = EXT_DATA_ITEM_MAP_OK_NO_DEST,
    };

    using UcesbMapUType = std::underlying_type<UcesbMap>::type;

    inline auto operator|(const UcesbMap& left, const UcesbMap& right) -> UcesbMap
    {
        auto res = static_cast<UcesbMapUType>(left) | static_cast<UcesbMapUType>(right);
        return static_cast<UcesbMap>(res);
    }
    inline void operator|=(UcesbMap& left, const UcesbMap& right) { left = left | right; }

    inline auto operator&(const UcesbMap& left, const UcesbMap& right) -> UcesbMap
    {
        auto res = static_cast<UcesbMapUType>(left) & static_cast<UcesbMapUType>(right);
        return static_cast<UcesbMap>(res);
    }

    inline void operator&=(UcesbMap& left, const UcesbMap& right) { left = left & right; }

    inline auto operator~(const UcesbMap& map) -> UcesbMap
    {
        return static_cast<UcesbMap>(~(static_cast<UcesbMapUType>(map)));
    }

    inline auto operator<<(std::ostream& ostream, const UcesbMap& map) -> std::ostream&
    {
        ostream << std::bitset<UCESB_MAP_BITSIZE>(static_cast<UcesbMapUType>(map));
        return ostream;
    }

    inline auto UcesbMap2String(UcesbMap map) -> std::string
    {
        switch (map)
        {
            case UcesbMap::zero:
                return std::string{ "zero" };
            case UcesbMap::match:
                return std::string{ "match" };
            case UcesbMap::no_dest:
                return std::string{ "no_dest" };
            case UcesbMap::not_found:
                return std::string{ "not_found" };
            case UcesbMap::type_mismatch:
                return std::string{ "type_mismatch" };
            case UcesbMap::ctrl_mismatch:
                return std::string{ "ctrl_mismatch" };
            case UcesbMap::array_fewer:
                return std::string{ "array_fewer" };
            case UcesbMap::array_more:
                return std::string{ "array_more" };
            case UcesbMap::not_done:
                return std::string{ "not_done" };
            case UcesbMap::ok:
                return std::string{ "ok" };
            case UcesbMap::ok_no_dest:
                return std::string{ "ok_no_dest" };
            default:
                return fmt::format("", std::bitset<UCESB_MAP_BITSIZE>(static_cast<UcesbMapUType>(map)).to_string());
        }
    }
} // namespace R3B

template <>
class fmt::formatter<R3B::UcesbMap>
{
  public:
    constexpr auto parse(format_parse_context& ctx)
    {
        if (ctx.begin() == ctx.end())
        {
            return ctx.end();
        }
        const auto* specifier_iter = ctx.begin();
        if (*specifier_iter == 's' or *specifier_iter == 'b')
        {
            presentation_ = *specifier_iter;
            ++specifier_iter; // NOLINT
        }

        if (specifier_iter == ctx.end())
        {
            throw format_error("UcesbMap format failed: missing right curly bracket!");
        }

        if (*specifier_iter != '}')
        {
            throw format_error("UcesbMap format failed: only one speicifier is allowed");
        }
        return specifier_iter;
    }

    template <typename FmtContent>
    constexpr auto format(const R3B::UcesbMap& flag, FmtContent& ctn) const
    {
        switch (presentation_)
        {
            default:
            case 's':
                return format_to(ctn.out(), "{}", R3B::UcesbMap2String(flag));
            case 'b':
                return format_to(ctn.out(), "0x{:b}", static_cast<R3B::UcesbMapUType>(flag));
        }
    }

  private:
    char presentation_ = 's';
};
