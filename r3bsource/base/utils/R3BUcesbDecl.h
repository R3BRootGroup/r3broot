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

// This header file can totally be avoided if, for some reasons, Ucesb client code base didn't hide these declarations
// in SOURCE file.

#include <cstddef>
#include <cstdint>

constexpr auto EXT_DATA_STATE_INIT = 0;
constexpr auto EXT_DATA_STATE_OPEN = 1;
constexpr auto EXT_DATA_STATE_OPEN_OUT = 2;
constexpr auto EXT_DATA_STATE_PARSED_HEADERS = 3;
constexpr auto EXT_DATA_STATE_SETUP_READ = 4;
constexpr auto EXT_DATA_STATE_SETUP_WRITE = 5;

struct ext_data_structure_item
{
    uint32_t _offset;   /* Not used within STRUCT_WRITER itself. */
    uint32_t _length;   /* not needed, info only */
    const char* _block; /* not needed, info only */

    const char* _var_name;
    uint32_t _var_array_len;
    const char* _var_ctrl_name;
    uint32_t _var_type;
    uint32_t _limit_min;
    uint32_t _limit_max;

    uint32_t _map_success;

#if STRUCT_WRITER
    uint32_t _ctrl_offset;
#endif
    /* Used for remapping. */
    struct ext_data_structure_item* _ctrl_item;
    struct ext_data_structure_item* _next_off_item;
    /* Temporary used while creating remap list. */
    struct ext_data_structure_item* _match_item;
    struct ext_data_structure_item* _child_item;
};

struct ext_data_structure_info
{
    struct ext_data_structure_item* _items;

    /* Used while returning items for ext_data_struct_info_map_success(). */
    struct ext_data_structure_item* _ret_item;
    int _ret_for_server;

    struct ext_data_structure_info* _server_struct_info;

    uint32_t _map_success;

    const char* _last_error;
};

struct ext_data_client_struct
{
    const char* _id;

    /* Todo: the following controls _raw_ptr in ext_data_client.  With
     * several structures...?
     */
    uint32_t _max_raw_words;

    uint32_t _orig_xor_sum_msg;
    size_t _orig_struct_size;
    void* _orig_array; /* for mapping */

    uint32_t _orig_max_pack_items;
    uint32_t _orig_static_pack_items;

    uint32_t* _orig_pack_list;
    uint32_t* _orig_pack_list_end;

    size_t _dest_struct_size;

    uint32_t _dest_max_pack_items;
    uint32_t _dest_static_pack_items;

    uint32_t* _dest_pack_list;
    uint32_t* _dest_pack_list_end;

    uint32_t* _dest_reverse_pack;

    uint32_t* _map_list;
    uint32_t* _map_list_end;

    struct ext_data_structure_info* _struct_info_msg;
};

struct ext_data_client
{
    int _fd;
    int _fd_close; /* If using STDIN, we are not to close it, so -1. */

    char* _buf;
    size_t _buf_alloc;
    size_t _buf_used;
    size_t _buf_filled;

    uint32_t* _raw_ptr; /* This is not allocated; just used. */
    uint32_t _raw_words;
    uint32_t* _raw_swapped;

    struct ext_data_client_struct* _structures;
    int _num_structures;

    uint32_t _sort_u32_words;

    const char* _last_error;

    int _state;

    int _fetched_event;
};

/* Layout of the structure information generated.
 */

struct ext_data_structure_layout_item
{
    uint32_t _offset;
    uint32_t _size;
    uint32_t _xor;
    const char* _name;
};

struct ext_data_structure_layout
{
    uint32_t _magic;
    uint32_t _size_info;
    uint32_t _size_struct;
    uint32_t _size_struct_onion;
    uint32_t _pack_list_items;

    uint32_t _num_items;
    // NOLINTNEXTLINE
    struct ext_data_structure_layout_item _items[1]; /* more than 1 for parts */
};
