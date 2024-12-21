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
#include <boost/asio.hpp>
#include <boost/process.hpp>
#include <ext_data_clnt.hh>

constexpr auto UCESB_NULL_STR_MSG = "Can't retrieve error message as last_error returns nullptr!";

namespace R3B
{
    class UcesbServerLauncher
    {
      public:
        explicit UcesbServerLauncher(ext_data_clnt* client)
            : client_{ client }
        {
        }
        void Launch(std::string command_string);
        void Setup(ext_data_struct_info& struct_info, size_t event_struct_size);
        void Close();

      private:
        ext_data_clnt* client_ = nullptr;
        std::unique_ptr<boost::process::child> ucesb_server_;
        boost::asio::io_service ios_;
        boost::process::async_pipe server_pipe_{ ios_ };
    };

} // namespace R3B
