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

#include "R3BSyncCheckData.h"

R3BSyncCheckData::R3BSyncCheckData()
    : fMaster(0)
    , fMasterRef(0)
    , fMusic(0)
    , fRpc(0)
    , fS2(0)
{
}

R3BSyncCheckData::R3BSyncCheckData(uint32_t master,
                                   uint32_t masterref,
                                   uint32_t music,
                                   uint32_t rpc,
                                   uint32_t s2,
                                   std::vector<uint32_t>& foot_input)
    : fMaster(master)
    , fMasterRef(masterref)
    , fMusic(music)
    , fRpc(rpc)
    , fS2(s2)
{
    for (auto& val : foot_input)
    {
        fFoot.push_back(val);
    }
}

ClassImp(R3BSyncCheckData);
