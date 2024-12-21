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

#ifndef R3BSyncCheckData_H
#define R3BSyncCheckData_H 1

#include "R3BLogger.h"
#include "TObject.h"
#include <Rtypes.h>

class R3BSyncCheckData : public TObject
{
  public:
    R3BSyncCheckData();
    R3BSyncCheckData(uint32_t master,
                     uint32_t masterref,
                     uint32_t music,
                     uint32_t rpc,
                     uint32_t s2,
                     std::vector<uint32_t>& foot_input);

    virtual ~R3BSyncCheckData() {}

    inline const uint32_t& GetMaster() const { return fMaster; }
    inline const uint32_t& GetMasterRef() const { return fMasterRef; }
    inline const uint32_t& GetMusic() const { return fMusic; }
    inline const uint32_t& GetRpc() const { return fRpc; }
    inline const uint32_t& GetS2() const { return fS2; }
    const uint32_t& GetFoot(int num) const
    {
        R3BLOG_IF(fatal, (num <= 0 || num > NFOOTS), "Bad FOOT number!");
        return fFoot.at(num - 1);
    }

  protected:
    // following the order from ext file
    uint32_t fMaster;
    uint32_t fMasterRef;
    uint32_t fMusic;
    uint32_t fRpc;
    uint32_t fS2;
    static constexpr int NFOOTS = 16;
    std::vector<uint32_t> fFoot;

  public:
    ClassDef(R3BSyncCheckData, 1)
};

#endif /* R3BSyncCheckData_H */
