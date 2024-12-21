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

#ifndef _R3BTOF_UNPACK_
#define _R3BTOF_UNPACK_

#include <map>

#include "FairUnpack.h"

class TH1F;

class R3BTofUnpack : public FairUnpack
{

  public:
    R3BTofUnpack(Short_t type, Short_t subType, Short_t procId, Short_t subCrate, Short_t control);
    virtual ~R3BTofUnpack();

    virtual Bool_t Init();
    virtual Bool_t DoUnpack(Int_t* data, Int_t size);
    virtual void Reset();

  private:
    std::map<Int_t, Int_t> fmap_leading_t;
    std::map<Int_t, Int_t> fmap_leading_c;
    std::map<Int_t, Int_t> fmap_trailing_t;
    std::map<Int_t, Int_t> fmap_trailing_c;

  protected:
    virtual void Register();

  public:
    ClassDef(R3BTofUnpack, 1)
};

#endif // _R3BTOF_UNPACK_