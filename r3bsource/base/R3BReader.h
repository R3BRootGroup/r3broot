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

#include <R3BUcesbMappingFlag.h>
#include <TNamed.h>
#include <TString.h>

#include <ext_data_struct_info.hh>

class TClonesArray;
class R3BEventHeader;

class R3BReader : public TNamed
{
  public:
    R3BReader(TString const&);
    virtual ~R3BReader();
    [[nodiscard]] virtual bool MismappedItemRequired(std::string_view /*item_name*/) const { return false; }
    bool AllowExtraCondition(R3B::UcesbMap map_flag, R3B::UcesbMap success_condition)
    {
        return (map_flag & ~(success_condition | extra_conditions_)) == R3B::UcesbMap::zero;
    }

    void SetExtraConditions(R3B::UcesbMap conditions) { extra_conditions_ = conditions; }
    void AddExtraConditions(R3B::UcesbMap conditions) { extra_conditions_ |= conditions; }

    /* Setup structure information */
    virtual Bool_t Init(ext_data_struct_info*) = 0;
    virtual void SetParContainers() {}
    virtual Bool_t ReInit() { return kTRUE; }
    /* Read data from full event structure */
    virtual Bool_t R3BRead() = 0;
    /* Reset */
    virtual void Reset() = 0;
    /* Return actual name of the reader */

    // actions when closed
    virtual void Close(){};

  private:
    R3B::UcesbMap extra_conditions_ = R3B::UcesbMap::zero;

  public:
    ClassDef(R3BReader, 0);
};
