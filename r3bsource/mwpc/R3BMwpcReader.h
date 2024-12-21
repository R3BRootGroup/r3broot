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

#pragma once

#include "R3BReader.h"
#include <Rtypes.h>
#include <TClonesArray.h>

struct EXT_STR_h101_SOFMWPC_t;
typedef struct EXT_STR_h101_SOFMWPC_t EXT_STR_h101_SOFMWPC;

class R3BMwpcReader : public R3BReader
{
  public:
    // Standard constructor
    R3BMwpcReader(EXT_STR_h101_SOFMWPC*, size_t);

    // Destructor
    virtual ~R3BMwpcReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select the max. number of mwpcs
    inline void SetMaxNbDet(UInt_t d) { fMaxDet = d; }

    // Accessor to select online mode
    inline void SetOnline(Bool_t option) { fOnline = option; }

  private:
    // Reader specific data structure from ucesb
    EXT_STR_h101_SOFMWPC* fData;
    // Max. number of mwpcs
    UInt_t fMaxDet;
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    // R3BMwpcMapped Item
    TClonesArray* fArrayMwpc0; /**< Output array. */
    TClonesArray* fArrayMwpc1; /**< Output array. */
    TClonesArray* fArrayMwpc2; /**< Output array. */
    TClonesArray* fArrayMwpc3; /**< Output array. */

  public:
    ClassDefOverride(R3BMwpcReader, 0);
};
