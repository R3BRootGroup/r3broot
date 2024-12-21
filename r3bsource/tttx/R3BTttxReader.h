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

struct EXT_STR_h101_TTTX_t;
typedef struct EXT_STR_h101_TTTX_t EXT_STR_h101_TTTX;
typedef struct EXT_STR_h101_TTTX_onion_t EXT_STR_h101_TTTX_onion;
class ext_data_struct_info;

class R3BTttxReader : public R3BReader
{
  public:
    // Standard constructor
    R3BTttxReader(EXT_STR_h101_TTTX*, size_t);

    // Destructor
    virtual ~R3BTttxReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }

    // Accessor to select events without pileup
    void SetSkipPileup() { fPileup = kTRUE; }

    void SetNumDets(Int_t num) { fDets = num; }
    void SetNumChannels(Int_t num) { fChannels = num; }
    void SetNumTref(Int_t num) { fTref = num; }
    void SetNumTtrig(Int_t num) { fTtrig = num; }

  private:
    Bool_t ReadData(EXT_STR_h101_TTTX_onion*, UShort_t);

    /* Reader specific data structure from ucesb */
    EXT_STR_h101_TTTX* fData;
    /* Data offset */
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    /* Geometry */
    Int_t fDets;
    Int_t fChannels;
    Int_t fTref;
    Int_t fTtrig;
    Bool_t fPileup;
    TClonesArray* fArray; /**< Output array. */

  public:
    ClassDefOverride(R3BTttxReader, 0);
};
