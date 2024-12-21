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

class TClonesArray;

/* The structure containing the data layout of RPC detector */
struct EXT_STR_h101_RPC_t;
typedef struct EXT_STR_h101_RPC_t EXT_STR_h101_RPC;
typedef struct EXT_STR_h101_RPC_onion_t EXT_STR_h101_RPC_onion;
class ext_data_struct_info;

class R3BRpcReader : public R3BReader
{
  public:
    /* Default constructor / destructor */
    R3BRpcReader(EXT_STR_h101_RPC*, size_t);
    virtual ~R3BRpcReader();

    /*
     * Initialise the detector specific data items
     * in the ext_data_struct_info list of items.
     * This is called by the R3BUcesbSource.
     */
    virtual Bool_t Init(ext_data_struct_info*) override;

    /*
     * Copy data from the ucesb data stream to
     * native data containers, called by R3BUcesbSource.
     */
    virtual Bool_t R3BRead() override;

    /*
     * Reset the internal output array.
     */
    virtual void Reset() override;

    // Accessor to select online mode
    inline void SetOnline(Bool_t option) { fOnline = option; }

  private:
    /* Reader specific data structure from ucesb */
    EXT_STR_h101_RPC* fData;
    // Data offset
    size_t fOffset;
    // Don't store data for online
    Bool_t fOnline;
    // An event counter
    unsigned int fNEvent;
    /*
     * An array where the read R3BPaddleTamexMappedData
     * items are stored.
     */
    TClonesArray* fArrayRpc; /**< Output array. */

    Int_t fCoarseReference; // earliest coarse counter to detect coarse counter wrap

  public:
    ClassDefOverride(R3BRpcReader, 0);
};
