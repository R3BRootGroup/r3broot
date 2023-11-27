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

#ifndef R3BPSPXREADER_H
#define R3BPSPXREADER_H

#include "R3BReader.h"
#include <Rtypes.h>

struct EXT_STR_h101_PSP_t;
typedef struct EXT_STR_h101_PSP_t EXT_STR_h101_PSP;
typedef struct EXT_STR_h101_PSP_onion_t EXT_STR_h101_PSP_onion;
class TClonesArray;

/**
 * Class to unpack (with ucesb) to Mapped data for PSPX detector data.
 * This includes: Checking for error messages.
 * @author Ralf Plag (?), Bastian Loeher(?), Ina Syndikus
 * Modified by M. Holl, Dec 2019
 */

class R3BPspxReader : public R3BReader
{
  public:
    // Standard Constructor
    R3BPspxReader(EXT_STR_h101_PSP*, size_t);

    // Destructor
    virtual ~R3BPspxReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    EXT_STR_h101_PSP* fData; /**< Reader specific data structure from ucesb */
    size_t fOffset;          /**< Data Offset */
    // Don't store data for online
    Bool_t fOnline;
    std::vector<TClonesArray*> fMappedItems; /**< Array holding output (Mapped) data */

  public:
    ClassDefOverride(R3BPspxReader, 4);
};

#endif
