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

#ifndef R3BWhiterabbitFootReader_H
#define R3BWhiterabbitFootReader_H 1

#include "R3BReader.h"
#include <Rtypes.h>
#include <vector>

struct EXT_STR_h101_WRFOOT_t;
typedef struct EXT_STR_h101_WRFOOT_t EXT_STR_h101_WRFOOT;
typedef struct EXT_STR_h101_WRFOOT_onion_t EXT_STR_h101_WRFOOT_onion;

class TClonesArray;
class R3BEventHeader;
class ext_data_struct_info;

/**
 * A reader of FOOT white rabbit data with UCESB.
 * Receives mapped raw data and converts it to R3BRoot objects.
 * @author J.L. Rodriguez
 * @since Feb 24, 2022
 */

class R3BWhiterabbitFootReader : public R3BReader
{
  public:
    // Standard constructor
    R3BWhiterabbitFootReader(EXT_STR_h101_WRFOOT_onion*, size_t, std::vector<UInt_t>);

    // Destructor
    virtual ~R3BWhiterabbitFootReader();

    // Setup structure information
    virtual Bool_t Init(ext_data_struct_info*) override;

    // Read data from full event structure
    virtual Bool_t R3BRead() override;

    // Reset
    virtual void Reset() override;

    // Accessor to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    // An event counter
    UInt_t fNEvent;
    // Reader specific data structure from ucesb
    EXT_STR_h101_WRFOOT_onion* fData;
    // Offset of detector specific data in full data structure
    size_t fOffset;
    // The whiterabbit subsystem ID
    std::vector<UInt_t> fWhiterabbitId;
    // A pointer to the R3BEventHeader structure
    R3BEventHeader* fEventHeader;
    // Don't store data for online
    Bool_t fOnline;
    // Output array
    TClonesArray* fArray;

  public:
    ClassDefOverride(R3BWhiterabbitFootReader, 0);
};
#endif // R3BWhiterabbitFootReader_H
