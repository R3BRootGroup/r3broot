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

// -----------------------------------------------------------------
// -----                   R3BFootMappedData                   -----
// -----          Created 16/07/2021 by J.L. Rodriguez         -----
// -----------------------------------------------------------------

#pragma once

#include <TObject.h>

class R3BFootMappedData : public TObject
{
  public:
    // Default Constructor
    R3BFootMappedData() = default;

    /** Standard Constructor
     *@param detId       Detector unique identifier
     *@param stripId     Strip unique identifier
     *@param energy      Total energy deposited in the strip
     **/
    explicit R3BFootMappedData(uint8_t detId, uint16_t stripId, int energy);

    // Destructor virtual
    virtual ~R3BFootMappedData() = default;

    // Getters
    [[nodiscard]] inline const uint8_t& GetDetId() const { return fDetId; }
    [[nodiscard]] inline const uint16_t& GetStripId() const { return fStripId; }
    [[nodiscard]] inline const int& GetEnergy() const { return fEnergy; }

  protected:
    uint8_t fDetId = 0;    // Detector unique identifiers
    uint16_t fStripId = 0; // Strip unique identifiers
    int fEnergy = 0;       // Energy per strip

  public:
    ClassDefOverride(R3BFootMappedData, 1)
};
