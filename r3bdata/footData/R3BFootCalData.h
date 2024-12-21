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

// --------------------------------------------------------------------
// -----              R3BFootCalData header file                  -----
// -----       Created 16/07/21 by J.L. Rodriguez-Sanchez         -----
// --------------------------------------------------------------------

#pragma once

#include <TObject.h>
#include <cmath>
#include <cstdint>

class R3BFootCalData : public TObject
{
  public:
    /** Default constructor **/
    R3BFootCalData() = default;

    /** Constructor with arguments
     *@param fDetId    Detector unique identifier
     *@param fStripId  Strip unique identifier
     *@param fEnergy   Total energy deposited on the strip ([GeV] in sim)
     **/
    explicit R3BFootCalData(uint8_t detid, uint16_t stripid, double energy);

    /** Destructor **/
    virtual ~R3BFootCalData() = default;

    // Accessors
    [[nodiscard]] inline const uint8_t& GetDetId() const { return fDetId; }
    [[nodiscard]] inline const uint16_t& GetStripId() const { return fStripId; }
    [[nodiscard]] inline const double& GetEnergy() const { return fEnergy; }

    // Modifiers
    inline void SetDetId(uint8_t detid) { fDetId = detid; }
    inline void SetStripId(uint16_t stripid) { fStripId = stripid; }
    inline void SetEnergy(double energy) { fEnergy = energy; }

  protected:
    uint8_t fDetId = 0;            // detector unique identifier
    uint16_t fStripId = 0;         // strip unique identifier
    double fEnergy = std::nan(""); // total energy in the strip

  public:
    ClassDefOverride(R3BFootCalData, 1)
};
