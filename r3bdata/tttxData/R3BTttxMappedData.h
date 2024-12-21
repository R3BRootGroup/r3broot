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

// -------------------------------------------------------------------
// -----           R3BTttxMappedData header file                 -----
// -----       Created 01/12/2023 by M. Whitehead                -----
// -------------------------------------------------------------------

#pragma once

#include <TObject.h>
#include <cstdint>
#include <iostream>
#include <string>

class R3BTttxMappedData : public TObject
{
  public:
    // Default constructor
    R3BTttxMappedData() = default;

    // Constructor with arguments (explicit)
    explicit R3BTttxMappedData(uint8_t detID,
                               uint8_t stripID,
                               int32_t time,
                               int32_t energy,
                               bool pileup,
                               bool overflow);

    // Destructor virtual
    virtual ~R3BTttxMappedData() = default;

    // Accessors with [[nodiscard]]
    [[nodiscard]] inline const uint8_t& GetDetID() const { return fDetID; }
    [[nodiscard]] inline const uint8_t& GetStripID() const { return fStripID; }
    [[nodiscard]] inline const int32_t& GetTime() const { return fTime; }
    [[nodiscard]] inline const int32_t& GetEnergy() const { return fEnergy; }
    [[nodiscard]] inline const bool& GetPileupStatus() const { return fPileup; }
    [[nodiscard]] inline const bool& GetOverflowStatus() const { return fOverflow; }

    // Modifiers
    inline void SetDetID(uint8_t id) { fDetID = id; }
    inline void SetStripID(uint8_t id) { fStripID = id; }
    inline void SetTime(int32_t time) { fTime = time; }
    inline void SetEnergy(int32_t energy) { fEnergy = energy; }
    inline void SetPileup(bool pu) { fPileup = pu; }
    inline void SetOverflow(bool ov) { fOverflow = ov; }

    // Support for printing
    [[nodiscard]] std::string toString() const;
    void Print(const Option_t*) const override;

  protected:
    uint8_t fDetID = 0;
    uint8_t fStripID = 0;
    int32_t fTime = 0;
    int32_t fEnergy = 0;
    bool fPileup = false;
    bool fOverflow = false;

  public:
    ClassDefOverride(R3BTttxMappedData, 1)
};

// Operator overloading for printing R3BTttxMappedData
std::ostream& operator<<(std::ostream& os, const R3BTttxMappedData& data);
