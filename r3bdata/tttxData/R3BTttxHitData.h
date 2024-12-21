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
// -----               R3BTttxHitData header file                -----
// -----      Created 11/01/2024 by J. L. Rodriguez-Sanchez      -----
// -------------------------------------------------------------------

#pragma once

#include <TObject.h>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

class R3BTttxHitData : public TObject
{
  public:
    // Default constructor
    R3BTttxHitData() = default;

    // Constructor with arguments (explicit)
    explicit R3BTttxHitData(int8_t idet, double xpos, double energy, double ang, double charge, double time = 0);

    // Destructor virtual
    virtual ~R3BTttxHitData() = default;

    // Accessors with [[nodiscard]]
    [[nodiscard]] inline const int8_t& GetDetID() const { return fDetID; }
    [[nodiscard]] inline const double& GetX() const { return fXpos; }
    [[nodiscard]] inline const double& GetEnergy() const { return fEnergy; }
    [[nodiscard]] inline const double& GetAngle() const { return fAng; }
    [[nodiscard]] inline const double& GetChargeZ() const { return fCharge; }
    [[nodiscard]] inline const double& GetTime() const { return fTime; }

    // Modifiers
    inline void SetDetID(int8_t idet) { fDetID = idet; }
    inline void SetXpos(double xpos) { fXpos = xpos; }
    inline void SetEnergy(double energy) { fEnergy = energy; }
    inline void SetAngle(double angle) { fAng = angle; }
    inline void SetChargeZ(double charge) { fCharge = charge; }
    inline void SetTime(double time) { fTime = time; }

    // Support for printing
    [[nodiscard]] std::string toString() const;
    void Print(const Option_t*) const override;

  protected:
    int8_t fDetID = -1; // 0: Reconstructed with all detectors, 1-: ID
    double fXpos = std::nan("");
    double fEnergy = std::nan("");
    double fAng = std::nan("");
    double fCharge = std::nan("");
    double fTime = std::nan("");

  public:
    ClassDefOverride(R3BTttxHitData, 1)
};

// Operator overloading for printing R3BTttxHitData
std::ostream& operator<<(std::ostream& os, const R3BTttxHitData& data);
