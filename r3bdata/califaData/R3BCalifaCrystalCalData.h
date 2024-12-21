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

#include <TObject.h>
#include <cstdint>
#include <iostream>
#include <string>

class R3BCalifaCrystalCalData : public TObject
{
  public:
    // Default constructor
    R3BCalifaCrystalCalData() = default;

    /** Constructor with arguments (explicit)
     *@param fCrystalId   Crystal unique identifier
     *@param fEnergy      Total energy deposited on the crystal ([MeV] in sim)
     *@param fNf  	  Total Nf (fast) ([MeV] in sim)
     *@param fNs	  Total Ns (slow) ([MeV] in sim)
     *@param fTime        Time since event start [ns]
     *@param fToT_Energy  Total energy deposited on the crystal from ToT ([MeV] in sim)
     **/
    explicit R3BCalifaCrystalCalData(uint16_t cryid,
                                     double energy,
                                     double Nf,
                                     double Ns,
                                     ULong64_t time = 0,
                                     double tot_energy = 0);

    // Destructor virtual
    virtual ~R3BCalifaCrystalCalData() = default;

    // Copy constructor
    R3BCalifaCrystalCalData(const R3BCalifaCrystalCalData&);
    R3BCalifaCrystalCalData& operator=(const R3BCalifaCrystalCalData&) { return *this; }

    // Accessors with [[nodiscard]]
    [[nodiscard]] inline const uint16_t& GetCrystalId() const { return fCrystalId; }
    [[nodiscard]] inline const double& GetEnergy() const { return fEnergy; }
    [[nodiscard]] inline const double& GetNf() const { return fNf; }
    [[nodiscard]] inline const double& GetNs() const { return fNs; }
    [[nodiscard]] inline const ULong64_t& GetTime() const { return fTime; }
    [[nodiscard]] inline const double& GetToTEnergy() const { return fToT_Energy; }

    // Modifiers
    inline void SetCrystalId(uint16_t cryid) { fCrystalId = cryid; }
    inline void SetEnergy(double energy) { fEnergy = energy; }
    inline void SetNf(double Nf) { fNf = Nf; }
    inline void SetNs(double Ns) { fNs = Ns; }
    inline void SetTime(ULong64_t time) { fTime = time; }
    inline void SetToTEnergy(double energy) { fToT_Energy = energy; }

    // Adding more values
    inline void AddMoreEnergy(double moreEnergy) { fEnergy += moreEnergy; }
    inline void AddMoreNf(double moreNf) { fNf += moreNf; }
    inline void AddMoreNs(double moreNs) { fNs += moreNs; }
    inline void AddMoreTot(double moreTot) { fToT_Energy += moreTot; }

    // Support for printing
    [[nodiscard]] std::string toString() const;
    void Print(const Option_t*) const override;

  protected:
    uint16_t fCrystalId = 0; // crystal unique identifier
    double fEnergy = 0;      // total energy in the crystal
    double fNf = 0;          // total Nf in the crystal
    double fNs = 0;          // total Nf in the crystal
    ULong64_t fTime = 0;     // time of the interaction
    double fToT_Energy = 0;  // total energy in the crystal from ToT

  public:
    ClassDefOverride(R3BCalifaCrystalCalData, 2)
};

// Operator overloading for printing R3BCalifaCrystalCalData
std::ostream& operator<<(std::ostream& os, const R3BCalifaCrystalCalData& data);
