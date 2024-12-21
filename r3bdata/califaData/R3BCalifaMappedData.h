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

class R3BCalifaMappedData : public TObject
{
  public:
    // Default Constructor
    R3BCalifaMappedData() = default;

    /** Standard Constructor
     *@param crystalId   Crystal unique identifier
     *@param energy      Total energy deposited in the crystal [GeV]
     *@param nf          Total fast amplitude deposited in the crystal [a.u.]
     *@param ns          Total slow amplitude deposited on the crystal [a.u.]
     *@param febextime   Internal febex time per crystal [ns]
     *@param wrts        Timestamp per crystal or time since event start in simulation [ns]
     *@param ov          Overflow bits
     *@param pu          Pileup bits
     *@param dc          Discard bits
     *@param tot         Time-over-threshold
     **/
    explicit R3BCalifaMappedData(UShort_t crystalId,
                                 int16_t energy,
                                 int16_t nf,
                                 int16_t ns,
                                 uint64_t febextime,
                                 uint64_t wrts,
                                 uint32_t ov,
                                 uint16_t pu,
                                 uint16_t dc,
                                 uint16_t tot);

    // Destructor
    virtual ~R3BCalifaMappedData() = default;

    // Accessors with [[nodiscard]]
    [[nodiscard]] inline const UShort_t& GetCrystalId() const { return fCrystalId; }
    [[nodiscard]] inline const int16_t& GetEnergy() const { return fEnergy; }
    [[nodiscard]] inline const int16_t& GetNf() const { return fNf; }
    [[nodiscard]] inline const int16_t& GetNs() const { return fNs; }
    [[nodiscard]] inline const uint64_t& GetFebexTime() const { return fFebexTime; }
    [[nodiscard]] inline const uint64_t& GetWrts() const { return fWrts; }
    [[nodiscard]] inline const uint32_t& GetOverFlow() const { return fOverFlow; }
    [[nodiscard]] inline const uint16_t& GetPileup() const { return fPileup; }
    [[nodiscard]] inline const uint16_t& GetDiscard() const { return fDiscard; }
    [[nodiscard]] inline const uint16_t& GetTot() const { return fTot; }

    // Support for printing
    [[nodiscard]] std::string toString() const;
    void Print(const Option_t*) const override;

  protected:
    UShort_t fCrystalId; // Crystal unique identifier
    int16_t fEnergy;     // Total energy in the crystal
    int16_t fNf;         // Total fast amplitude in the crystal
    int16_t fNs;         // Total slow amplitude in the crystal
    uint64_t fFebexTime; // Internal febex time
    uint64_t fWrts;      // Timestamp per crystal
    uint32_t fOverFlow;  // Overflow bits
    uint16_t fPileup;    // Pileup bits
    uint16_t fDiscard;   // Discard bits
    uint16_t fTot;       // Time-over-treshold

  public:
    ClassDefOverride(R3BCalifaMappedData, 4)
};

// Operator overloading for printing R3BCalifaMappedData
std::ostream& operator<<(std::ostream& os, const R3BCalifaMappedData& data);
