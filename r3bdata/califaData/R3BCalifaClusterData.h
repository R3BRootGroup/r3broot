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

class R3BCalifaClusterData : public TObject
{
  public:
    // Default constructor
    R3BCalifaClusterData() = default;

    // Constructor with arguments  (explicit)
    explicit R3BCalifaClusterData(std::vector<uint16_t> crystalList,
                                  double ene,
                                  double nf,
                                  double ns,
                                  double theta,
                                  double phi,
                                  ULong64_t time,
                                  uint8_t clusterType);

    // Copy constructor
    R3BCalifaClusterData(const R3BCalifaClusterData&);

    R3BCalifaClusterData& operator=(const R3BCalifaClusterData&) { return *this; }

    // Destructor
    virtual ~R3BCalifaClusterData() = default;

    // Accessors with [[nodiscard]]
    [[nodiscard]] inline const std::vector<uint16_t> GetCrystalList() const { return fCrystalList; }
    [[nodiscard]] inline const double& GetEnergy() const { return fEnergy; }
    [[nodiscard]] inline const double& GetNf() const { return fNf; }
    [[nodiscard]] inline const double& GetNs() const { return fNs; }
    [[nodiscard]] inline const double& GetTheta() const { return fTheta; }
    [[nodiscard]] inline const double& GetPhi() const { return fPhi; }
    [[nodiscard]] inline const ULong64_t& GetTime() const { return fTime; }
    [[nodiscard]] inline const uint8_t& GetClusterType() const { return fClusterType; }
    [[nodiscard]] inline const int GetNbOfCrystalHits() const { return fCrystalList.size(); }
    [[nodiscard]] inline const int GetMotherCrystal() const { return fCrystalList.at(0); }

    // Modifiers
    inline void SetCrystalList(std::vector<uint16_t> list) { fCrystalList = std::move(list); }
    inline void SetEnergy(double ene) { fEnergy = ene; }
    inline void SetNf(double nf) { fNf = nf; }
    inline void SetNs(double ns) { fNs = ns; }
    inline void SetTheta(double theta) { fTheta = theta; }
    inline void SetPhi(double phi) { fPhi = phi; }
    inline void SetTime(ULong64_t time) { fTime = time; }
    inline void SetClusterType(uint8_t id) { fClusterType = id; }

    // Support for printing
    [[nodiscard]] std::string toString() const;
    void Print(const Option_t*) const override;

  protected:
    std::vector<uint16_t> fCrystalList; // List with crystals inside the cluster. First one is the mother crystal
    double fEnergy = 0;                 // total energy deposited
    double fNf = 0;                     // total Nf deposited
    double fNs = 0;                     // total Ns deposited
    double fTheta = 0;                  // reconstructed theta
    double fPhi = 0;                    // reconstructed phi
    ULong64_t fTime = 0;                // WR time stamp
    uint8_t fClusterType = 0;           // Cluster type: 0 proton, 1 gamma, 2 saturation

  public:
    ClassDefOverride(R3BCalifaClusterData, 4)
};

// Operator overloading for printing R3BCalifaClusterData
std::ostream& operator<<(std::ostream& os, const R3BCalifaClusterData& data);
