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

// ------------------------------------------------------------------
// -----            R3BFootHitData header file                  -----
// -----       Created 16/07/21  by J.L. Rodriguez-Sanchez      -----
// ------------------------------------------------------------------

#pragma once

#include <TObject.h>
#include <TVector3.h>
#include <cmath>

class R3BFootHitData : public TObject
{
  public:
    /** Default constructor **/
    R3BFootHitData() = default;

    /** Constructor with arguments
     *@param fDetId     Detector unique identifier
     *@param fNbHit     Number of hits
     *@param fMulStrip  Number of strips with signals
     *@param fPos         Position in mm, local coordinate (detector frame)
     *@param fTheta     Master: Angle theta [rad] (lab frame)
     *@param fPhi       Master: Angle Phi [rad] (lab frame)
     *@param fEnergy    Total energy deposited by the hit ([GeV] in sim)
     **/
    explicit R3BFootHitData(uint8_t detid,
                            uint16_t nbhit,
                            double pos,
                            TVector3 master,
                            double energy,
                            uint16_t mulstrip = 0);

    /** Copy constructor **/
    R3BFootHitData(const R3BFootHitData&);

    R3BFootHitData& operator=(const R3BFootHitData&) { return *this; }

    /** Destructor **/
    virtual ~R3BFootHitData() = default;

    // Accessors
    [[nodiscard]] inline const uint8_t& GetDetId() const { return fDetId; }
    [[nodiscard]] inline const uint16_t& GetNbHit() const { return fNbHit; }
    [[nodiscard]] inline const uint16_t& GetMulStrip() const { return fMulStrip; }
    [[nodiscard]] inline const double& GetPos() const { return fPos; }
    [[nodiscard]] inline const double& GetTheta() const { return fTheta; }
    [[nodiscard]] inline const double& GetPhi() const { return fPhi; }
    [[nodiscard]] inline const TVector3 GetPosLab() const { return fmaster; }
    [[nodiscard]] inline const double& GetEnergy() const { return fEnergy; }

  protected:
    uint8_t fDetId = 0;
    uint16_t fNbHit = 0;
    double fPos = std::nan("");
    double fTheta = std::nan("");
    double fPhi = std::nan("");
    TVector3 fmaster;
    double fEnergy = std::nan("");
    uint16_t fMulStrip = 0;

  public:
    ClassDefOverride(R3BFootHitData, 1)
};
