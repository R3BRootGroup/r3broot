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

#ifndef R3BNEULANDPOINT_H
#define R3BNEULANDPOINT_H

#include "FairMCPoint.h"
#include "TObject.h"
#include "TVector3.h"
#include <iostream>

class R3BNeulandPoint : public FairMCPoint
{

  public:
    R3BNeulandPoint()
        : fLightYield{ 0 }
        , particle_id_{ 0 }
        , parent_particle_id_{ 0 }
    {
    }
    // NOLINTBEGIN
    R3BNeulandPoint(const Int_t trackID,
                    const Int_t detID,
                    const TVector3& pos,
                    const TVector3& mom,
                    const Double_t tof,
                    const Double_t length,
                    const Double_t eLoss,
                    const UInt_t EventId,
                    const Double_t lightYield,
                    const int particle_id,
                    const int parent_particle_id) // NOLINTEND
        : FairMCPoint{ trackID, detID, pos, mom, tof, length, eLoss, EventId }
        , fLightYield{ lightYield }
        , particle_id_{ particle_id }
        , parent_particle_id_{ parent_particle_id }
    {
    }

    R3BNeulandPoint(const FairMCPoint& point,
                    const Double_t lightYield,
                    const int particle_id,
                    const int parent_particle_id)
        : FairMCPoint{ point }
        , fLightYield{ lightYield }
        , particle_id_{ particle_id }
        , parent_particle_id_{ parent_particle_id }
    {
    }

    [[nodiscard]] auto GetMomentum() const -> TVector3;
    [[nodiscard]] auto GetPosition() const -> TVector3;
    [[nodiscard]] auto GetPaddle() const -> Int_t { return GetDetectorID(); }
    [[nodiscard]] auto GetLightYield() const -> Double_t { return fLightYield; }
    [[nodiscard]] auto GetPID() const -> int { return particle_id_; }
    ClassDefOverride(R3BNeulandPoint, 2);

    void Print(const Option_t* /*opt*/) const override;

  protected:
    Double_t fLightYield;

  private:
    int particle_id_;
    int parent_particle_id_;
};

auto operator<<(std::ostream&, const R3BNeulandPoint&) -> std::ostream&; // Support easy printing

#endif // R3BNEULANDPOINT_H
