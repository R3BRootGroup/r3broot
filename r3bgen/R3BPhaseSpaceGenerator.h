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

// Wrapper for TGenPhaseSpace

#include "FairGenerator.h"
#include "R3BBeamProperties.h"
#include "R3BDistribution.h"
#include "R3BParticleSelector.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PxPyPzM4D.h>
#include <Math/Vector4Dfwd.h>
#include <R3BIOConnector.h>

#include <set>
#include <vector>

struct R3BPhaseSpaceGenParticleInfo
{
    R3BPhaseSpaceGenParticleInfo() = default;
    int pdg_code = 0;
    double mass = 0.;                   // GeV
    double kinetic_energy = 0.;         // GeV. Defintion: E-m
    ROOT::Math::PxPyPzMVector momentum; // GeV
    ROOT::Math::XYZTVector position;    // cm

    ClassDefNV(R3BPhaseSpaceGenParticleInfo, 1);
};

class R3BPhaseSpaceGenerator : public FairGenerator, public R3BParticleSelector
{
  public:
    explicit R3BPhaseSpaceGenerator(unsigned int seed = 0U);

    // Modifiers:
    void SetParticleWhitelist(const std::set<int>& whitelist) { particle_whitelist_ = whitelist; }
    void AddParticleToWhitelist(int pdg_code) { particle_whitelist_.emplace(pdg_code); }
    void EnableWhitelist(bool is_enabled = true) { is_whitelist_enabled_ = is_enabled; }
    void EnableWrite(bool is_enabled = true);

    // realtive energy distribution in keV
    void SetErelDistribution(const R3BDistribution<1>& ErelDistribution) { fErel_keV = ErelDistribution; }

    // Getters:
    [[nodiscard]] auto GetBeam() const -> const R3BBeamProperties& { return beam_properties_; }
    [[nodiscard]] auto GetErelDistribution() const -> const R3BDistribution<1>& { return fErel_keV; }

    auto GetBeam() -> R3BBeamProperties& { return beam_properties_; }

  private:
    bool is_whitelist_enabled_ = false;
    bool is_written_enabled_ = false;
    double total_mass_ = 0.;
    R3BDistribution<1> fErel_keV; //!
    TRandom3 rnd_gen_;
    TGenPhaseSpace phase_space_gen_;
    R3BBeamProperties beam_properties_;
    std::vector<int> pdg_codes_;
    std::vector<double> masses;
    std::set<int> particle_whitelist_;

    R3B::OutputVectorConnector<R3BPhaseSpaceGenParticleInfo> particle_output_{
        "PhaseGenParticle"
    }; //! This won't be stored

    void addParticle(int pdgCode, double mass) override;
    auto Init() -> bool override;
    auto ReadEvent(FairPrimaryGenerator* primGen) -> bool override;

    ClassDefOverride(R3BPhaseSpaceGenerator, 3);
};
