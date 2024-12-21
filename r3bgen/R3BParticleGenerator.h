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

#ifndef R3BPARTICLEGENERATOR_H
#define R3BPARTICLEGENERATOR_H

#include "R3BBeamProperties.h"
#include "R3BDistribution.h"
#include "R3BDistribution1D.h"
#include "R3BParticleSelector.h"

#include "FairGenerator.h"

#include "TRandom3.h"

#include <iostream>

class FairPrimaryGenerator;

class R3BParticleGenerator : public FairGenerator, public R3BParticleSelector
{

  public:
    R3BParticleGenerator(int pdgCode, unsigned int seed = 0U);
    R3BParticleGenerator(unsigned int seed = 0U);

    R3BBeamProperties Beam; //!
    bool Init() override;
    bool ReadEvent(FairPrimaryGenerator* primGen) override;

  protected:
    void addParticle(const int pdgCode, const double mass) override;

  private:
    R3BDistribution<3> fVertex_cm;   //!
    R3BDistribution<2> fSpread_mRad; //!
    R3BDistribution<1> fEnergy_AMeV; //!

    int fPDGCode;     // The PDG Code of the Particle
    double fMass_GeV; // The Mass of the Particle
    TRandom3 fRNG;    // the RNG

    ClassDefOverride(R3BParticleGenerator, 1)
};

#endif
