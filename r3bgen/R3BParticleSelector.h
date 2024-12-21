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

#ifndef R3BPARTICLESELECTOR_H
#define R3BPARTICLESELECTOR_H

#include "Rtypes.h"

#include <array>

class TParticle;
class FairParticle;

class R3BParticleSelector
{
  public:
    virtual ~R3BParticleSelector() = default;

    void AddProton();
    void AddNeutron();
    void AddParticle(FairParticle* particle);
    void AddParticle(const char* particleName);
    void AddParticle(const char* particleName,
                     int z,
                     int a,
                     int s,
                     double mass_GeV,
                     int charge,
                     bool stable,
                     double decaytime_ns);
    void AddParticle(int pdgCode);
    void AddParticle(int z, int a);

  protected:
    virtual void addParticle(const int pdgCode, const double mass) = 0;
};

#endif
