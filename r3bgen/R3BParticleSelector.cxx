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

#include "R3BParticleSelector.h"

#include "FairLogger.h"
#include "FairParticle.h"
#include "FairRunSim.h"

#include "TDatabasePDG.h"

#include "G4NistManager.hh"

constexpr auto PDG_Proton = 2212;
constexpr auto PDG_Neutron = 2112;

// C++17 will make the following lines much better, but for now we have to stick with this...
// maybe not very good to read, but its at least constexpr
constexpr int ChargeMass2PDGCode(int z, int a)
{
    return (z > 0 ? (1000000000 + z * 10000 + a * 10) : -(1000000000 - z * 10000 + a * 10));
}
constexpr int PDGCode2Mass(int pdg)
{
    return (pdg > 1000000000 ? ((pdg / 10) % 1000) : (pdg < 1000000000 ? ((-pdg / 10) % 1000) : 0));
}
constexpr int PDGCode2Charge(int pdg)
{
    return (pdg > 1000000000 ? ((pdg / 10000) % 1000) : (pdg < 1000000000 ? -((-pdg / 10000) % 1000) : 0));
}
constexpr std::array<int, 2> PDGCode2ChargeMass(int pdg) { return { PDGCode2Charge(pdg), PDGCode2Mass(pdg) }; }

void R3BParticleSelector::AddProton() { AddParticle(PDG_Proton); }
void R3BParticleSelector::AddNeutron() { AddParticle(PDG_Neutron); }

void R3BParticleSelector::AddParticle(FairParticle* particle) { addParticle(particle->GetPDG(), particle->GetMass()); }

void R3BParticleSelector::AddParticle(const char* particleName)
{
    auto particle =
        dynamic_cast<FairParticle*>(FairRunSim::Instance()->GetUserDefParticles()->FindObject(particleName));
    if (!particle)
        LOG(fatal) << "R3BParticleSelector:  Particle is not defined!";

    addParticle(particle->GetPDG(), particle->GetMass());
}

void R3BParticleSelector::AddParticle(const char* particleName,
                                      int z,
                                      int a,
                                      int s,
                                      double mass_GeV,
                                      int charge,
                                      bool stable,
                                      double decaytime_ns)
{
    auto particle =
        dynamic_cast<FairParticle*>(FairRunSim::Instance()->GetUserDefParticles()->FindObject(particleName));
    if (!particle)
    {
        particle = new FairParticle(particleName, z, a, s, mass_GeV, charge, stable, decaytime_ns);
        FairRunSim::Instance()->AddNewParticle(particle);
    }

    addParticle(particle->GetPDG(), mass_GeV);
}

void R3BParticleSelector::AddParticle(int pdgCode)
{
    if (pdgCode > 1000000000)
    {
        // we will use G4NistManager because it knows more exotic isotopes
        addParticle(
            pdgCode,
            G4NistManager::Instance()->GetIsotopeMass(PDGCode2Charge(pdgCode), PDGCode2Mass(pdgCode)) / CLHEP::GeV);
        return;
    }

    // for other particles TDatabasePDG should be fine
    const auto particle = TDatabasePDG::Instance()->GetParticle(pdgCode);
    if (!particle)
        LOG(fatal) << "R3BParticleSelector: Unknown PDG Code!";

    addParticle(pdgCode, particle->Mass());
}

void R3BParticleSelector::AddParticle(int z, int a)
{
    if (a <= 0)
        LOG(fatal) << "R3BParticleSelector::AddParticle: A <= 0 given.";

    if (a == 1)
    {
        switch (z)
        {
            case (1):
                addParticle(PDG_Proton, G4NistManager::Instance()->GetIsotopeMass(z, a) / CLHEP::GeV);
                return;
            case (0):
                addParticle(PDG_Neutron, TDatabasePDG::Instance()->GetParticle(PDG_Neutron)->Mass());
                return;
            case (-1):
                addParticle(-PDG_Proton, G4NistManager::Instance()->GetIsotopeMass(-z, a) / CLHEP::GeV);
                return;

            default:
                break;
        }
    }

    const auto mass = G4NistManager::Instance()->GetIsotopeMass(z, a) / CLHEP::GeV;

    if (mass == 0.)
    {
        LOG(fatal) << "R3BParticleSelector::AddParticle: Could not get the mass of particles with A = " << a
                   << " and Z = " << z << ". Define the particle yourself! ";
    }

    const auto pdgCode = ChargeMass2PDGCode(z, a);
    const auto particleName = TString::Format("particle_%d", pdgCode);

    auto particle =
        dynamic_cast<FairParticle*>(FairRunSim::Instance()->GetUserDefParticles()->FindObject(particleName));

    if (!particle)
    {
        particle = new FairParticle(pdgCode, particleName, kPTIon, mass, z, -1);
        FairRunSim::Instance()->AddNewParticle(particle);
    }

    addParticle(pdgCode, mass);
}
