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

#include "ElasticScattering.h"
#include "TVector3.h"
#include <cmath>

using std::pow;
using std::sqrt;

namespace Neuland
{
    Double_t RecoilProtonEnergy(const R3BNeulandCluster* cluster)
    {
        // Range of the proton in material is proportional to its energy E(R) = aR^b
        // Here, the "energy moment" is used as a more robust representation of the neutron track length
        // Values for a and b fitted to simulated data by shooting protons at NeuLAND
        const Double_t a = 55.629;
        const Double_t b = 0.652103;
        return a * std::pow(cluster->GetEnergyMoment(), b);
    }

    Double_t RecoilScatteringAngle(const R3BNeulandCluster* cluster)
    {
        const TVector3 pNUnit = cluster->GetFirstHit().GetPosition().Unit();
        const TVector3 pp_Unit = (cluster->GetEnergyCentroid() - cluster->GetFirstHit().GetPosition()).Unit();
        const Double_t cosTheta = pNUnit.Dot(pp_Unit);
        return cosTheta;
    }

    Double_t ScatteredNeutronEnergy(const R3BNeulandCluster* first, const R3BNeulandCluster* second)
    {
        const Double_t E0n = 938.; // Rest Mass Neutron [MeV]
        const TVector3 pN_ = second->GetFirstHit().GetPosition() - first->GetFirstHit().GetPosition();
        const Double_t t = second->GetT() - first->GetT();

        const Double_t v2 = pN_.Mag2() / std::pow(t, 2); // cm²/ns²
        const Double_t c2 = 898.75517873681758374;       // cm²/ns²
        if (v2 > c2 || v2 > first->GetFirstHit().GetPosition().Mag2() / std::pow(first->GetT(), 2))
        {
            return 0;
        }

        const Double_t gamma = 1. / std::sqrt(1. - (v2 / c2));
        const Double_t En_ = (gamma - 1.) * E0n;
        return En_;
    }

    Double_t ScatteredNeutronAngle(const R3BNeulandCluster* first, const R3BNeulandCluster* second)
    {
        const TVector3 pNUnit = first->GetFirstHit().GetPosition().Unit();
        const TVector3 pN_Unit = (second->GetFirstHit().GetPosition() - first->GetFirstHit().GetPosition()).Unit();
        const Double_t cosTheta = pNUnit.Dot(pN_Unit); // \cos(\theta_{NN'})
        return cosTheta;
    }

    Double_t NeutronEnergyFromElasticProtonScattering(const R3BNeulandCluster* cluster)
    {
        const Double_t minimalProtonKineticEnergy = 0.; // MeV
        // E0n:       Rest mass proton/neutron
        // En:       Total energy of incoming neutron
        // Ep_:      Total energy of scattered proton
        // Ep_Kin:   Kinetic energy of scattered proton
        // pp_:      Momentum Vector of scattered proton
        // pn:       Momentum Vector of incoming neutron
        // cosTheta: Angle between incoming neutron and scattered proton

        // Use the Energy Moment of the cluster to determine the proton kinetic energy.
        const Double_t Ep_Kin = RecoilProtonEnergy(cluster);
        if (Ep_Kin < minimalProtonKineticEnergy)
        {
            // If the kinetic energy of the proton does not match, this method doesn't work
            return 0;
        }

        // Rest Mass of Proton/Neutron [MeV]
        const Double_t E0n = 938.; // MeV
        const Double_t Ep_ = E0n + Ep_Kin;

        const Double_t cosTheta = RecoilScatteringAngle(cluster);

        /*if (pp_Unit.Z() < 0. || cosTheta < 0.)
        {
            // No scattering in backward directions
            // No scattering angle > 90°
            return 0;
        }*/

        // From Four-Momentum conservation one can deduce:
        // \frac{E_p' + E0n}{E_p' - E0n} \cos^2(\theta_{np'}) = \frac{E_N + E0n}{E_N - E0n}
        const Double_t a = cosTheta * cosTheta * (Ep_ + E0n) / (Ep_ - E0n);
        const Double_t En = E0n * ((a + 1.) / (a - 1.));

        // return EKin
        return En - E0n;
    }

    Double_t NeutronEnergyFromElasticScattering(const R3BNeulandCluster* first,
                                                const R3BNeulandCluster* second,
                                                const Double_t targetMass)
    {
        const Double_t E0n = 938.;       // Rest Mass Neutron [MeV]
        const Double_t E0k = targetMass; // Rest Mass heavy nucleus [MeV], e.g. 11000 for Carbon

        const Double_t En_ = ScatteredNeutronEnergy(first, second) + E0n;
        // \cos(\theta_{NN'})
        const Double_t C = ScatteredNeutronAngle(first, second);

        // Formula for En. Yes its long, no the pows do not make it slow
        const Double_t En =
            (sqrt(pow(C, 4) * pow(En_, 4) * pow(E0n, 2) - 2 * pow(C, 4) * pow(En_, 2) * pow(E0n, 4) +
                  pow(C, 4) * pow(E0n, 6) + pow(C, 2) * pow(En_, 4) * pow(E0k, 2) -
                  pow(C, 2) * pow(En_, 4) * pow(E0n, 2) - 2 * pow(C, 2) * pow(En_, 2) * pow(E0k, 2) * pow(E0n, 2) +
                  2 * pow(C, 2) * pow(En_, 2) * pow(E0n, 4) + pow(C, 2) * pow(E0k, 2) * pow(E0n, 4) -
                  pow(C, 2) * pow(E0n, 6)) +
             pow(En_, 2) * E0k - En_ * pow(E0k, 2) - En_ * pow(E0n, 2) + E0k * pow(E0n, 2)) /
            (pow(C, 2) * pow(En_, 2) - pow(C, 2) * pow(E0n, 2) - pow(En_, 2) + 2 * En_ * E0k - pow(E0k, 2));

        return En - E0n;
    }

    Double_t MaybeElasticScattering(const R3BNeulandCluster* first,
                                    const R3BNeulandCluster* second,
                                    const Double_t targetMass)
    {
        const Double_t E0n = 938.;
        const Double_t E0k = targetMass;

        const Double_t En = first->GetFirstHit().GetEToF() + E0n;
        const Double_t En_ = ScatteredNeutronEnergy(first, second) + E0n;
        const Double_t Ek_ = first->GetE() + targetMass;
        const Double_t cosTheta = ScatteredNeutronAngle(first, second);

        const Double_t zero = En_ * En_ + En_ * Ek_ - En * E0k - E0n * E0n -
                              sqrt((En_ * En_ - E0n * E0n) * (En * En - E0n * E0n)) * cosTheta;
        return zero;
    }

    Double_t ElasticScatteringTargetMass(const R3BNeulandCluster* first, const R3BNeulandCluster* second)
    {
        const Double_t E0n = 938.;

        const Double_t En = first->GetFirstHit().GetEToF() + E0n;
        const Double_t En_ = ScatteredNeutronEnergy(first, second) + E0n;
        const Double_t Ek_kin = first->GetE();
        const Double_t cosTheta = ScatteredNeutronAngle(first, second);

        const Double_t E0k =
            (En_ * En_ - E0n * E0n - sqrt((En_ * En_ - E0n * E0n) * (En * En - E0n * E0n)) * cosTheta + En_ * Ek_kin) /
            (En - En_);

        return E0k;
    }
} // namespace Neuland
