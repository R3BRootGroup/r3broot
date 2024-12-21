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

#ifndef NEULAND_ELASTICPROTONENERGY
#define NEULAND_ELASTICPROTONENERGY

#include "R3BNeulandCluster.h"

namespace Neuland
{
    Double_t RecoilProtonEnergy(const R3BNeulandCluster*);
    Double_t RecoilScatteringAngle(const R3BNeulandCluster*);
    Double_t ScatteredNeutronEnergy(const R3BNeulandCluster* first, const R3BNeulandCluster* second);
    Double_t ScatteredNeutronAngle(const R3BNeulandCluster* first, const R3BNeulandCluster* second);
    Double_t NeutronEnergyFromElasticProtonScattering(const R3BNeulandCluster*);
    Double_t NeutronEnergyFromElasticScattering(const R3BNeulandCluster* first,
                                                const R3BNeulandCluster* second,
                                                const Double_t targetMass);
    Double_t MaybeElasticScattering(const R3BNeulandCluster* first,
                                    const R3BNeulandCluster* second,
                                    const Double_t targetMass);
    Double_t ElasticScatteringTargetMass(const R3BNeulandCluster* first, const R3BNeulandCluster* second);
} // namespace Neuland

#endif // NEULAND_ELASTICPROTONENERGY
