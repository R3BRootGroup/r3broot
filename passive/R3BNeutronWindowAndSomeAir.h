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

#ifndef R3BROOT_R3BNEUTRONWINDOWANDSOMEAIR_H
#define R3BROOT_R3BNEUTRONWINDOWANDSOMEAIR_H

/* Place a neutron window (4mm steel) at xstart and air between start and stop
 * This can be used for NeuLAND simulations as there is currently no vacuum chamber implemented.
 * Instead, invert the problem: Use vacuum in the Cave and place a bunch o' air. Presto!
 * Bonus: No need to adjust the magnetic field for different beam energies.
 * The heavy ion will always fly through vacuum, and the neutrons through the right amount of material
 */

#include "FairModule.h"

class TGeoMedium;

class R3BNeutronWindowAndSomeAir : public FairModule
{
  public:
    explicit R3BNeutronWindowAndSomeAir(double start = 700, double stop = 1300);

    void ConstructGeometry() override { ConstructRootGeometry(); }
    void ConstructRootGeometry(TGeoMatrix* _ = nullptr) override;

  private:
    TGeoMedium* FindMaterial(const std::string& mat) const;

    double fStart;
    double fStop;

    ClassDefOverride(R3BNeutronWindowAndSomeAir, 1)
};

#endif // R3BROOT_R3BNEUTRONWINDOWANDSOMEAIR_H
