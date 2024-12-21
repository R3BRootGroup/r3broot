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

#ifndef R3BBEAMPROPERTIES
#define R3BBEAMPROPERTIES

#include "R3BDistribution.h"
#include "R3BDistribution1D.h"

class R3BBeamProperties
{
  public:
    R3BBeamProperties();

    // Vertex distribution in cm
    R3BDistribution<3>& GetVertexDistribution() { return fVertex_cm; }
    void SetVertexDistribution(const R3BDistribution<3>& vertexDistribution_x_y_z_cm)
    {
        fVertex_cm = vertexDistribution_x_y_z_cm;
    }

    // spread in mRad
    R3BDistribution<2>& GetSpreadDistribution() { return fSpread_mRad; }
    void SetSpreadDistribution(const R3BDistribution<2>& spreadDistribution_theta_phi_mRad)
    {
        fSpread_mRad = spreadDistribution_theta_phi_mRad;
    }
    void SetSpreadDistribution(const R3BDistribution<1>& spreadDistribution_theta_mRad)
    {
        fSpread_mRad =
            R3BDistribution<2>({ spreadDistribution_theta_mRad, R3BDistribution1D::Flat(0, 2 * TMath::Pi() * 1e3) });
    }

    // beta distribution
    R3BDistribution<1>& GetBetaDistribution() { return fBeta; }
    void SetBetaDistribution(const R3BDistribution<1>& betaDistribution) { fBeta = betaDistribution; }

    // This method will calculate the beam beta using the given nucleon mass or approximate it from the number of
    // nucleons
    void SetEnergyDistribution(R3BDistribution<1> energyDistribution_AMeV, int nNucleons = 1, double mass_MeV = 0.);

  private:
    R3BDistribution<3> fVertex_cm;   //!
    R3BDistribution<2> fSpread_mRad; //!
    R3BDistribution<1> fBeta;        //!
};

#endif
