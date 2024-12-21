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

/**
 * Gamma (cascade) event generator.
 *
 * @author Max Winkel <max.winkel@ph.tum.de>
 */

#ifndef _R3BGAMMAGENERATOR_H_
#define _R3BGAMMAGENERATOR_H_

#include <map>
#include <utility>

#include <TRandom3.h>
#include <TVector3.h>

#include <FairGenerator.h>

/**
 * Gamma (cacade) event generator.
 *
 * You may consider an instance of this class as an excited nucleus performing gamma decay.
 *
 * Use SetEnergyLevel(id, E) to define energy levels with arbitrary, integer IDs and energies in GeV.
 * Then, use SetBranchingRatio(idLevelFrom, idLevelTo, ratio) to define the transition propabilities between the levels.
 * Finally, using SetInitalLevel() define the one and only initial energy level to start the cascade,
 * or SetInitialLevelDistribution() to define multiple possible initial states with a given distribution.
 */
class R3BGammaGenerator : public FairGenerator
{
  public:
    R3BGammaGenerator()
        : FairGenerator("R3BGammaGenerator", "R3BGammaGenerator")
    {
    }
    R3BGammaGenerator(const char* name, const char* title)
        : FairGenerator(name, title)
    {
    }

    virtual Bool_t Init();
    virtual Bool_t ReadEvent(FairPrimaryGenerator* primGen);

    /**
     * Defines an energy level.
     *
     * @param idLevel Arbitrary ID for the level for later usage in SetBranchingRatio(), SetInitialLevel() or
     * SetInitialLevelDistribution().
     * @param E Energy (GeV) of the level
     */
    void SetEnergyLevel(int idLevel, double E);

    /**
     * Defines the branching ratio (transition probability) between to levels.
     *
     * @remark The sum of all weight/probabilities from each initial level is renormalized to 100%.
     *
     * @param idLevelFrom ID of source level (as defined by SetEnergyLevel() )
     * @param idLevelTo ID of destination level (as defined by SetEnergyLevel() )
     * @param ratio Weight for this transition (the sum of all weights for each source level will be renormalized to
     * 100%)
     * @param radiative If false, no gamma ray will be emitted for this transition.
     */
    void SetBranchingRatio(int idLevelFrom, int idLevelTo, double ratio, bool radiative = true);

    /**
     * Defines the only initially populated energy level.
     *
     * @param idLevel ID of initial energy level (as defined by SetEnergyLevel() )
     */
    void SetInitialLevel(int idLevel);

    /**
     * Defines a distribution of initially populated levels.
     *
     * @param idLevel ID of initial energy level (as defined by SetEnergyLevel() )
     * @param ratio Weight of this level (all weights will be renormalized to 100%)
     */
    void SetInitialLevelDistribution(int idLevel, double ratio);

    /**
     * Simulate in-flight emission of Gamma rays: Defines beta of emitting system
     * relative to lab frame.
     *
     * @param beta Direction and velocity (relativistiv beta, i.e. v/c) of emitting system (moving ion)
     */
    void SetLorentzBoost(TVector3& beta);

  protected:
    void GenerateGamma(double E, FairPrimaryGenerator* primGem);

    std::map<int, double> fEnergyLevels;
    std::map<int, std::map<int, std::pair<double, bool>>> fBranchingRatios;
    std::map<int, double> fInitialLevelDistribution;

    TVector3 fBeta;

    TRandom3 fRandom;

    ClassDef(R3BGammaGenerator, 1);
};

#endif
