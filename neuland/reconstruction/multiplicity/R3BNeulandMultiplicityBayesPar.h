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

#ifndef R3BNEULANDMULTIPLICITYBAYESPAR_H
#define R3BNEULANDMULTIPLICITYBAYESPAR_H

#include "FairParGenericSet.h"
#include "FairParamList.h"
#include "R3BNeulandMultiplicity.h"
#include "TArrayD.h"
#include <array>

class R3BNeulandMultiplicityBayesPar : public FairParGenericSet
{
  public:
    R3BNeulandMultiplicityBayesPar(const char* name = "R3BNeulandMultiplicityBayesPar",
                                   const char* title = "Neuland Multiplicity Bayes Parameters",
                                   const char* context = "TestDefaultContext");
    ~R3BNeulandMultiplicityBayesPar() override;

    void clear() override;
    void putParams(FairParamList*) override;
    Bool_t getParams(FairParamList*) override;
    // void printParams() override;

    void Fill(int n, int nHits, int nClusters, int Edep);
    bool CheckIfProperlyLoaded() const;
    R3BNeulandMultiplicity::MultiplicityProbabilities GetProbabilities(int nHits, int nClusters, int Edep) const;

    TArrayD GetHitDensity(int n) const { return fHits.at(n); }
    TArrayD GetClusterDensity(int n) const { return fClusters.at(n); }
    TArrayD GetEdepDensity(int n) const { return fEdep.at(n); }

    void Finish();

  private:
    std::array<TArrayD, NEULAND_MAX_MULT> fHits;
    std::array<TArrayD, NEULAND_MAX_MULT> fClusters;
    std::array<TArrayD, NEULAND_MAX_MULT> fEdep;
    mutable bool fIsProperlyLoaded;

    R3BNeulandMultiplicityBayesPar(const R3BNeulandMultiplicityBayesPar&);
    R3BNeulandMultiplicityBayesPar& operator=(const R3BNeulandMultiplicityBayesPar&);

    ClassDefOverride(R3BNeulandMultiplicityBayesPar, 1)
};

#endif // R3BNEULANDMULTIPLICITYBAYESPAR_H
