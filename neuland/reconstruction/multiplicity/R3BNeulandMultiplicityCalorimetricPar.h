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

#ifndef R3BNEULANDMULTIPLICITYCALORIMETRICPAR_H
#define R3BNEULANDMULTIPLICITYCALORIMETRICPAR_H

#include "FairParGenericSet.h"
#include "FairParamList.h"
#include "TCutG.h"
#include "TMap.h"
#include <map>

/**
 * NeuLAND number of clusters / energy - neutron multiplicity parameter storage
 * @author Jan Mayer
 *
 * Stores the cuts for the 2D Calibr method, can be asked about the neutron multiplicity
 */

class R3BNeulandMultiplicityCalorimetricPar : public FairParGenericSet
{
  public:
    // needs to be public?
    // Note: There MUST NOT be a //! after the variable
    TMap* fNeutronCuts;

    R3BNeulandMultiplicityCalorimetricPar(const char* name = "R3BNeulandMultiplicityCalorimetricPar",
                                          const char* title = "Neuland Neutron 2D Parameters",
                                          const char* context = "TestDefaultContext");
    ~R3BNeulandMultiplicityCalorimetricPar() override;

    void clear() override;
    void putParams(FairParamList*) override;
    Bool_t getParams(FairParamList*) override;
    void printParams() override;

    std::map<UInt_t, TCutG*> GetNeutronCuts() const;
    TCutG* GetNeutronCut(const Int_t n) const;
    void SetNeutronCuts(const std::map<UInt_t, TCutG*>& cuts);
    UInt_t GetNeutronMultiplicity(const Double_t energy, const Double_t nClusters) const;

  private:
    R3BNeulandMultiplicityCalorimetricPar(const R3BNeulandMultiplicityCalorimetricPar&);
    R3BNeulandMultiplicityCalorimetricPar& operator=(const R3BNeulandMultiplicityCalorimetricPar&);

    ClassDefOverride(R3BNeulandMultiplicityCalorimetricPar, 2)
};

#endif // R3BNEULANDMULTIPLICITYCALORIMETRICPAR_H
