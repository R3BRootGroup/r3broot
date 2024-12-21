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

#ifndef R3BROOT_R3BNEULANDNEUTRONRECONSTRUCTIONSTATISTICS_H
#define R3BROOT_R3BNEULANDNEUTRONRECONSTRUCTIONSTATISTICS_H

#include "FairTask.h"
#include "R3BNeulandCluster.h"
#include "R3BNeulandNeutron.h"
#include "TCAConnector.h"
#include "TH1D.h"
#include <iostream>

class R3BNeulandNeutronReconstructionStatistics : public FairTask
{
  public:
    R3BNeulandNeutronReconstructionStatistics(const TString primary = "NeulandNeutronClusters",
                                              const TString secondary = "NeulandPrimaryClusters",
                                              const TString predicted = "NeulandSecondaryClusters",
                                              std::ostream& out = std::cout);
    ~R3BNeulandNeutronReconstructionStatistics() override = default;

    // No copy and no move is allowed (Rule of three/five)
    R3BNeulandNeutronReconstructionStatistics(const R3BNeulandNeutronReconstructionStatistics&) =
        delete;                                                                                      // copy constructor
    R3BNeulandNeutronReconstructionStatistics(R3BNeulandNeutronReconstructionStatistics&&) = delete; // move constructor
    R3BNeulandNeutronReconstructionStatistics& operator=(const R3BNeulandNeutronReconstructionStatistics&) =
        delete; // copy assignment
    R3BNeulandNeutronReconstructionStatistics& operator=(R3BNeulandNeutronReconstructionStatistics&&) =
        delete; // move assignment

  protected:
    InitStatus Init() override;
    void Finish() override;

  public:
    void Exec(Option_t*) override;

  private:
    TCAInputConnector<R3BNeulandCluster> fPrimaryClusters;
    TCAInputConnector<R3BNeulandCluster> fSecondaryClusters;
    TCAInputConnector<R3BNeulandNeutron> fPredictedNeutrons;
    TString fPredictedName;
    TH1D* fhTP;
    TH1D* fhFP;
    TH1D* fhFN;
    TH1D* fhTN;
    TH1D* fhF1;
    int fTP;
    int fFP;
    int fFN;
    int fTN;
    std::ostream& fOut;
    std::vector<int> fMult;

  public:
    ClassDefOverride(R3BNeulandNeutronReconstructionStatistics, 0);
};

#endif // R3BROOT_R3BNEULANDNEUTRONRECONSTRUCTIONSTATISTICS_H
