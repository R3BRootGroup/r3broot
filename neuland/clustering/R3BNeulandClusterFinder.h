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

#ifndef R3BNEULANDCLUSTERFINDER_H_NEW
#define R3BNEULANDCLUSTERFINDER_H_NEW

/**
 * NeuLAND cluster finder task
 * @author Jan Mayer
 *
 * For each event, get the R3BNeulandHits and group them into R3BNeulandClusters using the Neuland Clustering Engine.
 *   Input:  Branch NeulandHits    = TClonesArray("R3BNeulandDigi")
 *   Output: Branch NeulandClusters = TClonesArray("R3BNeulandCluster")
 *
 */

#include "ClusteringEngine.h"
#include "FairTask.h"
#include "R3BNeulandCluster.h"
#include "R3BNeulandHit.h"
#include "TCAConnector.h"
#include <TClonesArray.h>

class R3BNeulandClusterFinder : public FairTask
{
  public:
    R3BNeulandClusterFinder(Double_t dx = 1. * 7.5,
                            Double_t dy = 1. * 7.5,
                            Double_t dz = 2. * 7.5,
                            Double_t dt = 1.,
                            TString input = "NeulandHits",
                            TString output = "NeulandClusters");

    ~R3BNeulandClusterFinder() override = default;

    // No copy and no move is allowed (Rule of three/five)
    R3BNeulandClusterFinder(const R3BNeulandClusterFinder&) = delete;            // copy constructor
    R3BNeulandClusterFinder(R3BNeulandClusterFinder&&) = delete;                 // move constructor
    R3BNeulandClusterFinder& operator=(const R3BNeulandClusterFinder&) = delete; // copy assignment
    R3BNeulandClusterFinder& operator=(R3BNeulandClusterFinder&&) = delete;      // move assignment

  protected:
    InitStatus Init() override;

  public:
    void Exec(Option_t*) override;

  private:
    Neuland::ClusteringEngine<R3BNeulandHit> fClusteringEngine;
    TCAInputConnector<R3BNeulandHit> fDigis;
    TCAOutputConnector<R3BNeulandCluster> fClusters;

    ClassDefOverride(R3BNeulandClusterFinder, 0);
};

#endif // R3BNEULANDCLUSTERFINDER_H_NEW
