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

#ifndef R3BROOT_R3BNEULANDPRIMARYCLUSTERFINDER_H
#define R3BROOT_R3BNEULANDPRIMARYCLUSTERFINDER_H

#include "FairTask.h"
#include "R3BMCTrack.h"
#include "R3BNeulandCluster.h"
#include "R3BNeulandHit.h"
#include "TCAConnector.h"

class R3BNeulandPrimaryClusterFinder : public FairTask
{
  public:
    explicit R3BNeulandPrimaryClusterFinder(TString primaryHits = "NeulandPrimaryHits",
                                            TString clusters = "NeulandClusters",
                                            TString primaryClusters = "NeulandPrimaryClusters",
                                            TString secondaryClusters = "NeulandSecondaryClusters");

    ~R3BNeulandPrimaryClusterFinder() override = default;

    // No copy and no move is allowed (Rule of three/five)
    R3BNeulandPrimaryClusterFinder(const R3BNeulandPrimaryClusterFinder&) = delete;
    R3BNeulandPrimaryClusterFinder(R3BNeulandPrimaryClusterFinder&&) = delete;
    R3BNeulandPrimaryClusterFinder& operator=(const R3BNeulandPrimaryClusterFinder&) = delete;
    R3BNeulandPrimaryClusterFinder& operator=(R3BNeulandPrimaryClusterFinder&&) = delete;

  protected:
    InitStatus Init() override;

  public:
    void Exec(Option_t*) override;

  private:
    TCAInputConnector<R3BNeulandHit> fPrimaryHits;
    TCAInputConnector<R3BNeulandCluster> fClusters;

    TCAOutputConnector<R3BNeulandCluster> fPrimaryClusters;
    TCAOutputConnector<R3BNeulandCluster> fSecondaryClusters;

    ClassDefOverride(R3BNeulandPrimaryClusterFinder, 0);
};

#endif // R3BROOT_R3BNEULANDPRIMARYCLUSTERFINDER_H
