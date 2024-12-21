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

#ifndef R3BROOT_R3BNEULANDONLINERECONSTRUCTION_H
#define R3BROOT_R3BNEULANDONLINERECONSTRUCTION_H

#include "FairTask.h"
#include "R3BNeulandCluster.h"
#include "R3BNeulandHit.h"
#include "TCAConnector.h"

class TH1D;
class TH2D;
class TH3D;
class R3BEventHeader;

class R3BNeulandOnlineReconstruction : public FairTask
{
  public:
    R3BNeulandOnlineReconstruction();

    InitStatus Init() override;
    void Exec(Option_t*) override;
    void FinishTask() override;

    void ResetHistos();
    void SetDistanceToTarget(double x) { fDistanceToTarget = x; }

  private:
    R3BEventHeader* fEventHeader;

    TCAInputConnector<R3BNeulandHit> fNeulandHits;
    TCAInputConnector<R3BNeulandCluster> fNeulandClusters;

    TH1D* hHitX;
    TH1D* hHitY;
    TH1D* hHitZ;
    TH1D* hHitT;
    TH1D* hHitTadj;
    TH1D* hHitE;

    TH1D* hHitMult;
    TH1D* hClusterMult;
    TH1D* hClusterSize;
    TH3D* hHits3D;
    TH3D* hClusters3D;

    TH1D* hEtot;
    TH2D* hEtotVSNClusters;

    double fDistanceToTarget;

    ClassDefOverride(R3BNeulandOnlineReconstruction, 0)
};

#endif // R3BROOT_R3BNEULANDONLINERECONSTRUCTION_H
