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

#ifndef R3BROOT_R3BNEULANDPRIMARYINTERACTIONFINDER_H
#define R3BROOT_R3BNEULANDPRIMARYINTERACTIONFINDER_H

#include "FairTask.h"
#include "R3BMCTrack.h"
#include "R3BNeulandHit.h"
#include "R3BNeulandPoint.h"
#include "TCAConnector.h"
#include "TH1D.h"
#include "TH2D.h"

class R3BNeulandPrimaryInteractionFinder : public FairTask
{
  public:
    explicit R3BNeulandPrimaryInteractionFinder(TString pointsIn = "NeulandPoints",
                                                TString hitsIn = "NeulandHits",
                                                TString pointsOut = "NeulandPrimaryPoints",
                                                TString hitsOut = "NeulandPrimaryHits",
                                                TString tracksOut = "NeulandPrimaryTracks");

    ~R3BNeulandPrimaryInteractionFinder() override = default;

    // No copy and no move is allowed (Rule of three/five)
    R3BNeulandPrimaryInteractionFinder(const R3BNeulandPrimaryInteractionFinder&) = delete;
    R3BNeulandPrimaryInteractionFinder(R3BNeulandPrimaryInteractionFinder&&) = delete;
    R3BNeulandPrimaryInteractionFinder& operator=(const R3BNeulandPrimaryInteractionFinder&) = delete;
    R3BNeulandPrimaryInteractionFinder& operator=(R3BNeulandPrimaryInteractionFinder&&) = delete;

  protected:
    InitStatus Init() override;
    void Finish() override;

  public:
    void Exec(Option_t*) override;

  private:
    TCAInputConnector<R3BMCTrack> fTracksIn;
    TCAInputConnector<R3BNeulandPoint> fPointsIn;
    TCAInputConnector<R3BNeulandHit> fHitsIn;

    TCAOutputConnector<R3BMCTrack> fTracksOut;
    TCAOutputConnector<R3BNeulandPoint> fPointsOut;
    TCAOutputConnector<R3BNeulandHit> fHitsOut;

    TH1D* fhDistance;
    TH2D* fhPointsVsHits;
    TH2D* fhPointVsHitPaddle;

    ClassDefOverride(R3BNeulandPrimaryInteractionFinder, 0);
};

#endif // R3BROOT_R3BNEULANDPRIMARYINTERACTIONFINDER_H
