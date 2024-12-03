/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BNEULANDHITMON_H
#define R3BNEULANDHITMON_H

/** Neuland Digi Monitor
 * @author Jan Mayer
 * @since  10.06.2015
 *
 *  Input:  Digis. Should work both on Digis from Monte Carlo simulations
 *          as well as experimental data on digi level.
 *  Output: Vairous diagrams, 3D scatterplot for each event if option 3DTRACK is given.
 */

#include "FairTask.h"
#include "R3BNeulandHit.h"
#include "TCAConnector.h"

class TH1D;
class TH2D;
class TH3D;
class TH1I;

class R3BNeulandHitMon : public FairTask
{
  public:
    explicit R3BNeulandHitMon(TString input = "NeulandHits",
                              TString output = "NeulandHitMon",
                              const Option_t* option = "");

    ~R3BNeulandHitMon() override = default;

    // No copy and no move is allowed (Rule of three/five)
    R3BNeulandHitMon(const R3BNeulandHitMon&) = delete;            // copy constructor
    R3BNeulandHitMon(R3BNeulandHitMon&&) = delete;                 // move constructor
    R3BNeulandHitMon& operator=(const R3BNeulandHitMon&) = delete; // copy assignment
    R3BNeulandHitMon& operator=(R3BNeulandHitMon&&) = delete;      // move assignment

  protected:
    InitStatus Init() override;
    void Finish() override;

  public:
    void Exec(Option_t* /*option*/) override;

    void SetDistanceToTarget(double distance) { fDistanceToTarget = distance; }

  private:
    TString fOutput;

    TCAInputConnector<R3BNeulandHit> fHits;

    double fDistanceToTarget = 0.;

    Bool_t fIs3DTrackEnabled = false;
    TH3D* fh3 = nullptr;

    TH1D* hTime = nullptr;
    TH1D* hTimeAdj = nullptr;
    TH1I* hMult = nullptr;
    TH1D* hDepth = nullptr;
    TH1D* hForemostEnergy = nullptr;
    TH1D* hSternmostEnergy = nullptr;
    TH2D* hDepthVSForemostEnergy = nullptr;
    TH2D* hDepthVSSternmostEnergy = nullptr;
    TH1D* hEtot = nullptr;
    TH1D* hE = nullptr;
    TH1D* hX = nullptr;
    TH1D* hY = nullptr;
    TH1D* hT = nullptr;
    TH1D* hTNeigh = nullptr;
    TH2D* hDepthVSEtot = nullptr;
    TH2D* hPosVSEnergy = nullptr;
    TH2D* hdeltaEE = nullptr;
    TH1D* hBeta = nullptr;

    ClassDefOverride(R3BNeulandHitMon, 0); // NOLINT
};

#endif // R3BNEULANDHITMON_H
