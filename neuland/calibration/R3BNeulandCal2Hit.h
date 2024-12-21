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

#ifndef R3BNEULANDCAL2HIT_H
#define R3BNEULANDCAL2HIT_H

#include "FairTask.h"
#include "R3BNeulandCalData.h"
#include "R3BNeulandHit.h"
#include "R3BNeulandHitModulePar.h"
#include "TCAConnector.h"
#include <map>
#include <vector>

class R3BNeulandHitPar;
class R3BEventHeader;

class R3BNeulandCal2Hit : public FairTask
{
  public:
    R3BNeulandCal2Hit(const char* name = "R3BNeulandCal2Hit", const Int_t iVerbose = 0);
    ~R3BNeulandCal2Hit() override = default;

    InitStatus Init() override;
    InitStatus ReInit() override;
    void SetParContainers() override;
    void Exec(Option_t*) override;
    void FinishTask() override;

    inline void SetFirstPlaneHorizontal() { fFirstPlaneHorizontal = kTRUE; }

    // Distance to target in cm
    inline void SetDistanceToTarget(Double_t d) { fDistanceToTarget = d; }

    // Global time offset in ns
    inline void SetGlobalTimeOffset(Double_t t0) { fGlobalTimeOffset = t0; }

    // Energy cutoff in MeV
    inline void SetEnergyCutoff(Double_t enecut) { fEnergyCutoff = enecut; }

    // Walk correction
    inline void EnableWalk(Bool_t walk = kFALSE) { fWalkEnabled = walk; }

  private:
    void SetParameter();
    Double_t GetUnsaturatedEnergy(const Int_t qdc, const Double_t gain, const Double_t saturation) const;

    R3BEventHeader* fEventHeader;

    TCAInputConnector<R3BNeulandCalData> fCalData;
    TCAOutputConnector<R3BNeulandHit> fHits;

    R3BNeulandHitPar* fPar;

    Bool_t fFirstPlaneHorizontal;
    Int_t fNumberOfPlanes;
    Double_t fDistanceToTarget;
    std::vector<Double_t> fDistancesToFirstPlane;
    std::vector<Double_t> fAttenuationValues;
    Double_t fGlobalTimeOffset;
    Double_t fEnergyCutoff;

    Bool_t fWalkEnabled;
    Double_t WalkCorrection(Double_t);

    std::map<Int_t, R3BNeulandHitModulePar> fParMap;

    std::map<Int_t, R3BNeulandCalData*> fHitMap;

    UInt_t fEventNumber = 0;

  public:
    ClassDefOverride(R3BNeulandCal2Hit, 0)
};

#endif
