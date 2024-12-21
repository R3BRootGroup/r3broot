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

#ifndef R3B_FRAGMENTTRACKER_H
#define R3B_FRAGMENTTRACKER_H

#include "FairTask.h"

#include <string>
#include <vector>

class TClonesArray;
class R3BFieldPar;
class R3BTPropagator;
class R3BTrackingDetector;
class R3BTrackingParticle;
class R3BTrackingSetup;
class R3BFragmentFitterGeneric;

class TH1F;

class R3BFragmentTracker : public FairTask
{
  public:
    R3BFragmentTracker(const char* name, Bool_t vis = kFALSE, Int_t verbose = 1);
    virtual ~R3BFragmentTracker();

    virtual InitStatus Init();
    virtual InitStatus ReInit();

    virtual void SetParContainers();

    virtual void Exec(const Option_t* = "");

    virtual void Finish();

    void SetFragmentFitter(R3BFragmentFitterGeneric* fitter) { fFitter = fitter; }
    void SetEnergyLoss(Bool_t energyLoss) { fEnergyLoss = energyLoss; }

  private:
    Bool_t InitPropagator();

    R3BFieldPar* fFieldPar;
    R3BTPropagator* fPropagator;
    TClonesArray* fArrayMCTracks; // simulation output??? To compare?
    R3BTrackingSetup* fDetectors; // array of R3BTrackingDetector
    std::vector<R3BTrackingParticle*> fFragments;
    TClonesArray* fArrayFragments;
    Int_t fNEvents;
    Bool_t fVis;
    //   	Double_t amu = 0.93149410242;
    Double_t amu = 0.938272;

    R3BFragmentFitterGeneric* fFitter;
    Bool_t fEnergyLoss;

    Double_t fAfterGladResolution;

    TH1F* fh_mult_psp;
    TH1F* fh_mult_fi4;
    TH1F* fh_mult_fi5;
    TH1F* fh_mult_tofd;
    TH1F* fh_eloss_psp_mc;
    TH1F* fh_eloss_psp;
    TH1F* fh_eloss_fi4_mc;
    TH1F* fh_eloss_fi4;
    TH1F* fh_eloss_fi5_mc;
    TH1F* fh_eloss_fi5;
    TH1F* fh_eloss_fi6_mc;
    TH1F* fh_eloss_fi6;
    TH1F* fh_ncand;
    TH1F* fh_x_res[6];
    TH1F* fh_x_pull[6];
    TH1F* fh_A_reco1;
    TH1F* fh_A_reco2;
    TH1F* fh_mom_res;
    TH1F* fh_mass_res;
    TH1F* fh_chi2;
    TH1F* fh_vz_res;
    TH1F* fh_beta_res;

    ClassDef(R3BFragmentTracker, 1)
};

#endif
