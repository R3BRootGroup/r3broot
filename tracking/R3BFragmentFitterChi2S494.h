/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BFRAGMENTFITTERCHI2S494
#define R3BFRAGMENTFITTERCHI2S494

#include "R3BFragmentFitterGeneric.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"
#include "Minuit2/Minuit2Minimizer.h"

class R3BFragmentFitterChi2S494 : public R3BFragmentFitterGeneric
{
  public:
    R3BFragmentFitterChi2S494();
    ~R3BFragmentFitterChi2S494();

    void Init(R3BTPropagator* prop = nullptr, Bool_t energyLoss = kTRUE);

    Int_t FitTrack(R3BTrackingParticle*, R3BTrackingSetup*);
    
    Int_t FitTrackBeta(R3BTrackingParticle*, R3BTrackingSetup*);

    Int_t FitTrackMomentumForward(R3BTrackingParticle*, R3BTrackingSetup*);
    
    Int_t FitTrackMomentumBackward(R3BTrackingParticle*, R3BTrackingSetup*);

    Int_t FitTrackBackward(R3BTrackingParticle*, R3BTrackingSetup*);
    
    Int_t FitTrackBackward2D(R3BTrackingParticle*, R3BTrackingSetup*);

    Double_t TrackFragment(R3BTrackingParticle* particle,
                           Bool_t energyLoss,
                           Double_t& devTof,
                           Double_t& time,
                           Double_t& chi2);
    Double_t DbetaDx(R3BTrackingParticle* candidate);
    Double_t DbetaChi2(R3BTrackingParticle* candidate);
    Double_t DbetaDt(R3BTrackingParticle* candidate);
    Double_t DmDx(R3BTrackingParticle* candidate, Bool_t energy_loss);
    Double_t DmDxTof(R3BTrackingParticle* candidate, Bool_t energy_loss);
    Double_t DmDt(R3BTrackingParticle* candidate, Bool_t energy_loss);
    Double_t Velocity(R3BTrackingParticle* candidate);

  private:
    ROOT::Math::Minimizer* fMinimum;
    ROOT::Math::Minimizer* minimum_m;
    R3BTPropagator* fPropagator;
   	Double_t amu = 0.931494028;   // Gev/c**2
    TVector3 pinv;

    ClassDef(R3BFragmentFitterChi2S494, 1)
};

#endif
