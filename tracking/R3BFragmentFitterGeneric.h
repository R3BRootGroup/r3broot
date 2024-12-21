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

#ifndef R3BFRAGMENTFITTERGENERIC
#define R3BFRAGMENTFITTERGENERIC

#include "Rtypes.h"

class R3BTrackingParticle;
class R3BTrackingSetup;
class R3BTPropagator;

class R3BFragmentFitterGeneric
{
  public:
    R3BFragmentFitterGeneric();
    virtual ~R3BFragmentFitterGeneric();

    virtual void Init(R3BTPropagator* prop = nullptr, Bool_t energyLoss = kTRUE) = 0;

    virtual Int_t FitTrack(R3BTrackingParticle*, R3BTrackingSetup*) = 0;

    virtual Int_t FitTrackBeta(R3BTrackingParticle*, R3BTrackingSetup*) = 0;

    virtual Int_t FitTrackBackward(R3BTrackingParticle*, R3BTrackingSetup*) = 0;

    virtual Int_t FitTrackBackward2D(R3BTrackingParticle*, R3BTrackingSetup*) = 0;

    ClassDef(R3BFragmentFitterGeneric, 1)
};

#endif
