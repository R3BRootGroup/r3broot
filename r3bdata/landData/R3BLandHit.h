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

// -------------------------------------------------------------------------
// -----                      R3BLandHit header file                  -----
// -----                  Created 26/03/09  by D.Bertini               -----
// -------------------------------------------------------------------------

/**  R3BLandHit.h
 **/

#ifndef R3BLANDHIT_H
#define R3BLANDHIT_H

#include "TObject.h"
#include "TVector3.h"

#include "FairHit.h"

class R3BLandHit : public FairHit
{

  public:
    /** Default constructor **/
    R3BLandHit();

    /** Copy constructor **/
    R3BLandHit(const R3BLandHit&);

    /** Destructor **/
    virtual ~R3BLandHit();

    R3BLandHit& operator=(const R3BLandHit&) { return *this; }

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void SetEnergy(Double_t en) { fEnergy = en; }
    Double_t GetEnergy() { return fEnergy; }

    void SetToF(Double_t time) { fTof = time; }
    Double_t GetToF() { return fTof; }

    void SetDetID(Int_t id) { fDetectorID = id; }
    Int_t GetDetID() { return fDetectorID; }

    void SetSaturnPaddle(Int_t paddleNb) { fSaturnPaddleNb = paddleNb; }
    Int_t GetSaturnPaddle() { return fSaturnPaddleNb; }

  protected:
    Double32_t fEnergy;
    Double32_t fTof;
    Int_t fSaturnPaddleNb;

    ClassDef(R3BLandHit, 1)
};

#endif
