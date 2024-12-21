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
// -----                      R3BTraHit header file                  -----
// -----                  Created 07/07/09  by M. Labiche               -----
// -------------------------------------------------------------------------

#ifndef R3BTRAHIT_H
#define R3BTRAHIT_H

#include "TObject.h"
#include "TVector3.h"

#include "FairHit.h"

class R3BTraHit : public FairHit
{

  public:
    /** Default constructor **/
    R3BTraHit();

    /** Copy constructor **/
    R3BTraHit(const R3BTraHit& point) { *this = point; };

    /** Destructor **/
    virtual ~R3BTraHit();

    /** Output to screen **/
    virtual void Print(const Option_t* opt) const;

    void SetEnergy(Double_t en) { fEnergy = en; }
    Double_t GetEnergy() { return fEnergy; }

    void SetToF(Double_t time) { fTof = time; }
    Double_t GetToF() { return fTof; }

    void SetDetID(Int_t id) { fDetectorID = id; }
    Int_t GetDetID() { return fDetectorID; }

  protected:
    Double32_t fEnergy;
    Double32_t fTof;

    ClassDef(R3BTraHit, 1)
};

#endif
