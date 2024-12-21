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

// ----------------------------------------------------------------
// -----           R3BTwimDigitizer source file               -----
// -----         Created 06/11/19 by JL Rodriguez             -----
// ----------------------------------------------------------------

#ifndef R3BTwimDigitizer_H
#define R3BTwimDigitizer_H 1

#include "FairTask.h"
#include "R3BTwimHitData.h"
#include <map>
#include <string>

class TClonesArray;

class R3BTwimDigitizer : public FairTask
{
  public:
    /** Default constructor **/
    R3BTwimDigitizer();

    /** Standard constructor **/
    R3BTwimDigitizer(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BTwimDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Reset();

    /** Setters for sigmas **/
    void SetSigma_x(Float_t sig) { fsigma_x = sig; }
    void SetSigCharge(Float_t sig) { fZsig = sig; }

  private:
    TClonesArray* fMCTrack;
    TClonesArray* fTwimPoints;
    TClonesArray* fTwimHits;
    Int_t fDetId;
    Float_t fsigma_x;
    Float_t fZsig;
    TString fName;

    // Adds a R3BTwimHitData to the HitCollection
    R3BTwimHitData* AddR3BHitData(Int_t detId, Double_t theta, Double_t z, Double_t x, Double_t offset);

  public:
    // Class definition
    ClassDef(R3BTwimDigitizer, 1);
};

#endif /* R3BTwimDigitizer_H */
