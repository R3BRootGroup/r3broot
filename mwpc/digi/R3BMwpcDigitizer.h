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

// --------------------------------------------------------------
// -----          R3BMwpcDigitizer source file              -----
// -----       Created 03/11/19 by JL Rodriguez             -----
// --------------------------------------------------------------

#ifndef R3BMwpcDigitizer_H
#define R3BMwpcDigitizer_H 1

#include "FairTask.h"
#include "R3BMwpcHitData.h"
#include "TRotation.h"
#include "TVector3.h"
#include <map>
#include <string>

class TClonesArray;
class R3BTGeoPar;

class R3BMwpcDigitizer : public FairTask
{
  public:
    /** Default constructor **/
    R3BMwpcDigitizer();

    /** Standard constructor **/
    R3BMwpcDigitizer(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BMwpcDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    // Fair specific
    virtual void SetParContainers();

    virtual void Reset();

    /** Setters for sigmas **/
    void SetSigma_x(Float_t sigma_x) { fsigma_x = sigma_x; }
    void SetSigma_y(Float_t sigma_y) { fsigma_y = sigma_y; }

  private:
    void SetParameter();

    TClonesArray* fMCTrack;
    TClonesArray* fMwpcPoints;
    TClonesArray* fMwpcHits;
    R3BTGeoPar* fMwpcGeoPar;
    TString fName;
    Float_t fsigma_x;
    Float_t fsigma_y;
    TVector3 fTrans;
    TRotation fRot;

    /** Private method AddHitData **/
    R3BMwpcHitData* AddHitData(Double_t x, Double_t y);

  public:
    // Class definition
    ClassDef(R3BMwpcDigitizer, 1);
};

#endif
