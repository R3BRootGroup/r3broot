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

#ifndef R3BPSPDIGITIZER_H
#define R3BPSPDIGITISER_H 1

#include "FairTask.h"
#include "R3BPspDigi.h"
#include "R3BPspDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BPspDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BPspDigitizer();

    /** Destructor **/
    ~R3BPspDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BPspDigi* AddHit(Int_t psp3mul, Double_t psp3x, Double_t psp3y, Double_t psp3e);

  protected:
    TClonesArray* fPspPoints;
    TClonesArray* fPspMCTrack;
    TClonesArray* fPspDigi;

    // Parameter class
    R3BPspDigiPar* fPspDigiPar;

    //- Control Hitograms

    //   TH1F* PSPXhis;

    Int_t eventNoPsp;

  private:
    virtual void SetParContainers();

    ClassDef(R3BPspDigitizer, 1);
};

#endif
