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

// R3BsfiDigitizer.h 25. 08. 2017.
// Created by: Aron Kripko

#ifndef R3BSFIDIGITISER_H
#define R3BSFIDIGITISER_H 1

#include "FairTask.h"
#include "R3BFi4HitItem.h"
#include <TRandom3.h>
#include <map>
#include <string>

class TClonesArray;
class TH1F;
class TH2F;

class R3BsfiDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BsfiDigitizer();

    R3BsfiDigitizer(Double_t esigma, Double_t tsigma);

    /** Destructor **/
    ~R3BsfiDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    void SetEnergyResolution(Double_t e);
    void SetTimeResolution(Double_t t);

  protected:
    TClonesArray* fPoints;

    TClonesArray* fMCTrack;

    TClonesArray* Hits;

  private:
    TRandom3* prnd;
    Double_t esigma;
    Double_t tsigma;
    Int_t NumOfFibers;

    ClassDef(R3BsfiDigitizer, 1);
};

#endif
