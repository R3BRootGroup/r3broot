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

#ifndef R3BGFIDIGITIZER_H
#define R3BGFIDIGITIZER_H 1

#include "FairTask.h"
#include "R3BGfiDigi.h"
#include "R3BGfiDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BGfiDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BGfiDigitizer();

    /** Destructor **/
    ~R3BGfiDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BGfiDigi* AddHit(Int_t gf1mul, Double_t gf1_01x, Int_t gf2mul, Double_t gf2_01x);

  protected:
    TClonesArray* fGfiPoints;
    TClonesArray* fGfiMCTrack;
    TClonesArray* fGfiDigi;

    // Parameter class
    R3BGfiDigiPar* fGfiDigiPar;

    //- Control Hitograms

    Int_t eventNoGfi;

  private:
    virtual void SetParContainers();

    ClassDef(R3BGfiDigitizer, 1);
};

#endif
