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

#ifndef R3BTOFDIGITIZER_H
#define R3BTOFDDIGITISER_H 1

#include "FairTask.h"
#include "R3BTofDigi.h"
#include "R3BTofDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BTofDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BTofDigitizer();

    /** Destructor **/
    ~R3BTofDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BTofDigi* AddHit(Int_t tfmul,
                       Double_t tfwx_p1,
                       Double_t tfwy_p1,
                       Double_t tfwt_p1,
                       Double_t tfwpath_p1,
                       // Double_t tfwpx_p1,Double_t tfwpy_p1,Double_t tfwpz_p1);
                       Double_t tfwpx_p1,
                       Double_t tfwpy_p1,
                       Double_t tfwpz_p1,
                       Double_t tfwe_p1);

  protected:
    TClonesArray* fTofDigi;
    TClonesArray* fTofPoints;
    TClonesArray* fTofMCTrack;

    // Parameter class
    R3BTofDigiPar* fTofDigiPar;

    //- Control Hitograms

    TH1F* TOFelosshis;

    TH1F* TfwXhis;
    TH1F* TfwYhis;
    TH1F* TfwThis;

    Int_t eventNoTof;

  private:
    virtual void SetParContainers();

    ClassDef(R3BTofDigitizer, 1);
};

#endif
