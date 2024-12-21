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

#ifndef R3BTOF2PDIGITIZER_H
#define R3BTOF2PDIGITIZER_H 1

#include "FairTask.h"
#include "R3BTof2pDigi.h"
#include "R3BTof2pDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BTof2pDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BTof2pDigitizer();

    /** Destructor **/
    ~R3BTof2pDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BTof2pDigi* AddHit(
        Int_t tfmul,
        Double_t tfwx_p1,
        Double_t tfwy_p1,
        Double_t tfwt_p1,
        Double_t tfwx_p2,
        Double_t tfwy_p2,
        Double_t tfwt_p2,
        Double_t tfwpath_p1,
        Double_t tfwpath_p2,
        // Double_t tfwpx_p1,Double_t tfwpy_p1,Double_t tfwpz_p1,Double_t tfwpx_p2,Double_t tfwpy_p2,Double_t tfwpz_p2);
        Double_t tfwpx_p1,
        Double_t tfwpy_p1,
        Double_t tfwpz_p1,
        Double_t tfwpx_p2,
        Double_t tfwpy_p2,
        Double_t tfwpz_p2,
        Double_t tfwe_p1,
        Double_t tfwe_p2);

  protected:
    TClonesArray* fTof2pDigi;
    TClonesArray* fTofPoints;
    TClonesArray* fTofMCTrack;

    // Parameter class
    R3BTof2pDigiPar* fTof2pDigiPar;

    //- Control Hitograms

    TH1F* TOFelosshis;

    TH1F* TfwXhis;
    TH1F* TfwYhis;
    TH1F* TfwThis;

    Int_t eventNoTof;

  private:
    virtual void SetParContainers();

    ClassDef(R3BTof2pDigitizer, 1);
};

#endif
