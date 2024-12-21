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

#ifndef R3BMTOFDIGITIZER_H
#define R3BMTOFDIGITIZER_H 1

#include "FairTask.h"
#include "R3BmTofDigi.h"
#include "R3BmTofDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BmTofDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BmTofDigitizer();

    /** Destructor **/
    ~R3BmTofDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    // R3BmTofDigi* AddHit(Int_t ntmul,Double_t ntfx,Double_t ntfy,Double_t ntft,Double_t ntfpath,Double_t
    // ntfpx,Double_t ntfpy, Double_t ntfpz);
    R3BmTofDigi* AddHit(Int_t ntmul,
                        Double_t ntfx,
                        Double_t ntfy,
                        Double_t ntft,
                        Double_t ntfpath,
                        Double_t ntfpx,
                        Double_t ntfpy,
                        Double_t ntfpz,
                        Double_t ntfe);

  protected:
    TClonesArray* fmTofPoints;
    TClonesArray* fmTofMCTrack;
    TClonesArray* fmTofDigi;

    // Parameter class
    R3BmTofDigiPar* fmTofDigiPar;

    //- Control Hitograms

    TH1F* NtfXhis;
    TH1F* fMult;
    TH1F* fEloss;

    Int_t eventNomTof;

  private:
    virtual void SetParContainers();

    ClassDef(R3BmTofDigitizer, 1);
};

#endif
