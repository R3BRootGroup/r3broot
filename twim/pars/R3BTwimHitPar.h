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
// -----            R3BTwimHitPar source file               -----
// -----     Created 24/11/19 by J.L. Rodriguez-Sanchez     -----
// --------------------------------------------------------------

#ifndef R3BTwimHitPar_H
#define R3BTwimHitPar_H 1

#include "FairParGenericSet.h"
#include "R3BTsplinePar.h"
#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TSpline.h"
#include <TObjString.h>
#include <vector>

#include "FairDetParIo.h"
#include "FairLogger.h"
#include "FairParamList.h"
#include "TArrayF.h"
#include "TMath.h"
#include "TString.h"

#include <iostream>

class FairParamList;

class R3BTwimHitPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BTwimHitPar(const char* name = "twimHitPar",
                  const char* title = "TWIM Hit Parameters",
                  const char* context = "TwimHitParContext");

    /** Destructor **/
    virtual ~R3BTwimHitPar();

    /** Method to reset all parameters **/
    virtual void clear();

    /** Method to store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Method to retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list);

    /** Method to print values of parameters to the standard output **/
    virtual void print();
    void printParams();

    /** Accessor functions **/
    const Int_t GetNumSec() { return fNumSec; }
    const Int_t GetNumParZFit() { return fNumParamsZFit; }
    const Int_t GetNumAnodes() { return fNumAnodes; }
    const Int_t GetInUse(Int_t sec, Int_t anode) { return fIn_use->GetAt((sec - 1) * 16 + anode - 1); }
    const Float_t GetAnodePos(Int_t anode) { return fAnode_pos->GetAt(anode - 1); }
    Float_t GetEmean_tof(Int_t sec) { return fEmean_tof[sec - 1]; }
    Float_t GetEmean_dt(Int_t sec) { return fEmean_dt[sec - 1]; }
    TArrayF* GetZHitPar() { return fDetZHitParams; }
    TArrayF* GetZTofHitPar() { return fTofHitParams; }
    TSpline3* GetSpline(Int_t sec) { return fSpline[sec - 1]->GetSpline(); }

    TArrayF* GetTofCorrectionPar(Int_t sectwim)
    {
        // sectwim 1 for left and 2 for right
        TArrayF* fPar = new TArrayF(3);
        fPar->AddAt(fTofHitParams->GetAt((sectwim - 1) * 3), 0);
        fPar->AddAt(fTofHitParams->GetAt((sectwim - 1) * 3 + 1), 1);
        fPar->AddAt(fTofHitParams->GetAt((sectwim - 1) * 3 + 2), 2);
        return fPar;
    }

    void SetNumSec(Int_t nbsec) { fNumSec = nbsec; }
    void SetNumAnodes(Int_t nbAnodes) { fNumAnodes = nbAnodes; }
    void SetNumParZFit(Int_t nbParams) { fNumParamsZFit = nbParams; }
    void SetInUse(Int_t value, Int_t sec, Int_t anode) { fIn_use->AddAt(value, (sec - 1) * 16 + anode - 1); }
    void SetZHitPar(Double_t cc, Int_t ii) { fDetZHitParams->AddAt(cc, ii); }
    void SetZTofHitPar(Double_t cc, Int_t ii) { fTofHitParams->AddAt(cc, ii); }
    void SetAnodePos(Float_t value, Int_t anode) { fAnode_pos->AddAt(value, anode - 1); }
    void SetEmean_tofParam(Float_t val, Int_t sec) { fEmean_tof[sec - 1] = val; }
    void SetEmean_dtParam(Float_t val, Int_t sec) { fEmean_dt[sec - 1] = val; }
    void SetSpline(TSpline3* spline, Int_t nbsec) { fSpline[nbsec - 1]->SetSpline(spline); }

    void SetTofCorrectionPar(Float_t value[3], Int_t sectwim)
    {
        // sectwim 1 for left and 2 for right
        fTofHitParams->AddAt(value[0], (sectwim - 1) * 3);
        fTofHitParams->AddAt(value[1], (sectwim - 1) * 3 + 1);
        fTofHitParams->AddAt(value[2], (sectwim - 1) * 3 + 2);
    }

    // Create more Methods if you need them!

  private:
    TArrayF* fDetZHitParams; // Calibration Parameters for charge Z
    TArrayF* fTofHitParams;  // Tof calibration parameters for charge Z
    TArrayI* fIn_use;        // 1: anode ready, 0:otherwise
    TArrayF* fAnode_pos;     // Position of each anode along the beam direction
    Int_t fNumSec;
    Int_t fNumAnodes;     // Number of anodes
    Int_t fNumParamsZFit; // number of hit parameters in the fit for charge Z

    std::vector<Float_t> fEmean_tof;
    std::vector<Float_t> fEmean_dt;

    std::vector<R3BTsplinePar*> fSpline;

    const R3BTwimHitPar& operator=(const R3BTwimHitPar&); /*< an assignment operator>*/
    R3BTwimHitPar(const R3BTwimHitPar&);                  /*< a copy constructor >*/

    ClassDef(R3BTwimHitPar, 1);
};

#endif
