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

// ------------------------------------------------------------------
// -----             R3BMusicCalPar source file                 -----
// -----         Created 24/11/19  by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------------

#ifndef R3BMusicCalPar_H
#define R3BMusicCalPar_H 1

#include "FairParGenericSet.h" // for FairParGenericSet

#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TObject.h"
#include <TObjString.h>

class FairParamList;

class R3BMusicCalPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BMusicCalPar(const char* name = "musicCalPar",
                   const char* title = "MUSIC Cal Parameters",
                   const char* context = "MusicCalParContext");

    /** Destructor **/
    virtual ~R3BMusicCalPar();

    /** Method to reset all parameters **/
    virtual void clear();

    /** Method to store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Method to retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list);

    /** Method to print values of parameters to the standard output **/
    void printParams();

    /** Accessor functions **/
    Int_t GetNumAnodes() const { return fNumAnodes; }
    Int_t GetNumParamsEFit() const { return fNumParamsEFit; }
    Int_t GetNumParamsPosFit() const { return fNumParamsPosFit; }
    Int_t GetInUse(Int_t anode) const { return fIn_use->GetAt(anode - 1); }
    TArrayF* GetAnodeCalParams() const { return fAnodeCalParams; }
    TArrayF* GetPosParams() const { return fPosParams; }

    void SetNumAnodes(Int_t numberAnodes) { fNumAnodes = numberAnodes; }
    void SetNumParamsEFit(Int_t numberParams) { fNumParamsEFit = numberParams; }
    void SetNumParamsPosFit(Int_t numberParams) { fNumParamsPosFit = numberParams; }
    void SetInUse(Int_t value, Int_t anode) { fIn_use->AddAt(value, anode - 1); }
    void SetAnodeCalParams(Float_t cc, Int_t ii) { fAnodeCalParams->AddAt(cc, ii); }
    void SetPosParams(Float_t cc, Int_t ii) { fPosParams->AddAt(cc, ii); }

    // Create more Methods if you need them!

  private:
    TArrayF* fAnodeCalParams; // Energy calibration parameters for anodes
    TArrayF* fPosParams;      // Position calibration parameters for anodes
    TArrayI* fIn_use;         // 1: anode ready, 0:otherwise
    Int_t fNumAnodes;         // number of anodes
    Int_t fNumParamsEFit;     /* number of cal parameters in the fit
                                gaus: A_fit & B_fit & C_fit*/
    Int_t fNumParamsPosFit;   // number of parameters for position calibration

    const R3BMusicCalPar& operator=(const R3BMusicCalPar&); /*< an assignment operator>*/

    R3BMusicCalPar(const R3BMusicCalPar&); /*< a copy constructor >*/

  public:
    ClassDef(R3BMusicCalPar, 1);
};

#endif /* R3BMusicCalPar_H */
