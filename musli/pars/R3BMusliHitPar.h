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

#ifndef R3BMusliHitPar_H
#define R3BMusliHitPar_H 1

#include "FairParGenericSet.h"

#include "TArrayD.h"
#include <Rtypes.h>
#include <stdint.h>
#include <vector>

class FairParamList;

class R3BMusliHitPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BMusliHitPar(const char* name = "musliHitPar",
                   const char* title = "MUSLI Hit Parameters",
                   const char* context = "MusliHitParContext");

    /** Destructor **/
    virtual ~R3BMusliHitPar();

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
    const Int_t GetNumGroupsAnodes() { return fNumGroups; }
    const Int_t GetNumTypes() { return fNumTypes; }
    const Int_t GetMaxMult() { return fMaxMult; }
    const Int_t GetNumParamsEcorrBetaFit() { return fNumParamsEcorrBetaFit; }
    const Int_t GetNumParamsZFit() { return fNumParamsZFit; }
    TArrayD* GetEaveVsBetaHitParams() { return fEaveVsBetaHitParams; }
    TArrayD* GetEcorrBetaHitParams() { return fEcorrBetaHitParams; }
    TArrayD* GetZHitParams() { return fZHitParams; }

    void SetNumGroupsAnodes(Int_t n) { fNumGroups = n; }
    void SetNumTypes(Int_t n) { fNumTypes = n; }
    void SetMaxMult(Int_t max) { fMaxMult = max; }
    void SetNumParamsEcorrBetaFit(Int_t n) { fNumParamsEcorrBetaFit = n; }
    void SetNumParamsZFit(Int_t n) { fNumParamsZFit = n; }
    void SetEaveVsBetaHitParams(Double_t val, Int_t type) { fEaveVsBetaHitParams->AddAt(val, type - 1); }
    void SetEcorrBetaHitParams(Double_t val, Int_t type, Int_t deg)
    {
        // 0 - based
        int index = (type - 1) * fNumParamsEcorrBetaFit + deg - 1;
        fEcorrBetaHitParams->AddAt(val, index);
    }
    void SetZHitParams(Double_t val, Int_t type, Int_t deg)
    {
        // 0 - based
        int index = (type - 1) * fNumParamsZFit + deg - 1;
        fZHitParams->AddAt(val, index);
    }

  private:
    Int_t fNumGroups; // number of anodes groups at CAL level
    Int_t fNumTypes;  // number of types (1 type per group of anodes)
    Int_t fMaxMult;   // maximum number of hits per group of anodes
    Int_t fNumParamsEcorrBetaFit;
    Int_t fNumParamsZFit;

    TArrayD* fEaveVsBetaHitParams;
    TArrayD* fEcorrBetaHitParams;
    TArrayD* fZHitParams;

    const R3BMusliHitPar& operator=(const R3BMusliHitPar&);
    R3BMusliHitPar(const R3BMusliHitPar&);

  public:
    ClassDef(R3BMusliHitPar, 1);
};

#endif /* R3BMusliHitPar_H */
