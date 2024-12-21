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

// ---------------------------------------------------------------
// -----             R3BTwimCalPar source file               -----
// -----       Created 24/11/19 by J.L. Rodriguez-Sanchez    -----
// ---------------------------------------------------------------

#ifndef R3BTwimCalPar_H
#define R3BTwimCalPar_H 1

#include "FairParGenericSet.h"

#include "TArrayF.h"
#include "TArrayI.h"
#include <Rtypes.h>
#include <stdint.h>
#include <vector>

class FairParamList;

class R3BTwimCalPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BTwimCalPar(const char* name = "twimCalPar",
                  const char* title = "TWIM Cal Parameters",
                  const char* context = "TwimCalParContext");

    /** Destructor **/
    virtual ~R3BTwimCalPar();

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
    const Int_t GetNumSec() { return fNumSections; }
    const Int_t GetNumAnodes() { return fNumAnodes; }
    const Int_t GetNumParamsEFit() { return fNumParamsEFit; }
    const Int_t GetNumParamsPosFit() { return fNumParamsPosFit; }
    const Int_t GetNumAnodesTRef() { return fNumAnodesTRef; }
    const Int_t GetNumAnodesTrig() { return fNumAnodesTrig; }
    const Int_t GetMaxMult() { return fMaxMult; }
    const Int_t GetMinDT() { return fMinDT; }
    const Int_t GetMaxDT() { return fMaxDT; }
    const Int_t GetInUse(Int_t sec, Int_t anode) { return fIn_use[sec - 1]->GetAt(anode - 1); }
    TArrayF* GetAnodeCalEParams(Int_t sec) { return fAnodeECalParams[sec - 1]; }
    TArrayF* GetAnodeCalPosParams(Int_t sec) { return fAnodePosCalParams[sec - 1]; }

    void SetNumSec(Int_t numberSections) { fNumSections = numberSections; }
    void SetNumAnodes(Int_t numberAnodes) { fNumAnodes = numberAnodes; }
    void SetNumParamsEFit(Int_t numberParams) { fNumParamsEFit = numberParams; }
    void SetNumParamsPosFit(Int_t numberParams) { fNumParamsPosFit = numberParams; }
    void SetNumAnodesTRef(Int_t nb) { fNumAnodesTRef = nb; }
    void SetNumAnodesTrig(Int_t nb) { fNumAnodesTrig = nb; }
    void SetMaxMult(Int_t max) { fMaxMult = max; }
    void SetMinDT(Int_t minDT) { fMinDT = minDT; }
    void SetMaxDT(Int_t maxDT) { fMaxDT = maxDT; }
    void SetInUse(Int_t value, Int_t sec, Int_t anode) { fIn_use[sec - 1]->AddAt(value, anode - 1); }
    void SetAnodeCalParams(Float_t val, Int_t sec, Int_t anode, Int_t indexpar)
    {
        // 1-base
        int index = (anode - 1) * fNumParamsEFit + indexpar - 1;
        fAnodeECalParams[sec - 1]->AddAt(val, index);
    }
    void SetPosParams(Float_t val, Int_t sec, Int_t anode, Int_t indexpar)
    {
        // 1-base
        int index = (anode - 1) * fNumParamsPosFit + indexpar - 1;
        fAnodePosCalParams[sec - 1]->AddAt(val, index);
    }

  private:
    std::vector<TArrayF*> fAnodeECalParams;   // Energy calibration parameters for anodes
    std::vector<TArrayF*> fAnodePosCalParams; // Position calibration parameters for anodes
    std::vector<TArrayI*> fIn_use;            // 1: anode ready, 0:otherwise
    Int_t fNumSections;                       // number of sections
    Int_t fNumAnodes;                         // number of anodes
    Int_t fNumParamsEFit;                     // number of cal parameters in the fit
    Int_t fNumParamsPosFit;                   // number of parameters for position calibration
    Int_t fNumAnodesTRef;                     // number of anodes for TREF
    Int_t fNumAnodesTrig;                     // number of anodes for Trigger
    Int_t fMaxMult;                           // Max. multiplicity per anode
    Int_t fMinDT;                             // Min. drift time accepted
    Int_t fMaxDT;                             // Max. drift time accepted

    const R3BTwimCalPar& operator=(const R3BTwimCalPar&);
    R3BTwimCalPar(const R3BTwimCalPar&);

  public:
    ClassDef(R3BTwimCalPar, 1);
};

#endif /* R3BTwimCalPar_H */
