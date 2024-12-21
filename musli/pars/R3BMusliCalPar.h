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

#ifndef R3BMusliCalPar_H
#define R3BMusliCalPar_H 1

#include "FairParGenericSet.h"

#include "TArrayF.h"
#include "TArrayI.h"
#include <Rtypes.h>
#include <stdint.h>
#include <vector>

class FairParamList;

class R3BMusliCalPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BMusliCalPar(const char* name = "musliCalPar",
                   const char* title = "MUSLI Cal Parameters",
                   const char* context = "MusliCalParContext");

    /** Destructor **/
    virtual ~R3BMusliCalPar();

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
    const Int_t GetNumSignals() { return fNumSignals; }
    const Int_t GetNumGroupsAnodes() { return fNumGroupsAnodes; }
    const Int_t GetNumParamsEneFit() { return fNumParamsEneFit; }
    const Int_t GetNumParamsPosFit() { return fNumParamsPosFit; }
    const Int_t GetNumParamsMulHit() { return fNumParamsMultHit; }
    const Int_t GetMaxMult() { return fMaxMult; }
    const Int_t GetInUse(Int_t signal) { return fIn_use->GetAt(signal - 1); }
    TArrayF* GetEneCalParams() { return fEneCalParams; }
    TArrayF* GetPosCalParams() { return fPosCalParams; }
    TArrayF* GetMultHitCalParams() { return fMultHitCalParams; }

    void SetNumSignals(Int_t n) { fNumSignals = n; }
    void SetNumGroupsAnodes(Int_t n) { fNumGroupsAnodes = n; }
    void SetNumParamsEneFit(Int_t n) { fNumParamsEneFit = n; }
    void SetNumParamsPosFit(Int_t n) { fNumParamsPosFit = n; }
    void SetNumParamsMultHit(Int_t n) { fNumParamsMultHit = n; }
    void SetMaxMult(Int_t max) { fMaxMult = max; }
    void SetInUse(Int_t value, Int_t signal) { fIn_use->AddAt(value, signal - 1); }
    void SetEneParams(Float_t val, Int_t signal, Int_t indexpar)
    {
        // 0-base
        int index = (signal - 1) * fNumParamsEneFit + indexpar;
        fEneCalParams->AddAt(val, index);
    }
    void SetPosParams(Float_t val, Int_t signal, Int_t indexpar)
    {
        // 0-base
        int index = (signal - 1) * fNumParamsPosFit + indexpar;
        fPosCalParams->AddAt(val, index);
    }
    void SetMultHitParams(Float_t val, Int_t signal, Int_t indexpar)
    {
        // 0-base
        int index = (signal - 1) * fNumParamsMultHit + indexpar;
        fMultHitCalParams->AddAt(val, index);
    }

  private:
    TArrayF* fEneCalParams;     // Energy calibration parameters for anodes
    TArrayF* fPosCalParams;     // Position calibration parameters for anodes
    TArrayF* fMultHitCalParams; // Calibration parameters for using multihits
    TArrayI* fIn_use;           // 1: anode ready, 0:otherwise
    Int_t fNumSignals;          // number of signals at the Mapped level (group of anodes + Tref + Ttrig)
    Int_t fNumGroupsAnodes;     // number of groups of anodes signals (each of them have energy and time)
    Int_t fNumParamsEneFit;     // number of cal parameters in the fit
    Int_t fNumParamsPosFit;     // number of parameters for position calibration
    Int_t fNumParamsMultHit;    // number of parameters for Multhit Calibration
    Int_t fMaxMult;             // Max. multiplicity per anode

    const R3BMusliCalPar& operator=(const R3BMusliCalPar&);
    R3BMusliCalPar(const R3BMusliCalPar&);

  public:
    ClassDef(R3BMusliCalPar, 1);
};

#endif /* R3BMusliCalPar_H */
