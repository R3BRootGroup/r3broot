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
// -----          R3BMusicHitPar source file                -----
// -----     Created 24/11/19  by J.L. Rodriguez-Sanchez    -----
// --------------------------------------------------------------

#ifndef R3BMusicHitPar_H
#define R3BMusicHitPar_H 1

#include "FairParGenericSet.h" // for FairParGenericSet

#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TObject.h"
#include <TObjString.h>

class FairParamList;

class R3BMusicHitPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BMusicHitPar(const char* name = "musicHitPar",
                   const char* title = "MUSIC Hit Parameters",
                   const char* context = "MusicHitParContext");

    /** Destructor **/
    virtual ~R3BMusicHitPar();

    /** Method to reset all parameters **/
    virtual void clear();

    /** Method to store all parameters using FairRuntimeDB **/
    virtual void putParams(FairParamList* list);

    /** Method to retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list);

    /** Method to print values of parameters to the standard output **/
    void printParams();

    /** Accessor functions **/
    Int_t GetNumParZFit() const { return fNumParamsZFit; }
    Int_t GetNumAnodes() const { return fNumAnodes; }
    Int_t GetInUse(Int_t anode) const { return fIn_use->GetAt(anode - 1); }
    Float_t GetAnodePos(Int_t anode) const { return fAnode_pos->GetAt(anode - 1); }
    TArrayF* GetZHitPar() const { return fDetZHitParams; }
    TArrayF* GetAngCorPar() const { return fAngCorParams; }

    void SetNumParZFit(Int_t nbParams) { fNumParamsZFit = nbParams; }
    void SetNumAnodes(Int_t nbAnodes) { fNumAnodes = nbAnodes; }
    void SetInUse(Int_t value, Int_t anode) { fIn_use->AddAt(value, anode - 1); }
    void SetZHitPar(Double_t cc, Int_t ii) { fDetZHitParams->AddAt(cc, ii); }
    void SetAngCorPar(Double_t cc, Int_t ii) { fAngCorParams->AddAt(cc, ii); }
    void SetAnodePos(Float_t value, Int_t anode) { fAnode_pos->AddAt(value, anode - 1); }

    // Create more Methods if you need them!

  private:
    TArrayF* fDetZHitParams; // Calibration Parameters for charge Z
    TArrayI* fIn_use;        // 1: anode ready, 0:otherwise
    TArrayF* fAnode_pos;     // Position of each anode along the beam direction
    Int_t fNumAnodes;        // Number of anodes
    Int_t fNumParamsZFit;    // number of hit parameters in the fit for charge Z
    TArrayF* fAngCorParams;

    const R3BMusicHitPar& operator=(const R3BMusicHitPar&); /*< an assignment operator>*/

    R3BMusicHitPar(const R3BMusicHitPar&); /*< a copy constructor >*/

  public:
    ClassDef(R3BMusicHitPar, 1);
};

#endif /* R3BMusicHitPar_H */
