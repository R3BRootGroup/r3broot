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

#ifndef R3BRPCHITPAR_H
#define R3BRPCHITPAR_H
// 41 Strips +  5 Bars
#define N_STRIP_NB (41 + 5)

#include <TObjString.h>
#include <TVector3.h>

#include "FairParGenericSet.h"
#include "FairParamList.h"
#include "TArrayF.h"

class R3BRpcHitPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    R3BRpcHitPar(const char* name = "RpcHitPar",
                 const char* title = "RPC Hit Finder Parameters",
                 const char* context = "TestDefaultContext");
    /** Destructor **/
    virtual ~R3BRpcHitPar();

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
    const Double_t GetNumChannels() { return fNumChannels; }
    TArrayF* GetCalParams1() { return fHitCalParams1; }
    TArrayF* GetCalParams2() { return fHitCalParams2; }
    TArrayF* GetCalParams3() { return fHitCalParams3; }
    TArrayF* GetCalParams4() { return fHitCalParams4; }

    void SetNumChannels(Int_t numberCha) { fNumChannels = numberCha; }
    void SetCalParams1(Float_t cc, Int_t ii) { fHitCalParams1->AddAt(cc, ii); }
    void SetCalParams2(Float_t cc, Int_t ii) { fHitCalParams2->AddAt(cc, ii); }
    void SetCalParams3(Float_t cc, Int_t ii) { fHitCalParams3->AddAt(cc, ii); }
    void SetCalParams4(Float_t cc, Int_t ii) { fHitCalParams4->AddAt(cc, ii); }

  private:
    TArrayF* fHitCalParams1; /*< Calibration Parameters of Pos 1>*/
    TArrayF* fHitCalParams2; /*< Calibration Parameters of Pos 2>*/
    TArrayF* fHitCalParams3; /*< Calibration Parameters of Pos 3>*/
    TArrayF* fHitCalParams4; /*< Calibration Parameters of Pos 4>*/
    Int_t fNumChannels;      /*< number of channels>*/

    const R3BRpcHitPar& operator=(const R3BRpcHitPar&);
    R3BRpcHitPar(const R3BRpcHitPar&);

    ClassDef(R3BRpcHitPar, 1); //
};

#endif /* !R3BRPCHITPAR_H*/
