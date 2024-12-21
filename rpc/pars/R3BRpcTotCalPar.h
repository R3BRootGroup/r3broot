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

#ifndef R3BRPCTOTCALPAR_H
#define R3BRPCTOTCALPAR_H
// 41 Strips +  5 Bars
#define N_NUM (41 * 2 + 5 * 2)

#include "FairParGenericSet.h"
#include "TObject.h"

#include "TArrayF.h"
#include "TObjArray.h"
#include <TObjString.h>

class FairParamList;

class R3BRpcTotCalPar : public FairParGenericSet
{

  public:
    /** Standard constructor **/
    R3BRpcTotCalPar(const char* name = "RpcTotCalPar",
                    const char* title = "RPC Cal Parameters",
                    const char* context = "R3BRpcTotCalParContext");

    /** Destructor **/
    virtual ~R3BRpcTotCalPar();

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
    TArrayF* GetCalParams() { return fCalParams; }

    void SetNumChannels(Int_t numberCha) { fNumChannels = numberCha; }
    void SetCalParams(Float_t cc, Int_t ii) { fCalParams->AddAt(cc, ii); }

    /** Create more Methods if you need them! **/

  private:
    TArrayF* fCalParams; /*< Calibration Parameters of Channels>*/
    Int_t fNumChannels;  /*< number of channels>*/

    const R3BRpcTotCalPar& operator=(const R3BRpcTotCalPar&);
    R3BRpcTotCalPar(const R3BRpcTotCalPar&);

    ClassDef(R3BRpcTotCalPar, 1);
};

#endif /* R3BRPCTotCALPAR_H */
