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

#ifndef R3BRPCPRECAL2CALPAR_H
#define R3BRPCPRECAL2CALPAR_H

#include "FairTask.h"
#include "R3BRpcTotCalPar.h"
#include "TH1F.h"

class TClonesArray;
class R3BRpcTotCalPar;
class R3BEventHeader;

class R3BRpcPreCal2CalPar : public FairTask
{

  public:
    /** Default constructor **/
    R3BRpcPreCal2CalPar();

    /** Standard constructor **/
    R3BRpcPreCal2CalPar(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BRpcPreCal2CalPar();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Accessor functions **/
    const Int_t GetNumChannels() { return fNumChannels; }

    void SetNumChannels(Int_t numberCha) { fNumChannels = numberCha; }

    void SetDebugMode(Bool_t debug) { fDebugMode = debug; }

  private:
    Bool_t fDebugMode;
    Int_t fNumChannels;

    R3BEventHeader* fR3BEventHeader;
    R3BRpcTotCalPar* fTotCalPar;

    TH1F* fhTot[N_NUM];

    TClonesArray* fPreCalDataCA; /**< Array with RPC PreCal-input data. >*/

  public:
    ClassDef(R3BRpcPreCal2CalPar, 1);
};

#endif
