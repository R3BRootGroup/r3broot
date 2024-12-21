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

#ifndef R3BRPCCAL2HITPAR_H
#define R3BRPCCAL2HITPAR_H

#include "FairTask.h"
#include "R3BRpcHitPar.h"
#include "TGraph.h"
#include "TH1F.h"

#define CSCINT (140.)

class TClonesArray;
class R3BEventHeader;

class R3BRpcCal2HitPar : public FairTask
{

  public:
    /** Default constructor **/
    R3BRpcCal2HitPar();

    /** Standard constructor **/
    R3BRpcCal2HitPar(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BRpcCal2HitPar();

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

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /** Accessor functions **/
    const Int_t GetNumChannels() { return fNumChannels; }

    void SetNumChannels(Int_t numberCha) { fNumChannels = numberCha; }

    void SetDebugMode(Bool_t debug) { fDebugMode = debug; }

    void SetNBCalib() { fRpcCalib = false; }

  private:
    void CalculateParsStrip();
    void CalculateParsPmt();
    void SetParameter();
    Bool_t fDebugMode;
    Bool_t fRpcCalib;
    Int_t fNumChannels;

    TH1F* fhPos[N_STRIP_NB];
    TH1F* fhTime[N_STRIP_NB];

    R3BRpcHitPar* fHitPar;    /**< Container for Hit parameters. >*/
    TClonesArray* fCalDataCA; /**< Array with Cal RPC Strip - input data. >*/

  public:
    ClassDef(R3BRpcCal2HitPar, 1);
};

#endif
