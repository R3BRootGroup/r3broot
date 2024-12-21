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

#ifndef R3BRPCMAPPED2PRECALPAR_H
#define R3BRPCMAPPED2PRECALPAR_H

#include "FairTask.h"
#include "TH1F.h"

class TClonesArray;
class R3BEventHeader;

class R3BTCalEngine;
class R3BTCalPar;

class R3BRpcMapped2PreCalPar : public FairTask
{

  public:
    /** Default constructor **/
    R3BRpcMapped2PreCalPar();

    /** Standard constructor **/
    R3BRpcMapped2PreCalPar(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BRpcMapped2PreCalPar();

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

    R3BTCalPar* fTCalPar;        /**< Parameter container. */
    R3BTCalEngine* fEngine;      /**< Instance of the TCAL engine. */
    TClonesArray* fMappedDataCA; /**< Array with RPC Mapped-input data. >*/

  public:
    ClassDef(R3BRpcMapped2PreCalPar, 1);
};

#endif
