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

#ifndef R3BRPCDIGITIZER_H
#define R3BRPCDIGITIZER_H

#include "FairTask.h"
#include "R3BRpc.h"
#include "R3BRpcCalData.h"
#include "R3BRpcPars4Sim.h"
#include "R3BRpcPoint.h"
#include "Rtypes.h"

class TClonesArray;

class R3BRpcDigitizer : public FairTask
{
  public:
    /** Standard contructor **/
    R3BRpcDigitizer();

    /** Destructor **/
    virtual ~R3BRpcDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method Reset **/
    virtual void Reset();

    virtual void SetParContainers();

  private:
    void SetParameter();

    TClonesArray* fRpcPointDataCA; //!  The RPC hit collection
    TClonesArray* fRpcCalDataCA;   /**< Array with RPC Cal- output data. >*/

    R3BRpcPars4Sim* fSim_Par; // Parameter Container for a Realistic Simulation

    UInt_t fNumberOfChannels;

    /** Private method AddCal
     **
     ** Adds a RpcCal data
     **/
    R3BRpcCalData* AddCal(Int_t ident, Double_t energy, ULong64_t time, Double_t tot_energy = 0.);

    inline void ResetParameters(){};

  public:
    ClassDef(R3BRpcDigitizer, 1);
};

#endif /* R3BRPCDIGITIZER_H */
