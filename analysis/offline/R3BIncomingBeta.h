/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#ifndef R3BIncomingBeta_H
#define R3BIncomingBeta_H 1

// ROOT header
#include "TMath.h"
#include <TArrayF.h>
#include <array>
#include <vector>

// R3B headers
#include "R3BEventHeader.h"
#include "R3BFrsData.h"
#include "R3BFrsSciPosCalData.h"
#include "R3BMusicHitData.h"

// FAIR headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "FairTask.h"

class R3BIncomingIDPar;
class TClonesArray;
class R3BEventHeader;
class R3BCoarseTimeStitch;

class R3BIncomingBeta : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BIncomingBeta();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BIncomingBeta(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BIncomingBeta();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    virtual InitStatus ReInit();
    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    virtual void SetParContainers();

    virtual void Reset();

    // Accessor to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }
    void SetUseTref() { fUseTref = kTRUE; }
    void SetUseMultHit() { fUseMultHit = kTRUE; }
    void SetNumDets(int val) { fNumDet = val; }
    void SetStaId(int val) { fStaId = val; }
    void SetStoId(int val) { fStoId = val; } // if id == 0 and bFrsSci = true, then LOS will be used
    void SetLosRefCh(int val) { fLosRefCh = val; }
    void UseFrsSci(bool val) { fUseFrsSci = val; } // If false, sci2 will be used
    void SetLosCalRange(Double_t low, Double_t high)
    {
        fLosCalTrig_Low = low;
        fLosCalTrig_High = high;
    } // Only for S118/S091

  protected:
    R3BEventHeader* fHeader{}; // Event header

  private:
    void SetParameter();
    R3BCoarseTimeStitch* fTimeStitch;
    R3BIncomingIDPar* fIncomingID_Par; // Parameter container
    TClonesArray* fFrsDataCA;          /**< Array with FRS-output data. >*/

    TClonesArray* fHitSci2;
    TClonesArray* fPosCalFrsSci;
    TClonesArray* fHitLos;
    TClonesArray* fCalLos;
    TClonesArray* fCalLosTrig;
    TClonesArray* fTcalSci2; /**< Array with Tcal items. */

    int fStaId;
    int fStoId;
    int fLosRefCh;
    Double_t fLosCalTrig_Low, fLosCalTrig_High;

    Bool_t fOnline; // Don't store data for online

    UInt_t fNumDet;
    TArrayF* fToFoffset;
    TArrayF *fPosS2Left, *fPosS2Right;
    TArrayF *fTof2InvV_p0, *fTof2InvV_p1;
    Float_t fBeta_max, fBeta_min;
    Bool_t fUseTref;
    Bool_t fUseMultHit;
    bool fUseFrsSci;

    /** Private method FrsData **/
    //** Adds a FrsData to the analysis
    R3BFrsData* AddData(Int_t StaId,
                        Int_t StoId,
                        Double_t z,
                        Double_t aq,
                        Double_t betaval,
                        Double_t brhoval,
                        Double_t xs2,
                        Double_t xc,
                        Double_t tof);

  public:
    ClassDef(R3BIncomingBeta, 1)
};

#endif /* R3BIncomingBeta_H */
