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

#ifndef R3BAnalysisIncomingID_H
#define R3BAnalysisIncomingID_H 1

// ROOT headers
#include "TCutG.h"
#include "TMath.h"
#include <TArrayF.h>

// FAIR headers
#include "FairTask.h"

class R3BIncomingIDPar;
class TClonesArray;
class R3BEventHeader;
class R3BTcutPar;
class R3BCoarseTimeStitch;

/**
 * This taks reads all detector data items for the analysis of incoming
 * projectiles from FRS.
 */
class R3BAnalysisIncomingID : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BAnalysisIncomingID();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BAnalysisIncomingID(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BAnalysisIncomingID();

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

    void SetBetaCorrectionForZ(Double_t p0, Double_t p1, Double_t p2, Double_t Zprimary, Double_t Zoffset)
    {
        fP0 = p0;
        fP1 = p1;
        fP2 = p2;
        fZprimary = Zprimary;
        fZoffset = Zoffset;
    }

    // Accessor to select online mode
    void SetOnline(Bool_t option) { fOnline = option; }

    // Accessor to select the MUSIC for the incoming ID
    void SetMusicForPID() { fUseLOS = kFALSE, fUsePspx1 = kFALSE; }

    // Accessor to select the LOS for the incoming ID
    void SetLosForPID() { fUseLOS = kTRUE, fUsePspx1 = kFALSE; }

    // Acsessor to set use of trigger corrected times
    void SetUseTref() { fUseTref = kTRUE; }

    void SetNumDet(int val) { fNumDet = val; }

  protected:
    R3BEventHeader* fHeader{}; // Event header

  private:
    void SetParameter();
    R3BCoarseTimeStitch* fTimeStitch;
    R3BIncomingIDPar* fIncomingID_Par; // Parameter container
    TClonesArray* fHitItemsMus;
    TClonesArray* fHitItemsMusli;
    TClonesArray* fFrsDataCA; /**< Array with FRS-output data. >*/
    TClonesArray* fHitLos;
    TClonesArray* fTriggerLos;
    TClonesArray* fHitPspx1_x;
    TClonesArray* fHitPspx1_y;

    Bool_t fOnline;            // Don't store data for online
    Bool_t fUseLOS, fUsePspx1; // Use LOS or PSPX1 charge (otherwise MUSIC charge)
    Double_t fP0, fP1, fP2, fZprimary, fZoffset;
    Bool_t fUseTref; // Use trigger corrected times

    Double_t fPos_p0;
    Double_t fPos_p1;

    Int_t fNumDet;
    TArrayF* fToFoffset;
    TArrayF *fPosS2Left, *fPosS2Right;
    TArrayF *fTof2InvV_p0, *fTof2InvV_p1;
    TArrayF* fBrho0_S2toCC;
    TArrayF* fDispersionS2;
    Float_t fx0_point, fy0_point, frot_ang;
    Float_t fx0_Aq, fy0_Aq, fang_Aq;
    Float_t fBeta_max, fBeta_min;
    TCutG *fCutS2, *fCutCave;

  public:
    ClassDef(R3BAnalysisIncomingID, 1)
};

#endif /* R3BAnalysisIncomingID_H */
