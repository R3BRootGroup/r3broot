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

#ifndef R3BPTOFCAL2HITPAR_H
#define R3BPTOFCAL2HITPAR_H

#include "FairTask.h"

#define PtofPaddlesPerPlane 8

class TClonesArray;
class R3BPtofHitPar;
class TH1D;

class R3BPtofCal2HitPar : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BPtofCal2HitPar();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BPtofCal2HitPar(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BPtofCal2HitPar();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

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

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    /**
     * Method for setting the beam charge
     * @param charge the charge in terms of e
     */
    void SetBeamCharge(Int_t charge) { fBeamCharge = charge; }

    /**
     * Method for setting the y-position of the sweep.
     * If y != 0 this run will determine veff.
     * @param yPositon the y-position of the sweep
     */
    void SetYPosition(Double_t yPosition) { fYPosition = yPosition; }

  private:
    Int_t fBeamCharge;
    Int_t fYPosition;

    TH1D* fZScaleHistos[2 * PtofPaddlesPerPlane];
    TH1D* fTOffsetHistos[2 * PtofPaddlesPerPlane];

    R3BPtofHitPar* fPar;

    Int_t fNEvents;         /**< Event counter. */
    TClonesArray* fCalData; /**< Array with cal data - input data. */

  public:
    ClassDef(R3BPtofCal2HitPar, 1)
};

#endif
