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

// ------------------------------------------------------------------
// -----                  R3BSci8Cal2Hit                         -----
// -----            Created August 7th 2018 by A. Kelic-Heil     -----
// ----- Convert time calibrated data to hit level (single time) ----
// ----- Following R3BLosCal2Hit
// ------------------------------------------------------------------

#ifndef R3BSCI8CAL2HIT
#define R3BSCI8CAL2HIT

#include <map>

#include "FairTask.h"

class TClonesArray;
class TH1F;
class TH2F;
/*
 * TODO: This explanation is humbug.
 * An analysis task to apply TCAL calibration for NeuLAND.
 * This class reads NeuLAND mapped items with TDC values and
 * produces time items with time in [ns]. It requires TCAL
 * calibration parameters, which are produced in a separate
 * analysis run containing R3BSci8Cal2HitFill task.
 */
class R3BSci8Cal2Hit : public FairTask
{

  public:
    /*
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BSci8Cal2Hit();

    /*
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BSci8Cal2Hit(const char* name, Int_t iVerbose = 1);

    /*
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BSci8Cal2Hit();

    /*
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /*
     * Method for re-initialization of parameter containers
     * in case the Run ID has changed.
     */
    virtual InitStatus ReInit();

    /*
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    /*
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    virtual void FinishEvent();

    /*
     * Methods for setting position offset and effective velocity of light
     */
    inline void SetSci8ParamMCFD(Double_t offsetX, Double_t veffX)
    {
        fsci8OffsetX = offsetX;
        fsci8VeffX = veffX;
    }

    inline void SetSci8ParamTAMEX(Double_t offsetXT, Double_t veffXT)
    {
        fsci8OffsetXT = offsetXT;
        fsci8VeffXT = veffXT;
    }

    /*
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

  private:
    TClonesArray* fCalItems; /* < Array with Cal items - input data. */
    TClonesArray* fHitItems; /* < Array with Hit items - output data. */
    UInt_t fNofHitItems;     /* < Number of hit items for cur event. */
    Double_t fClockFreq;     /* < Clock cycle in [ns]. */
    Double_t fsci8VeffX;
    Double_t fsci8OffsetX;
    Double_t fsci8OffsetXT;
    Double_t fsci8VeffXT;

    TClonesArray* fMapped; /* < Array with mapped data - input data. */

    TH1F* fhTres_M;
    TH1F* fhTres_T;
    TH1F* fhQ;

    TH2F* fhQ_L;
    TH2F* fhQ_R;
    TH2F* fhQ_vs_X;
    TH2F* fhQ1_vs_Q2;

    Int_t Icount = 0;

  public:
    ClassDef(R3BSci8Cal2Hit, 1)
};

#endif
