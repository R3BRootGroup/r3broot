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
// -----                      R3BRoluCal2Hit                    -----
// -----            Created July 2019 by A. Kelic-Heil          -----
// ----- Convert calibrated data to hit level (single time)     -----
// -----                  Following R3BLosCal2Hit               -----
// ------------------------------------------------------------------

#ifndef R3BRoluCAL2HIT_H
#define R3BRoluCAL2HIT_H

#include <map>

#include "FairTask.h"

#define VFTX_CLOCK_MHZ 200

class TClonesArray;
class TH1F;
class TH2F;
/*
 * TODO: This explanation is humbug.
 * An analysis task to apply TCAL calibration for NeuLAND.
 * This class reads NeuLAND mapped items with TDC values and
 * produces time items with time in [ns]. It requires TCAL
 * calibration parameters, which are produced in a separate
 * analysis run containing R3BRoluCal2HitFill task.
 */
class R3BRoluCal2Hit : public FairTask
{

  public:
    /*
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BRoluCal2Hit();

    /*
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BRoluCal2Hit(const char* name, int iVerbose = 1);

    /*
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BRoluCal2Hit();

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
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    // Accessor to select online mode
    void SetOnline(bool option) { fOnline = option; }

    // Accessor to select online mode
    void SetSkipTrigger(bool skip) { fSkipTrigger = skip; }

    inline void SetNofModules(int nDets, int nChs)
    {
        fNofDetectors = nDets;
        fNofChannels = nChs; //=4
    }

  private:
    TClonesArray* fCalItems{}; /* < Array with Cal items - input data. */
    TClonesArray* fCalTriggerItems{};
    TClonesArray* fHitItems;                         /* < Array with Hit items - output data. */
    double fClockFreq = 1. / VFTX_CLOCK_MHZ * 1000.; /* < Clock cycle in [ns]. */
    unsigned int fNofDetectors = 1;                  /**< Number of detectors. */
    unsigned int fNofChannels = 4;                   /**< Number of channels per detector. */
    int fnEvents = 0;

    TH1F* fhQ_R[2]{};
    TH1F* fhQ_O[2]{};
    TH1F* fhQ_L[2]{};
    TH1F* fhQ_U[2]{};

    int maxevent = 100000;
    // Don't store data for online
    bool fOnline = false;
    bool fSkipTrigger = false;

  public:
    ClassDef(R3BRoluCal2Hit, 2)
};

#endif // R3BRoluCAL2HIT_H
