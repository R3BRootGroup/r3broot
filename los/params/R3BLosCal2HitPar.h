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

// ------------------------------------------------------------
// -----                  R3BLosCal2HitPar                -----
// -----         Created February 2020 by A. Kelic-Heil   -----
// ------------------------------------------------------------

#ifndef R3BLOSCAL2HITPAR
#define R3BLOSCAL2HITPAR

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "TClonesArray.h"
#include "TMath.h"
#include <cstdlib>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BLosCal2HitPar : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BLosCal2HitPar();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BLosCal2HitPar(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BLosCal2HitPar();

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
    virtual void Fit();

    /**
     * Method fitting the data.
     * Is called by the framework after events are read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* option);

    virtual Double_t walk_correction(Int_t PMT, Double_t TOT, const Double_t* par);

    virtual Double_t calc_time_residual(const Double_t* par);
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
     * Method for setting required statistics per module.
     */

    inline void SetStats(Int_t nStats)
    {
        fStats = nStats;
        tvftx.resize(fStats, std::vector<Double_t>(8, 0));
        tlead.resize(fStats, std::vector<Double_t>(8, 0));
        ttrai.resize(fStats, std::vector<Double_t>(8, 0));
        energ.resize(fStats, std::vector<Double_t>(8, 0));
    }

    /* Method for setting which LOS detector to calibrate */
    inline void SetWhichLos(Int_t nDet) { fidentDet = nDet; }

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

  private:
    TClonesArray* fCalItems; /**< Array with Cal items - input data. */

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fTpat;
    Int_t fidentDet; /**< Number of LOS detectors. */
    Int_t fStats;
    Int_t NPM = 8;
    Int_t icount = 2;
    Int_t NPAR = NPM * icount;
    Double_t fClockFreq;
    TH1F* hfit;

    // creating the storage for energy and time for each PMTs
    std::vector<std::vector<Double_t>> tvftx; //(fStats, vector<Double_t> (8, 0));
    std::vector<std::vector<Double_t>> tlead; //(fStats, vector<Double_t> (8, 0));
    std::vector<std::vector<Double_t>> ttrai; //(fStats, vector<Double_t> (8, 0));
    std::vector<std::vector<Double_t>> energ; //(fStats, vector<Double_t> (8, 0));

    unsigned long fNEvents = 0, fNEvents_start = 0; /**< Event counter. */

  public:
    ClassDef(R3BLosCal2HitPar, 2)
};

#endif
