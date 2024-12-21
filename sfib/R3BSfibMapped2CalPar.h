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

//
// Creates time calib parameters for bunched fiber detector.
//

#ifndef R3BSFIBMAPPED2CALPAR_H
#define R3BSFIBMAPPED2CALPAR_H

#include "FairTask.h"
#include <R3BTCalEngine.h>

class TClonesArray;
class R3BTCalPar;
class R3BEventHeader;

class R3BSfibMapped2CalPar : public FairTask
{
  public:
    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param a_name a name of the task.
     * @param a_verbose a FairTask verbosity level.
     * @param a_update_rate a update rate for online histograms.
     * @param a_min_stats a minimum statistics for calibration.
     */
    R3BSfibMapped2CalPar(Int_t,
                         enum R3BTCalEngine::CTDCVariant = R3BTCalEngine::CTDC_8_12_FWD_250,
                         Int_t = 1e6,
                         Int_t = 1e5);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BSfibMapped2CalPar();

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
    virtual void Exec(Option_t*);

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
     * Method for setting the update rate for control histograms
     * @param rate an update rate value (events).
     */
    virtual void SetUpdateRate(Int_t);

    /**
     * Method for setting minimum required statistics per module.
     * Only detector modules with number of entries in TDC
     * distribution greater than minimum statistics will be
     * calibrated.
     * @param minStats a value of minimum statistics required.
     */
    virtual void SetMinStats(Int_t);

  private:
    enum R3BTCalEngine::CTDCVariant fCTDCVariant;
    TClonesArray* fMapped;
    R3BTCalPar* fTCalPar;
    R3BTCalEngine* fEngine;
    Int_t fUpdateRate;
    Int_t fMinStats;

  public:
    ClassDef(R3BSfibMapped2CalPar, 1)
};

#endif
