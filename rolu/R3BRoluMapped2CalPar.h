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

// ----------------------------------------------------------------
// -----     Create time calib parameters for Rolu via        -----
// ------         Created July 2019 by A. Kelic-Heil          -----
// -----            Following R3BLosMapped2CalPar             -----
// ----------------------------------------------------------------

#pragma once

#include <FairTask.h>

class R3BTCalPar;
class TClonesArray;
class R3BEventHeader;
class R3BTCalEngine;

/**
 * An analysis task for TCAL calibration of NeuLAND Tamex data.
 * This class fills TDC distribution for each Photomultiplier
 * of the NeuLAND detector and calculates the calibration
 * parameters using the R3BTCalEngine.
 * @author D. Kresan
 * @since September 7, 2015
 */
class R3BRoluMapped2CalPar : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BRoluMapped2CalPar();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BRoluMapped2CalPar(const char* name, int iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BRoluMapped2CalPar();

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
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    virtual void FinishTask();

    /**
     * Method for setting the update rate for control histograms
     * @param rate an update rate value (events).
     */
    inline void SetUpdateRate(int rate) { fUpdateRate = rate; }

    /**
     * Method for setting minimum required statistics per module.
     * Only detector modules with number of entries in TDC
     * distribution greater than minimum statistics will be
     * calibrated.
     * @param minStats a value of minimum statistics required.
     */
    inline void SetMinStats(int minStats) { fMinStats = minStats; }

    /**
     * Method for selecting events with certain trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(int trigger) { fTrigger = trigger; }

    /**
     * Method for setting number of Rolu detectors and channels.
     * @param nDets number of detectors.
     * @param nCh number of channels per detector (4+master trigger?)
     */
    inline void SetNofModules(int nDets) { fNofDetectors = nDets; }

    // Accessor to select online mode
    void SetSkipTrigger(bool skip) { fSkipTrigger = skip; }

  private:
    int fUpdateRate = 1000000; /**< An update rate. */
    int fMinStats = 100000;    /**< Minimum statistics required per module. */
    int fTrigger = -1;         /**< Trigger value. */

    unsigned int fNofDetectors = 1; /**< Number of detectors. */
    unsigned int fNofChannels = 4;  /**< Number of channels per detector. */
    unsigned int fNofTypes = 2;     /**< Number of time-types per channel ( TAMEX leading/trailing). */

    bool fSkipTrigger = false;

    int counter = 0;
    int fNEventsi = 0;              /**< Event counter. */
    R3BTCalPar* fCal_Par{};         /**< Parameter container. */
    TClonesArray* fMapped{};        /**< Array with mapped data - input data. */
    TClonesArray* fMappedTrigger{}; /**< Array with mapped trigger data - input data. */
    R3BEventHeader* header{};       /**< Event header - input data. */

    R3BTCalEngine* fEngine; /**< Instance of the TCAL engine. */

  public:
    ClassDef(R3BRoluMapped2CalPar, 1)
};
