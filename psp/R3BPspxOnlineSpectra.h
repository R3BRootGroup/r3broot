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
// -----                  R3BPspxOnlineSpectra            -----
// -----          Created Feb 23th 2019 by J.L. Rodriguez -----
// -----              Modified Dec 2019 by M. Holl        -----
// ------------------------------------------------------------

#ifndef R3BPspxOnlineSpectra_H
#define R3BPspxOnlineSpectra_H

#include "FairTask.h"
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>

#include "TClonesArray.h"
#include "TMath.h"
#include <cstdlib>
class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

#define PSPX 3

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BPspxOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BPspxOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BPspxOnlineSpectra(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BPspxOnlineSpectra();

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
     * Method for setting the trigger value.
     * @param trigger 1 - onspill, 2 - offspill, -1 - all events.
     */
    inline void SetTrigger(Int_t trigger) { fTrigger = trigger; }

    void Reset_PSPX_Histo();

  private:
    std::vector<TClonesArray*> fMappedItemsPspx; /**< Array with Mapped items for Pspx. */
    std::vector<TClonesArray*> fPrecalItemsPspx; /**< Array with Precal items for Pspx. */
    std::vector<TClonesArray*> fCalItemsPspx;    /**< Array with Cal items for Pspx. */
    std::vector<TClonesArray*> fHitItemsPspx;    /**< Array with Hit items for Pspx. */

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fNEvents;         /**< Event counter. */

    std::vector<TH1F*> fh_pspx_multiplicity;   /**< PSPX multiplicity on mapped level */
    std::vector<TH1F*> fh_pspx_strip_1;        /**< PSPX channel profile on mapped level */
    std::vector<TH1F*> fh_pspx_strip_2;        /**< PSPX channel profile on mapped level */
    std::vector<TH2F*> fh_pspx_energy_strip_1; /**< PSPX energy vs strip on mapped level */
    std::vector<TH2F*> fh_pspx_energy_strip_2; /**< PSPX energy vs strip on mapped level */

    std::vector<TH2F*> fh_pspx_cal_strip_frontback;  /**< PSPX strip front vs back on cal level */
    std::vector<TH2F*> fh_pspx_cal_pos_frontback;    /**< PSPX position front vs back on cal level */
    std::vector<TH2F*> fh_pspx_cal_energy_frontback; /**< PSPX energy front vs back on cal level */
    std::vector<TH1F*> fh_pspx_cal_energyfront;      /**< PSPX energy front vs back on cal level */
    std::vector<TH1F*> fh_pspx_cal_energyback;       /**< PSPX energy front vs back on cal level */

    std::vector<TH1F*> fh_pspx_cal_strip; /**< PSPX precal strip*/
  public:
    ClassDef(R3BPspxOnlineSpectra, 2)
};

#endif
