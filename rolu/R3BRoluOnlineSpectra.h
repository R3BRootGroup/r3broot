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
// -----                  R3BRoluOnlineSpectra             -----
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <array>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <TClonesArray.h>
#include <TFolder.h>
#include <TMath.h>
#include <cstdlib>

#define VFTX_CLOCK_MHZ 200

class TClonesArray;
class R3BCoarseTimeStitch;
class TH1F;
class TH2F;
class R3BEventHeader;

/**
 * This taks reads ROLU data items and plots histograms
 * for online checks.
 */
class R3BRoluOnlineSpectra : public FairTask
{
  public:
    // Default constructor.
    R3BRoluOnlineSpectra();

    // Standard constructor.
    R3BRoluOnlineSpectra(const char* name, int iVerbose = 1);

    // Destructor.
    virtual ~R3BRoluOnlineSpectra();

    // Method for task initialization.
    virtual InitStatus Init();

    // Method for event loop implementation.
    virtual void Exec(Option_t* option);

    // A method for finish of processing of an event.
    virtual void FinishEvent();

    // Method for finish of the task execution.
    virtual void FinishTask();

    // Method for setting the trigger value.
    //  @param trigger 1 - physics, 2 - offspill, -1 - all events.
    inline void SetTrigger(int trigger) { fTrigger = trigger; }
    inline void SetTpat(int tpat) { fTpat = tpat; }

    void Reset_ROLU_Histo();

    /* Method for setting number of ROLU detectors */
    void SetNofRoluModules(int nDets) { fNofRoluDetectors = nDets; }

  private:
    TClonesArray* fMappedItems;
    TClonesArray* fCalItems;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header{}; /**< Event header. */
    int fTrigger = -1;        /**< Trigger value. */
    int fTpat = -1;
    double fClockFreq = 1. / VFTX_CLOCK_MHZ * 1000.; /**< Clock cycle in [ns]. */
    int nRoluEvents = 0;
    int fNofRoluDetectors = 1;         /**< Number of ROLU detectors. */
    static constexpr int fNofRolu = 1; /**< Number of ROLU detectors. */

    int counter[4] = { 0 };

    const char* fDetectorNames[fNofRolu] = { "Rolu" };

    TH1F* fh_rolu_channels[fNofRolu]{};
    TH2F* fh_rolu_tot[fNofRolu]{};
    TH1F* fh_rolu_tot_1D[fNofRolu][4]{};
    TH1F* fh_rolu_LE_raw[fNofRolu][4]{};

  public:
    ClassDef(R3BRoluOnlineSpectra, 0)
};
