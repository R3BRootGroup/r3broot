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
// -----                  R3BGlobalAnalysis                -----
// -----            Created 13-04-2016 by M.Heil          -----
// -----               Fill online histograms             -----
// ------------------------------------------------------------

#ifndef R3BGLOBALANALYSIS
#define R3BGLOBALANALYSIS

#define N_PLANE_MAX_TOFD 4
#define N_PADDLE_MAX_TOFD 50
#define N_PADDLE_MAX_PTOF 100
#define N_FIBER_PLOT_2019 520 // range to plot
#define N_FIBER_PLOT 1050     // range to plot

#define N_PSPX 2
#define N_STRIPS_PSPX 32
#define PSPX_SPLITTER 0 // sstorck: set this variable to 1 if signals are split into energy and position readout

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

/**
 * This taks reads all detector data items and plots histograms
 * for online checks.
 */
class R3BGlobalAnalysis : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BGlobalAnalysis();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BGlobalAnalysis(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BGlobalAnalysis();

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

    /**
     * Methods for setting number of planes and paddles
     */
    inline void SetNofModules(Int_t planes, Int_t ppp)
    {
        fNofPlanes = planes;
        fPaddlesPerPlane = ppp;
    }

    //   virtual void SetParContainers();

  private:
    std::vector<TClonesArray*> fMappedItems;
    std::vector<TClonesArray*> fCalItems;
    std::vector<TClonesArray*> fHitItems;

    enum DetectorInstances
    {
        DET_AMS,
        DET_CALIFA,
        DET_FI_FIRST,
        DET_FI1A = DET_FI_FIRST,
        DET_FI1B,
        DET_FI2A,
        DET_FI2B,
        DET_FI3A,
        DET_FI3B,
        DET_FI4,
        DET_FI5,
        DET_FI6,
        DET_FI7,
        DET_FI8,
        DET_FI9,
        DET_FI10,
        DET_FI11,
        DET_FI12,
        DET_FI13,
        DET_FI_LAST = DET_FI13,
        DET_L3T,
        DET_LOS,
        DET_MUSIC,
        DET_NEULAND,
        DET_PSPX,
        DET_PTOF,
        DET_ROLU,
        DET_SCI8,
        DET_STRAW,
        DET_TOFD,
        DET_MAX
    };

#define NOF_FIB_DET (DET_FI_LAST - DET_FI_FIRST + 1)

    const char* fDetectorNames[DET_MAX + 1] = { "Ams",  "Califa", "Fi1a",  "Fi1b", "Fi2a",  "Fi2b",    "Fi3a", "Fi3b",
                                                "Fi4",  "Fi5",    "Fi6",   "Fi7",  "Fi8",   "Fi9",     "Fi10", "Fi11",
                                                "Fi12", "Fi13",   "L3t",   "Los",  "Music", "Neuland", "Pspx", "Ptof",
                                                "Rolu", "Sci8",   "Straw", "Tofd", NULL };

    // If FiberI is present or not:
    Int_t ifibdet;
    // Number of fibers per detector
    Double_t n_fiber[NOF_FIB_DET] = { 256.,  256.,  256., 256., 512., 512.,  2048.,
                                      2048., 1024., 512., 512., 512., 1024., 1024. };

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Double_t fClockFreq;    /**< Clock cycle in [ns]. */
    UInt_t fNofPlanes;
    UInt_t fPaddlesPerPlane; /**< Number of paddles per plane. */
    unsigned long long t0_prev = 0;
    Double_t time_previous_event = 0;

    Int_t fNEvents = 0;   /**< Event counter. */
    Int_t fFibEvents = 0; /**< Event counter. */
    Double_t max_values[NOF_FIB_DET][2048];
    Int_t FibMax[NOF_FIB_DET];

    TH2F* fhChargeLosTofD;
    TH2F* fh_los_pos;

    TH1F* fh_channels_Fib[NOF_FIB_DET];
    TH1F* fh_fibers_Fib[NOF_FIB_DET];
    TH1F* fh_fiber_Fib[NOF_FIB_DET];
    TH1F* fh_mult_Fib[NOF_FIB_DET];
    TH2F* fh_Fib_ToF[NOF_FIB_DET];
    TH1F* fh_xpos_Fib[NOF_FIB_DET];
    TH1F* fh_ypos_Fib[NOF_FIB_DET];

    TH2F* fh_time_Fib[NOF_FIB_DET];
    TH2F* fh_multihit_m_Fib[NOF_FIB_DET];
    TH2F* fh_multihit_s_Fib[NOF_FIB_DET];
    TH2F* fh_ToT_m_Fib[NOF_FIB_DET];
    TH2F* fh_ToT_s_Fib[NOF_FIB_DET];
    TH2F* fh_Fib_vs_Events[NOF_FIB_DET];
    TH2F* fh_Fibs_vs_Events[NOF_FIB_DET];
    TH2F* fh_Fibs_vs_Tofd[NOF_FIB_DET];
    TH2F* fh_Fibs_vs_Tofd_y[NOF_FIB_DET];

    TH2F* fh_Fib_vs_Fib[NOF_FIB_DET][NOF_FIB_DET];
    TH2F* fh_Fib_dx[NOF_FIB_DET][NOF_FIB_DET];

    TH2F* fh_Cave_position;

    TH2F* fh_tofd_pos;
    TH1F* fh_tofd_charge;
    TH1F* fh_TimePreviousEvent;
    TH1F* fh_tofd_mult;

    TH1F* fh_ptof_channels;
    TH1F* fh_ptof_channels_cut;
    TH1F* fh_ptof_test1;
    TH1F* fh_ptof_test2;
    TH1F* fh_ptof_TotPm1[N_PADDLE_MAX_PTOF];
    TH1F* fh_ptof_TotPm2[N_PADDLE_MAX_PTOF];

    TH1F* fh_pspx_strips_psp[N_PSPX];
    TH1F* fh_pspx_energy_psp[N_PSPX];
    TH1F* fh_pspx_multiplicity_psp[N_PSPX];

    TH2F* fh_pspx_pos1_strips;
    TH2F* fh_pspx_pos2_strips;
    TH2F* fh_pspx_pos1_energy;
    TH2F* fh_pspx_pos2_energy;

    TH2F* fh_pspx_cor_x_strips;
    TH2F* fh_pspx_cor_y_strips;
    TH2F* fh_pspx_cor_x_energy;
    TH2F* fh_pspx_cor_y_energy;

    TH2F* fh_ToF_vs_events[NOF_FIB_DET];

  public:
    ClassDef(R3BGlobalAnalysis, 1)
};

#endif
