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
// -----               R3BTofDOnlineSpectra               -----
// -----    Created 18/03/22 by J.L. Rodriguez-Sanchez    -----
// -----           Fill Tofd online histograms            -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <memory>
#include <vector>

class TClonesArray;
class TH1F;
class TH2F;
class R3BCoarseTimeStitch;
class R3BEventHeader;
class R3BTofDMappingPar;

class R3BTofDOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTofDOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BTofDOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTofDOnlineSpectra();

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

    virtual void SetParContainers();

    /**
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

    /**
     * Method for setting the trigger value.
     * @param trigger 1 - physics, 2 - offspill, -1 - all events.
     */
    void SetTrigger(Int_t trigger) { fTrigger = trigger; }
    void SetTpat(Int_t tpat1, Int_t tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }

    /**
     * Methods for setting number of planes and paddles
     */
    void SetNofModules(Int_t planes, Int_t p)
    {
        fNofPlanes = planes;
        fPaddlesPerPlane = p;
    }

    void SetMaxMult(Int_t m) { fMaxmul = m; }

    void SetClockRange(Double_t r) { fC_range_ns = r; }

    void SetClockCoincidenceRange(Double_t r) { fC_bar_coincidence_ns = r; }

    void SetToFCorrection(std::vector<double> tofcor) { fTofcor = tofcor; }

  private:
    void SetParameter();
    std::unique_ptr<R3BCoarseTimeStitch> fTimeStitch;
    TClonesArray* fMappedItems = nullptr;
    TClonesArray* fCalItems = nullptr;
    TClonesArray* fHitItems = nullptr;
    TClonesArray* fCalTriggerItems = nullptr; /**< Array with trigger Cal items - input data. */

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header = nullptr; /**< Event header. */
    R3BTofDMappingPar* fMapPar = nullptr;
    Int_t fTrigger = -1; /**< Trigger value. */
    Int_t fTpat1 = -1, fTpat2 = -1;
    UInt_t fNofPlanes = 4;
    Int_t fMaxmul = 100;
    UInt_t fPaddlesPerPlane = 44;         /**< Number of paddles per plane. */
    Double_t fC_range_ns = 2048 * 5.;     // ns
    Double_t fC_bar_coincidence_ns = 20.; // ns
    unsigned long fNEvents = 0;           /**< Event counter. */
    std::vector<double> fTofcor;

    std::vector<TH1F*> fh_tofd_channels;
    std::vector<TH2F*> fh_tofd_multihit;
    std::vector<TH2F*> fh_tofd_TotPm;
    std::vector<TH2F*> fh_tofd_dt;
    std::vector<TH2F*> fh_tofd_multihit_coinc;
    std::vector<TH2F*> fh_tofd_TotPm_coinc;
    std::vector<TH2F*> fh_num_side;
    std::vector<TH1F*> fh_tofd_bars;
    std::vector<TH2F*> fh_tofd_Tot_hit;
    std::vector<TH2F*> fh_tofd_time_hit;
    std::vector<TH1F*> fh_tofd_multihit_hit;
    std::vector<TH2F*> fh_tofd_dt_hit;
    std::vector<TH2F*> fh2_tofd_ypos_cal;
    std::vector<TH2F*> fh2_tofd_timedif_cal;

    std::vector<TH2F*> fh_tofd_time_los_h2;
    std::vector<TH2F*> fh2_tofd_time_los_cal;
    std::vector<TH1F*> fh_tofd_time_los[44];

    TH2F* fh2_tofd_time_vs_charge;
    TH2F* fh2_tofd_time_wouttrig_vs_charge;
    TH2F* fh2_tofd_charges13;
    TH2F* fh2_tofd_charges24;
    TH2F* fh2_tofd_charges12;
    TH2F* fh2_tofd_charges23;
    TH2F* fh2_tofd_charges34;
    TH1F* fh1_Zcharge;
    TH2F* fh2_zcharge_tpat;

  public:
    ClassDef(R3BTofDOnlineSpectra, 1)
};
