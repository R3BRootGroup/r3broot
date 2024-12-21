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

#ifndef R3BCALIFADEMOONLINESPECTRA
#define R3BCALIFADEMOONLINESPECTRA
#define N_MAX_CRY 64
#define N_MAX_PETALS 8

#include "FairTask.h"
#include "TCanvas.h"
#include "TMath.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

class TClonesArray;
class TH1F;
class TH2F;
class R3BEventHeader;

/**
 * This taks reads CALIFA data and plots online histograms
 */
class R3BCalifaDemoOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BCalifaDemoOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BCalifaDemoOnlineSpectra(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BCalifaDemoOnlineSpectra();

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
     * Method for setting number of petals
     */
    inline void SetPetals(Int_t petals) { fCalifaNumPetals = petals; }

    /**
     * Method for setting the configuration parameters file
     */
    inline void SetCalifaConfigFile(TString file) { fCalifaFile = file; }

    /**
     * Method to select binning and max range
     */
    void SetRange_bins(Int_t Histos_bins) { fMapHistos_bins = Histos_bins; }
    void SetRange_max(Int_t Histos_max) { fMapHistos_max = Histos_max; }

    /**
     * Method to reset histograms
     */
    void Reset_CALIFA_DEMO_Histo();

  private:
    /**
     * Method to change histogram scales
     */
    void Log_CALIFA_Histo();

    /**
     * Method for setting the Data Level (Mapped or Cal)
     */
    void Map2Cal_CALIFA_Histo();

    /**
     * Method for setting histogram sequence (Febex or Preamp. channels)
     */
    void Febex2Preamp_CALIFA_Histo();

    /**
     * Mapping s444 crystals
     */
    Int_t Map_For_s444(Int_t val);

    Int_t fMapHistos_max;
    Int_t fMapHistos_bins;

    TClonesArray* fMappedItemsCalifa; /**< Array with mapped items.    */
    TClonesArray* fCalItemsCalifa;    /**< Array with cal items.       */
    TClonesArray* fHitItemsCalifa;    /**< Array with hit items.       */
    TClonesArray* fWRItemsCalifa;     /**< Array with WR-Califa items. */
    TClonesArray* fWRItemsMaster;     /**< Array with WR-Master items. */

    // Check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.  */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fNEvents;         /**< Event counter. */

    Int_t fCalifaNumPetals;      /**< Number of Petals.   */
    Int_t fNumCrystalPetal;      /**< Crystals per Petal. */
    Int_t fOrderFebexPreamp[16]; /**< Selector for febex or preamp sequence. */

    // Multiplicities
    TH1F* fh_Califa_Mult[N_MAX_PETALS + 1];
    TH1F* fh_Califa_MultHit[N_MAX_PETALS + 1];

    // WR data
    TH1F* fh_Califa_wr;

    // Raw data
    TH2F* fh_Califa_cryId_petal;
    TH1F* fh_Califa_energy_per_petal[N_MAX_PETALS + 1]; //+1 for proton range:s444
    TH1F* fh_Califa_crystals[N_MAX_PETALS][N_MAX_CRY];
    TH2F* fh_Califa_cryId_energy;
    TH2F* fh_Califa_coinc_petal1;
    TH2F* fh_Califa_coinc_petal2;
    TH2F* fh_Califa_coinc_petal3;

    // Cal data
    TH2F* fh_Califa_cryId_energy_cal;
    TH1F* fh_Califa_energy_per_petal_cal[N_MAX_PETALS + 1]; //+1 for proton range:s444
    TH1F* fh_Califa_crystals_cal[N_MAX_PETALS][N_MAX_CRY];

    // Hit data
    TH2F* fh_Califa_theta_phi;
    TH2F* fh_Califa_theta_energy[N_MAX_PETALS];
    TH1F* fh_Califa_total_energy;

    TString fCalifaFile;  /**< Config file name. */
    Bool_t fCalON;        /**< Cal selector. */
    Bool_t fLogScale;     /**< Selecting scale. */
    Bool_t fRaw2Cal;      /**< Mapped or Cal selector. */
    Bool_t fFebex2Preamp; /**< Febex or Preamp selector. */

    TCanvas* cMap;
    TCanvas* cCalifa_wr;
    TCanvas* cCalifa1;
    TCanvas* cCalifa2;
    TCanvas* cCalifa3;
    TCanvas* cCalifa3b; // proton range
    TCanvas* cCalifa4[N_MAX_PETALS][4];
    TCanvas* cCalifa5;
    TCanvas* cCalifa6;
    TCanvas* cCalifa7;
    TCanvas* cCalifa8;
    TCanvas* cCalifa_hitpetal[N_MAX_PETALS];
    TCanvas* cCalifa10;
    TCanvas* cCalifa_hitenergy;

  public:
    ClassDef(R3BCalifaDemoOnlineSpectra, 1)
};

#endif
