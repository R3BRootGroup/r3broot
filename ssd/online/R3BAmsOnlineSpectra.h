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
// -----                  R3BAmsOnlineSpectra             -----
// -----    Created 22/07/18  by J.L. Rodriguez-Sanchez   -----
// -----           Fill AMS online histograms             -----
// ------------------------------------------------------------

#ifndef R3BAmsOnlineSpectra_H
#define R3BAmsOnlineSpectra_H

#include "FairTask.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TMath.h"
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

class TClonesArray;
class R3BEventHeader;
class R3BAmsMappingPar;

/**
 * This taks reads AMS data and plots online histograms
 */
class R3BAmsOnlineSpectra : public FairTask
{

  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BAmsOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BAmsOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BAmsOnlineSpectra();

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

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /**
     * Method to reset histograms
     */
    void Reset_AMS_Histo();

  private:
    void SetParameter();

    TClonesArray* fMappedItemsAms; /**< Array with mapped items. */
    TClonesArray* fCalItemsAms;    /**< Array with cal items. */
    TClonesArray* fHitItemsAms;    /**< Array with hit items. */

    R3BAmsMappingPar* fMap_Par; /**< Container with mapping parameters. >*/

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header. */
    Int_t fTrigger;         /**< Trigger value. */
    Int_t fNEvents;         /**< Event counter. */

    TCanvas* cMap;
    TCanvas *cCalL, *cCalR;
    TCanvas* cHit[6];
    TCanvas* cHitAngles;
    TCanvas* cHitEnergyCor;
    TCanvas* cHitThetaCor;
    TCanvas* cHitPhiCor;

    TH2F* fh_Ams_energy_allStrips[6];
    TH2F* fh_Ams_energy_allCalStrips[12];
    TH2F* fh_Ams_hit_Pos[6];
    TH2F* fh_Ams_hit_E[6];
    TH2F* fh_Ams_hit_E_theta[6];
    TH1F* fh_Ams_hit_Mul[6][2];
    TH2F* fh2_ams_theta_phi[2];
    TH2F* fh2_ams_e1_e2[2];
    TH2F* fh2_ams_theta1_theta2[4];
    TH2F* fh2_ams_phi1_phi2[2];

    // TString fAmsFile;        	      /**< Config file name. */
    Int_t fNbDet; /**< Number of AMS detectors. */

  public:
    ClassDef(R3BAmsOnlineSpectra, 1)
};

#endif
