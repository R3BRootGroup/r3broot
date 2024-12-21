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
// -----            R3BCalifavsTofDOnlineSpectra          -----
// -----    Created 21/05/22 by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <TCanvas.h>
#include <memory>
#include <vector>

class TClonesArray;
class TH2F;
class R3BEventHeader;

class R3BCalifavsTofDOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BCalifavsTofDOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BCalifavsTofDOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BCalifavsTofDOnlineSpectra() = default;

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    virtual void Exec(Option_t* /*option*/);

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
     * Method for setting min proton energy (in keV) for opening angle histogram
     */
    inline void SetMinProtonEnergyForOpening(Float_t min) { fMinProtonE = min; }

    /**
     * Method to reset histograms
     */
    void Reset_Histo();

    /**
     * Method for setting the fTpat
     */
    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    /**
     * Method for setting the charge of the nuclear residue
     */
    inline void SetZCharge(Float_t z) { fZselection = z; }

  private:
    TClonesArray* fHitItemsCalifa = nullptr;
    TClonesArray* fHitItemsTofd = nullptr;

    R3BEventHeader* header = nullptr;
    Int_t fNEvents = 0;
    Int_t fTpat = -1;
    Float_t fZselection = 5.;
    Float_t fMinProtonE = 50000.; /* Min proton energy (in keV) to calculate the opening angle */

    std::vector<TH2F*> fh2_Califa_theta_phi; // 0: all, 1: with TofD
    TH2F* fh2_Califa_coinPhi;
    TH2F* fh2_Califa_coinTheta;
    TCanvas* cCalifa_angles;

  public:
    ClassDef(R3BCalifavsTofDOnlineSpectra, 1)
};
