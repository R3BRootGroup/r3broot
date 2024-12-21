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
// -----            R3BTofDvsTttxOnlineSpectra            -----
// -----    Created 10/02/24 by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <TCanvas.h>
#include <memory>
#include <vector>

class TH2F;
class TClonesArray;
class R3BEventHeader;

constexpr const int8_t fNbHist_tofd_tttx = 3;

class R3BTofDvsTttxOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BTofDvsTttxOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    explicit R3BTofDvsTttxOnlineSpectra(const TString& name, int iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BTofDvsTttxOnlineSpectra() = default;

    /** Virtual method SetParContainers **/
    void SetParContainers() override;

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    InitStatus Init() override;

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    /**
     * Method for event loop implementation.
     * Is called by the framework every time a new event is read.
     * @param option an execution option.
     */
    void Exec(Option_t* /*option*/) override;

    /**
     * A method for finish of processing of an event.
     * Is called by the framework for each event after executing
     * the tasks.
     */
    void FinishEvent() override;

    /**
     * Method for finish of the task execution.
     * Is called by the framework after processing the event loop.
     */
    void FinishTask() override;

    /**
     * Method for setting the fTpat
     */
    inline void SetTpat(int tpat) { fTpat = tpat; }

    /**
     * Method for setting the charge selection for tofd
     */
    inline void SetZGateTofD(double minz, double maxz)
    {
        fMinZTofd = minz;
        fMaxZTofd = maxz;
    }

    /**
     * Method for setting the charge selection for tofd
     */
    inline void SetPlaneTofD(int plane) { fPlaneTofd = plane; }

    /**
     * Method to reset histograms
     */
    void Reset_Histo();

  private:
    TClonesArray* fHitItemsTofd = nullptr;
    TClonesArray* fHitItemsTttx = nullptr;
    R3BEventHeader* fHeader = nullptr;

    int fNEvents = 0;
    int fTpat = 0;
    double fMinZTofd = 0.;   // Atomic number
    double fMaxZTofd = 100.; // Atomic number
    int fPlaneTofd = 1;

    std::vector<TH2F*> fh2_charges;
    std::unique_ptr<TCanvas> cCTofDTx;

  public:
    ClassDefOverride(R3BTofDvsTttxOnlineSpectra, 1)
};
