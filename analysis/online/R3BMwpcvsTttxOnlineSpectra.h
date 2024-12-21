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
// -----            R3BMwpcvsTttxOnlineSpectra            -----
// -----    Created 26/01/24 by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <TCanvas.h>
#include <memory>
#include <vector>

class TH2F;
class TClonesArray;
class R3BEventHeader;

constexpr const int8_t fNbHist = 3;

class R3BMwpcvsTttxOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BMwpcvsTttxOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    explicit R3BMwpcvsTttxOnlineSpectra(const TString& name, int iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BMwpcvsTttxOnlineSpectra() = default;

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
     * Method for setting min difference in X (in mm) for MW-Tttx correlations
     */
    inline void SetMinDist(double min) { fMinDist = min; }

    /**
     * Method to reset histograms
     */
    void Reset_Histo();

  private:
    TClonesArray* fHitItemsMw = nullptr;
    TClonesArray* fHitItemsTttx = nullptr;
    R3BEventHeader* fHeader = nullptr;

    int fNEvents = 0;
    int fTpat = 0;
    double fMinDist = 10.; // mm

    std::vector<TH2F*> fh2_xcor;
    std::vector<TH2F*> fh2_multcor;
    std::unique_ptr<TCanvas> cCMwTx;
    std::unique_ptr<TCanvas> cCMwTx_mul;

  public:
    ClassDefOverride(R3BMwpcvsTttxOnlineSpectra, 1)
};
