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
// -----          R3BIncomingTrackingOnlineSpectra        -----
// -----    Created 25/02/22 by J.L. Rodriguez-Sanchez    -----
// -----         Fill tracking online histograms          -----
// ------------------------------------------------------------

#pragma once

#include <FairTask.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TH1.h>
#include <TH2F.h>
#include <memory>

class TClonesArray;
class R3BEventHeader;
class R3BTGeoPar;

class R3BIncomingTrackingOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BIncomingTrackingOnlineSpectra();

    R3BIncomingTrackingOnlineSpectra(const TString& name, float, float, float, float);

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BIncomingTrackingOnlineSpectra(const TString& name, int iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BIncomingTrackingOnlineSpectra() = default;

    /**
     * Method for task initialization.
     * This function is called by the framework before
     * the event loop.
     * @return Initialization status. kSUCCESS, kERROR or kFATAL.
     */
    InitStatus Init() override;

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    void SetParContainers() override;

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
     * Methods to set the year configuration.
     */
    inline void SetYearConfiguration(uint16_t year = 2024) { fYearConf = year; }

    /**
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

  private:
    R3BEventHeader* fHeader = nullptr;
    TClonesArray* fMwpc0HitDataCA = nullptr;
    TClonesArray* fMwpc1HitDataCA = nullptr;
    TClonesArray* fMwpc2HitDataCA = nullptr;
    TClonesArray* fFrsHitDataCA = nullptr;

    // Parameters
    R3BTGeoPar* fTargetGeoPar = nullptr;
    R3BTGeoPar* fMw0GeoPar = nullptr;
    R3BTGeoPar* fMw1GeoPar = nullptr;
    R3BTGeoPar* fMw2GeoPar = nullptr;

    int fNEvents = 0;
    uint16_t fYearConf = 2024;
    float fPosTarget = 2576.0;
    float fWidthTarget = 40.0;
    float fDist_acelerator_glad = 3855.0;
    TString fNameCut = "";

    // Canvas
    TCanvas* cTrackingXZ;
    TCanvas* cTrackingYZ;
    TCanvas* cBeamProfileTarget;
    TCanvas* cZvsBeta;

    // Histograms for Hit data
    TH2F* fh2_tracking_planeXZ;
    TH2F* fh2_tracking_planeYZ;
    TH2F* fh2_target_PosXY;
    TH2F* fh2_ZvsBeta;
    TH2F* fh2_angvsposx;
    TH2F* fh2_angvsposy;
    std::unique_ptr<TCutG> fCutIncoming = nullptr;

  public:
    ClassDefOverride(R3BIncomingTrackingOnlineSpectra, 1)
};
