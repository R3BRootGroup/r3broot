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

#pragma once

#include "FairTask.h"
#include "TCanvas.h"
#include "TH2.h"

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads SCI data and plots online histograms
 */
class R3BOnlineSpectraFrsSciVsMusli : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectraFrsSciVsMusli();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectraFrsSciVsMusli(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectraFrsSciVsMusli();

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
     * Methods to clean histograms.
     */
    virtual void Reset_Histo();

    /** Virtual method Reset **/
    virtual void Reset() {}

    void SetIdS2(UShort_t id) { fIdS2 = id; }
    void SetIdCaveC(UShort_t id) { fIdCaveC = id; }

  private:
    R3BEventHeader* fEventHeader; /**< Event header.      */

    TClonesArray* fFrsSci_PosCal; /**< Array with R3BFrsSciPosCalData */
    TClonesArray* fFrsSci_TofCal; /**< Array with R3BFrsSciTofCalData */
    TClonesArray* fMusli_Map;     /**< Array with R3BMusliMappedData */
    TClonesArray* fMusli_Cal;     /**< Array with R3BMusliCalData */
    TClonesArray* fMusli_Hit;     /**< Array with R3BMusliHitData */
    Int_t fNEvents;               /**< Event counter.     */
    UShort_t fIdS2;
    UShort_t fIdCaveC;

    // Canvas FrsSciCal vs MusliMap
    TCanvas* cMap_Pos;

    // Canvas FrsSciCal vs MusliCal
    TCanvas* cCal_Pos;
    TCanvas* cCal_EvsAoQ;

    // Histograms FrsSciCal vs MusliCal
    TH2D** fh2_Map_PosRawVsDT; // at Cave C only, one histo per pair of anodes [8]
    TH2D** fh2_Cal_PosCalVsDT; // at Cave C only, one histo per pair of anodes [8]
    TH2D** fh2_Cal_EvsAoQ;     // AoQ from S2 to cave C, one histo per signal [15]

    // Canvas FrsSciCal vs MusliHit
    TCanvas* cHit_EvsTofRaw;
    TCanvas* cHit_EvsBeta;
    TCanvas* cHit_EvsAoQ;
    TCanvas* cHit_ZvsAoQ;

    //// Histograms FrsSciCal vs MusliHit
    TH2D** fh2_Hit_EvsTofRaw; // TofRaw calculated from S2 to Cave C, one histo per type of data
    TH2D** fh2_Hit_EvsBeta;   // Beta calculated from S2 to Cave C, one histo per type of data
    TH2D** fh2_Hit_EvsAoQ;    // AoQ calculated from S2 to Cave C, one histo per type of data
    TH2D** fh2_Hit_ZvsAoQ;    // AoQ calculated from S2 to Cave C, one histo per type of data

  public:
    ClassDef(R3BOnlineSpectraFrsSciVsMusli, 1)
};
