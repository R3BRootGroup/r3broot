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
// -----              R3BGeneralOnlineSpectra             -----
// -----    Created 05/03/22 by J.L. Rodriguez-Sanchez    -----
// -----           Fill R3B online histograms             -----
// ------------------------------------------------------------

#ifndef R3BGeneralOnlineSpectra_H
#define R3BGeneralOnlineSpectra_H 1

#include "FairTask.h"
#include "TCanvas.h"
#include "TMath.h"

class TH1F;
class TH2F;
class TClonesArray;
class R3BEventHeader;
class R3BSofAtOnlineSpectra;
class R3BMwpcOnlineSpectra;
class R3BMwpcCorrelationOnlineSpectra;
class R3BSofTrimOnlineSpectra;
class R3BTwimOnlineSpectra;
class R3BSofTwimvsMusicOnlineSpectra;
class R3BSofTwimvsTrimOnlineSpectra;
class R3BSofSciOnlineSpectra;
class R3BSofSciVsMusicOnlineSpectra;
class R3BSofSciVsTrimOnlineSpectra;
class R3BSofSciVsMwpc0OnlineSpectra;
class R3BSofTofWOnlineSpectra;
class R3BSofTrimVsTofwOnlineSpectra;
class R3BSofScalersOnlineSpectra;
class R3BMusicOnlineSpectra;
class R3BAmsOnlineSpectra;
class R3BLosOnlineSpectra;
class R3BFootOnlineSpectra;
class R3BCalifaOnlineSpectra;
class R3BIncomingTrackingOnlineSpectra;
class R3BTwimvsFootOnlineSpectra;
class R3BTofDOnlineSpectra;

class R3BGeneralOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BGeneralOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BGeneralOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BGeneralOnlineSpectra();

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
    virtual void Reset_GENERAL_Histo();

  private:
    TClonesArray* fWRItemsMaster;  /**< Array with WR-Master items. */
    TClonesArray* fWRItemsSofia;   /**< Array with WR-Sofia items. */
    TClonesArray* fWRItemsCalifa;  /**< Array with WR-Califa items. */
    TClonesArray* fWRItemsNeuland; /**< Array with WR-Neuland items. */
    TClonesArray* fWRItemsS2;      /**< Array with WR-S2 items. */
    TClonesArray* fWRItemsS8;      /**< Array with WR-S8 items. */

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* fEventHeader; /**< Event header.      */
    Int_t fTrigger;               /**< Trigger value.     */
    Int_t fNEvents;               /**< Event counter.     */

    R3BMwpcOnlineSpectra* fMwpc0Online;
    R3BMwpcCorrelationOnlineSpectra* fMwpc01Online;
    R3BMwpcCorrelationOnlineSpectra* fMwpc02Online;
    R3BMwpcCorrelationOnlineSpectra* fMwpc12Online;
    R3BMwpcOnlineSpectra* fMwpc1Online;
    R3BMwpcOnlineSpectra* fMwpc2Online;
    R3BTwimOnlineSpectra* fTwimOnline;
    R3BMusicOnlineSpectra* fMusicOnline;
    R3BAmsOnlineSpectra* fAmsOnline;
    R3BCalifaOnlineSpectra* fCalifaOnline;
    R3BFootOnlineSpectra* fFootOnline;
    R3BLosOnlineSpectra* fLosOnline;
    R3BIncomingTrackingOnlineSpectra* fIncomingTrackingOnline;
    R3BTwimvsFootOnlineSpectra* fTwimvsFootOnlineSpectra;
    R3BTofDOnlineSpectra* fTofdOnlineSpectra;

    // Canvas
    TCanvas *cTpat, *cTrigger, *cWr, *cWrs;

    // Unpack
    TH1F *fh1_trigger, *fh1_tpat, *fh1_wr[2];
    TH1F* fh1_wrs[5];
    TH2F* fh2_tpat;
    std::vector<int> tpatbin;

  public:
    ClassDef(R3BGeneralOnlineSpectra, 0)
};

#endif
