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

#ifndef R3BOnlineSpectraSci2_H
#define R3BOnlineSpectraSci2_H

#include "FairTask.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads SCI data and plots online histograms
 */
class R3BOnlineSpectraSci2 : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectraSci2();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectraSci2(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectraSci2();

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
    virtual void Reset();

    void SetNbDetectors(Int_t ndets) { fNbDetectors = ndets; }
    void SetNbChannels(Int_t nchs) { fNbChannels = nchs; }
    // Int_t GetNbDetectors() { return fNbDetectors; }
    // Int_t GetNbChannels() { return fNbChannels; }

    void SetTpat(Int_t tpat1, Int_t tpat2)
    {
        fTpat1 = tpat1;
        fTpat2 = tpat2;
    }

  private:
    TClonesArray* fMapped; /**< Array with R3BSofSciMappedData */
    TClonesArray* fTcal;   /**< Array with R3BSofSciTcalData */

    Int_t fNbDetectors; // fNbDetectors is also equal to fIdCaveC
    Int_t fNbChannels;
    Int_t fTpat1, fTpat2;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* fEventHeader; /**< Event header.      */
    Int_t fNEvents;               /**< Event counter.     */

    // Canvas
    TCanvas** cMapped; // [fNbDetectors]
    TCanvas** cPos;    // [fNbDetectors]

    TCanvas* cMultMap;
    TCanvas* cMultTcal;
    TCanvas* cMultMap2D;
    TCanvas* cMultMap2D_RvsL;

    // Histograms - 1D multiplicity
    TH1I** fh1_multMap;  // [fNbDetectors * fNbChannels];
    TH1I** fh1_multTcal; // [fNbDetectors * fNbChannels];
    // Histograms - 1D multiplicity with condition on Tpat
    TH1I** fh1_multMap_condTpat;  // [fNbDetectors * fNbChannels];
    TH1I** fh1_multTcal_condTpat; // [fNbDetectors * fNbChannels];

    // Histograms - 2D multiplicity
    TH2I** fh2_mult_RvsL;      // [fNbDetectors];
    TH2I** fh2_mult_TrefVsPmt; //[fNbDetectors * (NbChannels-1)]

    // Histograms - 2D multiplicity with condition on Tpat
    TH2I** fh2_mult_RvsL_condTpat;      // [fNbDetectors];
    TH2I** fh2_mult_TrefVsPmt_condTpat; //[fNbDetectors * (NbChannels-1)]

    // Histograms for Mapped data : Fine Time and Mult
    TH1I** fh1_finetime; // [fNbDetectors * NbChannels];
    TH2I** fh2_mult;     // [fNbDetectors];
    TH1D** fh1_DeltaTref;

    // Histograms for X position at Tcal
    TH1F** fh1_RawPos_TcalMult1; // [fNbDetectors];

  public:
    ClassDef(R3BOnlineSpectraSci2, 1)
};

#endif
