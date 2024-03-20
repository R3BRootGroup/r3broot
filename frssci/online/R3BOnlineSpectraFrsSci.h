/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
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
#include "TH1.h"
#include "TH2.h"

class TClonesArray;
class R3BEventHeader;

/**
 * This taks reads SCI data and plots online histograms
 */
class R3BOnlineSpectraFrsSci : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BOnlineSpectraFrsSci();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BOnlineSpectraFrsSci(const char* name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BOnlineSpectraFrsSci();

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

    void SetNbDets(UShort_t ndets) { fNbDets = ndets; }
    void SetNbPmts(UShort_t npmts) { fNbPmts = npmts; }
    void SetNbTofs(UShort_t ndets)
    {
        switch (ndets)
        {
            case 1:
                fNbTofs = 0;
                break;
            case 2:
                fNbTofs = 1;
                break;
            default:
                UShort_t fact1 = 1;
                UShort_t fact2 = 1;
                for (UShort_t i = 2; i <= ndets; i++)
                    fact1 *= i;
                for (UShort_t i = 2; i <= (ndets - 2); i++)
                    fact2 *= i;
                fNbTofs = fact1 / (2 * fact2);
        }
    }

    void SetPosRange(float min, float max)
    {
        fpos_range_min = min;
        fpos_range_max = max;
    }
    void SetTofRange(float min, float max)
    {
        ftof_range_min = min;
        ftof_range_max = max;
    }

  private:
    R3BEventHeader* fEventHeader; /**< Event header.      */

    TClonesArray* fMapped; /**< Array with R3BFrsSciMappedData */
    TClonesArray* fTcal;   /**< Array with R3BFrsSciTcalData */
    TClonesArray* fPosCal; /**< Array with R3BFrsSciPosCalData */
    TClonesArray* fTofCal; /**< Array with R3BFrsSciTofCalData */
    Int_t fNEvents;        /**< Event counter.     */
    UShort_t fNbDets;
    UShort_t fNbPmts;
    UShort_t fNbTofs;

    // Canvas at Mapped level
    TCanvas* cMap_FT;
    TCanvas* cMap_Mult1D;
    TCanvas* cMap_Mult2D;

    // Histograms at Mapped level
    TH1I** fh1_Map_finetime; // [fNbDets * NbPmts]
    TH1I** fh1_Map_mult;     // [fNbDets * fNbPmts]
    TH2I** fh2_Map_mult;     // [fNbDets]
    TH2I** fh2_Map_multRvsL; // [fNbDets! / (2! * (fNbDets-2)!]

    // Canvas at Tcal level
    TCanvas* cTcal_Pos;
    TCanvas* cTcal_Tof;

    // Histograms at Tcal level
    TH1D** fh1_Tcal1Hit_PosRaw; // [fNbDets]
    TH1D** fh1_Tcal1Hit_TofRaw; // [fNbTofs]

    // Canvas at Cal level
    TCanvas* cCal_PosRaw;
    TCanvas* cCal_PosCal;
    TCanvas* cCal_TofRaw;
    TCanvas* cCal_TofCal;
    TCanvas* cCal_Beta;
    TCanvas* cCal_AoQ;
    TCanvas* cCal_BRhoVsPos;
    TCanvas* cCal_AoQVsPos;

    // Histograms at Tcal level
    TH1D** fh1_Cal_PosRaw;       // [fNbDets]
    TH1D** fh1_Cal_PosCal;       // [fNbDets]
    TH1D** fh1_Cal_TofRaw;       // [fNbTofs]
    TH1D** fh1_Cal_TofCal;       // [fNbTofs]
    TH1D** fh1_Cal_Beta;         // [fNbTofs]
    TH1D** fh1_Cal_BRho;         // [fNbTofs]
    TH1D** fh1_Cal_AoQ;          // [fNbTofs]
    TH2D** fh2_Cal_BRhoVsPosSta; // [fNbTofs]
    TH2D** fh2_Cal_BRhoVsPosSto; // [fNbTofs]
    TH2D** fh2_Cal_AoQVsPosSta;  // [fNbTofs]
    TH2D** fh2_Cal_AoQVsPosSto;  // [fNbTofs]

    // Define ranges
    float fpos_range_min;
    float fpos_range_max;
    float ftof_range_min;
    float ftof_range_max;

  public:
    ClassDef(R3BOnlineSpectraFrsSci, 1)
};
