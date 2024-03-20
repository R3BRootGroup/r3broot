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

// ------------------------------------------------------------
// -----           R3BIncomingIDOnlineSpectra             -----
// -----    Created 05/05/22 by J.L. Rodriguez-Sanchez    -----
// -----           Fill PID online histograms             -----
// ------------------------------------------------------------

#ifndef R3BIncomingIDOnlineSpectra_H
#define R3BIncomingIDOnlineSpectra_H 1

#include "FairTask.h"

class TClonesArray;
class R3BEventHeader;
class R3BTGeoPar;
class TCanvas;
class TH1F;
class TH2F;

class R3BIncomingIDOnlineSpectra : public FairTask
{
  public:
    /**
     * Default constructor.
     * Creates an instance of the task with default parameters.
     */
    R3BIncomingIDOnlineSpectra();

    /**
     * Standard constructor.
     * Creates an instance of the task.
     * @param name a name of the task.
     * @param iVerbose a verbosity level.
     */
    R3BIncomingIDOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /**
     * Destructor.
     * Frees the memory used by the object.
     */
    virtual ~R3BIncomingIDOnlineSpectra();

    virtual void SetParContainers();

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

    /**
     * Method for setting the fTpat
     */
    void SetTpat(Int_t tpat) { fTpat = tpat; }

    // Setting parameters
    void SetStartPlaId(Int_t id) { fStaId = id; };
    void SetRangeRawTof(Double_t min, Double_t max)
    {
        fMin_RawTof = min;
        fMax_RawTof = max;
    }

    void SetRangeBeta(Double_t min, Double_t max)
    {
        fMin_Beta = min;
        fMax_Beta = max;
    }

    void SetRangeZ(Double_t min, Double_t max)
    {
        fMin_Z = min;
        fMax_Z = max;
    }

    void SetRangeAq(Double_t min, Double_t max)
    {
        fMin_Aq = min;
        fMax_Aq = max;
    }

    void SetRangeBrho(Double_t min, Double_t max)
    {
        fMin_Brho = min;
        fMax_Brho = max;
    }

    void SetIsoGate(Double_t minz, Double_t maxz, Double_t minaq, Double_t maxaq)
    {
        fMin_Z_gate = minz;
        fMax_Z_gate = maxz;
        fMin_Aq_gate = minaq;
        fMax_Aq_gate = maxaq;
    }

  private:
    TClonesArray* fHitFrs; /**< Array with hit items. */
    TClonesArray* fHitLos; /**< Array with hit items. */
    TClonesArray* fMwpc0HitDataCA;
    TClonesArray* fMwpc1HitDataCA;

    // Parameters
    R3BTGeoPar* fMw0GeoPar;
    R3BTGeoPar* fMw1GeoPar;

    // check for trigger should be done globablly (somewhere else)
    R3BEventHeader* header; /**< Event header.      */
    Int_t fNEvents;         /**< Event counter.     */
    Int_t fTpat;

    // Start Plastic ID
    Int_t fStaId;

    // Ranges for the histograms;
    Double_t fMin_Z, fMax_Z, fMin_Aq, fMax_Aq, fMin_Brho, fMax_Brho, fMin_RawTof, fMax_RawTof, fMin_Beta, fMax_Beta;
    Double_t fMin_Z_gate, fMax_Z_gate, fMin_Aq_gate, fMax_Aq_gate;

    // Canvas
    TCanvas* cBeta;
    TCanvas* cTof;
    TCanvas* cBrho;
    TCanvas* cAqvsq;
    TCanvas* cXs2vsBeta;
    TCanvas* cAoQvsPosS2;
    TCanvas* cIsoGated;
    TCanvas* cLosE_Tof;
    TCanvas* cLosE_Tof2;

    // Histograms for Hit data
    TH1F* fh1_beta;
    TH1F* fh1_tof;
    TH1F* fh1_brho;
    TH2F* fh2_Aqvsq;
    TH2F* fh2_Xs2vsbeta;
    TH2F* fh2_Pos2vsAoQ_m1;
    TH2F* fh2_Z_xc;
    TH2F* fh2_IsoGated_Z_xc;
    TH2F* fh2_IsoGated_xs2_xc;
    TH2F* fh2_IsoGated_xc_anglec;
    TH2F* fh2_LosE_Tof;
    TH1F* fh1_LosE_withTof;
    TH1F* fh1_LosE_withoutTof;

  public:
    ClassDef(R3BIncomingIDOnlineSpectra, 1)
};

#endif /* R3BIncomingIDOnlineSpectra_H */
