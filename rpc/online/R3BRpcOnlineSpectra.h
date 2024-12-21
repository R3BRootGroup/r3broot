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

// -----------------------------------------------------------
// -----             R3BRpcOnlineSpectra                 -----
// -----   Created 23/01/22 by J.L. Rodriguez-Sanchez    -----
// -----         Fill RPC online histograms              -----
// -----------------------------------------------------------

#ifndef R3BRpcOnlineSpectra_H
#define R3BRpcOnlineSpectra_H

#include "FairTask.h"

#include "TCanvas.h"
#include "TLatex.h"
#include "TMath.h"
#include <Rtypes.h>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

class TClonesArray;
class R3BEventHeader;
class TH1F;
class TH2F;

class R3BRpcOnlineSpectra : public FairTask
{

  public:
    /* --------- Default Constructor ---------- */
    R3BRpcOnlineSpectra();

    /* --------- Standard Constructor ---------- */
    R3BRpcOnlineSpectra(const TString& name, Int_t iVerbose = 1);

    /* ------------ Standard  Destructor ----------- */
    virtual ~R3BRpcOnlineSpectra();

    /* ------------ Init Method ----------- */
    virtual InitStatus Init() override;

    /* ------------ Exec Method ----------- */
    virtual void Exec(Option_t* option) override;

    /* ------------ Finish Event ----------- */
    virtual void FinishEvent() override;

    /*-------- Finish Task : After the eventloop -------*/
    virtual void FinishTask() override;

    /* ------- Reset Methods : Clear histograms and time counters -------*/
    void Reset_RPC_Histo();
    void Reset_Efficiencies();

    /* ------ Choose to fill online histograms only if onspill or offspill ------ */
    void GO_ONSPILL();
    void GO_OFFSPILL();

    void SetTofAxisRange(Int_t bins, Float_t leftLimit, Float_t rightLimit)
    {

        fLeftTofLim = leftLimit;
        fRightTofLim = rightLimit;
        fTofBins = bins;
    };

    void SetRpcTimeAxisRange(Int_t bins, Float_t leftLimit, Float_t rightLimit)
    {

        fLeftRpcTimeLim = leftLimit;
        fRightRpcTimeLim = rightLimit;
        fRpcTimeBins = bins;
    };

    void SetRpcToTAxisRange(Int_t bins, Float_t leftLimit, Float_t rightLimit)
    {

        fLeftRpcToTLim = leftLimit;
        fRightRpcToTLim = rightLimit;
        fRpcToTBins = bins;
    };

    void SetOnspillTPatRange(Int_t tpat1, Int_t tpat2)
    {

        fFirstTPat = tpat1;
        fLastTPat = tpat2;
    };

    void SetTrigger(Int_t trigg) { fTrigger = trigg; }

  private:
    R3BEventHeader* fEventHeader;
    int fcounter;
    TClonesArray* fMappedDataItems;
    TClonesArray* fPreCalDataItems;
    TClonesArray* fCalDataItems;
    TClonesArray* fHitDataItems;
    TClonesArray* fPmtHitDataItems;

    Float_t meanCharges[41][50] = { { 0.0 } };
    Float_t stripHitsVec[41] = { 0.0 };

    Int_t counts = 0;
    Int_t tpatbin, fTPat;
    int64_t fTimeStart, fTimeEnd;
    Bool_t fSpill;
    Bool_t execBool;

    Float_t fLeftTofLim;
    Float_t fRightTofLim;

    Float_t fLeftRpcTimeLim;
    Float_t fRightRpcTimeLim;

    Float_t fLeftRpcToTLim;
    Float_t fRightRpcToTLim;

    Int_t fTofBins;
    Int_t fRpcTimeBins;
    Int_t fRpcToTBins;
    Float_t fBarHits;
    Float_t fStrip21Hits;

    R3BEventHeader* header; // Event header.
    Int_t fTrigger;         // Trigger value.
    Int_t fNEvents;         // Event counter.

    /* ----- Map Histograms -----*/
    TH1F* stripCoarseRightHisto[41];
    TH1F* stripFineRightHisto[41];

    TH1F* stripCoarseLeftHisto[41];
    TH1F* stripFineLeftHisto[41];

    TH1F* pmtCoarseHistoTop[8];
    TH1F* pmtFineHistoTop[8];
    TH1F* pmtCoarseHistoBottom[8];
    TH1F* pmtFineHistoBottom[8];

    TH1F* refCoarseHisto[9];
    TH1F* refFineHisto[9];

    TH2F* stripCoarseLeftCorr;
    TH2F* stripCoarseRightCorr;

    TH2F* stripFineLeftCorr;
    TH2F* stripFineRightCorr;

    /* ----- Map Canvases*/
    TCanvas* leftStripCanvasCoarse;
    TCanvas* rightStripCanvasCoarse;

    TCanvas* leftStripCanvasFine;
    TCanvas* rightStripCanvasFine;

    TCanvas* stripCoarseLeftCorrCanvas;
    TCanvas* stripCoarseRightCorrCanvas;

    TCanvas* stripFineLeftCorrCanvas;
    TCanvas* stripFineRightCorrCanvas;

    TCanvas* pmtCoarseCanvas;
    TCanvas* pmtFineCanvas;
    TCanvas* refFineCanvas;
    TCanvas* refCoarseCanvas;

    /* ----- Pre Cal Histograms ----- */
    TH2F* stripLeftTotCorr;
    TH2F* stripRightTotCorr;
    TH2F* stripLeftTimeCorr;
    TH2F* stripRightTimeCorr;
    TH2F* stripLeftBanana;
    TH2F* stripRightBanana;
    TH1F* pmtPreCalTimeHistoTop[8];
    TH1F* pmtPreCalTotHistoTop[8];
    TH1F* pmtPreCalTimeHistoBottom[8];
    TH1F* pmtPreCalTotHistoBottom[8];

    /* ----- Pre Cal Canvases ----- */
    TCanvas* stripLeftTotCorrCanvas;
    TCanvas* stripRightTotCorrCanvas;
    TCanvas* stripLeftTimeCorrCanvas;
    TCanvas* stripRightTimeCorrCanvas;
    TCanvas* stripLeftBananaCanvas;
    TCanvas* stripRightBananaCanvas;
    TCanvas* pmtPreCalTimeCanvas;
    TCanvas* pmtPreCalTotCanvas;

    /* ----- Cal Histograms ----- */
    TH2F* stripCalTimeCorr;
    TH2F* stripCalToTCorr;

    /* ----- Cal Canvases ----- */
    TCanvas* stripCalTimeCorrCanvas;
    TCanvas* stripCalTotCorrCanvas;

    /* ----- Hit Histograms ----- */
    TH2F* stripPosHitCorr;
    TH2F* meanChargeCorr;
    TH1F* totalChargeHist;

    /* ----- Hit Canvases ----- */
    TCanvas* hitMapCanvas;

    /* ------ ToF Correlations Histograms ------*/
    TH2F* tofCorr;
    TH1F* stripTofHisto[41];
    TH2F* tofVsPosCorr;

    /* ------ ToF Correlations Canvases ------*/
    TCanvas* tofCorrCanvas;
    TCanvas* stripTofCanvas;
    TCanvas* tofPosCanvas;

    /* ------ Strip - Bar Correlation -----*/
    TCanvas* timeDiffStripPmtCanvas;
    TCanvas* timeDiffStripPmtCorrCanvas;
    TCanvas* posRpcCorrCanvas;
    TCanvas* timeDiffLosPmtCanvas;
    TCanvas* stripIdVsHitsCanvas;

    /* ------ Strip - Bar Correlation Histograms -----*/
    TH1F* timeDiffStripPmtHisto[41];
    TH2F* timeDiffStripPmtCorr;

    TH2F* deltaTVsPosCorr;
    TH2F* posRpcVsPosBar;

    TH2F* stripIdVsNHitsCorr;

    /*------------ Efficiency Canvases --------*/
    TCanvas* stripMultCanvas;
    TCanvas* efficiencyCanvas;

    /*--------- Efficiency Histograms ----------*/
    TH1F* stripMultHisto;
    TH1F* hEfficiencyHisto_H;
    TH1F* hStripEffHisto_H;
    TH1F* hBarEffHisto_H;

    TH1F* hEfficiencyHisto_V1;
    TH1F* hStripEffHisto_V1;
    TH1F* hBarEffHisto_V1;

    TH1F* hEfficiencyHisto_V2;
    TH1F* hStripEffHisto_V2;
    TH1F* hBarEffHisto_V2;

    Int_t fFirstTPat;
    Int_t fLastTPat;

    TLatex *tex1, *tex2;

  public:
    ClassDefOverride(R3BRpcOnlineSpectra, 1)
};

#endif /* R3BRpcOnlineSpectra_H */
