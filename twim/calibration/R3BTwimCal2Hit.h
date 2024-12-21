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

// -------------------------------------------------------------
// -----             R3BTwimCal2Hit                        -----
// -----    Created 30/11/19 by J.L. Rodriguez-Sanchez     -----
// -------------------------------------------------------------

#pragma once

#include "FairTask.h"
#include "R3BLogger.h"
#include "R3BTwimHitData.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVectorD.h"
#include <TRandom.h>

class TClonesArray;
class R3BTwimHitPar;
class R3BEventHeader;

class R3BTwimCal2Hit : public FairTask
{

  public:
    /** Default constructor **/
    R3BTwimCal2Hit();

    /** Standard constructor **/
    R3BTwimCal2Hit(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BTwimCal2Hit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method Reset **/
    virtual void Reset();

    virtual void SetParContainers();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Method to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }
    void SetExpId(Int_t exp)
    {
        R3BLOG(info, "fExpId is set locally. Original:" << fExpId << ", New value:" << exp);
        R3BLOG(info, "Using R3BEventHeader::SetExpId() is recommended instead.");
        fExpId = exp;
    }

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    inline void SetTpat(Int_t tpat) { fTpat = tpat; }

    inline void SetDebug() { fDebug = kTRUE; }

  private:
    void SetParameter();
    /** Private method Experiment s4551 => s455+fission **/
    virtual void S4551();
    /** Private method Experiment s455 **/
    virtual void S455();
    /** Private method Experiment s467 **/
    virtual void S467();

    R3BEventHeader* header; /**< Event header. */

    Bool_t fOnline; // Don't store data for online
    Int_t fExpId;
    Int_t fTpat;
    Bool_t fDebug;

    Int_t fNumSec; // Number of sections
    Int_t fNumAnodes;
    Int_t fNumAnodesAngleFit;
    // Int_t fNumAnodesAngleFit4551[4];
    Int_t fIndex[16];
    Int_t fNumParamsTof;
    Int_t fNumParams;
    Int_t fMaxEnergyperanode;
    Float_t fZ0[4], fZ1[4], fZ2[4], fZ3[4], fZ4[4];
    Int_t StatusAnodes[4][16]; // Sections and anodes
    Double_t fPosAnodes[16];   // Position-Z of each anode
    TArrayF* CalZTofParams;
    TArrayF* CalZParams;
    Float_t fEmean_tof[4];
    Float_t fEmean_dt[4];
    TVectorD fPosZ;

    R3BTwimHitPar* fCal_Par;      /**< Parameter container. >*/
    TClonesArray* fTwimCalDataCA; /**< Array with Twim Cal-input data. >*/
    TClonesArray* fTwimHitDataCA; /**< Array with Twim Hit-output data. >*/
    TClonesArray* fHitItemsTofW;  /**< Array with tofw-hit items. */

    /** Private method TwimHitData **/
    //** Adds a TwimHitData to the detector
    R3BTwimHitData* AddHitData(UShort_t secID,
                               Double_t theta,
                               Double_t charge_z,
                               Double_t xpos,
                               Double_t offset = 0.0,
                               Double_t ene_ave = 0.0);

    TH1F** dx;
    TH1F** dx_all;
    TH2F** dx_vs_pos;

    TCanvas** dx_canvas;
    TCanvas** dx_all_canvas;
    TCanvas** dx_vs_pos_canvas;

  public:
    // Class definition
    ClassDef(R3BTwimCal2Hit, 1)
};
