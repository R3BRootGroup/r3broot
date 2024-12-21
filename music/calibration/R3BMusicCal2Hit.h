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
// -----                                                   -----
// -----             R3BMusicCal2Hit                       -----
// -----    Created 30/11/19  by J.L. Rodriguez-Sanchez    -----
// -------------------------------------------------------------

#ifndef R3BMusicCal2Hit_H
#define R3BMusicCal2Hit_H

#include "FairTask.h"
#include "R3BMusicHitData.h"
#include "TH1F.h"
#include "TVectorD.h"
#include <TRandom.h>

class TClonesArray;
class R3BMusicHitPar;

class R3BMusicCal2Hit : public FairTask
{

  public:
    /** Default constructor **/
    R3BMusicCal2Hit();

    /** Standard constructor **/
    R3BMusicCal2Hit(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BMusicCal2Hit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t*);

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

    /** Method for simulation reconstruction **/
    void SetSimulation() { fSim = true; }

  private:
    void SetParameter();
    void ExecSim(int);

    Int_t fNumAnodes;
    Int_t fNumAnodesAngleFit;
    Int_t fNumParams;
    Float_t fZ0, fZ1, fZ2;
    Float_t fx0_point, fy0_point, frot_ang;
    TArrayF* CalZParams;
    TArrayF* CalAngCorParams;
    Int_t fStatusAnodes[8]; // Status anodes
    Double_t fPosAnodes[8]; // Position-Z of each anode
    TVectorD fPosZ;
    bool fSim;
    Bool_t fOnline; // Don't store data for online

    R3BMusicHitPar* fCal_Par;      /**< Parameter container. >*/
    TClonesArray* fMusicCalDataCA; /**< Array with Music Cal-input data. >*/
    TClonesArray* fMusicHitDataCA; /**< Array with Music Hit-output data. >*/

    /** Private method MusicHitData **/
    //** Adds a MusicHitData to the detector
    R3BMusicHitData* AddHitData(Double_t theta, Double_t charge_z, Double_t ene_ave = 0., Double_t good_dt = 0.);

  public:
    // Class definition
    ClassDef(R3BMusicCal2Hit, 1)
};

#endif
