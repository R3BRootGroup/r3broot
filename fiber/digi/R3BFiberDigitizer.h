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

// ---------------------------------------------------------------
// -----                    R3BFiberDigitizer                -----
// -----          Created 10/05/2022 by J.L. Rodriguez       -----
// ---------------------------------------------------------------

#ifndef R3BFiberDigitizer_H
#define R3BFiberDigitizer_H 1

#include "FairTask.h"
#include "R3BBunchedFiberHitData.h"
#include "TRotation.h"
#include "TVector3.h"
#include <TRandom3.h>
#include <string>

class TClonesArray;
class TH1F;
class TH2F;
class R3BTGeoPar;

class R3BFiberDigitizer : public FairTask
{

  public:
    /** Default constructor **/
    R3BFiberDigitizer(const TString& name = "Fi");

    R3BFiberDigitizer(const TString& name, Double_t esigma, Double_t tsigma, Double_t ysigma);

    /** Destructor **/
    virtual ~R3BFiberDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();
    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Reset();
    virtual void SetParContainers();

    void SetEnergyResolution(Double_t e);
    void SetTimeResolution(Double_t t);
    void SetYPositionResolution(Double_t y);
    void SetMinPID(Int_t min) { fMinPID = min; }

  private:
    void SetParameter();
    TString fName;
    TRandom3* rand;
    Int_t fMinPID;
    Double_t esigma;
    Double_t tsigma;
    Double_t ysigma;
    Double_t xsigma;
    R3BTGeoPar* fFiGeoPar;
    TRotation fRot;
    TVector3 fTrans;

    TClonesArray* fFiPoints;
    TClonesArray* fFiHits;
    TClonesArray* fMCTrack;

    Float_t fiber_thickness;
    Int_t fiber_nbr;
    Float_t air_layer;
    Float_t detector_width;

    /** Private method AddHitData **/
    // Adds a R3BFiberHitData
    R3BBunchedFiberHitData* AddHitData(Int_t DetId,
                                       Double_t x,
                                       Double_t y,
                                       Double_t Eloss,
                                       Double_t time,
                                       Int_t fiber,
                                       Double_t a_bottom_time_ns,
                                       Double_t a_top_time_ns,
                                       Double_t a_bottom_tot_ns,
                                       Double_t a_top_tot_ns);

    ClassDef(R3BFiberDigitizer, 1);
};

#endif /* R3BFiberDigitizer_H */
