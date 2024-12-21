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

#ifndef R3BTOFDDIGITIZERCAL_H
#define R3BTOFDDIGITIZERCAL_H 1

#include "FairTask.h"
#include "R3BTofdCalData.h"
#include "R3BTofdDigi.h"
#include "R3BTofdDigiPar.h"
#include "R3BTofdPoint.h"
#include <map>
#include <string>

class TClonesArray;
class TH1F;
class TH2F;
class TRandom3;

class R3BTofdDigitizerCal : public FairTask
{

  public:
    /** Default constructor **/
    R3BTofdDigitizerCal();

    /** Destructor **/
    ~R3BTofdDigitizerCal();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    /** Setters for sigmas **/
    void SetSigma_y(Float_t sigma_y) { ysigma = sigma_y; }
    void SetSigma_t(Float_t sigma_t) { tsigma = sigma_t; }
    void SetSigma_ELoss(Float_t sigma_ELoss) { esigma = sigma_ELoss; }

  protected:
    TClonesArray* fTofdPoints;
    TClonesArray* fMCTrack;
    TClonesArray* fTofdCals;
    TClonesArray* fCalTriggerItems;

  private:
    TRandom3* prnd;
    Float_t ysigma;
    Float_t tsigma;
    Float_t esigma;

    Int_t maxevent;
    Int_t counter = 0;
    // detector parameter
    Int_t number_layers = 2;
    Int_t number_paddles = 44; // per layer
    Int_t number_sides = 2;
    Int_t number_channels = number_layers * 100 + 50;
    Float_t paddle_width = 2.70000;
    Float_t paddle_thickness = 0.50000;
    Float_t air_gap_paddles = 0.04;
    Float_t air_gap_layer = 5.;

    Float_t detector_width = number_paddles * paddle_width + (number_paddles - 1) * air_gap_paddles + paddle_width;
    Float_t detector_thickness = (number_layers - 1) * air_gap_layer + number_layers * paddle_thickness;

    ClassDef(R3BTofdDigitizerCal, 1);
};

#endif
