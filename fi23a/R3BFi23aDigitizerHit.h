/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// R3BFi23aDigitizer.h

#ifndef R3BFi23ADIGITIZERHIT_H
#define R3BFi23ADIGITIZERHIT_H 1

#include "FairTask.h"
#include "R3BFiberMAPMTHitData.h"
#include <TRandom3.h>
#include <map>
#include <string>

class TClonesArray;
class TH1F;
class TH2F;

class R3BFi23aDigitizerHit : public FairTask
{

  public:
    /** Default constructor **/
    R3BFi23aDigitizerHit();

    R3BFi23aDigitizerHit(Double_t esigma, Double_t tsigma, Double_t xsigma, Double_t ysigma);

    /** Destructor **/
    ~R3BFi23aDigitizerHit();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    void SetEnergyResolution(Double_t e);
    void SetTimeResolution(Double_t t);
    void SetXPositionResolution(Double_t x);
    void SetYPositionResolution(Double_t y);

  protected:
    TClonesArray* fFi23aPoints;
    TClonesArray* fMCTrack;
    TClonesArray* fFi23aHits;

  private:
    TRandom3* prnd;
    Double_t esigma;
    Double_t tsigma;
    Double_t ysigma;
    Double_t xsigma;
	Float_t detector_height = 10.752;
    Float_t fiber_thickness = 0.028000;
    Int_t fiber_nbr = 384;
    Float_t dead_layer = 0.9;
    Float_t air_layer = 0.01*0.; // relative to fiber_thickness
    Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);

    ClassDef(R3BFi23aDigitizerHit, 1);
};

#endif
