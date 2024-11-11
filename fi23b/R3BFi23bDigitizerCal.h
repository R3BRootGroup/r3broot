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

// R3BFi23bDigitizerCal.h

#ifndef R3BFI23BDIGITIZERCAL_H
#define R3BFI23BDIGITIZERCAL_H 1

#include "FairTask.h"
#include "R3BFiberMAPMTCalData.h"
#include <TRandom3.h>
#include <map>
#include <string>

class TClonesArray;
class TH1F;
class TH2F;

class R3BFi23bDigitizerCal : public FairTask
{

  public:
    /** Default constructor **/
    R3BFi23bDigitizerCal();

    R3BFi23bDigitizerCal(const char* name, Int_t iVerbose = 1);
   
    /** Destructor **/
   virtual  ~R3BFi23bDigitizerCal();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    
    virtual void Reset();

    inline void SetEnergyResolution(Double_t e)
    {
		esigma = e;
    }    
    inline void SetTimeResolution(Double_t t)
    {
	    tsigma = t;	
	}	
    inline void SetYPositionResolution(Double_t y)
    {
		ysigma = y;
	}	

  protected:
    TClonesArray* fFi23bPoints;

    TClonesArray* fMCTrack;

    TClonesArray* fFi23bCals;
    TClonesArray* fFi23bTriggerCals;

  private:
    TRandom3* prnd;
    Double_t esigma;
    Double_t tsigma;
    Double_t ysigma;

   Float_t detector_height = 10.752;
    Float_t fiber_thickness = 0.028000;
    Int_t fiber_nbr = 384;
    Float_t dead_layer = 0.9;
    Float_t air_layer = 0.01*0.; // relative to fiber_thickness
    Float_t detector_width = fiber_nbr * fiber_thickness * (1 + air_layer);

    ClassDef(R3BFi23bDigitizerCal, 1);
};

#endif
