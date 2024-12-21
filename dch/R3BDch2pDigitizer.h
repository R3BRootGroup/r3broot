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

#ifndef R3BDCH2PDIGITIZER_H
#define R3BDCH2PDIGITIZER_H 1

#include "FairTask.h"
#include "R3BDch2pDigi.h"
#include "R3BDch2pDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BDch2pDigitizer : public FairTask
{

  public:
    /** Detector resolution **/
    int detectorResolution; // 0=all off, 1= all on, ...

    /** Default constructor **/
    R3BDch2pDigitizer();

    /** Detector resolution constructor **/
    R3BDch2pDigitizer(int detres);

    /** Destructor **/
    ~R3BDch2pDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BDch2pDigi* AddHit(Int_t pd1mul,
                         Double_t Pdx1_p1,
                         Double_t Pdy1_p1,
                         Double_t Pdx1_p2,
                         Double_t Pdy1_p2,
                         Int_t pd2mul,
                         Double_t Pdx2_p1,
                         Double_t Pdy2_p1,
                         Double_t Pdx2_p2,
                         Double_t Pdy2_p2);

  protected:
    TClonesArray* fDch2pPoints;
    TClonesArray* fDch2pMCTrack;
    TClonesArray* fDch2pDigi;

    // Parameter class
    R3BDch2pDigiPar* fDch2pDigiPar;

    //- Control Hitograms

    Int_t eventNoDch;

  private:
    virtual void SetParContainers();

    ClassDef(R3BDch2pDigitizer, 1);
};

#endif
