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

#ifndef R3BDCHDIGITIZER_H
#define R3BDCHDIGITIZER_H 1

#include "FairTask.h"
#include "R3BDchDigi.h"
#include "R3BDchDigiPar.h"
#include <map>
#include <string>

class TClonesArray;
class TObjectArray;
class TH1F;
class TH2F;

class R3BDchDigitizer : public FairTask
{

  public:
    /** Detector resolution **/
    int detectorResolution; // 0=all off, 1= all on, ...

    /** Default constructor **/
    R3BDchDigitizer();

    /** Detector resolution constructor **/
    R3BDchDigitizer(int detres);

    /** Destructor **/
    ~R3BDchDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Finish();
    virtual void Reset();

    R3BDchDigi* AddHit(Int_t pd1mul,
                       Double_t Pdx1_p1,
                       Double_t Pdy1_p1,
                       Int_t pd2mul,
                       Double_t Pdx2_p1,
                       Double_t Pdy2_p1);

  protected:
    TClonesArray* fDchPoints;
    TClonesArray* fDchMCTrack;
    TClonesArray* fDchDigi;

    // Parameter class
    R3BDchDigiPar* fDchDigiPar;

    //- Control Hitograms
    TH1F* DCH1Px;
    TH1F* DCH2Px;

    TH1F* DCH1Py;
    TH1F* DCH2Py;

    TH1F* DCH1Pz;
    TH1F* DCH2Pz;

    TH1F* DCH1X;
    TH1F* DCH2X;

    TH1F* DCH1Y;
    TH1F* DCH2Y;

    TH1F* DCH1elosshis;
    TH1F* DCH2elosshis;

    TH1F* TrackPx;
    TH1F* TrackPy;
    TH1F* TrackPz;

    TH2F* TrackPxVSDCH1Px;

    Int_t eventNoDch;

  private:
    virtual void SetParContainers();

    ClassDef(R3BDchDigitizer, 1);
};

#endif
