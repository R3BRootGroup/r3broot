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

#ifndef _R3B_LAND_DIGITIZER_QA_H_
#define _R3B_LAND_DIGITIZER_QA_H_

#include "FairTask.h"

class TClonesArray;
class TH1F;
class TH2F;

/**
 * \brief R3BLandDigitizerQA task.
 *
 * Task for filling controll histograms of R3BLandDigitizer
 * performance. Plots are stored in the output file.
 */

class R3BLandDigitizerQA : public FairTask
{

  public:
    /** Default constructor. */
    R3BLandDigitizerQA();

    /** Copy constructor. */
    R3BLandDigitizerQA(const R3BLandDigitizerQA&);

    /** Destructor. */
    virtual ~R3BLandDigitizerQA();

    /** Assignment operator. */
    R3BLandDigitizerQA& operator=(const R3BLandDigitizerQA&) { return *this; }

    /** Initialization of task. */
    virtual InitStatus Init();

    /** Execute one event. */
    virtual void Exec(Option_t* option);

    /** Finish event. */
    virtual void FinishEvent();

    /** Finish task. */
    virtual void Finish();

    /** Specify beam energy to be used. */
    inline void SetBeamEnergy(const Double_t& beamEnergy) { fBeamEnergy = beamEnergy; }

  private:
    Int_t fEvent;
    TClonesArray* fTracks;
    TClonesArray* fPoints;
    TClonesArray* fDigis;

    Double_t fBeamEnergy;

    TH1F* fhTotEnergy;
    TH1F* fhTotLight;
    TH1F* fhTotEnergyRel;
    TH1F* fhTotLightRel;
    TH2F* fhElossLight;
    TH1F* fhTotEnergyReco;
    TH1F* fhPaddleE;
    TH1F* fhFirstTvsBar;
    TH1F* fhFirstBvsBar;
    TH1F* fhFirstE;
    TH2F* fhElossTime;
    TH2F* fhElossPdg;
    TH2F* fhQdcXL;
    TH2F* fhQdcXR;
    TH2F* fhThresh;
    TH1F* fhQdcRatio;

    void CreateHisto();

    ClassDef(R3BLandDigitizerQA, 1)
};

#endif
