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
#include <map>
#include <string>

class TClonesArray;
class TRandom3;

class R3BTofDDigitizerCal : public FairTask
{
  public:
    /** Default constructor **/
    R3BTofDDigitizerCal();

    /** Destructor **/
    virtual ~R3BTofDDigitizerCal();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    virtual void Reset();

    /** Setters for sigmas **/
    void SetSigma_y(Float_t sigma_y) { ysigma = sigma_y; }
    void SetSigma_t(Float_t sigma_t) { tsigma = sigma_t; }
    void SetSigma_ELoss(Float_t sigma_ELoss) { esigma = sigma_ELoss; }

  private:
    TClonesArray* fTofdPoints;
    TClonesArray* fMCTrack;
    TClonesArray* fTofdCals;
    TClonesArray* fCalTriggerItems;

    TRandom3* prnd;
    Float_t ysigma;
    Float_t tsigma;
    Float_t esigma;

    Int_t maxevent;
    // detector parameter
    Int_t number_layers;
    Int_t number_paddles; // per layer
    Int_t number_sides;
    Int_t number_channels;

    ClassDef(R3BTofDDigitizerCal, 1);
};

#endif /* R3BTOFDDIGITIZERCAL_H */
