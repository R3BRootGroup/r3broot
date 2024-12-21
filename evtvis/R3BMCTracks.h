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

#pragma once

#include <FairTask.h>

#include <Rtypes.h>              // for Double_t, etc
#include <TEveTrackPropagator.h> // IWYU pragma: keep needed by cint
#include <TString.h>             // for TString

class FairEventManager;
class TClonesArray;
class TEveTrackList;
class TObjArray;
class TParticle;

class R3BMCTracks : public FairTask
{
  public:
    /** Default constructor **/
    R3BMCTracks();

    /** Standard constructor
     *@param name        Name of task
     *@param iVerbose    Verbosity level
     **/
    R3BMCTracks(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BMCTracks();

    /** Set verbosity level. For this task and all of the subtasks. **/
    void SetVerbose(Int_t iVerbose) { fVerbose = iVerbose; }
    /** Executed task **/
    virtual void Exec(Option_t* option);
    virtual InitStatus Init();
    virtual void SetParContainers();

    /** Action after each event**/
    virtual void Finish();
    void Reset();
    TEveTrackList* GetTrGroup(TParticle* P);

  protected:
    TClonesArray* fTrackList; //!
    TEveTrackPropagator* fTrPr;
    FairEventManager* fEventManager; //!
    TObjArray* fEveTrList;
    TString fEvent;         //!
    TEveTrackList* fTrList; //!
    // TEveElementList *fTrackCont;

    Double_t MinEnergyLimit;
    Double_t MaxEnergyLimit;
    Double_t PEnergy;

  private:
    R3BMCTracks(const R3BMCTracks&);
    R3BMCTracks& operator=(const R3BMCTracks&);

    ClassDef(R3BMCTracks, 1);
};
