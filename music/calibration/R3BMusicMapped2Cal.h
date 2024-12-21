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

// -----------------------------------------------------------------
// -----                                                       -----
// -----                R3BMusicMapped2Cal                     -----
// -----        Created 24/11/19  by J.L. Rodriguez-Sanchez    -----
// -----------------------------------------------------------------

#ifndef R3BMusicMapped2Cal_H
#define R3BMusicMapped2Cal_H

#include "FairTask.h"
#include "R3BMusicCalData.h"
#include "R3BMusicMappedData.h"
#include "TH1F.h"
#include <TRandom.h>

#define MAX_MULT_MUSIC_CAL 100
#define MAX_NB_MUSICANODE 8
#define MAX_NB_MUSICTREF 2

class TClonesArray;
class R3BCoarseTimeStitch;
class R3BMusicCalPar;

class R3BMusicMapped2Cal : public FairTask
{
  public:
    /** Default constructor **/
    R3BMusicMapped2Cal();

    /** Standard constructor **/
    R3BMusicMapped2Cal(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BMusicMapped2Cal();

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

    /** Accessor to select online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    void SetParameters();

    Int_t fNumAnodes;
    Int_t fMaxMult;
    Int_t fNumParams;
    Int_t fNumPosParams;
    Int_t fNumAnodesRef;
    Int_t fMaxSigma;
    TArrayF* CalParams;
    TArrayF* PosParams;

    Int_t mulanode[MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];
    Double_t energy[MAX_MULT_MUSIC_CAL][MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];
    Double_t dtime[MAX_MULT_MUSIC_CAL][MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];

    Bool_t fOnline; // Don't store data for online

    R3BCoarseTimeStitch* fTimeStitch;
    R3BMusicCalPar* fCal_Par;         /**< Parameter container. >*/
    TClonesArray* fMusicMappedDataCA; /**< Array with Music Mapped-input data. >*/
    TClonesArray* fMusicCalDataCA;    /**< Array with Music Cal-output data. >*/

    /** Private method AddCalData **/
    //** Adds a MusicCalData to the anodeCalCollection
    R3BMusicCalData* AddCalData(UShort_t anodeid, Double_t dtime, Double_t energy);

  public:
    // Class definition
    ClassDef(R3BMusicMapped2Cal, 1)
};

#endif
