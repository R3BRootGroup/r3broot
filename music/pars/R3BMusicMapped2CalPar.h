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
// -----                R3BMusicMapped2CalPar                  -----
// -----        Created 29/01/20 by J.L. Rodriguez-Sanchez     -----
// -----------------------------------------------------------------

#ifndef R3BMusicMapped2CalPar_H
#define R3BMusicMapped2CalPar_H 1

#include "FairTask.h"
#include "R3BMusicMapped2Cal.h"
#include "R3BMusicMappedData.h"
#include "TGraph.h"
#include "TH1F.h"

class TClonesArray;
class R3BCoarseTimeStitch;
class R3BMusicCalPar;

class R3BMusicMapped2CalPar : public FairTask
{
  public:
    /** Default constructor **/
    R3BMusicMapped2CalPar();

    /** Standard constructor **/
    R3BMusicMapped2CalPar(const TString& name,
                          Int_t iVerbose = 1,
                          const TString& namedeta = "Mwpc0",
                          const TString& namedetb = "Mwpc2");

    /** Destructor **/
    virtual ~R3BMusicMapped2CalPar();

    /** Virtual method Exec **/
    virtual void Exec(Option_t*);

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Method to set up the position of DetA **/
    void SetPosDetA(Float_t pos) { fPosDetA = pos; }

    /** Method to set up the position of DetB **/
    void SetPosDetB(Float_t pos) { fPosDetB = pos; }

    /** Method to set up the position of Music **/
    void SetPosMusic(Float_t pos) { fPosMusic = pos; }

    /** Method to set up the limits for fit **/
    void SetFitLimits(Int_t left, Int_t right)
    {
        fLimit_left = left;
        fLimit_right = right;
    }

  private:
    Int_t fNumAnodes;
    Int_t fMaxMult;
    Int_t fMinStadistics;
    Int_t fNumParams;
    Int_t fNumPosParams;
    Int_t fNumAnodesRef;
    Int_t fLimit_left;
    Int_t fLimit_right;
    Int_t fMaxSigma;
    TArrayF* CalParams;
    TArrayF* PosParams;

    TString fNameDetA;
    Float_t fPosDetA; // Position in the beam direction in mm
    TString fNameDetB;
    Float_t fPosDetB;  // Position in the beam direction in mm
    Float_t fPosMusic; // Position in the beam direction in mm

    Int_t mulanode[MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];
    Double_t energy[MAX_MULT_MUSIC_CAL][MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];
    Double_t dtime[MAX_MULT_MUSIC_CAL][MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];

    R3BCoarseTimeStitch* fTimeStitch;
    R3BMusicCalPar* fCal_Par;         /**< Parameter container. >*/
    TClonesArray* fMusicMappedDataCA; /**< Array with Music Mapped-input data. >*/
    TClonesArray* fHitItemsDetA;      /**< Array with hit items. */
    TClonesArray* fHitItemsDetB;      /**< Array with hit items. */

    TGraph** fg_anode;

  public:
    // Class definition
    ClassDef(R3BMusicMapped2CalPar, 1)
};

#endif
