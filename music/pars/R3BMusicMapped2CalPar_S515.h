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
// -----                R3BMusicMapped2CalPar_S515             -----
// -----            Created 23/03/2022 M. Feijoo Fontán        -----
// -----------------------------------------------------------------

#ifndef R3BMusicMapped2CalPar_S515_H
#define R3BMusicMapped2CalPar_S515_H 1

#include "FairTask.h"
#include "R3BMusicMapped2Cal.h"
#include "R3BMusicMappedData.h"
#include "TGraph.h"
#include "TH1F.h"

class TClonesArray;
class R3BMusicCalPar;

class R3BMusicMapped2CalPar_S515 : public FairTask
{

  public:
    /** Default constructor **/
    R3BMusicMapped2CalPar_S515();

    /** Standard constructor **/
    R3BMusicMapped2CalPar_S515(const TString& name,
                               Int_t iVerbose = 1,
                               const TString& namedeta = "Mwpc0",
                               const TString& namedetb = "Los");

    /** Destructor **/
    virtual ~R3BMusicMapped2CalPar_S515();

    /** Virtual method Exec **/
    virtual void Exec(Option_t*);

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Method to set up the position of MwpcA **/
    void SetPosMwpcA(Float_t pos) { fPosMwpcA = pos; }

    /** Method to set up the position of Los **/
    void SetPosLos(Float_t pos) { fPosLos = pos; }

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
    Float_t fPosMwpcA; // Position in the beam direction in mm
    TString fNameDetB;
    Float_t fPosLos;   // Position in the beam direction in mm
    Float_t fPosMusic; // Position in the beam direction in mm

    Int_t mulanode[MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];
    Double_t energy[MAX_MULT_MUSIC_CAL][MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];
    Double_t dtime[MAX_MULT_MUSIC_CAL][MAX_NB_MUSICANODE + MAX_NB_MUSICTREF];

    R3BMusicCalPar* fCal_Par;         /**< Parameter container. >*/
    TClonesArray* fMusicMappedDataCA; /**< Array with Music Mapped-input data. >*/
    TClonesArray* fHitItemsMwpcA;     /**< Array with hit items. */
    TClonesArray* fHitItemsLos;       /**< Array with hit items. */

    TGraph** fg_anode;

  public:
    // Class definition
    ClassDef(R3BMusicMapped2CalPar_S515, 1)
};

#endif
