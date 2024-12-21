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

#ifndef R3BMusliCal2Hit_H
#define R3BMusliCal2Hit_H 1

#include "TArrayD.h"

#include "FairTask.h"
#include "R3BFrsData.h"
#include "R3BFrsSciTofCalData.h"
#include "R3BMusliCalData.h"
#include "R3BMusliHitData.h"

#define MAX_MULT_MUSLI 100
#define MAX_NUM_GROUPS_CAL 15
#define MAX_NUM_TYPES_HIT 4

class TClonesArray;
class R3BMusliHitPar;
class R3BTGeoPar;
// class R3BEventHeader;

class R3BMusliCal2Hit : public FairTask
{

  public:
    /** Default constructor **/
    R3BMusliCal2Hit();

    /** Standard constructor **/
    R3BMusliCal2Hit(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BMusliCal2Hit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

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
    void SetExpId(Int_t exp)
    {
        fExpId = exp;
    } // Mutator to set fExpId manually. It should be globally defined by EventHeader.
    void SetDirectBeta(Double_t beta) { fDirectBeta = beta; }
    void SetIdS2(UShort_t id) { fIdS2 = id; }
    void SetIdCaveC(UShort_t id) { fIdCaveC = id; }

  private:
    void SetParameters();
    Double_t BetaCorr_pol1(Double_t beta, Double_t p0, Double_t p1) { return p0 + p1 * beta; }
    Double_t BetaCorr_std(Double_t beta, Double_t p0, Double_t p1) { return p0 + p1 * pow(beta, -5. / 3.); }
    Double_t SqrtE2Z(Double_t SqrtE, Double_t p0, Double_t p1, Double_t p2)
    {
        return p0 + p1 * SqrtE + p2 * SqrtE * SqrtE;
    }

    //    R3BEventHeader* header; /**< Event header. */

    Int_t fNumGroups;
    Int_t fNumTypes;
    Int_t fMaxMult;
    Int_t fExpId;
    Int_t fNumEcorrBeta;
    Int_t fNumZ;
    Double_t fDirectBeta; // if direct beam from SIS
    UShort_t fIdS2;
    UShort_t fIdCaveC; // for beta we want to use the S2 -> Cave C Beta
    TArrayD* fEaveVsBetaHitParams;
    TArrayD* fEcorrBetaHitParams;
    TArrayD* fZHitParams;

    UInt_t mult_cal[MAX_NUM_GROUPS_CAL];

    // variables to calculate Zcharge
    Double_t e_cal[MAX_MULT_MUSLI][MAX_NUM_GROUPS_CAL];
    Double_t nba[MAX_NUM_TYPES_HIT];
    Double_t Esum[MAX_NUM_TYPES_HIT];
    Double_t eave[MAX_NUM_TYPES_HIT];
    Double_t e_hit[MAX_NUM_TYPES_HIT];
    Double_t z_hit[MAX_NUM_TYPES_HIT];

    // variables to calculate theta
    Double_t dt_cal[MAX_MULT_MUSLI][MAX_NUM_GROUPS_CAL];
    Double_t x_hit[MAX_NUM_TYPES_HIT];
    Double_t theta_hit[MAX_NUM_TYPES_HIT];

    R3BMusliHitPar* fHit_Par;      /**< Parameter container. > */
    R3BTGeoPar* fMusliGeo_Par;     /**< Parameter container. > */
    TClonesArray* fMusliCalDataCA; /**< Array with Musli Cal-input data. >*/
    TClonesArray* fMusliHitDataCA; /**< Array with Musli Hit-output data. >*/
    TClonesArray* fFrsDataCA;      /**< Array with Frs input data (for beta). >*/
    Bool_t fFrsSciTofCal;          // need to know if FrsSciTofCal or FrsData is used

    Bool_t fOnline; // Don't store data for online

    /** Private method AddHitData **/
    //** Adds a MusliHitData to the HitCollection
    R3BMusliHitData* AddHitData(UInt_t t, Double_t e, Double_t z, Double_t x, Double_t theta);

  public:
    // Class definition
    ClassDef(R3BMusliCal2Hit, 1)
};

#endif /* R3BMusliCal2Hit_H */
