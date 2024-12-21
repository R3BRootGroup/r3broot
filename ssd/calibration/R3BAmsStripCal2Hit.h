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

// ----------------------------------------------------------------------
// -----                                                            -----
// -----                     R3BAmsStripCal2Hit                     -----
// -----             Created 01/06/18 by J.L. Rodriguez-Sanchez    -----
// ----------------------------------------------------------------------

#ifndef R3BAmsStripCal2Hit_H
#define R3BAmsStripCal2Hit_H 1

#include "FairTask.h"

#include "R3BAmsHitData.h"
#include "R3BAmsStripCal2Hit.h"
#include "R3BAmsStripCalData.h"

#include "TVector3.h"
#include <Rtypes.h>

class TClonesArray;
class R3BAmsMappingPar;

class R3BAmsStripCal2Hit : public FairTask
{

  public:
    /** Default constructor **/
    R3BAmsStripCal2Hit();

    /** Standard constructor **/
    R3BAmsStripCal2Hit(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BAmsStripCal2Hit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method Reset **/
    virtual void Reset();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method Finish **/
    virtual void Finish();

    /** Virtual method SetParContainers **/
    virtual void SetParContainers();

    /** Accessor for selecting online mode **/
    void SetOnline(Bool_t option) { fOnline = option; }

    /** Accessor for selecting max. number of clusters per ams detector **/
    void SetMaxNumClusters(Int_t max) { fMaxNumClusters = max; }

    /** Accessor to set up the threshold for the cluster energy sum **/
    void SetClusterEnergy(Float_t thsum) { fThSum = thsum; }

  private:
    void SetParameter();
    virtual void DefineClusters(Int_t* nfound, Double_t fPitch, Double_t* fChannels, TH1F* hsst, Double_t cluster[][2]);

    Double_t fPitchK, fPitchS;
    Double_t fScen, fKcen;
    Float_t fThSum;
    Int_t fMaxNumDet, fMaxNumClusters;
    TH1F* hams[16];

    R3BAmsMappingPar* fMap_Par;       /**< Parameter container with mapping. >*/
    TClonesArray* fAmsStripCalDataCA; /**< Array with AMS Cal-input data. >*/
    TClonesArray* fAmsHitDataCA;      /**< Array with AMS Hit-output data. >*/

    Bool_t fOnline; // Don't store data for online
    Double_t* fChannelPeaks;

    /** Private method AddHitData **/
    //** Adds a AmsHitData to the HitCollection
    R3BAmsHitData* AddHitData(Int_t detid,
                              Int_t numhit,
                              Double_t s,
                              Double_t k,
                              TVector3 master,
                              Double_t energy_s,
                              Double_t energy_k,
                              Int_t mulS,
                              Int_t mulK);

  public:
    // Class definition
    ClassDef(R3BAmsStripCal2Hit, 1)
};

#endif
