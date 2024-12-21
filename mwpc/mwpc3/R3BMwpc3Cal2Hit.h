/******************************************************************************
 **   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 **   Copyright (C) 2019-2025 Members of R3B Collaboration                     *
 **                                                                            *
 **             This software is distributed under the terms of the            *
 **                 GNU General Public Licence (GPL) version 3,                *
 **                    copied verbatim in the file "LICENSE".                  *
 **                                                                            *
 ** In applying this license GSI does not waive the privileges and immunities  *
 ** granted to it by virtue of its status as an Intergovernmental Organization *
 ** or submit itself to any jurisdiction.                                      *
 *******************************************************************************/

//----------------------------------------------------------------------
//-----                        R3BMwpc3Cal2Hit                     -----
//-----             Created 14/10/19 by G. García Jiménez          -----
//-----             by modifying J.L classes for MWPC0             -----
//-----               s455 method 17/06/22 by Antia GG             -----
//----------------------------------------------------------------------

#ifndef R3BMWPC3CAL2HIT_H
#define R3BMWPC3CAL2HIT_H 1

#include "FairTask.h"
#include "R3BMwpcCalData.h"
#include "R3BMwpcHitData.h"
#include "R3BSofTofWHitData.h"
#include "TH1F.h"
#include <TRandom.h>

using namespace std;

#define Mw3PadsX 288
#define Mw3PadsY 120

class TClonesArray;
class R3BEventHeader;

class R3BMwpc3Cal2Hit : public FairTask
{

  public:
    /* Default Constructor */
    R3BMwpc3Cal2Hit();

    /** Standard constructor **/
    R3BMwpc3Cal2Hit(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BMwpc3Cal2Hit();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method Reset **/
    virtual void Reset();

    // Fair specific
    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    void SetOnline(Bool_t option) { fOnline = option; }
    void SetExpId(Int_t exp)
    {
        fExpId = exp;
    } // Mutator to set fExpId manually. It should be globally defined by EventHeader.

  private:
    /** Private method Experiment s455 **/
    virtual void S455();
    /** Private method Experiment s467 **/
    virtual void S467();

    R3BEventHeader* header; /**< Event header. */

    Double_t fSizeX; // Detector size in X and Y
    Double_t fSizeY; // Detector size in X and Y
    Double_t fwx;    // Pad width in X
    Double_t fwy;    // Pad width in Y
    Double_t fx[Mw3PadsX], fy[Mw3PadsY];
    vector<Int_t> fQX;
    vector<Int_t> fQY;
    vector<Int_t> fPadX;
    vector<Int_t> fPadY;
    vector<pair<Int_t, Int_t>> fPairX;
    vector<pair<Int_t, Int_t>> fPairY;
    vector<vector<pair<Int_t, Int_t>>> fClusterX;
    vector<vector<pair<Int_t, Int_t>>> fClusterY;

    double fThresholdX;
    double fThresholdY;
    Bool_t fOnline;          // Don't store data for online
    Bool_t fTofWallMatching; // boolean to know if position matching with the tofwall is required
    Int_t fExpId;

    TClonesArray* fMwpcCalDataCA;    /**< Array with Cal input data. >*/
    TClonesArray* fMwpcHitDataCA;    /**< Array with Hit output data. >*/
    TClonesArray* fTofWallHitDataCA; /**< Array with Hit output data. >*/

    /** Private method AddHitData **/

    // Adds a MwpcHitData to the MwpcHitCollection
    R3BMwpcHitData* AddHitData(Double_t x, Double_t y);
    /** **/
    void TofWallMatching(vector<double>, vector<double>, vector<double>, vector<double>);
    /** **/
    double GetChargeMax(vector<pair<int, int>> p1);
    /** **/
    vector<pair<int, int>> FindCluster(vector<pair<int, int>>& p1);
    /** **/
    void ReconstructHitWithTofWallMatching();
    /** Private method to obtain the position X **/
    Double_t GetPositionX(Double_t qmax, Int_t padmax, Double_t qleft, Double_t qright);
    /** Private method to obtain the position Y **/
    Double_t GetPositionY(Double_t qmax, Int_t padmax, Double_t qdown, Double_t qup);

    Double_t FittedHyperbolicSecant(string XorY, vector<double> vQ, vector<int> vStrip, int Qmax, int StripMax);

  public:
    // Class definition
    ClassDef(R3BMwpc3Cal2Hit, 1)
};

#endif
