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

// -----------------------------------------------------------------------
// -----                       R3BMwpc1 header file                  -----
// -----                 Created 06/10/19 by JL Rodriguez            -----
// -----                  s455 method 17/06/22 by Antia GG           -----
// -----------------------------------------------------------------------

#ifndef R3BMwpc1Cal2Hit_H
#define R3BMwpc1Cal2Hit_H 1

#include "FairTask.h"
#include "R3BMwpcCalData.h"
#include "R3BMwpcHitData.h"
#include "TH1F.h"
#include <TRandom.h>

#define Mw1PadsX 64
#define Mw1PadsY 40

// Mw1PadsX 64

using namespace std;

class TClonesArray;
class R3BEventHeader;

class R3BMwpc1Cal2Hit : public FairTask
{

  public:
    /** Default constructor **/
    R3BMwpc1Cal2Hit();

    /** Standard constructor **/
    R3BMwpc1Cal2Hit(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BMwpc1Cal2Hit();

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

    Double_t fSize; // Detector size in X and Y
    Double_t fwx;   // Pad width in X
    Double_t fwy;   // Pad width in Y
    Int_t fx[Mw1PadsX], fx_p1[Mw1PadsX], fx_p2[Mw1PadsX], fy[Mw1PadsY];

    Bool_t fOnline; // Don't store data for online
    Int_t fExpId;

    TClonesArray* fMwpcCalDataCA; /**< Array with Cal input data. >*/
    TClonesArray* fMwpcHitDataCA; /**< Array with Hit output data. >*/

    /** Private method AddHitData **/
    // Adds a MwpcHitData to the MwpcHitCollection
    R3BMwpcHitData* AddHitData(Double_t x, Double_t y, Int_t plane = 1);

    /** Private method to obtain the position X **/
    Double_t GetPositionX(Double_t qmax, Int_t padmax, Double_t qleft, Double_t qright);
    /** Private method to obtain the position Y **/
    Double_t GetPositionY(Double_t qmax, Int_t padmax, Double_t qdown, Double_t qup);
    /** Private method to obtain the position Y **/
    [[nodiscard]] Double_t GetPositionYCoG(Double_t qmax, Int_t padmax, Double_t qdown, Double_t qup) const;

  public:
    // Class definition
    ClassDef(R3BMwpc1Cal2Hit, 1)
};

#endif
