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

// -------------------------------------------------------------
// -----         R3BTwimGainMatching source file           -----
// -----    Created 16/02/22 by Antia Grana Gonzalez       -----
// -------------------------------------------------------------

#ifndef R3BTwimGainMatching_H
#define R3BTwimGainMatching_H 1

#include "FairTask.h"
#include <Rtypes.h>

class TClonesArray;
class R3BTwimCalPar;
class TCanvas;
class TGraph;

class R3BTwimGainMatching : public FairTask
{
  public:
    /** Default constructor **/
    R3BTwimGainMatching();

    /** Standard constructor **/
    R3BTwimGainMatching(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BTwimGainMatching();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* option);

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    // Method to set up the limits for fit
    void SetFitLimits(Int_t left, Int_t right)
    {
        fLimit_left = left;
        fLimit_right = right;
    }

    // Method to set up the reference anode
    // a is given in 1-base
    void SetRefAnode(Int_t a) { fRefAnode = a - 1; }

  private:
    Int_t fNumSec;
    Int_t fNumAnodes;
    Int_t fMinStadistics;
    Int_t fNumParams;
    Int_t fNumAnodesRef;
    Int_t fRefAnode;
    Int_t fLimit_left;
    Int_t fLimit_right;
    TArrayF* CalParams;
    TGraph** Anode8_vs_anodes;
    TCanvas** canvas;

    R3BTwimCalPar* fCal_Par;         /**< Parameter container. >*/
    TClonesArray* fTwimMappedDataCA; /**< Array with Twim Mapped-input data. >*/

  public:
    // Class definition
    ClassDef(R3BTwimGainMatching, 1)
};

#endif /* R3BTwimGainMatching_H */
