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

// --------------------------------------------------------------
// -----             R3BFootMapped2StripCal                 -----
// -----     Created 08/09/21 by J.L. Rodriguez-Sanchez     -----
// --------------------------------------------------------------

#pragma once

#include "R3BFootCalData.h"

#include <FairTask.h>

#include <Rtypes.h>
#include <TArrayF.h>

class TClonesArray;
class R3BFootCalPar;

class R3BFootMapped2StripCal : public FairTask
{
  public:
    /** Default constructor **/
    R3BFootMapped2StripCal();

    /** Standard constructor **/
    R3BFootMapped2StripCal(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BFootMapped2StripCal();

    /** Virtual method Exec **/
    void Exec(Option_t* /*option*/) override;

    /** Virtual method Reset **/
    virtual void Reset();

    void SetParContainers() override;

    // Fair specific
    /** Virtual method Init **/
    InitStatus Init() override;

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    // Method to setup online mode
    void SetOnline(bool option) { fOnline = option; }

    /**
     * Method for setting the thresholds: Signal>sigma_strip*fTimesSigma
     */
    void SetThresholdSigma(Double_t th) { fTimesSigma = th; }

    /**
     * Method for setting the N of strip hit per foot: fNStrip
     */
    void SetNStrip(Double_t nstrip) { fNStrip = nstrip; }

  private:
    void SetParameter();

    int fNDets = 16;
    int fNStrip = 640;
    int MaxSigma = 5;
    int fNumParams = 2;
    double fTimesSigma = 2.;
    TArrayF* CalParams = nullptr;

    bool fOnline = false; // Don't store data for online

    R3BFootCalPar* fCal_Par = nullptr;       // Parameter container
    TClonesArray* fFootMappedData = nullptr; // Array with FOOT Mapped input data
    TClonesArray* fFootCalData = nullptr;    // Array with FOOT Cal output data

    // Private method AddCalData
    R3BFootCalData* AddCalData(int8_t detid, int16_t stripid, double energy);

  public:
    // Class definition
    ClassDefOverride(R3BFootMapped2StripCal, 1)
};
