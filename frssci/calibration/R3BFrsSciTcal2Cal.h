/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// *** *************************************************************** *** //
// ***                  R3BFrsSciTcal2Cal                              *** //
// ***    convert tcal data to cal data : PosCal and TofCal            *** //
// ***     1/ select the proper hit in the multihit tcal data          *** //
// ***     2/ calculate PosRaw and apply calibration parameters        *** //
// ***     2/ calculate TofRaw and apply calibration parameters        *** //
// *** *************************************************************** *** //

#pragma once

#include "FairTask.h"
#include "R3BFrsSciCalPar.h"
#include "R3BFrsSciPosCalData.h"
#include "R3BFrsSciTcalData.h"
#include "R3BFrsSciTofCalData.h"
#include "TClonesArray.h"
#include "TRandom3.h"

class TRandom3;
class R3BFrsSciTcalData;

class R3BFrsSciTcal2Cal : public FairTask
{
  public:
    // --- Default constructor --- //
    R3BFrsSciTcal2Cal();

    // --- Standard constructor --- //
    R3BFrsSciTcal2Cal(const char* name, Int_t iVerbose = 1);

    // --- Destructor --- //
    virtual ~R3BFrsSciTcal2Cal();

    virtual void Exec(Option_t* option);

    virtual void SetParContainers();

    virtual InitStatus Init();

    virtual InitStatus ReInit();

    virtual void Reset();

    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    UInt_t fNevent;
    TClonesArray* fTcal;      // input data - FrsSci
    R3BFrsSciCalPar* fCalPar; // cal parameters container - FrsSci
    TClonesArray* fPosCal;    // output data for position per detector
    TClonesArray* fTofCal;    // output data for Tof
    Bool_t fOnline;           // Don't store data for online

    TRandom3 rand;

    R3BFrsSciPosCalData* AddPosCalData(UShort_t det,
                                       Double_t traw,
                                       Double_t traw_wtref,
                                       Float_t rawpos,
                                       Float_t calpos);
    R3BFrsSciTofCalData* AddTofCalData(UShort_t rank,
                                       UShort_t detsta,
                                       UShort_t detsto,
                                       Float_t calpossta,
                                       Float_t calpossto,
                                       Double_t rawtof,
                                       Double_t caltof,
                                       Double_t beta,
                                       Double_t gamma,
                                       Double_t aoq);

  public:
    ClassDef(R3BFrsSciTcal2Cal, 1)
};
