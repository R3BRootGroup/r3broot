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

#pragma once

#include "FairTask.h"
#include "R3BMSOffsetPar.h"
class TClonesArray;
class R3BMSOffsetPar;
class R3BEventHeader;
class R3BSamplerMappedData;
class R3BTprevTnext : public FairTask
{
  public:
    /** Default constructor **/
    R3BTprevTnext();

    /** Standard constructor **/
    explicit R3BTprevTnext(const TString& name, Int_t iVerbose = 1);

    /** Virtual method Init **/
    InitStatus Init() override;

    /** Virtual method Exec **/
    void Exec(Option_t* opt) override;

    void SetNoRecordTimeGap(Double_t No_Record_Time_Gap)
    {
        fNo_record_time = No_Record_Time_Gap;
    } // Setter for a value of Tprev/Tnext if there is no recorded hit before and/or after the MS
    void SetDelta_clk(Double_t clock) { fDelta_clk = clock; }

    static constexpr auto CLOCK_PERIOD = 10;    // ns
    static constexpr auto ERROR_NO_MS = -30;    // Assigned value for events where there is no MS
    static constexpr auto INVALID_TPTN = -10;   // Assigned value for invalid or nonexistent TPrev or TNext
    static constexpr auto ERROR_MULTI_MS = -20; // Assigned value for events where there are multiple MS
    static constexpr auto INVALID_EVENT = -40;  // Assigned value for events where the MS was not correctly written
  private:
    R3BMSOffsetPar* fMSOffsetPar = nullptr;
    TClonesArray* fSamplerMapped = nullptr;
    TClonesArray* fSamplerMSMapped = nullptr;
    R3BEventHeader* fR3BEventHeader = nullptr; /**< Event header - input data. */
    Double_t fDelta_clk = 1.0;
    Double_t fNo_record_time = 2e3;

    ClassDefOverride(R3BTprevTnext, 1); // NOLINT
};
