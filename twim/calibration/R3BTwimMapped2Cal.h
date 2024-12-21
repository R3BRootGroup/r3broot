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
// -----                 R3BTwimMapped2Cal                     -----
// -----        Created 24/11/19 by J.L. Rodriguez-Sanchez     -----
// -----------------------------------------------------------------

#pragma once

#include "FairTask.h"
#include "R3BLogger.h"
#include "R3BTwimCalData.h"
#include "R3BTwimMappedData.h"
#include "TH1F.h"
#include <TRandom.h>
#include <stdint.h>
#include <vector>

constexpr Int_t MAXNUMSEC = 4;
constexpr Int_t MAXNUMANODE = 16;
constexpr Int_t MAXNUMTREF = 2;
constexpr Int_t MAXNUMTRIG = 2;
constexpr Int_t MAXNUMCH = MAXNUMANODE + MAXNUMTREF + MAXNUMTRIG;
constexpr Int_t MAXMULT = 20;

class TClonesArray;
class R3BEventHeader;
class R3BTwimCalPar;

class R3BTwimMapped2Cal : public FairTask
{
  public:
    /** Default constructor **/
    R3BTwimMapped2Cal();

    /** Standard constructor **/
    R3BTwimMapped2Cal(const TString& name, Int_t iVerbose = 1);

    /** Destructor **/
    ~R3BTwimMapped2Cal() override = default;

    /** Virtual method Exec **/
    void Exec(Option_t*) override;

    void SetParContainers() override;

    /** Method Init **/
    InitStatus Init() override;

    /** Method ReInit **/
    InitStatus ReInit() override;

    /** Method to set up the experiment ID **/
    void SetExpId(Int_t exp)
    {
        R3BLOG(info, "fExpId is set locally. Original:" << fExpId << ", New value:" << exp);
        R3BLOG(info, "Using R3BEventHeader::SetExpId() is recommended instead.");
        fExpId = exp;
    }

    void SetOnline(Bool_t option) { fOnline = option; }

  private:
    /** Method to set parameters **/
    void SetParameter();

    /** Method Reset **/
    void Reset();

    R3BEventHeader* header;

    Int_t fNumSec;
    Int_t fNumAnodes;
    Int_t fNumEParams;
    Int_t fNumPosParams;
    Int_t fNumAnodesRef;
    Int_t fNumAnodesTrig;
    Int_t fMaxMult;
    Int_t fMinDT;
    Int_t fMaxDT;
    Int_t fExpId;
    std::vector<TArrayF*> CalEParams;
    std::vector<TArrayF*> PosParams;

  public:
    class CalAnode
    {
      public:
        CalAnode() = default;
        void Init();
        void AddMult() { fmult++; };
        void SetE(Double_t val) { fE.at(fmult) = val; };
        void SetDt(Int_t val) { fDT.at(fmult) = val; };
        void SetVal(Double_t ene, Int_t dtime)
        {
            fE.at(fmult) = ene;
            fDT.at(fmult) = dtime;
            fmult++;
        };
        [[nodiscard]] auto GetMult() const { return fmult; };
        auto GetE(Int_t mult) { return fE.at(mult); };
        auto GetDT(Int_t mult) { return fDT.at(mult); };

      private:
        Int_t fmult{ 0 };
        std::vector<Double_t> fE;
        std::vector<Int_t> fDT;
    };

    class CalSection
    {
      public:
        CalSection();
        void Init();
        auto& GetAnode(Int_t val) { return fAnode.at(val); };
        auto& GetTrig(Int_t val = 0) { return fTrig.at(val); };
        auto& GetTref(Int_t val = 0) { return fTref.at(val); };

      private:
        std::vector<CalAnode> fAnode;
        std::vector<CalAnode> fTref;
        std::vector<CalAnode> fTrig;
    };

  private:
    std::vector<CalSection> fTwimCal;

    Bool_t fOnline;

    R3BTwimCalPar* fCal_Par;         /**< Parameter container. >*/
    TClonesArray* fTwimMappedDataCA; /**< Array with Mapped-input data. >*/
    TClonesArray* fTwimCalDataCA;    /**< Array with Cal-output data. >*/

    /** Private method AddCalData **/
    R3BTwimCalData* AddCalData(Int_t secID, Int_t anodeID, Double_t dtime, Double_t energy);

  public:
    // Class definition
    ClassDefOverride(R3BTwimMapped2Cal, 1)
};
