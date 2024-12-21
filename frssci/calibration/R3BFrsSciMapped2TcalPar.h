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

#pragma once

#include "FairTask.h"
#include "TH1F.h"

class TClonesArray;
class R3BFrsSciTcalPar;

class R3BFrsSciMapped2TcalPar : public FairTask
{
  public:
    /** Default constructor **/
    R3BFrsSciMapped2TcalPar();

    /** Standard constructor **/
    R3BFrsSciMapped2TcalPar(const char* name, Int_t iVerbose = 1);

    /** Destructor **/
    virtual ~R3BFrsSciMapped2TcalPar();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method FinishEvent **/
    virtual void FinishEvent();

    /** Virtual method FinishTask **/
    virtual void FinishTask();

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method ReInit **/
    virtual InitStatus ReInit();

    /** Virtual method calculate the Vftx Tcal Parameters **/
    virtual void CalculateVftxTcalParams();

    void SetOutputFile(const char* outFile);

    const UShort_t GetNumDets() { return fNumDets; }
    const UShort_t GetNumPmts() { return fNumPmts; }
    const UShort_t GetNumPars() { return fNumPars; }
    const UInt_t GetMinStat() { return fMinStat; }

    void SetNumDets(UShort_t nDets) { fNumDets = nDets; }
    void SetNumPmts(UShort_t nPmts) { fNumPmts = nPmts; }
    void SetNumPars(UShort_t nPars) { fNumPars = nPars; }
    void SetMinStat(UInt_t minstat) { fMinStat = minstat; }

  private:
    Int_t fMinStatistics; // Minimum statistics to proceed to the calibration

    TClonesArray* fMapped; // Array with input mapped data
    UShort_t fNumDets;
    UShort_t fNumPmts;
    UShort_t fNumPars;
    UInt_t fMinStat;
    R3BFrsSciTcalPar* fTcalPar; // Tcal Parameters

    // histograms
    TH1F** fh_TimeFineBin;
    TH1D** fh_TimeFineNs;
    char* fOutputFile;

  public:
    ClassDef(R3BFrsSciMapped2TcalPar, 0);
};
