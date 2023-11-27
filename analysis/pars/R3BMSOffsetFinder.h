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
#include "TH1F.h"

class TClonesArray;
class R3BMSOffsetPar;
class R3BEventHeader;
class R3BSamplerMappedData;
class R3BMSOffsetFinder : public FairTask // NOLINT
{

  public:
    /** Default constructor **/
    R3BMSOffsetFinder();

    R3BMSOffsetFinder(const R3BMSOffsetFinder&) = delete;
    R3BMSOffsetFinder(R3BMSOffsetFinder&&) = delete;
    R3BMSOffsetFinder& operator=(const R3BMSOffsetFinder&) = delete;
    R3BMSOffsetFinder& operator=(R3BMSOffsetFinder&&) = delete;
    /** Standard constructor **/
    explicit R3BMSOffsetFinder(const TString& name, Int_t iVerbose = 1);

  private:
    /** Virtual method Init **/
    InitStatus Init() override;

    /** Virtual method Exec **/
    void Exec(Option_t* opt) override;

    /** Virtual method FinishEvent **/
    void FinishEvent() override;

    /** Virtual method FinishTask **/
    void FinishTask() override;

    /** Virtual method Reset **/
    virtual void Reset();

    /** Virtual method ReInit **/
    InitStatus ReInit() override;

    /** Virtual method Search MS Offset **/
    virtual void SearchMSOffset();

    /** Virtual method SetParContainers **/
    void SetParContainers() override;

  public:
    /** Accessor functions **/
    [[nodiscard]] Double_t GetMSOffset() const { return fMSOffset; }

    void SetMSOffset(Double_t MSOffset) { fMSOffset = MSOffset; }

    void SetMinStatistics(Int_t MinStatistics) { fMinStatistics = MinStatistics; }

  private:
    // Number of histograms, limits and bining
    Double_t fMSOffset = 0.0;

    // Minimum statistics and parameters
    Int_t fMinStatistics = 1;

    R3BMSOffsetPar* fMSOffsetPar = nullptr;   /**< Parameter container. >*/
    TClonesArray* fSamplerMapped = nullptr;   /**< Array with SAMP Mapped input data. >*/
    TClonesArray* fSamplerMSMapped = nullptr; /**< Array with SAMPMS MApped input data. >*/

    TH1F* fh_Offset_Finder = nullptr;

  public:
    ClassDefOverride(R3BMSOffsetFinder, 1); // NOLINT
};
