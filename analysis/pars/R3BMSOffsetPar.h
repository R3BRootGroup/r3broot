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

#include "FairParGenericSet.h"

#include "TArrayD.h"
#include <Rtypes.h>
#include <cstdint>
#include <vector>

class FairParamList;

class R3BMSOffsetPar : public FairParGenericSet
{
  public:
    /** Standard constructor **/
    explicit R3BMSOffsetPar(std::string_view name = "MSOffsetPar",
                            std::string_view title = "MS Offset Parameter",
                            std::string_view context = "MSOffsetContext");

    /** Method to reset all parameters **/
    void clear() override;

    /** Method to store all parameters using FairRuntimeDB **/
    void putParams(FairParamList* list) override;

    /** Method to retrieve all parameters using FairRuntimeDB**/
    Bool_t getParams(FairParamList* list) override;

    /** Method to print values of parameters to the standard output **/
    void print() override;
    void printParams() override;

    /** Accessor functions **/
    [[nodiscard]] Double_t GetMSOffset() const { return fMSOffset; }

    void SetMSOffset(Double_t n) { fMSOffset = n; }

  private:
    Double_t fMSOffset; // Master Start Offset from SAMPMSv and SAMPv

  public:
    ClassDefOverride(R3BMSOffsetPar, 1); // NOLINT
};
