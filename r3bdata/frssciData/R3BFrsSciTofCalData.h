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

#include "TObject.h"

class R3BFrsSciTofCalData : public TObject
{
  public:
    // Default Constructor
    R3BFrsSciTofCalData();

    // Standard Constructor
    R3BFrsSciTofCalData(UShort_t,
                        UShort_t,
                        UShort_t,
                        Float_t,
                        Float_t,
                        Double_t,
                        Double_t,
                        Double_t,
                        Double_t,
                        Double_t);

    // Destructor
    virtual ~R3BFrsSciTofCalData() {}

    // Getters
    inline const UShort_t& GetRank() const { return fRank; }
    inline const UShort_t& GetDetIdSta() const { return fDetIdSta; }
    inline const UShort_t& GetDetIdSto() const { return fDetIdSto; }
    inline const Float_t& GetCalPosStaMm() const { return fCalPosStaMm; }
    inline const Float_t& GetCalPosStoMm() const { return fCalPosStoMm; }
    inline const Double_t& GetRawTofNs() const { return fRawTofNs; }
    inline const Double_t& GetCalTofNs() const { return fCalTofNs; }
    inline const Double_t& GetBeta() const { return fBeta; }
    inline const Double_t& GetBRho() const { return fBRho; }
    inline const Double_t& GetAoQ() const { return fAoQ; }

  private:
    UShort_t fRank;
    UShort_t fDetIdSta;
    UShort_t fDetIdSto;
    Float_t fCalPosStaMm;
    Float_t fCalPosStoMm;
    Double_t fRawTofNs;
    Double_t fCalTofNs;
    Double_t fBeta;
    Double_t fBRho;
    Double_t fAoQ;

  public:
    ClassDef(R3BFrsSciTofCalData, 2)
};
