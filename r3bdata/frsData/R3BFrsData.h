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

// ---------------------------------------------------------------
// -----                  R3BFrsData                         -----
// -----     Created 30/04/2021 by J.L. Rodriguez            -----
// ---------------------------------------------------------------

#ifndef R3BFrsData_H
#define R3BFrsData_H
#include "TObject.h"

class R3BFrsData : public TObject
{

  public:
    // Default Constructor
    R3BFrsData();

    /** Standard Constructor
     *@param fStaId  ID of Start Sci starting from 1
     *@param fStoId  ID of Stop Sci starting from 1
     *@param fZ      Z of fragments
     *@param fAq     A/q of fragments
     *@param fBeta   Beta of fragments
     *@param fBrho   Brho of fragments
     *@param cs2     Position at S2
     *@param xc      Position at Cave
     **/
    R3BFrsData(Int_t StaId,
               Int_t StoId,
               Double_t z,
               Double_t aq,
               Double_t beta,
               Double_t brho,
               Double_t xs2 = 0.,
               Double_t xc = 0.,
               Double_t tof = 0.);

    // Destructor
    virtual ~R3BFrsData() {}

    // Getters
    inline const Double_t GetStaId() const { return fStaId; }
    inline const Double_t GetStoId() const { return fStoId; }
    inline const Double_t GetZ() const { return fZ; }
    inline const Double_t GetAq() const { return fAq; }
    inline const Double_t GetBeta() const { return fBeta; }
    inline const Double_t GetTof() const { return fTof; }
    inline const Double_t GetBrho() const { return fBrho; }
    inline const Double_t GetXS2() const { return fxs2; }
    inline const Double_t GetXCave() const { return fxc; }

    void SetStaId(Int_t v) { fStaId = v; }
    void SetStoId(Int_t v) { fStoId = v; }
    void SetZ(Double_t v) { fZ = v; }
    void SetAq(Double_t v) { fAq = v; }
    void SetBeta(Double_t v) { fBeta = v; }
    void SetBrho(Double_t v) { fBrho = v; }
    void SetXS2(Double_t v) { fxs2 = v; }
    void SetXCave(Double_t v) { fxc = v; }

  protected:
    Int_t fStaId, fStoId;
    Double_t fZ, fAq;
    Double_t fBeta, fTof, fBrho, fxs2, fxc;

  public:
    ClassDef(R3BFrsData, 1)
};

#endif
