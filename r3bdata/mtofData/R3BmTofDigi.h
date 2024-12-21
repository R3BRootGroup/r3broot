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

// -------------------------------------------------------------------------
// -----                      R3BmTofDigi header file                  -----
// -----                  Created 28/06/11  by D.Bertini/Justyna               -----
// -------------------------------------------------------------------------

/**  R3BmTofDigi.h
 **/

#ifndef R3BMTOFDIGI_H
#define R3BMTOFDIGI_H

#include "TObject.h"
#include "TVector3.h"

class R3BmTofDigi : public TObject
{

  public:
    /** Default constructor **/
    R3BmTofDigi();
    R3BmTofDigi(Int_t ntmul,
                Double_t ntfx,
                Double_t ntfy,
                Double_t ntft,
                Double_t ntfpath,
                Double_t ntfpx,
                // Double_t ntfpy,Double_t ntfpz);
                Double_t ntfpy,
                Double_t ntfpz,
                Double_t ntfe);

    /** Destructor **/
    virtual ~R3BmTofDigi();

    /** Output to screen **/
    virtual void Print(const Option_t* opt = "") const;

    void Setntmul(Int_t mul) { Ntmul = mul; }
    Double_t Getntmul() { return Ntmul; }

    void Setntfx(Double_t pos) { Ntx = pos; }
    Double_t Getntfx() { return Ntx; }

    void Setntfy(Double_t pos) { Nty = pos; }
    Double_t Getntfy() { return Nty; }

    void Setntft(Double_t time) { Ntt = time; }
    Double_t Getntft() { return Ntt; }

    void Setntfpath(Double_t length) { Ntpath = length; }
    Double_t Getntfpath() { return Ntpath; }

    void Setntfpx(Int_t mom) { Ntfpx = mom; }
    Double_t Getntfpx() { return Ntfpx; }

    void Setntfpy(Int_t mom) { Ntfpy = mom; }
    Double_t Getntfpy() { return Ntfpy; }

    void Setntfpz(Int_t mom) { Ntfpz = mom; }
    Double_t Getntfpz() { return Ntfpz; }

    void Setntfe(Double_t en) { Nte = en; }
    Double_t Getntfe() { return Nte; }

  protected:
    Int_t Ntmul;
    Double32_t Ntx;
    Double32_t Nty;
    Double32_t Ntt;
    Double32_t Ntpath;
    Double32_t Ntfpx;
    Double32_t Ntfpy;
    Double32_t Ntfpz;
    Double32_t Nte;

    ClassDef(R3BmTofDigi, 1)
};

#endif
