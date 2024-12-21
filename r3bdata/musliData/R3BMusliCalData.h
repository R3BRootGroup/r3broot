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

// ----------------------------------------------------------------------------
// MUSLI stands for MUSIC LIGHT ION
// ----------------------------------------------------------------------------
//   The Cal Data has the following structure
//     UInt_t fSignal = [1..8], [9..12], [13..14], [15]
//     Double_t fDT; is the Drift Time Tstop - Tstart(=Tref)
//     Double_t fE;  is the energy
// ----------------------------------------------------------------------------

#ifndef R3BMusliCalData_H
#define R3BMusliCalData_H 1

#include "TObject.h"

class R3BMusliCalData : public TObject
{
  public:
    /** Default constructor **/
    R3BMusliCalData();

    R3BMusliCalData(UInt_t s, Double_t dt, Double_t e);

    /** Destructor **/
    virtual ~R3BMusliCalData() {}

    /** Accessors **/
    inline const UInt_t& GetSignal() const { return fSignal; }
    inline const Double_t& GetDT() const { return fDT; }
    inline const Double_t& GetEnergy() const { return fE; }

    /** Modifiers **/
    void SetSignal(UInt_t signal) { fSignal = signal; };
    void SetDT(Double_t dt) { fDT = dt; };
    void SetEnergy(Double_t energy) { fE = energy; };

  protected:
    UInt_t fSignal;
    Double_t fDT;
    Double_t fE;

  public:
    ClassDef(R3BMusliCalData, 1)
};

#endif
