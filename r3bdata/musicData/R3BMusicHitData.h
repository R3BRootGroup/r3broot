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
// -----                      R3BMusicHitData header file              -----
// -----                  Created 15/10/19 by JL Rodriguez             -----
// -------------------------------------------------------------------------

#ifndef R3BMusicHitData_H
#define R3BMusicHitData_H 1

#include "TObject.h"

class R3BMusicHitData : public TObject
{
  public:
    /** Default constructor **/
    R3BMusicHitData();

    /** Constructor with arguments
     *@param theta    Angle theta in [mrad]
     *@param z        Atomic number Z in charge units
     **/
    R3BMusicHitData(Double_t theta, Double_t z);
    R3BMusicHitData(Double_t theta, Double_t z, Double_t ene);
    R3BMusicHitData(Double_t theta, Double_t z, Double_t ene, Double_t good_dt);

    /** Destructor **/
    virtual ~R3BMusicHitData() {}

    /** Accessors **/
    Double_t GetTheta() const { return fTheta; }
    Double_t GetZcharge() const { return fZ; }
    Double_t GetEave() const { return fE; }
    Double_t GetGoodDt() const { return fDt; }

    /** Modifiers **/
    void SetTheta(Double_t theta) { fTheta = theta; };
    void SetZcharge(Double_t z) { fZ = z; };
    void SetEave(Double_t ene) { fE = ene; };
    void SetDt(Double_t dt) { fDt = dt; };

  private:
    Double_t fTheta, fZ, fE, fDt;

    ClassDef(R3BMusicHitData, 1)
};

#endif
