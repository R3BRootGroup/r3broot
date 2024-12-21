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
// -----                      R3BMusicCalData header file              -----
// -----                  Created 15/10/19 by JL Rodriguez             -----
// -------------------------------------------------------------------------

#ifndef R3BMusicCalData_H
#define R3BMusicCalData_H 1

#include "TObject.h"

class R3BMusicCalData : public TObject
{
  public:
    /** Default constructor **/
    R3BMusicCalData();

    /** Constructor with arguments
     *@param anodeID  Anode ID
     *@param dtime    Drift time calibrated in [mm] for the position
     *@param energy   Energy deposit [channels], for simulation in [MeV]
     **/
    R3BMusicCalData(UShort_t anodeID, Double_t dtime, Double_t energy);

    /** Destructor **/
    virtual ~R3BMusicCalData() {}

    /** Accessors **/
    UShort_t GetAnodeID() const { return fAnodeID; }
    Double_t GetDTime() const { return fDT; }
    Double_t GetEnergy() const { return fEnergy; }

    /** Modifiers **/
    void SetAnodeID(UShort_t id) { fAnodeID = id; };
    void SetDTime(Double_t dt) { fDT = dt; };
    void SetEnergy(Double_t energy) { fEnergy = energy; };

  private:
    UShort_t fAnodeID;
    Double_t fDT, fEnergy;

    ClassDef(R3BMusicCalData, 1)
};

#endif
