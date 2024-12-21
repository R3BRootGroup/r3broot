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

#ifndef R3BAmsDigitizer_H
#define R3BAmsDigitizer_H 1

#include "FairTask.h"
#include "R3BAmsStripCalData.h"
#include "R3BTraPoint.h"
#include "Rtypes.h"

class TClonesArray;

class R3BAmsDigitizer : public FairTask
{
  public:
    /** Standard contructor **/
    R3BAmsDigitizer();

    /** Destructor **/
    virtual ~R3BAmsDigitizer();

    /** Virtual method Init **/
    virtual InitStatus Init();

    /** Virtual method Exec **/
    virtual void Exec(Option_t* opt);

    /** Virtual method Reset **/
    virtual void Reset();

  private:
    TClonesArray* fPointData;
    TClonesArray* fCalData;

    /** Private method AddCal
     **
     ** Adds a R3BAmsStripCalData data
     **/
    R3BAmsStripCalData* AddCal(Int_t detid, Int_t sideid, Int_t stripid, Double_t energy);

  public:
    ClassDef(R3BAmsDigitizer, 1);
};

#endif /* R3BAmsDigitizer_H */
