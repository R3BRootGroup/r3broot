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

#ifndef R3BPSPXMAPPEDDATA_H
#define R3BPSPXMAPPEDDATA_H

#include "TObject.h"

/**
 * Class containing PSPX detector data on Mapped level.
 * @author Ralf Plag, Ina Syndikus
 * @since January 2016
 * Modified Dec 2019 by M. Holl
 */

class R3BPspxMappedData : public TObject
{
  public:
    /** Default Constructor **/
    R3BPspxMappedData();

    /** Standard Constructor **/
    R3BPspxMappedData(Int_t strip1, Int_t strip2, Int_t energy1, Int_t energy2);

    /** Destructor **/
    virtual ~R3BPspxMappedData() {}

    void SetValue(Int_t, Int_t, Int_t);
    // Getters
    inline const Int_t& GetStrip1() const { return fStrip[0]; }
    inline const Int_t& GetStrip2() const { return fStrip[1]; }
    inline const Int_t& GetEnergy1() const { return fEnergy[0]; }
    inline const Int_t& GetEnergy2() const { return fEnergy[1]; }

  private:
    Int_t fStrip[2];  // Channel number, counting from 1, one entry for each side of each strip
    Int_t fEnergy[2]; // Energy/Collected charge, one entry for each side of each strip

  public:
    ClassDef(R3BPspxMappedData, 5)
};

#endif
