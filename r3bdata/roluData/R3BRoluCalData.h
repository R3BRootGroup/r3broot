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

#ifndef R3BROLUCALITEM_H
#define R3BROLUCALITEM_H

#include "TObject.h"

// for the data analysis of the ROLU following LOS.

class R3BRoluCalData : public TObject
{
  public:
    // Default Constructor
    R3BRoluCalData();

    // Standard Constructor
    R3BRoluCalData(UInt_t);

    // Destructor
    virtual ~R3BRoluCalData() {}

    UInt_t GetDetector() const;
    UInt_t GetTAMEXLNcha();
    UInt_t GetTAMEXTNcha();

    Double_t GetTimeL_ns(UInt_t) const;
    Double_t GetTimeT_ns(UInt_t) const;

    Double_t GetMeanTime();
    Double_t GetMeanTimeTAMEXL();
    Double_t GetMeanTimeTAMEXT();

  private:
    UInt_t fDetector;
    void Reset();

  public:
    Double_t fTimeL_ns[4]; // TAMEX leading
    Double_t fTimeT_ns[4]; // TAMEX trailing

    ClassDef(R3BRoluCalData, 1)
};

#endif // R3BROLUCALITEM_H
