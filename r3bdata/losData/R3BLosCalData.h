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

#ifndef R3BLOSCALITEM_H
#define R3BLOSCALITEM_H

#include "TObject.h"

// for the data analysis of the LOS detectors.
// Introduced by Ralf, FEB 2016
// September 2016 - TAMEX readout of MCFD amplitudes added by Aleksandra.
class R3BLosCalData : public TObject
{
  public:
    // Default Constructor
    R3BLosCalData();

    // Standard Constructor
    R3BLosCalData(UInt_t);

    // Destructor
    virtual ~R3BLosCalData() {}

    UInt_t GetDetector() const;
    UInt_t GetVFTXNcha();
    UInt_t GetTAMEXLNcha();
    UInt_t GetTAMEXTNcha();
    UInt_t GetMTDC32Ncha();

    Double_t GetTimeV_ns(UInt_t) const;
    Double_t GetTimeL_ns(UInt_t) const;
    Double_t GetTimeT_ns(UInt_t) const;
    Double_t GetTimeM_ns(UInt_t) const;

    Double_t GetMeanTime();
    Double_t GetMeanTimeTAMEXL();
    Double_t GetMeanTimeTAMEXT();
    Double_t GetMeanTimeVFTX();
    Double_t GetMeanTimeMTDC32();

  private:
    UInt_t fDetector;
    void Reset();

  public:
    Double_t fTimeV_ns[8]; // VFTX
    Double_t fTimeL_ns[8]; // TAMEX leading
    Double_t fTimeT_ns[8]; // TAMEX trailing
    Double_t fTimeM_ns[8]; // MTDC32

    ClassDef(R3BLosCalData, 2)
};

#endif
