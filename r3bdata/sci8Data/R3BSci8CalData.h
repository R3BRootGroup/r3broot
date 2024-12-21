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

#ifndef R3BSCI8CALITEM_H
#define R3BSCI8CALITEM_H

#include "TObject.h"

// for the data analysis of the Sci8 detectors.

class R3BSci8CalData : public TObject
{
  public:
    // Default Constructor
    R3BSci8CalData();

    // Standard Constructor
    R3BSci8CalData(UInt_t detector);

    // Destructor
    virtual ~R3BSci8CalData() {}

    // Getters
    inline const UInt_t& GetDetector() const { return fDetector; }

    Double_t GetMeanTime();
    Double_t GetMeanTimeTAMEXL();
    Double_t GetMeanTimeTAMEXT();
    Double_t GetMeanTimeVFTX();
    UInt_t GetVFTXNcha();
    UInt_t GetTAMEXLNcha();
    UInt_t GetTAMEXTNcha();

  private:
    UInt_t fDetector;

  public:
    Double_t fTimeV_r_ns; // right VFTX
    Double_t fTimeV_l_ns; // left VFTX
    Double_t fTimeL_r_ns; // right TAMEX leading
    Double_t fTimeL_l_ns; // left TAMEX leading
    Double_t fTimeT_r_ns; // right  TAMEX trailing
    Double_t fTimeT_l_ns; // left  TAMEX trailing

    ClassDef(R3BSci8CalData, 1)
};

#endif
