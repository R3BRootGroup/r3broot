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
// -----------------------------------------------------------------
// -----           R3BPspxMappedPar header file                -----
// -----           Created 16/05/12  by I.Syndikus             -----
// -----           Modified Dec 2019 by M. Holl                -----
// -----------------------------------------------------------------

#ifndef R3BPSPXCALPAR_H
#define R3BPSPXCALPAR_H

#include "FairParGenericSet.h"

#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TObject.h"

class FairParIo;
class FairParamList;

/**
 * Class for Parameters for Precal2Cal Conversion for PSPX detector data.
 * @author Ina Syndikus
 * @since May 12, 2016
 * Modified Dec 2019 by M.Holl
 */

class R3BPspxCalPar : public FairParGenericSet
{

  public:
    /** Standard constructor **/
    R3BPspxCalPar(const char* name = "R3BPspxCalPar",
                  const char* title = "Pspx CAL parameters",
                  const char* context = "Default");

    /** Destructor **/
    virtual ~R3BPspxCalPar();

    // Getter & Setter
    inline const Int_t& GetNumDetectors() const { return fNumDetectors; }
    inline const TArrayI& GetNumStrips() const { return fNumStrips; }
    inline const TArrayF& GetCalPar() const { return fCalPar; }

    /** Print parameters **/
    virtual void printParams();

    /** Reset all parameters **/
    virtual void clear();

    void putParams(FairParamList*);
    Bool_t getParams(FairParamList*);

  private:
    Int_t fNumDetectors; // number of detectors
    TArrayI fNumStrips;  // number of strips per detector
    TArrayF fCalPar;     // calibration parameters for each strip

    R3BPspxCalPar(const R3BPspxCalPar&);
    R3BPspxCalPar& operator=(const R3BPspxCalPar&);

    ClassDef(R3BPspxCalPar, 5);
};

#endif
