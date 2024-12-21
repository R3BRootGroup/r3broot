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
// ------------------------------------------------------------
// -----           R3BPspxPrecalPar header file           -----
// -----           Created 17/03/13  by I.Syndikus        -----
// -----           Modified Dec 2019 by M. Holl           -----
// ------------------------------------------------------------

#ifndef R3BPSPXPRECALPAR_H
#define R3BPSPXPRECALPAR_H

#include "FairParGenericSet.h"

#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TObject.h"

class FairParIo;
class FairParamList;

/**
 * Class for Parameters for Mapped2Precal Conversion for PSPX detector data.
 * @author Ina Syndikus
 * @since March 13, 2017
 * Modified Dec 2019 by M. Holl
 */

class R3BPspxPrecalPar : public FairParGenericSet
{

  public:
    /** Standard constructor **/
    R3BPspxPrecalPar(const char* name = "R3BPspxPrecalPar",
                     const char* title = "Pspx PRECAL parameters",
                     const char* context = "Default");

    /** Destructor **/
    virtual ~R3BPspxPrecalPar();

    // Getter & Setter
    inline const Int_t& GetNumDetectors() const { return fNumDetectors; }
    inline const TArrayI& GetNumStrips() const { return fNumStrips; }
    inline const TArrayF& GetPrecalPar() const { return fPrecalPar; }

    /** Print parameters **/
    virtual void printParams();

    /** Reset all parameters **/
    virtual void clear();

    void putParams(FairParamList*);
    Bool_t getParams(FairParamList*);

  private:
    Int_t fNumDetectors; // number of detectors
    TArrayI fNumStrips;  // number of strips per detector
    TArrayF fPrecalPar;  // calibration parameters for each strip

    R3BPspxPrecalPar(const R3BPspxPrecalPar&);
    R3BPspxPrecalPar& operator=(const R3BPspxPrecalPar&);

    ClassDef(R3BPspxPrecalPar, 3);
};

#endif
