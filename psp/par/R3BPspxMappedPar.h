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

#ifndef R3BPSPXMAPPEDPAR_H
#define R3BPSPXMAPPEDPAR_H

#include "FairParGenericSet.h"

#include "TArrayF.h"
#include "TArrayI.h"
#include "TObjArray.h"
#include "TObject.h"

class FairParIo;
class FairParamList;

/**
 * Class for Parameters for unpacking with PspxReader for PSPX detector data.
 * @author Ina Syndikus
 * @since May 12, 2016
 * Modified Dec 2019 by M. Holl
 */

class R3BPspxMappedPar : public FairParGenericSet
{

  public:
    /** Standard constructor **/
    R3BPspxMappedPar(const char* name = "R3BPspxMappedPar",
                     const char* title = "Pspx MAPPED parameters",
                     const char* context = "Default");

    /** Destructor **/
    virtual ~R3BPspxMappedPar();

    // Getter & Setter
    inline const Int_t& GetPspxParDetector() const { return pspxmappedpardetector; }
    inline const TArrayI& GetPspxParStrip() const { return pspxmappedparstrip; }

    /** Print parameters **/
    virtual void printparams();

    /** Reset all parameters **/
    virtual void clear();

    void putParams(FairParamList*);
    Bool_t getParams(FairParamList*);

  private:
    TArrayI pspxmappedparstrip;  // Number of strips per detector
    Int_t pspxmappedpardetector; // Number of detectors

    R3BPspxMappedPar(const R3BPspxMappedPar&);
    R3BPspxMappedPar& operator=(const R3BPspxMappedPar&);

    ClassDef(R3BPspxMappedPar, 1);
};

#endif
