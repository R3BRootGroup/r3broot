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
// -----           Created 16/06/22  by I.Syndikus             -----
// -----           Modified Dec 2019 by M. Holl                -----
// -----------------------------------------------------------------

#ifndef R3BPSPXHITPAR_H
#define R3BPSPXHITPAR_H

#include "FairParGenericSet.h"

#include "TArrayF.h"
#include "TObjArray.h"
#include "TObject.h"

class FairParIo;
class FairParamList;

/**
 * Class for Parameters for Cal2Hit Conversion for PSPX detector data.
 * @author Ina Syndikus
 * @since June 22, 2016
 * Modified Dec 2019 by M.Holl
 */

class R3BPspxHitPar : public FairParGenericSet
{

  public:
    /** Standard constructor **/
    R3BPspxHitPar(const char* name = "R3BPspxHitPar",
                  const char* title = "Pspx HIT parameters",
                  const char* context = "Default");

    /** Destructor **/
    virtual ~R3BPspxHitPar();

    // Getter & Setter
    inline const Int_t& GetNumDetectors() const { return fNumDetectors; }
    inline const TArrayF& GetHitPosPar() const { return fHitPosPar; }
    inline const TArrayF& GetHitEPar() const { return fHitEPar; }

    /** Print parameters **/
    virtual void printParams();

    /** Reset all parameters **/
    virtual void clear();

    void putParams(FairParamList*);
    Bool_t getParams(FairParamList*);

  private:
    Int_t fNumDetectors; // number of detectors
    TArrayF fHitPosPar;  // parameters for position
    TArrayF fHitEPar;    // parameters for energy

    R3BPspxHitPar(const R3BPspxHitPar&);
    R3BPspxHitPar& operator=(const R3BPspxHitPar&);

    ClassDef(R3BPspxHitPar, 5);
};

#endif
