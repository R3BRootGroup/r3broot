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

#ifndef R3BNEULANDCAL2HITPAR_H
#define R3BNEULANDCAL2HITPAR_H

#include <memory>

#include "TGraph.h"

#include "FairTask.h"

#include "R3BNeulandCommon.h"

class TClonesArray;

class R3BEventHeader;
class R3BNeulandHitPar;
class R3BNeulandHitModulePar;

namespace Neuland
{
    namespace Calibration
    {
        class HitCalibrationEngine;
        class CosmicTracker;
    } // namespace Calibration
} // namespace Neuland

class R3BNeulandCal2HitPar : public FairTask
{
  public:
    R3BNeulandCal2HitPar()
        : R3BNeulandCal2HitPar("R3BNeulandCal2HitPar", 1){};
    R3BNeulandCal2HitPar(const char* name, const Int_t iVerbose = 1);
    virtual ~R3BNeulandCal2HitPar();

    virtual InitStatus Init();

    virtual void Exec(Option_t* option);

    virtual void FinishTask();

    void SetTpat(UInt_t CosmicTpat = 0) { fCosmicTpat = CosmicTpat; }

    void SavePlots(Bool_t savePlots = true) { fSavePlots = savePlots; }

  private:
    bool IsCosmicEvent() const;

    std::unique_ptr<Neuland::Calibration::HitCalibrationEngine> fHitCalEngine; //!
    std::unique_ptr<Neuland::Calibration::CosmicTracker> fCosmicTracker;       //!

    TClonesArray* fCalNeuland;        //!
    TClonesArray* fMappedLos;         //!
    R3BEventHeader* fEventHeader;     //!
    R3BNeulandHitPar* fNeulandHitPar; //!

    UInt_t fEventNumber = 0;
    UInt_t fAcceptedEventNumber = 0;

    UInt_t fIgorcnt0 = 0;
    UInt_t fIgorcnt1 = 0;
    UInt_t fIgorcnt2 = 0;

    UInt_t fCosmicTpat = 0;

    Bool_t fSavePlots;

  public:
    ClassDef(R3BNeulandCal2HitPar, 1)
};

#endif
