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

#ifndef R3BNEULANDHITCALIBRATIONENGINE_H
#define R3BNEULANDHITCALIBRATIONENGINE_H

#include <vector>

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

#include "R3BNeulandCosmicTrack.h"
#include "R3BNeulandHitCalibrationBar.h"
#include "R3BNeulandHitModulePar.h"
#include "R3BNeulandTSyncer.h"

class R3BNeulandHitPar;
class TDirectory;

namespace Neuland
{
    namespace Calibration
    {
        class HitCalibrationEngine
        {
          public:
            HitCalibrationEngine();

            void Init(const R3BNeulandHitPar* hitpar);
            void Set(const Int_t id, const Int_t side, const Double_t time, const Int_t qdc);
            inline Bool_t IsValid(const Int_t id) { return fBars[id].IsValid(); }
            inline Double_t GetPosition(const Int_t id) { return fBars[id].GetPosition(); }
            void Add(const R3BNeulandCosmicTrack& track, const UInt_t eventNumber);
            void Reset();
            std::vector<R3BNeulandHitModulePar> Calibrate(TDirectory* histoDir = nullptr);

          private:
            void draw() const;

            TSyncer fTSyncer;
            std::vector<HitCalibrationBar> fBars;
            std::vector<ULong64_t> fHitMask;

            TH1F fBarDistribution;
            TH1F fStoppedDistribution;
            TH1F fInteractionsDistribution;
            TH1F fStoppedInteractionsDistribution;
            TH1F fTrackLengthDistribution;
            TH1F fTotalTrackLengthDistribution;
            TH1F fTotalStoppedTrackLengthDistribution;
            TH2F fCorrelationMatrix;
            TH3F fTrackEntryPointDistribution; // Z X Y
            TH3F fTrackDirectionDistribution;  // Z X Y
        };
    } // namespace Calibration
} // namespace Neuland
#endif
