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

#ifndef R3BNEULANDCOSMICTRACKER_H
#define R3BNEULANDCOSMICTRACKER_H

#include <array>
#include <vector>

#include "TF1.h"
#include "TGraphErrors.h"
#include "TVector3.h"

#include "R3BNeulandCommon.h"
#include "R3BNeulandCosmicTrack.h"

namespace Neuland
{
    namespace Calibration
    {
        class CosmicTracker
        {
            using DPair = std::array<Double_t, 2>;

          public:
            CosmicTracker();

            void SetDistances(const std::vector<Double_t>& distances) { fDistances = distances; }
            void AddPoint(const Int_t barID, const Double_t pos = NaN);
            const R3BNeulandCosmicTrack& GetTrack();
            void Reset();

          private:
            void filter(TGraphErrors& graph) const;
            DPair fit(TGraphErrors& graph);
            DPair linearRegression(const Double_t* x, const Double_t* y, const Int_t points) const;
            void fillInteractions(R3BNeulandCosmicTrack& track) const;
            Double_t getCrossPointTime(const TVector3& point,
                                       const TVector3& direction,
                                       const TVector3& invDirection,
                                       const DPair& xBounds,
                                       const DPair& yBounds,
                                       const DPair& zBounds) const;

            std::vector<Double_t> fDistances;
            std::vector<Int_t> fBarIDs;

            R3BNeulandCosmicTrack fTrack;

            TF1 fFit;
            TGraphErrors fXZ; // i.e. Vertical Bars
            TGraphErrors fYZ; // i.e. Horizontal Bars
        };
    } // namespace Calibration
} // namespace Neuland
#endif
