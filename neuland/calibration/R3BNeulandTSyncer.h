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

#ifndef R3BNEULANDTSYNCER_H
#define R3BNEULANDTSYNCER_H

#include <array>
#include <vector>

#include "TH1F.h"

#include "R3BNeulandCommon.h"

struct LSQR_DOUBLE_VECTOR;

namespace Neuland
{
    namespace Calibration
    {
        class TSyncer
        {
          public:
            struct ValueErrorPair
            {
                Double_t Value;
                Double_t Error;
            };

            TSyncer();

            void AddBarData(const Int_t barID, const Double_t time);
            void ClearBarData(const Int_t barID);
            void DoEvent();

            std::vector<ValueErrorPair> GetTSync(UInt_t nPlanes = Neuland::MaxNumberOfPlanes);

          private:
            void calcTSyncs();

            struct Bar
            {
                ValueErrorPair TSyncNextBar;
                std::array<ValueErrorPair, BarsPerPlane> TSyncNextPlane;

                TH1F NextBar;
                std::array<TH1F, BarsPerPlane> NextPlane;

                std::vector<Double_t> NextBarLog;
                std::array<std::vector<Double_t>, BarsPerPlane> NextPlaneLog;
            };

            std::array<ULong64_t, Neuland::MaxNumberOfPlanes> HitMask;
            std::array<Double_t, Neuland::MaxNumberOfBars> EventData;
            std::array<Bar, Neuland::MaxNumberOfBars> Data;

            TH1F SamplingHistogram;
        };
    } // namespace Calibration
} // namespace Neuland

#endif
