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

#ifndef NEULAND_DIGITIZING_TACQUILA_H
#define NEULAND_DIGITIZING_TACQUILA_H

#include "R3BDigitizingEngine.h"
#include "TRandom3.h"
#include "Validated.h"

namespace R3B::Digitizing::Neuland::TacQuila
{
    struct Params
    {
        double fPMTThresh;             // [MeV]
        double fSaturationCoefficient; //
        Bool_t fExperimentalDataIsCorrectedForSaturation;
        double fTimeRes;         // time + Gaus(0., fTimeRes) [ns]
        double fEResRel;         // Gaus(e, fEResRel * e) []
        double fIntegrationTime; // [ns]
        std::shared_ptr<TRandom3> fRnd;

        Params();
    };

    extern const Params TACQUILA_DEFAULT_PARAM;
    class Channel : public Digitizing::Channel
    {
      public:
        Channel(const Channel&) = delete;
        Channel(Channel&&) = delete;
        Channel& operator=(const Channel&) = delete;
        Channel& operator=(Channel&&) = delete;
        explicit Channel(ChannelSide, const TacQuila::Params& = TACQUILA_DEFAULT_PARAM);
        ~Channel() override = default;
        void AddHit(Hit newHit) override;
        bool HasFired() override;
        double GetQDC();
        double GetTDC();
        double GetEnergy();
        double GetTrigTime() override { return GetTDC(); }

      private:
        // NOTE: Some expensive calculations and random distributions are cached
        // so they do not need to be recomputed every time a Getter is called
        auto ConstructSignals() -> Signals override;
        std::vector<Hit> fPMTHits;
        mutable Validated<std::vector<Hit>::const_iterator> cachedFirstHitOverThresh;
        std::vector<Hit>::const_iterator FindThresholdExceedingHit() const;

        double BuildQDC();
        mutable Validated<double> cachedQDC;

        double BuildTDC();
        mutable Validated<double> cachedTDC;

        double BuildEnergy();
        mutable Validated<double> cachedEnergy;

        const TacQuila::Params& par;
    };

} // namespace R3B::Digitizing::Neuland::TacQuila

#endif // NEULAND_DIGITIZING_TACQUILA_H
