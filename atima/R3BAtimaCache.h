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

#ifndef R3BATIMACACHE_H
#define R3BATIMACACHE_H

#include "TGraph2D.h"
#include "TString.h"

#include "R3BAtima.h"

namespace R3BAtima
{
    class Cache
    {
      public:
        struct RangeSelector
        {
            // RangeSelector(const Double_t min, const Double_t max, const Int_t steps) : MinValue(min), MaxValue(max),
            // Steps(steps){}

            bool operator==(const RangeSelector& other) const
            {
                return (MinValue == other.MinValue && MaxValue == other.MaxValue && Steps == other.Steps);
            }
            bool operator!=(const RangeSelector& other) const { return !(*this == other); }

            Double_t MinValue;
            Double_t MaxValue;
            Int_t Steps;
        };

        Cache(const Double_t pMass_u,
              const Double_t pCharge_e,
              const RangeSelector& energies_MeV_per_u,
              const TargetMaterial& targetMaterial,
              const RangeSelector& distances_mm);
        Cache(const Double_t pMass_u,
              const Double_t pCharge_e,
              const RangeSelector& energies_MeV_per_u,
              const TargetMaterial& targetMaterial,
              const RangeSelector& distances_mm,
              const TString& path);

        TransportResult operator()(const Double_t energy_MeV_per_u, const Double_t distance_mm) const;

      private:
        Bool_t read(const TString& path);
        void write(const TString& path) const;
        void calculate();

        Double_t fProjMass;
        Double_t fProjCharge;
        RangeSelector fEnergies;
        TargetMaterial fTargetMaterial; //!
        RangeSelector fDistances;

        mutable TGraph2D fG_ELoss;
        mutable TGraph2D fG_EStrag;
        mutable TGraph2D fG_AngStrag;
        mutable TGraph2D fG_Range;
        mutable TGraph2D fG_RemainigRange;
        mutable TGraph2D fG_dEdXIn;
        mutable TGraph2D fG_dEdXOut;
        mutable TGraph2D fG_ToF;
        mutable TGraph2D fG_InterpolatedTargetThickness;
    };
} // namespace R3BAtima
#endif
