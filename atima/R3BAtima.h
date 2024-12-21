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

#ifndef R3BATIMA_H
#define R3BATIMA_H

#include "Rtypes.h"

#include <array>
#include <vector>

namespace R3BAtima
{

    struct MaterialCompound
    {
        MaterialCompound(const Double_t mass_u = 1, const Double_t charge_e = 1)
            : Mass_u(mass_u)
            , Charge_e(charge_e)
            , Ratio(1.)
        {
        }
        MaterialCompound(const Double_t mass_u, const Double_t charge_e, const Double_t ratio)
            : Mass_u(mass_u)
            , Charge_e(charge_e)
            , Ratio(ratio)
        {
        }
        MaterialCompound(const std::array<Double_t, 3>& arr)
            : Mass_u(arr[0])
            , Charge_e(arr[1])
            , Ratio(arr[2])
        {
        }

        bool operator==(const MaterialCompound& other) const
        {
            return (Mass_u == other.Mass_u && Charge_e == other.Charge_e && Ratio == other.Ratio);
        }
        bool operator!=(const MaterialCompound& other) const { return !(*this == other); }

        Double_t Mass_u;
        Double_t Charge_e;
        Double_t Ratio;
    };

    struct TargetMaterial
    {
        static const TargetMaterial Air;
        static const TargetMaterial LH2;
        static const TargetMaterial Si;
        static const TargetMaterial BC400;

        TargetMaterial(const std::vector<MaterialCompound>& compounds, const Double_t density, const Bool_t isGas)
            : Compounds(compounds)
            , Density(density)
            , IsGas(isGas)
        {
        }

        std::vector<MaterialCompound> Compounds;
        Double_t Density;
        Bool_t IsGas;
    };

    struct TransportResult
    {
        Double_t EnergyIn_MeV_per_u;
        Double_t EnergyOut_MeV_per_u;
        Double_t ELoss_MeV_per_u;
        Double_t EStrag_MeV_per_u;
        Double_t AngStrag_mRad;
        Double_t Range_mg_per_cm2;
        Double_t RangeStrag_mg_per_cm2;
        Double_t RemainingRange_mg_per_cm2;
        Double_t dEdXIn_MeVcm2_per_mg;
        Double_t dEdXOut_MeVcm2_per_mg;
        Double_t ToF_ns;
        Double_t InterpolatedTargetThickness;
    };

    TransportResult Calculate(Double_t projMass_u,
                              Double_t projCharge_e,
                              Double_t projEnergy_MeV_per_u,
                              const TargetMaterial& targetMaterial,
                              Double_t tarThickness_mg_per_cm2);

    TransportResult Calculate_mm(Double_t projMass_u,
                                 Double_t projCharge_e,
                                 Double_t projEnergy_per_u,
                                 const TargetMaterial& targetMaterial,
                                 Double_t tarThickness_mm);

}; // namespace R3BAtima

#endif
