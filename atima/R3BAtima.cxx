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

#include "R3BAtima.h"

#include <numeric>

extern "C"
{
    void calculate_(double* proj,
                    int* pn,
                    double* ein,
                    double* targ,
                    int* nnuca,
                    int* tm,
                    double* rhoa,
                    int* igasa,
                    double* thick,
                    double* eout,
                    double* range,
                    double* dedxi,
                    double* dedxo,
                    double* remra,
                    double* rstr,
                    double* estr,
                    double* astr,
                    double* tof,
                    double* intrthick);
}

namespace R3BAtima
{

    const TargetMaterial TargetMaterial::Air = TargetMaterial({ MaterialCompound(14, 7, 0.755267),
                                                                MaterialCompound(16, 8, 0.231781),
                                                                MaterialCompound(39.948, 18, 0.012827),
                                                                MaterialCompound(12, 6, 0.000124) },
                                                              0.0012,
                                                              kTRUE);
    const TargetMaterial TargetMaterial::LH2 = TargetMaterial({ MaterialCompound(1, 1) }, 0.0708, kFALSE);
    const TargetMaterial TargetMaterial::Si = TargetMaterial({ MaterialCompound(28, 14) }, 2.336, kFALSE);
    const TargetMaterial TargetMaterial::BC400 =
        TargetMaterial({ MaterialCompound(1, 1, 0.085292), MaterialCompound(12, 6, 0.914708) }, 1.032, kFALSE);

    TransportResult Calculate(Double_t projMass_u,
                              Double_t projCharge_e,
                              Double_t projEnergy_MeV_per_u,
                              const TargetMaterial& targetMaterial,
                              Double_t tarThickness_mg_per_cm2)
    {
        TransportResult res;

        const auto invTotRatio = std::accumulate(targetMaterial.Compounds.begin(),
                                                 targetMaterial.Compounds.end(),
                                                 0.,
                                                 [](const Double_t sum, const MaterialCompound& comp) -> Double_t
                                                 { return sum + comp.Ratio; });

        Int_t tcompsize = 3, tcompnum = targetMaterial.Compounds.size(), pn = 2;
        std::vector<Double_t> tarFortran;
        tarFortran.resize(tcompnum * tcompsize);
        for (int i = 0; i < tcompnum; ++i)
        {
            tarFortran[i + tcompnum * 0] = targetMaterial.Compounds[i].Mass_u;
            tarFortran[i + tcompnum * 1] = targetMaterial.Compounds[i].Charge_e;
            tarFortran[i + tcompnum * 2] = targetMaterial.Compounds[i].Ratio * invTotRatio;
        }

        Double_t proj[2] = { projMass_u, projCharge_e };
        Double_t density = targetMaterial.Density;
        Int_t tGas = 0;
        if (targetMaterial.IsGas)
            tGas = 1;

        calculate_(proj,
                   &pn,
                   &projEnergy_MeV_per_u,
                   &tarFortran[0],
                   &tcompnum,
                   &tcompsize,
                   &density,
                   &tGas,
                   &tarThickness_mg_per_cm2,
                   &res.EnergyOut_MeV_per_u,
                   &res.Range_mg_per_cm2,
                   &res.dEdXIn_MeVcm2_per_mg,
                   &res.dEdXOut_MeVcm2_per_mg,
                   &res.RemainingRange_mg_per_cm2,
                   &res.RangeStrag_mg_per_cm2,
                   &res.EStrag_MeV_per_u,
                   &res.AngStrag_mRad,
                   &res.ToF_ns,
                   &res.InterpolatedTargetThickness);

        res.EnergyIn_MeV_per_u = projEnergy_MeV_per_u;
        res.EStrag_MeV_per_u *= res.dEdXOut_MeVcm2_per_mg;
        res.ELoss_MeV_per_u = res.EnergyIn_MeV_per_u - res.EnergyOut_MeV_per_u;
        res.AngStrag_mRad *= 1e3;
        res.dEdXIn_MeVcm2_per_mg *= projMass_u;
        res.dEdXOut_MeVcm2_per_mg *= projMass_u;
        return res;
    }
    TransportResult Calculate_mm(Double_t projMass_u,
                                 Double_t projCharge_e,
                                 Double_t projEnergy_MeV_per_u,
                                 const TargetMaterial& targetMaterial,
                                 Double_t tarThickness_mm)
    {
        return Calculate(projMass_u,
                         projCharge_e,
                         projEnergy_MeV_per_u,
                         targetMaterial,
                         targetMaterial.Density * tarThickness_mm * 100);
    }
} // namespace R3BAtima
