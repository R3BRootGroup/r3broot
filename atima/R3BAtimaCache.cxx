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

#include "R3BAtimaCache.h"

#include "FairLogger.h"

#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>

constexpr UInt_t MagicWord = 0xF0F00000;

namespace R3BAtima
{

    Cache::Cache(const Double_t pMass_u,
                 const Double_t pCharge_e,
                 const RangeSelector& energies_MeV_per_u,
                 const TargetMaterial& targetMaterial,
                 const RangeSelector& distances_mm)
        : fProjMass(pMass_u)
        , fProjCharge(pCharge_e)
        , fEnergies(energies_MeV_per_u)
        , fTargetMaterial(targetMaterial)
        , fDistances(distances_mm)
    {
        calculate();
    }

    Cache::Cache(const Double_t pMass_u,
                 const Double_t pCharge_e,
                 const RangeSelector& energies_MeV_per_u,
                 const TargetMaterial& targetMaterial,
                 const RangeSelector& distances_mm,
                 const TString& path)
        : fProjMass(pMass_u)
        , fProjCharge(pCharge_e)
        , fEnergies(energies_MeV_per_u)
        , fTargetMaterial(targetMaterial)
        , fDistances(distances_mm)
    {
        if (!read(path))
        {
            LOG(info) << "AtimaCache-File '" << path
                      << "' does not exist or does not match the provided Data! Will (re)calculate and (over)write!";

            fProjMass = pMass_u;
            fProjCharge = pCharge_e;
            fEnergies = energies_MeV_per_u;
            fTargetMaterial = targetMaterial;
            fDistances = distances_mm;

            calculate();
            write(path);
        }
    }

    TransportResult Cache::operator()(const Double_t energy_MeV_per_u, const Double_t distance_mm) const
    {

        if (energy_MeV_per_u < fEnergies.MinValue || energy_MeV_per_u > fEnergies.MaxValue ||
            distance_mm < fDistances.MinValue || distance_mm > fDistances.MaxValue)
            LOG(fatal) << "R3BAtima::Cache: given value outside of calculated range!";

        TransportResult res;

        res.EnergyIn_MeV_per_u = energy_MeV_per_u;
        res.ELoss_MeV_per_u = fG_ELoss.Interpolate(energy_MeV_per_u, distance_mm);
        res.EnergyOut_MeV_per_u = res.EnergyIn_MeV_per_u - res.ELoss_MeV_per_u;
        res.EStrag_MeV_per_u = fG_EStrag.Interpolate(energy_MeV_per_u, distance_mm);
        res.AngStrag_mRad = fG_AngStrag.Interpolate(energy_MeV_per_u, distance_mm);
        res.Range_mg_per_cm2 = fG_Range.Interpolate(energy_MeV_per_u, distance_mm);
        res.RemainingRange_mg_per_cm2 = fG_RemainigRange.Interpolate(energy_MeV_per_u, distance_mm);
        res.dEdXIn_MeVcm2_per_mg = fG_dEdXIn.Interpolate(energy_MeV_per_u, distance_mm);
        res.dEdXOut_MeVcm2_per_mg = fG_dEdXOut.Interpolate(energy_MeV_per_u, distance_mm);
        res.ToF_ns = fG_ToF.Interpolate(energy_MeV_per_u, distance_mm);
        res.InterpolatedTargetThickness = fG_InterpolatedTargetThickness.Interpolate(energy_MeV_per_u, distance_mm);

        return res;
    }

    Bool_t Cache::read(const TString& path)
    {
        std::ifstream fstream(path.Data(), std::ifstream::binary);

        if (!fstream.good())
            return kFALSE;

        UInt_t Version;
        fstream.read(reinterpret_cast<char*>(&Version), sizeof(Version));
        if (Version != MagicWord)
            return kFALSE;

        Double_t projM, projCharge;

        fstream.read(reinterpret_cast<char*>(&projM), sizeof(projM));
        fstream.read(reinterpret_cast<char*>(&projCharge), sizeof(projCharge));

        if (projM != fProjMass || projCharge != fProjCharge)
            return kFALSE;

        Double_t mass, charge, ratio;

        Int_t nComps;
        fstream.read(reinterpret_cast<char*>(&nComps), sizeof(nComps));

        if (nComps != fTargetMaterial.Compounds.size())
            return kFALSE;

        for (auto i = 0; i < nComps; ++i)
        {
            fstream.read(reinterpret_cast<char*>(&mass), sizeof(mass));
            fstream.read(reinterpret_cast<char*>(&charge), sizeof(charge));
            fstream.read(reinterpret_cast<char*>(&ratio), sizeof(ratio));

            if (fTargetMaterial.Compounds[i].Mass_u != mass || fTargetMaterial.Compounds[i].Charge_e != charge ||
                fTargetMaterial.Compounds[i].Ratio != ratio)
                return kFALSE;
        }

        Double_t minVal, maxVal;
        Int_t steps;

        fstream.read(reinterpret_cast<char*>(&minVal), sizeof(minVal));
        fstream.read(reinterpret_cast<char*>(&maxVal), sizeof(maxVal));
        fstream.read(reinterpret_cast<char*>(&steps), sizeof(steps));

        if (minVal != fEnergies.MinValue || maxVal != fEnergies.MaxValue || steps != fEnergies.Steps)
            return kFALSE;

        fstream.read(reinterpret_cast<char*>(&minVal), sizeof(minVal));
        fstream.read(reinterpret_cast<char*>(&maxVal), sizeof(maxVal));
        fstream.read(reinterpret_cast<char*>(&steps), sizeof(steps));

        if (minVal != fDistances.MinValue || maxVal != fDistances.MaxValue || steps != fDistances.Steps)
            return kFALSE;

        Double_t density;
        Bool_t isGas;

        fstream.read(reinterpret_cast<char*>(&density), sizeof(density));
        fstream.read(reinterpret_cast<char*>(&isGas), sizeof(isGas));

        if (density != fTargetMaterial.Density || isGas != fTargetMaterial.IsGas)
            return kFALSE;

        const auto maxStep = (fEnergies.Steps + 1) * (fDistances.Steps + 1);

        fG_ELoss.Set(maxStep);
        fG_EStrag.Set(maxStep);
        fG_AngStrag.Set(maxStep);
        fG_Range.Set(maxStep);
        fG_RemainigRange.Set(maxStep);
        fG_dEdXIn.Set(maxStep);
        fG_dEdXOut.Set(maxStep);
        fG_ToF.Set(maxStep);
        fG_InterpolatedTargetThickness.Set(maxStep);

        Int_t nPoints = 0;
        Double_t energy = fEnergies.MinValue;
        for (int i = 0; i <= fEnergies.Steps; ++i)
        {
            Double_t distance = fDistances.MinValue;
            for (int j = 0; j <= fDistances.Steps; ++j)
            {

                TransportResult res;

                fstream.read(reinterpret_cast<char*>(&res.ELoss_MeV_per_u), sizeof(res.ELoss_MeV_per_u));
                fstream.read(reinterpret_cast<char*>(&res.EStrag_MeV_per_u), sizeof(res.EStrag_MeV_per_u));
                fstream.read(reinterpret_cast<char*>(&res.AngStrag_mRad), sizeof(res.AngStrag_mRad));
                fstream.read(reinterpret_cast<char*>(&res.Range_mg_per_cm2), sizeof(res.Range_mg_per_cm2));
                fstream.read(reinterpret_cast<char*>(&res.RemainingRange_mg_per_cm2),
                             sizeof(res.RemainingRange_mg_per_cm2));
                fstream.read(reinterpret_cast<char*>(&res.dEdXIn_MeVcm2_per_mg), sizeof(res.dEdXIn_MeVcm2_per_mg));
                fstream.read(reinterpret_cast<char*>(&res.dEdXOut_MeVcm2_per_mg), sizeof(res.dEdXOut_MeVcm2_per_mg));
                fstream.read(reinterpret_cast<char*>(&res.ToF_ns), sizeof(res.ToF_ns));
                fstream.read(reinterpret_cast<char*>(&res.InterpolatedTargetThickness),
                             sizeof(res.InterpolatedTargetThickness));

                fG_ELoss.SetPoint(nPoints, energy, distance, res.ELoss_MeV_per_u);
                fG_EStrag.SetPoint(nPoints, energy, distance, res.EStrag_MeV_per_u);
                fG_AngStrag.SetPoint(nPoints, energy, distance, res.AngStrag_mRad);
                fG_Range.SetPoint(nPoints, energy, distance, res.Range_mg_per_cm2);
                fG_RemainigRange.SetPoint(nPoints, energy, distance, res.RemainingRange_mg_per_cm2);
                fG_dEdXIn.SetPoint(nPoints, energy, distance, res.dEdXIn_MeVcm2_per_mg);
                fG_dEdXOut.SetPoint(nPoints, energy, distance, res.dEdXOut_MeVcm2_per_mg);
                fG_ToF.SetPoint(nPoints, energy, distance, res.ToF_ns);
                fG_InterpolatedTargetThickness.SetPoint(nPoints, energy, distance, res.InterpolatedTargetThickness);

                ++nPoints;
                distance += (fDistances.MaxValue - fDistances.MinValue) / fDistances.Steps;
            }
            energy += (fEnergies.MaxValue - fEnergies.MinValue) / fEnergies.Steps;
        }

        if (!fstream.good())
            return kFALSE;

        Version = 0;

        fstream.read(reinterpret_cast<char*>(&Version), sizeof(Version));

        if (Version != MagicWord)
            return kFALSE;

        return kTRUE;
    }

    void Cache::write(const TString& path) const
    {
        std::ofstream fstream(path.Data(), std::ofstream::binary | std::ofstream::trunc);

        UInt_t Version = MagicWord;
        fstream.write(reinterpret_cast<const char*>(&Version), sizeof(Version));

        fstream.write(reinterpret_cast<const char*>(&fProjMass), sizeof(fProjMass));
        fstream.write(reinterpret_cast<const char*>(&fProjCharge), sizeof(fProjCharge));

        Int_t nComps = fTargetMaterial.Compounds.size();
        fstream.write(reinterpret_cast<const char*>(&nComps), sizeof(nComps));

        for (auto i = 0; i < nComps; ++i)
        {
            fstream.write(reinterpret_cast<const char*>(&fTargetMaterial.Compounds[i].Mass_u),
                          sizeof(fTargetMaterial.Compounds[i].Mass_u));
            fstream.write(reinterpret_cast<const char*>(&fTargetMaterial.Compounds[i].Charge_e),
                          sizeof(fTargetMaterial.Compounds[i].Charge_e));
            fstream.write(reinterpret_cast<const char*>(&fTargetMaterial.Compounds[i].Ratio),
                          sizeof(fTargetMaterial.Compounds[i].Ratio));
        }

        fstream.write(reinterpret_cast<const char*>(&fEnergies.MinValue), sizeof(fEnergies.MinValue));
        fstream.write(reinterpret_cast<const char*>(&fEnergies.MaxValue), sizeof(fEnergies.MaxValue));
        fstream.write(reinterpret_cast<const char*>(&fEnergies.Steps), sizeof(fEnergies.Steps));

        fstream.write(reinterpret_cast<const char*>(&fDistances.MinValue), sizeof(fDistances.MinValue));
        fstream.write(reinterpret_cast<const char*>(&fDistances.MaxValue), sizeof(fDistances.MaxValue));
        fstream.write(reinterpret_cast<const char*>(&fDistances.Steps), sizeof(fDistances.Steps));

        fstream.write(reinterpret_cast<const char*>(&fTargetMaterial.Density), sizeof(fTargetMaterial.Density));
        fstream.write(reinterpret_cast<const char*>(&fTargetMaterial.IsGas), sizeof(fTargetMaterial.IsGas));

        Int_t nPoints = 0;
        for (int i = 0; i <= fEnergies.Steps; ++i)
        {
            for (int j = 0; j <= fDistances.Steps; ++j)
            {
                fstream.write(reinterpret_cast<char*>(&fG_ELoss.GetZ()[nPoints]), sizeof(Double_t));
                fstream.write(reinterpret_cast<char*>(&fG_EStrag.GetZ()[nPoints]), sizeof(Double_t));
                fstream.write(reinterpret_cast<char*>(&fG_AngStrag.GetZ()[nPoints]), sizeof(Double_t));
                fstream.write(reinterpret_cast<char*>(&fG_Range.GetZ()[nPoints]), sizeof(Double_t));
                fstream.write(reinterpret_cast<char*>(&fG_RemainigRange.GetZ()[nPoints]), sizeof(Double_t));
                fstream.write(reinterpret_cast<char*>(&fG_dEdXIn.GetZ()[nPoints]), sizeof(Double_t));
                fstream.write(reinterpret_cast<char*>(&fG_dEdXOut.GetZ()[nPoints]), sizeof(Double_t));
                fstream.write(reinterpret_cast<char*>(&fG_ToF.GetZ()[nPoints]), sizeof(Double_t));
                fstream.write(reinterpret_cast<char*>(&fG_InterpolatedTargetThickness.GetZ()[nPoints]),
                              sizeof(Double_t));

                ++nPoints;
            }
        }

        fstream.write(reinterpret_cast<const char*>(&Version), sizeof(Version));
    }

    void Cache::calculate()
    {
        UInt_t currStep = 0;
        const auto maxStep = (fEnergies.Steps + 1) * (fDistances.Steps + 1);

        fG_ELoss.Clear();
        fG_ELoss.Set(maxStep);
        fG_EStrag.Clear();
        fG_EStrag.Set(maxStep);
        fG_AngStrag.Clear();
        fG_AngStrag.Set(maxStep);
        fG_Range.Clear();
        fG_Range.Set(maxStep);
        fG_RemainigRange.Clear();
        fG_RemainigRange.Set(maxStep);
        fG_dEdXIn.Clear();
        fG_dEdXIn.Set(maxStep);
        fG_dEdXOut.Clear();
        fG_dEdXOut.Set(maxStep);
        fG_ToF.Clear();
        fG_ToF.Set(maxStep);
        fG_InterpolatedTargetThickness.Clear();
        fG_InterpolatedTargetThickness.Set(maxStep);

        auto energy = fEnergies.MinValue;
        for (int i = 0; i <= fEnergies.Steps; ++i)
        {
            auto distance = fDistances.MinValue;
            for (int j = 0; j <= fDistances.Steps; ++j)
            {
                auto res = Calculate_mm(fProjMass, fProjCharge, energy, fTargetMaterial, distance);

                fG_ELoss.SetPoint(currStep, res.EnergyIn_MeV_per_u, distance, res.ELoss_MeV_per_u);
                fG_EStrag.SetPoint(currStep, res.EnergyIn_MeV_per_u, distance, res.EStrag_MeV_per_u);
                fG_AngStrag.SetPoint(currStep, res.EnergyIn_MeV_per_u, distance, res.AngStrag_mRad);
                fG_Range.SetPoint(currStep, res.EnergyIn_MeV_per_u, distance, res.Range_mg_per_cm2);
                fG_RemainigRange.SetPoint(currStep, res.EnergyIn_MeV_per_u, distance, res.RemainingRange_mg_per_cm2);
                fG_dEdXIn.SetPoint(currStep, res.EnergyIn_MeV_per_u, distance, res.dEdXIn_MeVcm2_per_mg);
                fG_dEdXOut.SetPoint(currStep, res.EnergyIn_MeV_per_u, distance, res.dEdXOut_MeVcm2_per_mg);
                fG_ToF.SetPoint(currStep, res.EnergyIn_MeV_per_u, distance, res.ToF_ns);
                fG_InterpolatedTargetThickness.SetPoint(
                    currStep, res.EnergyIn_MeV_per_u, distance, res.InterpolatedTargetThickness);

                distance += (fDistances.MaxValue - fDistances.MinValue) / fDistances.Steps;

                std::cout << "\rR3BAtima::Cache: Creating Cache: " << ++currStep << "/" << maxStep << std::flush;
            }

            energy += (fEnergies.MaxValue - fEnergies.MinValue) / fEnergies.Steps;
        }
        std::cout << std::endl;
    }
} // namespace R3BAtima
