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

#include "R3BNeulandCal2Hit.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "R3BEventHeader.h"
#include "R3BNeulandCommon.h"
#include "R3BNeulandHitPar.h"
#include "R3BTCalPar.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TGraph2D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TVector.h"

#include <array>

R3BNeulandCal2Hit::R3BNeulandCal2Hit(const char* name, const Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fEventHeader(nullptr)
    , fCalData("NeulandCalData")
    , fHits("NeulandHits")
    , fFirstPlaneHorizontal(true)
    , fDistanceToTarget(Neuland::NaN)
    , fGlobalTimeOffset(Neuland::NaN)
    , fEnergyCutoff(Neuland::NaN)
    , fWalkEnabled(kFALSE)
{
}

InitStatus R3BNeulandCal2Hit::Init()
{
    auto ioman = FairRootManager::Instance();
    if (ioman == nullptr)
    {
        throw std::runtime_error("R3BNeulandCal2Hit: No FairRootManager");
    }

    fEventHeader = dynamic_cast<R3BEventHeader*>(ioman->GetObject("EventHeader."));
    if (fEventHeader == nullptr)
    {
        throw std::runtime_error("R3BNeulandCal2Hit: No R3BEventHeader");
    }

    fCalData.Init();
    fHits.Init();

    SetParameter();
    return kSUCCESS;
}

void R3BNeulandCal2Hit::SetParContainers()
{
    fPar = dynamic_cast<R3BNeulandHitPar*>(FairRuntimeDb::instance()->getContainer("NeulandHitPar"));
}

void R3BNeulandCal2Hit::SetParameter()
{
    fParMap.clear();

    if (std::isnan(fDistanceToTarget))
    {
        fDistanceToTarget = fPar->GetDistanceToTarget();
    }
    else
    {
        fPar->SetDistanceToTarget(fDistanceToTarget);
    }

    if (std::isnan(fGlobalTimeOffset))
    {
        fGlobalTimeOffset = fPar->GetGlobalTimeOffset();
    }
    else
    {
        fPar->SetGlobalTimeOffset(fGlobalTimeOffset);
    }

    if (std::isnan(fEnergyCutoff))
    {
        fEnergyCutoff = fPar->GetEnergyCutoff();
    }
    else
    {
        fPar->SetEnergyCutoff(fEnergyCutoff);
    }

    fNumberOfPlanes = fPar->GetNumberOfPlanes();
    fDistancesToFirstPlane = fPar->GetDistancesToFirstPlane();
    fAttenuationValues.resize(fNumberOfPlanes * Neuland::BarsPerPlane, 0.);
    const auto nPars = fPar->GetNumModulePar();

    for (auto i = 0; i < nPars; i++)
    {
        const auto id = fPar->GetModuleParAt(i)->GetModuleId() - 1;
        fParMap[id] = *fPar->GetModuleParAt(i);
        fAttenuationValues[id] = exp(Neuland::TotalBarLength / fParMap[id].GetLightAttenuationLength());
    }

    LOG(info) << "R3BNeulandCal2Hit::SetParameter : Number of Parameters: " << fPar->GetNumModulePar();
}

InitStatus R3BNeulandCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BNeulandCal2Hit::Exec(Option_t*)
{

    if (fVerbose > 0 && ++fEventNumber % 10000 == 0)
        std::cout << "\rR3BNeulandCal2Hit " << fEventNumber << " Events converted." << std::flush;

    fHits.Reset();
    fHitMap.clear();

    auto calData = fCalData.Retrieve();

    const auto start = fEventHeader->GetTStart();
    const bool beam = !std::isnan(start);

    for (auto calDataPtrPtr = calData.begin(); calDataPtrPtr != calData.end(); calDataPtrPtr++)
    {
        const auto barID = (*calDataPtrPtr)->GetBarId() - 1;

        if (fParMap.find(barID) == fParMap.end())
            continue; // We do not have parameters for this module

        if (fHitMap.find(barID) == fHitMap.end()) // TODO what about multi hit?
        {
            // This is the first PMT of this bar, store the pointer and go to next cal entry
            fHitMap[barID] = *calDataPtrPtr;
            continue;
        }

        const auto& parameter = fParMap[barID];

        std::array<const R3BNeulandCalData*, 2> cal;
        if ((*calDataPtrPtr)->GetSide() == 1)
            cal = { *calDataPtrPtr, fHitMap[barID] };
        else
            cal = { fHitMap[barID], *calDataPtrPtr };

        const std::array<int, 2> qdc = { std::max(cal[0]->GetQdc() - parameter.GetPedestal(1), 1),
                                         std::max(cal[1]->GetQdc() - parameter.GetPedestal(2), 1) };

        const std::array<Double_t, 2> unsatEnergy = {
            GetUnsaturatedEnergy(qdc[0], parameter.GetEnergyGain(1), parameter.GetPMTSaturation(1)),
            GetUnsaturatedEnergy(qdc[1], parameter.GetEnergyGain(2), parameter.GetPMTSaturation(2))
        };

        const auto energy = TMath::Sqrt(fAttenuationValues[barID] * unsatEnergy[0] * unsatEnergy[1]);

        if (energy < fEnergyCutoff)
            continue;

        std::array<Double_t, 2> tdc;

        if (std::isnan(cal[0]->GetTriggerTime()) || std::isnan(cal[0]->GetTriggerTime()))
        {
            tdc = { cal[0]->GetTime() + parameter.GetTimeOffset(1), cal[1]->GetTime() + parameter.GetTimeOffset(2) };
        }
        else
        {
            tdc = { cal[0]->GetTime() - cal[0]->GetTriggerTime() + parameter.GetTimeOffset(1),
                    cal[1]->GetTime() - cal[1]->GetTriggerTime() + parameter.GetTimeOffset(2) };
        }

        // FIXME this should be done in Mapped2Cal
        // In Cal2Hit the difference between all bars should be checked
        if (tdc[0] - tdc[1] < -0.5 * Neuland::MaxCalTime)
            tdc[1] -= Neuland::MaxCalTime;
        else if (tdc[0] - tdc[1] > 0.5 * Neuland::MaxCalTime)
            tdc[0] -= Neuland::MaxCalTime;

        auto time = (tdc[0] + tdc[1]) * 0.5 - fGlobalTimeOffset;

        if (fWalkEnabled)
            time = time + WalkCorrection(energy);

        if (beam)
        {
            // the shift is to get fmod to work as indented: 4 peaks -> 1 peak w/o stray data (e.g. at 5 * 2048)
            // fmod gives the time between 0 and 5*2048
            // remainder() gives the time between -0.5*5*2048 and +0.5*5*2048
            time = remainder(time - start, Neuland::MaxCalTime);
        }
        else
        {
            time = Neuland::NaN;
        }

        const auto plane = Neuland::GetPlaneNumber(barID);
        const auto bar = (barID) % 50;

        TVector3 pos;
        TVector3 pixel;

        if (Neuland::IsPlaneHorizontal(plane) == fFirstPlaneHorizontal)
        {
            pos[0] = parameter.GetEffectiveSpeed() * (tdc[1] - tdc[0]);
            pos[1] = (bar + 0.5 - Neuland::BarsPerPlane * 0.5) * Neuland::BarSize_XY;

            pixel[0] = std::min(std::max(0., pos[0] / 5. + 25), 49.);
            pixel[1] = bar;
        }
        else
        {
            pos[0] = (bar + 0.5 - Neuland::BarsPerPlane * 0.5) * Neuland::BarSize_XY;
            pos[1] = parameter.GetEffectiveSpeed() * (tdc[1] - tdc[0]);

            pixel[0] = bar;
            pixel[1] = std::min(std::max(0., pos[1] / 5. + 25), 49.);
        }

        pos[2] = fDistanceToTarget + fDistancesToFirstPlane[plane];
        pixel[2] = plane;

        fHits.Insert({ barID + 1, tdc[0], tdc[1], time, unsatEnergy[0], unsatEnergy[1], energy, pos, pixel });
    }
}

void R3BNeulandCal2Hit::FinishTask()
{
    if (fVerbose > 0)
        std::cout << "\rR3BNeulandCal2Hit " << fEventNumber << " Events converted." << std::endl;
}

Double_t R3BNeulandCal2Hit::GetUnsaturatedEnergy(const Int_t qdc, const Double_t gain, const Double_t saturation) const
{
    return qdc / (gain - saturation * qdc);
}

Double_t R3BNeulandCal2Hit::WalkCorrection(Double_t x)
{
    Double_t y = 0;

    // if (x<=12.) y=24.78;

    if (x <= 12.)
        y = 35. - 7.1897575 * log(x) + 1.2 * log(x) * log(x) + 0.23616;

    if (x > 12.)
        y = 25.5 - 0.055 * x;

    return 24.78 - y;
}

ClassImp(R3BNeulandCal2Hit)
