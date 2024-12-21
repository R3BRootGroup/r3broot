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

#include "R3BNeulandCal2HitPar.h"
#include <FairRootManager.h>

#include "R3BEventHeader.h"
#include "R3BNeulandCalData.h"
#include "R3BNeulandCosmicTracker.h"
#include "R3BNeulandHitCalibrationEngine.h"
#include "R3BNeulandHitModulePar.h"
#include "R3BNeulandHitPar.h"

#include "FairLogger.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TH2F.h"
#include "TROOT.h"

using namespace Neuland;

R3BNeulandCal2HitPar::R3BNeulandCal2HitPar(const char* name, const Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fHitCalEngine(new Calibration::HitCalibrationEngine())
    , fCosmicTracker(new Calibration::CosmicTracker())
    , fNeulandHitPar(nullptr)
    , fCalNeuland(nullptr)
    , fMappedLos(nullptr)
    , fCosmicTpat(0)
    , fSavePlots(kFALSE)
{
}

R3BNeulandCal2HitPar::~R3BNeulandCal2HitPar() {}

InitStatus R3BNeulandCal2HitPar::Init()
{
    FairRootManager* mgr = FairRootManager::Instance();

    if (!mgr)
    {
        LOG(fatal) << "R3BNeulandCal2HitPar::Init: "
                   << "FairRootManager not found";
        return kFATAL;
    }

    fEventHeader = static_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    if (!fEventHeader)
    {
        LOG(info) << "R3BNeulandCal2HitPar::Init: "
                  << "Branch EventHeader not found, its ok";
        // return kFATAL;
    }

    fCalNeuland = static_cast<TClonesArray*>(mgr->GetObject("NeulandCalData"));

    if (!fCalNeuland)
    {
        LOG(fatal) << "R3BNeulandCal2HitPar::Init: "
                   << "Branch NeulandCalData not found";
        return kFATAL;
    }

    fNeulandHitPar = static_cast<R3BNeulandHitPar*>(FairRuntimeDb::instance()->getContainer("NeulandHitPar"));
    LOG(info) << "R3BNeulandCal2HitPar::Init: "
              << "Number of Hit-Paramteres found: " << fNeulandHitPar->GetNumModulePar();

    LOG(info) << "R3BNeulandCal2HitPar::Init: "
              << "Initializing NeulandHitCalibrationEngine.";
    fHitCalEngine->Init(fNeulandHitPar);

    if (fNeulandHitPar)
        fCosmicTracker->SetDistances(fNeulandHitPar->GetDistancesToFirstPlane());

    return kSUCCESS;
}

void R3BNeulandCal2HitPar::Exec(Option_t* option)
{
    if (++fEventNumber % 100000 == 0 && fVerbose) // 10000
    {
        const auto msg = TString::Format(
            "R3BNeulandCal2HitPar::Exec : Event: %10d,    accepted Events: %10d", fEventNumber, fAcceptedEventNumber);
        std::cout << msg << "\r" << std::flush;

        std::cout << " cnts   " << fIgorcnt0 << " ---- " << fIgorcnt1 << " ---- " << fIgorcnt2 << "     " << std::endl;
    }

    // Skip if this is not a pure cosmic event
    if (!IsCosmicEvent())
    {
        return;
    }

    ++fIgorcnt0;

    LOG(debug) << "R3BNeulandCal2HitPar::Exec: Event " << fEventNumber - 1;
    const auto nItems = fCalNeuland->GetEntriesFast();

    if (nItems < 8)
    {
        LOG(debug) << "   Event cannot be used: too few signals : " << nItems << "!";
        return;
    }

    fHitCalEngine->Reset();
    fCosmicTracker->Reset();

    Int_t addedPoints = 0;

    for (Int_t i = 0; i < nItems; i++)
    {
        const auto pmt = static_cast<R3BNeulandCalData*>(fCalNeuland->At(i));

        const auto id = pmt->GetBarId() - 1;
        const auto plane = GetPlaneNumber(id);
        const auto side = pmt->GetSide() - 1;

        if (std::isnan(pmt->GetTriggerTime()))
        {
            fHitCalEngine->Set(id, side, pmt->GetTime(), pmt->GetQdc());
        }
        else
        {
            fHitCalEngine->Set(id, side, pmt->GetTime() - pmt->GetTriggerTime(), pmt->GetQdc());
        }

        if (fHitCalEngine->IsValid(id))
        {
            ++addedPoints;
            fCosmicTracker->AddPoint(id, fHitCalEngine->GetPosition(id));
        }
    }

    ++fIgorcnt1;

    if (addedPoints < 6)
    {
        LOG(debug) << "   Event cannot be used: too few Points : " << addedPoints << " (" << nItems << ")!";
        return;
    }

    ++fIgorcnt2;

    const auto& cosmicTrack = fCosmicTracker->GetTrack();

    if (cosmicTrack.Interactions.size() == 0)
    {
        LOG(debug) << "   Getting Cosmic Track :  Failure!";
        return;
    }

    LOG(debug) << "   Getting Cosmic Track:  Success!";
    ++fAcceptedEventNumber;

    LOG(debug) << "   Adding data to calibration";
    fHitCalEngine->Add(cosmicTrack, fEventNumber);
}

void R3BNeulandCal2HitPar::FinishTask()
{
    LOG(info) << "R3BNeulandCal2HitPar::FinishTask: "
              << "Saved " << fAcceptedEventNumber << " Events.                       ";

    LOG(info) << "R3BNeulandCal2HitPar::FinishTask: "
              << "Starting Neuland Hit Calibration with " << fAcceptedEventNumber << " Events.";

    const auto batchMode = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);
    const auto defdir = gDirectory;

    TDirectory* neulandDir = nullptr;
    if (fSavePlots)
        neulandDir = defdir->mkdir("NeuLAND");

    auto parameters = fHitCalEngine->Calibrate(neulandDir);

    defdir->cd();
    gROOT->SetBatch(batchMode);

    fNeulandHitPar->GetListOfModulePar()->Clear();

    auto maxThreshold = 0.;
    for (auto& parameter : parameters)
    {
        fNeulandHitPar->AddModulePar(new R3BNeulandHitModulePar(parameter));

        for (auto side = 1; side <= 2; ++side)
        {
            const auto threshold =
                parameter.GetPMTThreshold(side) * exp(TotalBarLength / parameter.GetLightAttenuationLength());
            if (threshold > maxThreshold)
                maxThreshold = threshold;
        }
    }

    fNeulandHitPar->SetEnergyCutoff(maxThreshold);
    LOG(info) << TString::Format("R3BNeulandCal2HitPar::FinishTask: Recommended Minimum Energy Cutoff : >%4.2f MeV\n",
                                 maxThreshold);

    fNeulandHitPar->setChanged();

    LOG(info) << "R3BNeulandCal2HitPar::FinishTask: "
              << "Number of calibrated Bars: " << fNeulandHitPar->GetNumModulePar();
}

bool R3BNeulandCal2HitPar::IsCosmicEvent() const { return ((fEventHeader->GetTpat() & fCosmicTpat) == fCosmicTpat); }

ClassImp(R3BNeulandCal2HitPar)
