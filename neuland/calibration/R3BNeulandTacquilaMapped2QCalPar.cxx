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

#include "R3BNeulandTacquilaMapped2QCalPar.h"
#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "R3BEventHeader.h"
#include "R3BNeulandQCalPar.h"
#include "R3BNeulandTacquilaMappedData.h"
#include "R3BTCalEngine.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include <FairRootManager.h>

#define nPMTs 2 * fPaddles* fPlanes

R3BNeulandTacquilaMapped2QCalPar::R3BNeulandTacquilaMapped2QCalPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fPar(NULL)
    , fHits(NULL)
{
}

R3BNeulandTacquilaMapped2QCalPar::R3BNeulandTacquilaMapped2QCalPar()
    : FairTask("R3BNeulandTacquilaMapped2QCalPar")
    , fPar(NULL)
    , fHits(NULL)
{
}

R3BNeulandTacquilaMapped2QCalPar::~R3BNeulandTacquilaMapped2QCalPar()
{
    if (fData.size() > 0)
        for (Int_t plane = fPlanes - 1; plane >= 0; plane--)
            for (Int_t bar = fPaddles - 1; bar >= 0; bar--)
                for (Int_t side = 1; side <= 2; side++)
                    fData.at(plane).at(bar).at(side)->Delete();

    if (fPar)
        delete fPar;
}

InitStatus R3BNeulandTacquilaMapped2QCalPar::Init()
{
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan)
    {
        LOG(fatal) << " FairRootManager not found.";
        return kFATAL;
    }
    fHits = dynamic_cast<TClonesArray*>(fMan->GetObject("NeulandTacquilaMappedData"));
    if (!fHits)
    {
        LOG(fatal) << " Branch: NeulandTacquilaMappedData not found in Tree.";
        return kFATAL;
    }
    header = dynamic_cast<R3BEventHeader*>(fMan->GetObject("R3BEventHeader"));
    if (!header)
    {
        LOG(fatal) << " Branch: R3BEventHeader not found in Tree.";
        return kFATAL;
    }

    TH1F* histo;

    for (Int_t plane = 1; plane <= fPlanes; plane++)
    {
        std::vector<std::vector<TH1F*>> v_plane;
        for (Int_t bar = 1; bar <= fPaddles; bar++)
        {
            std::vector<TH1F*> v_bar;
            for (Int_t side = 1; side <= 2; side++)
            {
                std::string s =
                    "h_pedpar_" + std::to_string(plane) + "_" + std::to_string(bar) + "_" + std::to_string(side);
                histo = new TH1F(s.c_str(), "PedPar", 2001, -0.5, 2000.5);
                v_bar.push_back(histo);
            }
            v_plane.push_back(v_bar);
        }
        fData.push_back(v_plane);
    }
    fPar = dynamic_cast<R3BNeulandQCalPar*>(FairRuntimeDb::instance()->getContainer("NeulandQCalPar"));

    return kSUCCESS;
}

void R3BNeulandTacquilaMapped2QCalPar::Exec(Option_t* option)
{
    R3BNeulandTacquilaMappedData* hit;
    Int_t nHits = fHits->GetEntries();

    if (header->GetTrigger() != 2)
        return;

    if (nHits < nPMTs)
        return;

    for (Int_t i = 0; i < nHits; i++)
    {
        hit = dynamic_cast<R3BNeulandTacquilaMappedData*>(fHits->At(i));
        if (!hit)
            continue;

        fData.at(hit->GetPlane() - 1).at(hit->GetPaddle() - 1).at(hit->GetSide() - 1)->Fill(hit->GetQdcData());
    }
}

void R3BNeulandTacquilaMapped2QCalPar::FinishTask()
{
    Int_t i = 0;
    TH1F* pars = new TH1F("QCalPar", "Pedestal Offset", nPMTs, 0.5, nPMTs + 0.5);
    TH1F* histo;
    fPar->SetSize(nPMTs);
    for (Int_t plane = 1; plane <= fPlanes; plane++)
        for (Int_t bar = 1; bar <= fPaddles; bar++)
            for (Int_t side = 1; side <= 2; side++)
            {
                histo = fData.at(plane - 1).at(bar - 1).at(side - 1);
                Int_t maxBinPos = histo->GetBinCenter(histo->GetMaximumBin());
                pars->SetBinContent(i++ + 1, maxBinPos);
                fPar->SetParAt(plane, bar, side, maxBinPos);
            }
    pars->Write();
    fPar->setChanged();
}

ClassImp(R3BNeulandTacquilaMapped2QCalPar)
