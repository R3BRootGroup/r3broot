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

// -------------------------------------------------------------
// -----         R3BTwimGainMatching source file           -----
// -----    Created 16/02/22 by Antia Grana Gonzalez       -----
// -------------------------------------------------------------

// ROOT headers
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TMath.h"

// FAIR headers
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRuntimeDb.h"

// TWIM headers
#include "R3BLogger.h"
#include "R3BTwimCalPar.h"
#include "R3BTwimGainMatching.h"
#include "R3BTwimMappedData.h"

// R3BTwimGainMatching::Default Constructor --------------------------
R3BTwimGainMatching::R3BTwimGainMatching()
    : R3BTwimGainMatching("R3BTwimGainMatching", 1)
{
}

// R3BTwimGainMatchingPar::Standard Constructor --------------------------
R3BTwimGainMatching::R3BTwimGainMatching(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNumSec(4)
    , fNumAnodes(16)
    , fNumAnodesRef(1)
    , fMinStadistics(1000)
    , fRefAnode(7)
    , fLimit_left(5000)
    , fLimit_right(90000)
    , fNumParams(2)
    , fTwimMappedDataCA(NULL)
    , fCal_Par(NULL)
{
}

// Virtual R3BTwimGainMatching::Destructor
R3BTwimGainMatching::~R3BTwimGainMatching() { R3BLOG(debug1, ""); }

// -----   Public method Init   --------------------------------------------
InitStatus R3BTwimGainMatching::Init()
{
    R3BLOG(info, "");
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    fTwimMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("TwimMappedData"));
    if (!fTwimMappedDataCA)
    {
        R3BLOG(fatal, "TwimMappedData not found");
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        R3BLOG(error, "FairRuntimeDb not found");
        return kFATAL;
    }

    fCal_Par = dynamic_cast<R3BTwimCalPar*>(rtdb->getContainer("twimCalPar"));
    if (!fCal_Par)
    {
        R3BLOG(fatal, "Couldn't get handle on twimCalPar container");
        return kFATAL;
    }

    Anode8_vs_anodes = new TGraph*[fNumSec * fNumAnodes];
    char Name1[255];
    for (Int_t s = 0; s < fNumSec; s++)
        for (Int_t i = 0; i < fNumAnodes; i++)
        {
            sprintf(Name1, "Section_%i_anode_%i_vs_anode_%i", s + 1, fRefAnode + 1, i + 1);
            Anode8_vs_anodes[s * fNumAnodes + i] = new TGraph(1);
            Anode8_vs_anodes[s * fNumAnodes + i]->SetName(Name1);
            Anode8_vs_anodes[s * fNumAnodes + i]->SetMarkerColor(4);
            Anode8_vs_anodes[s * fNumAnodes + i]->SetMarkerStyle(20);
            Anode8_vs_anodes[s * fNumAnodes + i]->SetMarkerSize(1.2);
            Anode8_vs_anodes[s * fNumAnodes + i]->SetLineWidth(0.);
        }

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BTwimGainMatching::ReInit() { return kSUCCESS; }

// -----   Public method Execution   --------------------------------------------
void R3BTwimGainMatching::Exec(Option_t* option)
{
    // Reading the Input -- Mapped Data --
    Int_t nHits = fTwimMappedDataCA->GetEntriesFast();
    R3BTwimMappedData** mappedData = new R3BTwimMappedData*[nHits];
    Int_t secId = 0;
    Int_t anodeId = 0;
    Double_t energyperanode[fNumSec][fNumAnodes];
    Int_t multanode[fNumSec][fNumAnodes];

    if (nHits > 0)
    {

        for (Int_t s = 0; s < fNumSec; s++)
        {
            for (Int_t i = 0; i < fNumAnodes; i++)
            {
                energyperanode[s][i] = 0.;
                multanode[s][i] = 0;
            }
        }
        // std::cout << "nHits = " << nHits << "\n";
        for (Int_t i = 0; i < nHits; i++)
        {
            mappedData[i] = dynamic_cast<R3BTwimMappedData*>(fTwimMappedDataCA->At(i));
            secId = mappedData[i]->GetSecID() - 1;
            anodeId = mappedData[i]->GetAnodeID() - 1;
            if (multanode[secId][anodeId] == 0 && mappedData[i]->GetPileupStatus() == 0 &&
                anodeId < fNumAnodes) // multi=1
            {
                energyperanode[secId][anodeId] = mappedData[i]->GetEnergy();
                multanode[secId][anodeId]++;
            }
        }

        for (Int_t s = 0; s < fNumSec; s++)
            for (Int_t i = 0; i < fNumAnodes; i++)
                if (energyperanode[s][i] > 0. && energyperanode[s][fRefAnode] > 0. &&
                    abs(energyperanode[s][i] - energyperanode[s][fRefAnode]) < 6000. &&
                    energyperanode[s][i] < fLimit_right && energyperanode[s][fRefAnode] < fLimit_right)
                {
                    Anode8_vs_anodes[s * fNumAnodes + i]->SetPoint(Anode8_vs_anodes[s * fNumAnodes + i]->GetN(),
                                                                   energyperanode[s][i],
                                                                   energyperanode[s][fRefAnode]);
                }
    }

    if (mappedData)
        delete[] mappedData;
    return;
}

void R3BTwimGainMatching::FinishTask()
{
    fCal_Par->SetNumSec(fNumSec);
    fCal_Par->SetNumAnodes(fNumAnodes);
    fCal_Par->SetNumParamsEFit(fNumParams);

    TF1* anodes_fit[fNumSec][fNumAnodes];
    char Name3[255];
    char Name4[255];

    canvas = new TCanvas*[fNumSec * fNumAnodes];

    for (Int_t s = 0; s < fNumSec; s++)
    {
        for (Int_t i = 0; i < fNumAnodes; i++)
            if (Anode8_vs_anodes[s * fNumAnodes + i]->GetN() > fMinStadistics)
            {
                sprintf(Name3, "FIT_Section_%i_anode_8_vs_anode_%i", s + 1, i + 1);
                anodes_fit[s][i] = new TF1(Name3, "pol1", fLimit_left, fLimit_right);
                sprintf(Name4, "Anode8_vs_anodes_%i_%i", s + 1, i + 1);
                canvas[s * fNumAnodes + i] = new TCanvas(Name4, Name4, 0, 400, 0, 400);
                canvas[s * fNumAnodes + i]->cd();
                Anode8_vs_anodes[s * fNumAnodes + i]->Fit(anodes_fit[s][i], "R");
                Double_t par1 = anodes_fit[s][i]->GetParameter(0);
                Double_t par2 = anodes_fit[s][i]->GetParameter(1);
                fCal_Par->SetAnodeCalParams(par1, s + 1, i + 1, 1);
                fCal_Par->SetAnodeCalParams(par2, s + 1, i + 1, 2);
                Anode8_vs_anodes[s * fNumAnodes + i]->Write();
            }
    }
    fCal_Par->setChanged();
}

ClassImp(R3BTwimGainMatching);
