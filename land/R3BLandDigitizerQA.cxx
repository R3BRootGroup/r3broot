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

#include <stdio.h>

#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

#include "FairLogger.h"
#include "FairRootManager.h"

#include "R3BLandDigi.h"
#include "R3BLandDigitizerQA.h"
#include "R3BLandPoint.h"
#include "R3BMCTrack.h"

R3BLandDigitizerQA::R3BLandDigitizerQA()
    : FairTask("R3BLandDigitizerQA")
    , fEvent(0)
    , fTracks(NULL)
    , fPoints(NULL)
    , fDigis(NULL)
{
}

R3BLandDigitizerQA::R3BLandDigitizerQA(const R3BLandDigitizerQA& right)
    : FairTask("R3BLandDigitizerQA")
    , fEvent(right.fEvent)
    , fTracks(NULL)
    , fPoints(NULL)
    , fDigis(NULL)
{
}

R3BLandDigitizerQA::~R3BLandDigitizerQA() {}

InitStatus R3BLandDigitizerQA::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (NULL == rm)
    {
        LOG(error) << "R3BLandDigitizerQA: no pointer to FairRootManager...";
        return kERROR;
    }

    fTracks = dynamic_cast<TClonesArray*>(rm->GetObject("MCTrack"));
    if (NULL == fTracks)
    {
        LOG(error) << "R3BLandDigitizerQA: no MCTrack array...";
        return kERROR;
    }

    fPoints = dynamic_cast<TClonesArray*>(rm->GetObject("LandPoint"));
    if (NULL == fPoints)
    {
        LOG(error) << "R3BLandDigitizerQA: no LandPoint array...";
        return kERROR;
    }

    fDigis = dynamic_cast<TClonesArray*>(rm->GetObject("LandDigi"));
    if (NULL == fDigis)
    {
        LOG(error) << "R3BLandDigitizerQA: no LandDigi array...";
        return kERROR;
    }

    CreateHisto();

    return kSUCCESS;
}

void R3BLandDigitizerQA::Exec(Option_t* option)
{
    if (NULL == fTracks || NULL == fPoints || NULL == fDigis)
    {
        return;
    }

    Double_t totEnergy = 0.;
    Double_t totEnergyLee = 0.;
    Double_t totEnergyReco = 0.;

    Int_t nTracks = fTracks->GetEntriesFast();

    Int_t nPoints = fPoints->GetEntriesFast();
    R3BLandPoint* point;
    Double_t eloss;
    Int_t media;
    Double_t light;
    Double_t time;
    Int_t pdg;
    for (Int_t i = 0; i < nPoints; i++)
    {
        point = dynamic_cast<R3BLandPoint*>(fPoints->At(i));
        eloss = point->GetEnergyLoss() * 1000.;
        media = Int_t(point->GetPaddleType());
        light = point->GetLightYield() * 1000.;
        if (eloss > 0. && 3 == media)
        {
            totEnergy += eloss;
            totEnergyLee += light;
            time = point->GetTime();
            pdg = (dynamic_cast<R3BMCTrack*>(fTracks->At(point->GetTrackID())))->GetPdgCode();
            fhElossPdg->Fill(pdg, eloss);
        }
    }

    Int_t nDigis = fDigis->GetEntriesFast();
    R3BLandDigi* digi;
    Double_t qdc;
    Double_t tdc = nan("");
    Double_t x, y, z;
    Double_t beta;
    Double_t tdc_first = 1e10;
    Int_t first = -1;
    Double_t beta_first = 1e10;
    Int_t bfirst = -1;
    Double_t qdc_first;
    for (Int_t i = 0; i < nDigis; i++)
    {
        digi = dynamic_cast<R3BLandDigi*>(fDigis->At(i));
        qdc = digi->GetQdc();
        fhPaddleE->Fill(qdc);
        fhElossTime->Fill(tdc, qdc);
        if (qdc > 0.)
        {
            totEnergyReco += qdc;
            tdc = digi->GetTdc();
            x = digi->GetXX();
            y = digi->GetYY();
            z = digi->GetZZ();
            fhQdcXL->Fill(y, digi->GetQdcL());
            fhQdcXR->Fill(y, digi->GetQdcR());
            fhThresh->Fill(digi->GetPaddleNr(), qdc);
            fhQdcRatio->Fill(digi->GetQdcL() / digi->GetQdcR());
            beta = TMath::Sqrt(x * x + y * y + z * z) / tdc / (29.9792458);
            if (tdc < tdc_first)
            {
                tdc_first = tdc;
                first = digi->GetPaddleNr();
                qdc_first = qdc;
            }
            if (beta < beta_first)
            {
                beta_first = beta;
                bfirst = digi->GetPaddleNr();
            }
        }
    }
    if (first > 0)
    {
        fhFirstTvsBar->Fill((Double_t)first);
        fhFirstE->Fill(qdc_first);
    }
    if (bfirst > 0)
    {
        fhFirstBvsBar->Fill(bfirst);
    }

    if (totEnergy > 0.)
    {
        fhTotEnergy->Fill(totEnergy);
        fhTotLight->Fill(totEnergyLee);
        fhTotEnergyRel->Fill(totEnergy / fBeamEnergy);
        fhTotLightRel->Fill(totEnergyLee / fBeamEnergy);
        fhElossLight->Fill(totEnergy, totEnergyLee);
        fhTotEnergyReco->Fill(totEnergyReco);
    }

    if (0 == (fEvent % 100))
    {
        char str1[100];
        char str2[100];
        char str3[100];
        char str4[100];
        sprintf(str1, "%8d", fEvent);
        sprintf(str2, "%4d", nTracks);
        sprintf(str3, "%4d", nPoints);
        sprintf(str4, "%4d", nDigis);
        LOG(info) << "R3BLandDigitizerQA: event # " << str1 << ", MCTrack: " << str2 << ", LandPoint: " << str3
                  << ", LandDigi: " << str4;
    }
}

void R3BLandDigitizerQA::FinishEvent() { fEvent += 1; }

void R3BLandDigitizerQA::Finish()
{
    fhTotEnergy->Write();
    fhTotLight->Write();
    fhTotEnergyRel->Write();
    fhTotLightRel->Write();
    fhElossLight->Write();
    fhTotEnergyReco->Write();
    fhPaddleE->Write();
    fhFirstTvsBar->Write();
    fhFirstBvsBar->Write();
    fhFirstE->Write();
    fhElossTime->Write();
    fhElossPdg->Write();
    fhQdcXL->Write();
    fhQdcXR->Write();
    fhThresh->Write();
    fhQdcRatio->Write();
}

void R3BLandDigitizerQA::CreateHisto()
{
    fhTotEnergy = new TH1F("Total_Energy", "Total energy detected", 5000, 0., 5000.);
    fhTotLight = new TH1F("Total_Light", "Total light detected (energy equivalent)", 5000, 0., 5000.);
    fhTotEnergyRel = new TH1F("Total_Energy_rel", "Relative Total energy detected", 55, -0.02, 1.08);
    fhTotLightRel = new TH1F("Total_Light_rel", "Relative Total light detected (energy equivalent)", 55, -0.02, 1.08);
    fhElossLight = new TH2F("ElossLight", "Light quenching for protons", 5000, 0., 5000., 5000, 0., 5000.);
    fhTotEnergyReco = new TH1F("Total_Energy_Reco", "Total energy reconstructed", 5000, 0., 5000.);
    //  fhPaddleE = new TH1F("PaddleEnergy", "", 506, 0., 101.2);
    fhPaddleE = new TH1F("PaddleEnergy", "", 1012, 0., 202.4);
    fhFirstTvsBar = new TH1F("FirstTvsBar", "", 155, 0., 155.);
    fhFirstBvsBar = new TH1F("FirstBvsBar", "", 155, 0., 155.);
    fhFirstE = new TH1F("FirstE", "", 200, 0., 200.);
    fhElossTime = new TH2F("ElossTime", "", 400, 0., 200., 2000, 0., 200.);
    fhElossPdg = new TH2F("ElossPdg", "", 30000, 0., 30000., 2000, 0., 200.);
    fhQdcXL = new TH2F("QdcXL", "", 300, -150., 150., 100, 0., 100.);
    fhQdcXR = new TH2F("QdcXR", "", 300, -150., 150., 100, 0., 100.);
    fhThresh = new TH2F("Thresh", "", 150, 0., 150., 500, 0., 100.);
    fhQdcRatio = new TH1F("QdcRatio", "", 100, -5., 5.);
}

ClassImp(R3BLandDigitizerQA)
