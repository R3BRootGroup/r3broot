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

// -----------------------------------------------------------------------------
// -----                            R3BNeulandHitHist                      -----
// -----                     Created 22-04-2014 by D.Kresan                -----
// -----------------------------------------------------------------------------

#include "R3BNeulandHitHist.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BLosHitData.h"
#include "R3BNeulandHit.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

const Double_t C_LIGHT = 29.9792458; // Speed of light [cm/ns]

R3BNeulandHitHist::R3BNeulandHitHist()
    : fNofBars(100)
    , fFirstPlaneHorisontal(kFALSE)
    , fnEvents(0)
    , fLandDigi(NULL)
    , fLosHit(NULL)
{
}

R3BNeulandHitHist::R3BNeulandHitHist(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNofBars(100)
    , fFirstPlaneHorisontal(kFALSE)
    , fnEvents(0)
    , fLandDigi(NULL)
    , fLosHit(NULL)
{
}

R3BNeulandHitHist::~R3BNeulandHitHist() {}

InitStatus R3BNeulandHitHist::Init()
{
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan)
    {
        LOG(fatal) << "FairRootManager not found";
    }
    fLandDigi = dynamic_cast<TClonesArray*>(fMan->GetObject("NeulandHits"));
    fLosHit = dynamic_cast<TClonesArray*>(fMan->GetObject("LosHit"));
    CreateHistos();

    return kSUCCESS;
}

void R3BNeulandHitHist::Exec(Option_t* option)
{
    if (!fLosHit || fLosHit->GetEntriesFast() == 0)
        return;

    R3BLosHitData* losHit = dynamic_cast<R3BLosHitData*>(fLosHit->At(0));
    // Double_t startTime = losHit->GetTime(); TODO
    Double_t startTime = 0;
    Int_t count = 0;
    Double_t tmean = 0.;

    if (fLandDigi && fLandDigi->GetEntriesFast() > 0)
    {
        Int_t nLandDigi = fLandDigi->GetEntriesFast();
        R3BNeulandHit* digi;
        Int_t barId;
        Double_t time, tof, norm_tof;
        Double_t qdc, qdcL, qdcR;
        Double_t x, y, z;
        Double_t l;
        Double_t beta;

        for (Int_t i = 0; i < nLandDigi; i++)
        {
            digi = dynamic_cast<R3BNeulandHit*>(fLandDigi->At(i));

            barId = digi->GetPaddle();
            time = digi->GetT();
            qdc = digi->GetE();
            qdcL = digi->GetQdcL();
            qdcR = digi->GetQdcR();
            x = digi->GetPosition().X();
            y = digi->GetPosition().Y();
            z = digi->GetPosition().Z();
            l = TMath::Sqrt(x * x + y * y + z * z);

            tmean += time;
            count += 1;

            fh_land_qdcbarid->Fill(barId, qdc);
            fh_land_yx->Fill(x, y);
            if (barId < 51)
            {
                fh_land_yx1->Fill(x, y);
            }
            else if (barId < 101)
            {
                fh_land_yx2->Fill(x, y);
            }

            fh_los_corr->Fill(time, startTime);

            fh_land_barid->Fill(barId);
            fh_land_timebarid->Fill(barId, time);
            tof = time - startTime - fTimeOffset;
            beta = l / tof;
            norm_tof = tof - (l - 1300.) / beta;
            fh_land_qdctof->Fill(norm_tof, qdc);
            fh_land_tofbarid->Fill(barId, tof);
            fh_land_lbarid->Fill(barId, l);
            fh_land_ltime->Fill(tof, l);
            fh_land_betabarid->Fill(barId, beta);
            fh_land_tof->Fill(tof);
            if (beta > 29.7 && beta < 30.2)
                fh_land_norm_tof->Fill(norm_tof);

            fh_land_beta->Fill(beta);
            fh_land_qdc->Fill(qdc);
            if (beta > 16 && beta < 26)
                fh_land_qdc_cut->Fill(qdc);
        }

        fnEvents += 1;
        if (0 == (fnEvents % 100000))
        {
            LOG(info) << "R3BNeulandHitHist : " << fnEvents << " events collected,  start time=" << startTime;
            if (count)
            {
                LOG(info) << ",  mean time=" << (tmean / (Double_t)count);
            }
            LOG(info);
        }
    }
}

void R3BNeulandHitHist::FinishTask() { WriteHistos(); }

void R3BNeulandHitHist::CreateHistos()
{
    fh_land_barid = new TH1F("h_land_barid", "Bar ID", fNofBars, (Double_t)fNofBars + 0.5, 100.5);
    fh_land_qdcbarid =
        new TH2F("h_land_qdcbarid", "QDC vs Bar ID", fNofBars, 0.5, (Double_t)fNofBars + 0.5, 2000, 0., 40.);
    fh_land_timebarid =
        new TH2F("h_land_timebarid", "Time vs Bar ID", fNofBars, 0.5, (Double_t)fNofBars + 0.5, 4000, 0., 2000.);
    fh_land_tofbarid =
        new TH2F("h_land_tofbarid", "ToF vs Bar ID", fNofBars, 0.5, (Double_t)fNofBars + 0.5, 1000, 0., 100.);
    fh_land_tof = new TH1F("h_land_tof", "ToF", 1000, 0., 100.);
    fh_land_qdctof = new TH2F("h_land_qdctof", "QDC vs ToF", 2000, -50., 150., 2000, 0., 40.);
    fh_land_betabarid =
        new TH2F("h_land_betabarid", "Velocity vs Bar ID", fNofBars, 0.5, (Double_t)fNofBars + 0.5, 500, 0., 50.);
    fh_land_lbarid =
        new TH2F("h_land_lbarid", "Length vs Bar ID", fNofBars, 0.5, (Double_t)fNofBars + 0.5, 1000, 700., 900.);
    fh_land_ltime = new TH2F("h_land_ltime", "Length vs Bar ID", 4000, 0., 2000., 1000, 700., 900.);
    fh_land_yx = new TH2F("h_land_yx", "Y vs X", 340, -170., 170., 340, -170., 170.);
    if (fFirstPlaneHorisontal)
    {
        fh_land_yx1 = new TH2F("h_land_yx1", "Y vs X", 340, -170., 170., 68, -170., 170.);
        fh_land_yx2 = new TH2F("h_land_yx2", "Y vs X", 68, -170., 170., 340, -170., 170.);
    }
    else
    {
        fh_land_yx1 = new TH2F("h_land_yx1", "Y vs X", 68, -170., 170., 340, -170., 170.);
        fh_land_yx2 = new TH2F("h_land_yx2", "Y vs X", 340, -170., 170., 68, -170., 170.);
    }
    fh_land_beta = new TH1F("h_land_beta", "Velocity", 4500, -10., 35.);
    fh_land_qdc = new TH1F("h_land_qdc", "QDC", 2000, 0., 40.);
    fh_land_qdc_cut = new TH1F("h_land_qdc_cut", "QDC cut", 2000, 0., 40.);
    fh_los_time = new TH1F("h_los_time", "LOS Time", 2000, -1000., -800.);
    fh_los_corr = new TH2F("h_los_corr", "LOS vs NeuLAND", 2000, 700., 900., 2000, -1000., -800.);

    fh_land_norm_tof = new TH1F("h_land_norm_tof", "norm. ToF", 1000, 40., 50.);
}

void R3BNeulandHitHist::WriteHistos()
{
    fh_land_barid->Write();
    fh_land_qdcbarid->Write();
    fh_land_tof->Write();
    fh_land_qdctof->Write();
    fh_land_timebarid->Write();
    fh_land_tofbarid->Write();
    fh_land_yx->Write();
    fh_land_yx1->Write();
    fh_land_yx2->Write();
    fh_land_lbarid->Write();
    fh_land_ltime->Write();
    fh_land_beta->Write();
    fh_land_qdc->Write();
    fh_land_qdc_cut->Write();
    fh_land_betabarid->Write();
    fh_los_time->Write();
    fh_los_corr->Write();
    fh_land_norm_tof->Write();
}

ClassImp(R3BNeulandHitHist)
