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

#include "R3BNeulandCalTest.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BLosHit.h"
#include "R3BNeulandCalData.h"
#include "R3BNeulandPmt.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"

Double_t walk(Double_t x)
{
    Double_t y = 0;

    Double_t par1 = 1500.;                                                // +-0.2238
    Double_t par2 = 0.00075;                                              //+-2.355e-05
    y = par1 * TMath::Power(x, par2) - (par1 * TMath::Power(400., par2)); // Michael's

    // y=2.29083*log(x)-0.0870157*log(x)*log(x)-4.57824;  // mine

    return y;
    // return 0.;
}

R3BNeulandCalTest::R3BNeulandCalTest()
    : fnEvents(0)
    , fLandPmt(NULL)
    , fNeulandPmt(NULL)
    , fLosHit(NULL)
{
}

R3BNeulandCalTest::R3BNeulandCalTest(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fnEvents(0)
    , fLandPmt(NULL)
    , fNeulandPmt(NULL)
    , fLosHit(NULL)
{
}

R3BNeulandCalTest::~R3BNeulandCalTest() {}

InitStatus R3BNeulandCalTest::Init()
{
    FairRootManager* fMan = FairRootManager::Instance();
    if (!fMan)
    {
        LOG(fatal) << "FairRootManager not found";
    }
    fLandPmt = (TClonesArray*)fMan->GetObject("NeulandCalData");
    fNeulandPmt = (TClonesArray*)fMan->GetObject("NeulandPmt");
    fLosHit = (TClonesArray*)fMan->GetObject("LosHit");
    CreateHistos();

    return kSUCCESS;
}

void R3BNeulandCalTest::Exec(Option_t* option)
{
    Bool_t startSeen = kTRUE;
    Double_t startTime = 0.;
    if (fLosHit)
    {
        Int_t nLosHit = fLosHit->GetEntriesFast();
        Double_t stime[4];
        R3BLosHit* losHit;
        Int_t ch;
        Bool_t seen[] = { kFALSE, kFALSE, kFALSE, kFALSE };
        for (Int_t i = 0; i < nLosHit; i++)
        {
            losHit = (R3BLosHit*)fLosHit->At(i);
            ch = losHit->GetChannel();
            if (ch > 3)
            {
                continue;
            }
            else
            {
                seen[ch] = kTRUE;
                stime[ch] = losHit->GetTime();
            }
        }

        for (Int_t i = 0; i < 4; i++)
        {
            if (1 == i)
            {
                continue;
            }
            if (!seen[i])
            {
                startSeen = kFALSE;
                break;
            }
        }
        if (startSeen)
        {
            startTime = (stime[0] + stime[2] + stime[3]) / 3.;
        }
    }

    if (fLandPmt)
    {
        Int_t nLandPmt = fLandPmt->GetEntriesFast();
        R3BNeulandCalData* pmt1;
        R3BNeulandCalData* pmt2;
        Int_t barId;
        Double_t time1, time2;
        for (Int_t i1 = 0; i1 < nLandPmt; i1++)
        {
            pmt1 = (R3BNeulandCalData*)fLandPmt->At(i1);
            barId = pmt1->GetBarId();
            time1 = pmt1->GetTime();
            time1 += walk(pmt1->GetQdc());

            for (Int_t i2 = 0; i2 < nLandPmt; i2++)
            {
                if (i2 == i1)
                {
                    continue;
                }
                pmt2 = (R3BNeulandCalData*)fLandPmt->At(i2);
                if (barId != pmt2->GetBarId())
                {
                    continue;
                }
                time2 = pmt2->GetTime();
                time2 += walk(pmt2->GetQdc());

                if (startSeen && barId == 125)
                {
                    fh_los_corr->Fill((time1 + time2) / 2., startTime);
                    fh_tof->Fill((time1 + time2) / 2. - startTime);
                    fh_qdctof->Fill((time1 + time2) / 2. - startTime - 1563. + 30.,
                                    TMath::Sqrt(pmt1->GetQdc() * pmt2->GetQdc()));
                }
                if (startSeen && barId == 225)
                {
                    fh_qdctof->Fill((time1 + time2) / 2. - startTime - 1571. + 30.,
                                    TMath::Sqrt(pmt1->GetQdc() * pmt2->GetQdc()));
                }
            }
        }
    }
    if (fNeulandPmt)
    {
        Int_t nNeulandPmt = fNeulandPmt->GetEntriesFast();
        R3BNeulandPmt* pmt1;
        R3BNeulandPmt* pmt2;
        Int_t planeId;
        Int_t barId;
        Double_t time1, time2;
        for (Int_t i1 = 0; i1 < nNeulandPmt; i1++)
        {
            pmt1 = (R3BNeulandPmt*)fNeulandPmt->At(i1);
            planeId = pmt1->GetPlaneId();
            barId = pmt1->GetBarId();
            time1 = pmt1->GetTime();
            time1 += walk(pmt1->GetCharge());

            for (Int_t i2 = i1 + 1; i2 < nNeulandPmt; i2++)
            {
                if (i2 == i1)
                {
                    continue;
                }
                pmt2 = (R3BNeulandPmt*)fNeulandPmt->At(i2);
                if (barId != pmt2->GetBarId())
                {
                    continue;
                }
                time2 = pmt2->GetTime();
                time2 += walk(pmt2->GetCharge());

                fh_tdiff->Fill(time1 - time2);
                /*
                                if (startSeen && barId == 125)
                                {
                                    fh_los_corr->Fill((time1 + time2) / 2., startTime);
                                    fh_tof->Fill((time1 + time2) / 2. - startTime);
                                    fh_qdctof->Fill((time1 + time2) / 2. - startTime - 1563. + 30.,
                   TMath::Sqrt(pmt1->GetQdc()*pmt2->GetQdc()));
                                }
                                if (startSeen && barId == 225)
                                {
                                    fh_qdctof->Fill((time1 + time2) / 2. - startTime - 1571. + 30.,
                   TMath::Sqrt(pmt1->GetQdc()*pmt2->GetQdc()));
                                }
                */
            }
        }
    }

    fnEvents += 1;
}

void R3BNeulandCalTest::FinishTask() { WriteHistos(); }

void R3BNeulandCalTest::CreateHistos()
{
    fh_los_corr = new TH2F("h_los_corr", "LOS vs NeuLAND", 1000, 100, 400., 1000, -1360., -1310.);
    fh_tof = new TH1F("h_tof", "ToF", 1000, 1500., 1700.);
    fh_qdctof = new TH2F("h_qdctof", "QDC vs ToF", 1000, 0., 100., 200, 0., 2000.);
    fh_qdctof_2 = new TH2F("h_qdctof_2", "QDC vs ToF", 1000, 0., 100., 200, 0., 2000.);
    fh_tdiff = new TH1F("h_tdiff", "Tdiff", 10000, -50, 50.);
}

void R3BNeulandCalTest::WriteHistos()
{
    fh_los_corr->Write();
    fh_tof->Write();
    fh_qdctof->Write();
    fh_qdctof_2->Write();
    fh_tdiff->Write();
}

ClassImp(R3BNeulandCalTest)
