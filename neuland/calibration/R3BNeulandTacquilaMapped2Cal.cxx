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

// ------------------------------------------------------------
// -----               R3BNeulandTacquilaMapped2Cal               -----
// -----          Created 22-04-2014 by D.Kresan          -----
// ------------------------------------------------------------

#include "R3BNeulandTacquilaMapped2Cal.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "R3BEventHeader.h"
#include "R3BNeulandCalData.h"
#include "R3BNeulandQCalPar.h"
#include "R3BNeulandTacquilaMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include "TClonesArray.h"
#include "TH1F.h"
#include "TMath.h"

#define planes fNofPMTs / 100
#define toID(x, y, z) (((x - 1) * 50 + (y - 1)) * 2 + (z - 1))

Double_t wlk(Double_t x)
{
    Double_t y = 0;

    Double_t par1 = 1500.;                                                // +-0.2238
    Double_t par2 = 0.00075;                                              //+-2.355e-05
    y = par1 * TMath::Power(x, par2) - (par1 * TMath::Power(400., par2)); // Michael's

    // y=2.29083*log(x)-0.0870157*log(x)*log(x)-4.57824;  // mine

    return y;
    // return 0.;
}

R3BNeulandTacquilaMapped2Cal::R3BNeulandTacquilaMapped2Cal()
    : FairTask("NeulandTacquilaMapped2Cal", 1)
    , fNEvents(0)
    , fPulserMode(kFALSE)
    , fWalkEnabled(kTRUE)
    , fRawHit(NULL)
    , fPmt(new TClonesArray("R3BNeulandCalData"))
    , fNPmt(0)
    , fTcalPar(NULL)
    , fQCalPar(NULL)
    , fTrigger(-1)
    , fMap17Seen()
    , fMapStopTime()
    , fMapStopClock()
    , fMapQdcOffset()
    , fClockFreq(1. / TACQUILA_CLOCK_MHZ * 1000.)
{
}

R3BNeulandTacquilaMapped2Cal::R3BNeulandTacquilaMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fNEvents(0)
    , fPulserMode(kFALSE)
    , fWalkEnabled(kTRUE)
    , fRawHit(NULL)
    , fPmt(new TClonesArray("R3BNeulandCalData"))
    , fNPmt(0)
    , fTcalPar(NULL)
    , fQCalPar(NULL)
    , fTrigger(-1)
    , fMap17Seen()
    , fMapStopTime()
    , fMapStopClock()
    , fMapQdcOffset()
    , fClockFreq(1. / TACQUILA_CLOCK_MHZ * 1000.)
{
}

R3BNeulandTacquilaMapped2Cal::~R3BNeulandTacquilaMapped2Cal()
{
    if (fPmt)
    {
        delete fPmt;
        fPmt = NULL;
        fNPmt = 0;
    }
}

InitStatus R3BNeulandTacquilaMapped2Cal::Init()
{
    LOG(info) << "R3BNeulandTacquilaMapped2Cal::Init : read " << fTcalPar->GetNumModulePar() << " calibrated modules";
    // fTcalPar->printParams();

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        LOG(fatal) << "FairRootManager not found";
    }

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (NULL == header)
    {
        LOG(fatal) << "Branch R3BEventHeader not found";
    }

    fRawHit = dynamic_cast<TClonesArray*>(mgr->GetObject("NeulandTacquilaMappedData"));
    if (NULL == fRawHit)
    {
        LOG(fatal) << "Branch NeulandTacquilaMappedData not found";
    }

    mgr->Register("NeulandCalData", "Neuland", fPmt, kTRUE);

    fh_pulser_5_2 = new TH1F("h_pulser_5_2", "Single PMT resolution Bar 5 vs 2", 40000, -200., 200.);
    fh_pulser_105_2 = new TH1F("h_pulser_105_2", "Single PMT resolution Bar 105 vs 2", 40000, -200., 200.);

    SetParameter();

    return kSUCCESS;
}

void R3BNeulandTacquilaMapped2Cal::SetParContainers()
{
    FairRunAna* ana = FairRunAna::Instance();
    FairRuntimeDb* rtdb = ana->GetRuntimeDb();
    fTcalPar = dynamic_cast<R3BTCalPar*>(rtdb->getContainer("LandTCalPar"));
    fQCalPar = dynamic_cast<R3BNeulandQCalPar*>(rtdb->getContainer("NeulandQCalPar"));
}

void R3BNeulandTacquilaMapped2Cal::SetParameter()
{

    std::map<Int_t, Double_t> tempMapQdcOffset;
    Int_t i = 0;
    for (Int_t plane = 1; i <= planes; plane++)
        for (Int_t bar = 1; bar <= 50; bar++)
            for (Int_t side = 1; side <= 2; side++)
            {
                tempMapQdcOffset[i] = fQCalPar->GetParAt(plane, bar, side);
                i++;
            }

    LOG(info) << "R3BNeulandTacquilaMapped2Cal::SetParameter : Number of Parameters: " << i;

    fMapQdcOffset = tempMapQdcOffset;
}

InitStatus R3BNeulandTacquilaMapped2Cal::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BNeulandTacquilaMapped2Cal::Exec(Option_t* option)
{
    if (fTrigger >= 0)
    {
        if (header->GetTrigger() != fTrigger)
        {
            return;
        }
    }

    Int_t nHits = fRawHit->GetEntriesFast();
    if (fPulserMode)
    {
        if (nHits < fNofPMTs)
        {
            return;
        }
    }
    else
    {
        if (nHits > (fNofPMTs / 2))
        {
            return;
        }
    }

    if (nHits > 0)
    {
        MakeCal();
    }

    if (fPulserMode)
    {
        R3BNeulandCalData* pmt1;
        Double_t time1 = nan("");
        for (Int_t i = 0; i < fNPmt; i++)
        {
            pmt1 = dynamic_cast<R3BNeulandCalData*>(fPmt->At(i));
            if (pmt1->GetBarId() == 2 && pmt1->GetSide() == 1)
            {
                time1 = pmt1->GetTime();
                break;
            }
        }
        for (Int_t i = 0; i < fNPmt; i++)
        {
            pmt1 = dynamic_cast<R3BNeulandCalData*>(fPmt->At(i));
            if (pmt1->GetBarId() == 5 && pmt1->GetSide() == 1)
            {
                fh_pulser_5_2->Fill(pmt1->GetTime() - time1);
            }
            if (pmt1->GetBarId() == 105 && pmt1->GetSide() == 1)
            {
                fh_pulser_105_2->Fill(pmt1->GetTime() - time1);
            }
        }
    }
}

void R3BNeulandTacquilaMapped2Cal::MakeCal()
{
    Int_t nHits = fRawHit->GetEntriesFast();
    R3BNeulandTacquilaMappedData* hit;
    R3BNeulandTacquilaMappedData* hit2;
    Int_t iPlane;
    Int_t iPaddle;
    Int_t iSide;
    Int_t channel;
    Int_t tdc;
    R3BTCalModulePar* par;
    Double_t time;
    Double_t time2;
    Int_t qdc;

    for (Int_t khit = 0; khit < nHits; khit++)
    {
        hit2 = dynamic_cast<R3BNeulandTacquilaMappedData*>(fRawHit->At(khit));
        if (NULL == hit2)
        {
            continue;
        }

        iPlane = hit2->GetPlane();
        iPaddle = hit2->GetPaddle();
        iSide = hit2->GetSide();

        if (!(par = fTcalPar->GetModuleParAt(iPlane, iPaddle, iSide)))
        {
            LOG(debug) << "R3BNeulandTacquilaMapped2Cal::Exec : Tcal par not found, channel: " << iPlane << " / "
                       << iPaddle << " / " << iSide;
            continue;
        }

        tdc = hit2->GetTacData();
        time = par->GetTimeTacquila(tdc);
        if (time < 0. || time > fClockFreq)
        {
            LOG(error) << "R3BNeulandTacquilaMapped2Cal::Exec : error in time calibration: ch=" << channel
                       << ", tdc=" << tdc << ", time=" << time;
            continue;
        }

        if (!(par = fTcalPar->GetModuleParAt(iPlane, iPaddle, iSide + 2)))
        {
            LOG(debug) << "R3BNeulandTacquilaMapped2Cal::Exec : Tcal par not found, channel: " << iPlane << " / "
                       << iPaddle << " / " << (iSide + 2);
            continue;
        }

        tdc = hit2->GetStopT();
        time2 = par->GetTimeTacquila(tdc);
        if (time2 < 0. || time2 > fClockFreq)
        {
            LOG(error) << "R3BNeulandTacquilaMapped2Cal::Exec : error in time calibration: ch=" << channel
                       << ", tdc=" << tdc << ", time=" << time2;
            continue;
        }

        qdc = hit2->GetQdcData() - fMapQdcOffset[toID(iPlane, iPaddle, iSide)];
        qdc = std::max(qdc, 0);

        time = time - time2 + hit2->GetClock() * fClockFreq;
        if (fWalkEnabled)
        {
            time += wlk(qdc);
        }
        new ((*fPmt)[fNPmt]) R3BNeulandCalData((iPlane - 1) * 50 + iPaddle, iSide, time, NAN, qdc);
        fNPmt += 1;
    }
}

void R3BNeulandTacquilaMapped2Cal::FinishEvent()
{
    if (fVerbose && 0 == (fNEvents % 1000))
    {
        LOG(info) << "R3BNeulandTacquilaMapped2Cal::Exec : event=" << fNEvents << " nPMTs=" << fNPmt;
    }

    if (fPmt)
    {
        fPmt->Clear();
        fNPmt = 0;
    }
    fMap17Seen.clear();
    fMapStopTime.clear();
    fMapStopClock.clear();

    fNEvents += 1;
}

void R3BNeulandTacquilaMapped2Cal::FinishTask()
{
    fh_pulser_5_2->Write();
    fh_pulser_105_2->Write();
}

ClassImp(R3BNeulandTacquilaMapped2Cal)
