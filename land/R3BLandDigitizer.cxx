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

#include <fstream>
#include <iostream>
#include <string>

#include "TClonesArray.h"
#include "TF1.h"
#include "TGeoManager.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TVector3.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "R3BLandDigitizer.h"
#include "R3BLandPoint.h"

const Double_t R3BLandDigitizer::DEFAULT_SATURATION_COEFFICIENT = 0.012;

double fun1(double* x, double* par)
{
    double sigma = par[1];
    if (x[0] >= 14.)
    {
        sigma *= 4.;
    }
    return TMath::Exp(-0.5 * TMath::Power((x[0] - par[0]) / sigma, 2));
}

Double_t R3BLandDigitizer::BuildTOFRangeFromBeamEnergy(const Double_t& e) // [ns]
{
    if (e > 1100)
    {
        return 15.;
    }
    if (e < 350)
    {
        return 13.;
    }
    return 11. + (e / 498. - 1.) * 0.9;
}

// ----------------------------------------------------------------------------
R3BLandDigitizer::R3BLandDigitizer()
    : FairTask("R3B Land Digitization scheme")
    , f1(new TF1("f1", fun1, 0., 1000., 2))
    , fRnd(new TRandom3())
    , fThreshFileName("")
    , fNChannels(0)
    , fSaturationCoefficient(R3BLandDigitizer::DEFAULT_SATURATION_COEFFICIENT)
{
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
R3BLandDigitizer::R3BLandDigitizer(Int_t verbose)
    : FairTask("R3B Land Digitization scheme ", verbose)
    , f1(new TF1("f1", fun1, 0., 1000., 2))
    , fRnd(new TRandom3())
    , fThreshFileName("")
    , fNChannels(0)
    , fSaturationCoefficient(R3BLandDigitizer::DEFAULT_SATURATION_COEFFICIENT)
{
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
R3BLandDigitizer::~R3BLandDigitizer() {}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void R3BLandDigitizer::SetParContainers()
{
    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fLandDigiPar = dynamic_cast<R3BLandDigiPar*>(rtdb->getContainer("R3BLandDigiPar"));

    if (fVerbose && fLandDigiPar)
    {
        LOG(info) << "R3BLandDigitizer::SetParContainers() ";
        LOG(info) << "Container R3BLandDigiPar loaded ";
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
InitStatus R3BLandDigitizer::Init()
{
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fLandPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("LandPoint"));

    // Register output array LandDigi
    fLandDigi = new TClonesArray("R3BLandDigi", 1000);
    ioman->Register("LandDigi", "Digital response in Land", fLandDigi, kTRUE);

    // Parameter retrieval
    // Only after Init one retrieve the Digitization Parameters!
    npaddles = fLandDigiPar->GetMaxPaddle() + 1;
    nplanes = fLandDigiPar->GetMaxPlane();
    paddle_per_plane = (npaddles - 1) / nplanes;
    LOG(info) << "R3BLandDigitizer: # paddles: " << npaddles - 1 << ", # planes: " << nplanes;
    LOG(info) << "R3BLandDigitizer: paddles per plane: " << paddle_per_plane;
    c = 2.99792458E8;
    cMedia = 14.; // speed of light in material in [cm/ns]
                  //  cMedia = c / 1.58 * 1e9; // speed of light in material in [cm/ns]
    eventNo = 0;

    PM_res = new PM_RES*[npaddles];
    for (int i = 0; i < npaddles; i++)
    {
        PM_res[i] = new PM_RES[1500];
    }

    // Initialise control histograms
    hPMl = new TH1F("PM_left", "Arrival times of left PM", 1000, 0., 1000.);
    hPMr = new TH1F("PM_right", "Arrival times of right PM", 1000, 0., 1000.);
    Int_t nbins = 100;
    Float_t nmax = 100.;
    if (!fLandDigiPar->GetGeometryFileName().Contains("proto") &&
        !fLandDigiPar->GetGeometryFileName().Contains("jap") && npaddles < 202)
    {
        // LAND detector
        nbins = 20;
        nmax = 20.;
    }
    hMult1 = new TH1F("Multiplicity1", "Paddle multiplicity", nbins, 0., nmax);
    hMult2 = new TH1F("Multiplicity2", "Paddle multiplicity", nbins, 0., nmax);

    if (fNChannels > 0 && fThreshFileName.Length() > 0)
    {
        std::ifstream* tf = new std::ifstream(fThreshFileName.Data());
        for (Int_t i = 0; i < fNChannels; i++)
        {
            Float_t tl;
            Float_t tr;
            (*tf) >> tl >> tr;
            if (tl < 1e-6)
            {
                tl = 1e10;
            }
            if (tr < 1e-6)
            {
                tr = 1e10;
            }
            threshL[i] = tl;
            threshR[i] = tr;
        }
        tf->close();
    }
    else
    {
        for (Int_t i = 0; i < 5000; i++)
        {
            threshL[i] = 0.01;
            threshR[i] = 0.01;
        }
    }

    // If integration time has not been set otherwise, set it from beam energy - mimicks previous behavior
    if (fTOFRange < 0.01)
    {
        fTOFRange = BuildTOFRangeFromBeamEnergy(fBeamEnergy);
    }
    hRLTimeToTrig = new TH1F("hRLTimeToTrig", "R/Ltime-triggerTime", 200, -100, 100);

    return kSUCCESS;
}
// -------------------------------------------------------------------------

// -----   Public method Exec   --------------------------------------------
void R3BLandDigitizer::Exec(Option_t* opt)
{
    // Reset entries in output arrays
    Reset();

    eventNo += 1;

    Double_t timeRes = 0.15; // ns
    if (!fLandDigiPar->GetGeometryFileName().Contains("proto") &&
        !fLandDigiPar->GetGeometryFileName().Contains("jap") && npaddles < 202)
    {
        // LAND detector
        timeRes = 0.25; // ns
    }
    if (fLandDigiPar->GetGeometryFileName().Contains("jap"))
    {
        timeRes = 0.010;
    }

    // light attenuation of plastic scintillator
    att = 0.008; // [1/cm]
    // att = 0.0047619047619; // [1/cm]
    // decay constant of BC408 scintillator
    Double_t lambda = 1. / 2.1; // [1/ns]

    // half of the length of a scintillator
    plength = fLandDigiPar->GetPaddleLength(); // [cm]

    Double_t Thresh = 0.16;
    //  Double_t Thresh = 0.05;

    Int_t nentries = fLandPoints->GetEntries();
    Int_t PMmult[npaddles];
    Double_t QDC_temp[npaddles], TDC_temp[npaddles];
    Int_t mult1 = 0;
    Int_t mult2 = 0;
    Double_t xpos_temp[npaddles], ypos_temp[npaddles], zpos_temp[npaddles];
    Double_t xpaddle[npaddles], ypaddle[npaddles], zpaddle[npaddles];

    // reset
    for (Int_t j = 0; j < npaddles; j++)
    {
        PMmult[j] = 0;
        memset(PM_res[j], 0, 1500 * sizeof(PM_RES));
    }

    Double_t edep = 0.;

    for (Int_t l = 0; l < nentries; l++)
    {
        R3BLandPoint* land_obj = dynamic_cast<R3BLandPoint*>(fLandPoints->At(l));

        Int_t paddle = int(land_obj->GetSector()) - 1; // note that paddle starts at 1
        Int_t scint = int(land_obj->GetPaddleNb());
        Double_t eloss = land_obj->GetEnergyLoss() * 1000.;
        Double_t light = land_obj->GetLightYield() * 1000.;
        Double_t x = land_obj->GetXIn();
        Double_t y = land_obj->GetYIn();
        Double_t z = land_obj->GetZIn();
        Double_t time = land_obj->GetTime();
        Int_t media = int(land_obj->GetPaddleType());

        if (eloss > 0. && media == 3)
        {
            PMmult[paddle] = PMmult[paddle] + 1;

            edep += eloss;

            Int_t m = PMmult[paddle];

            if (m > 1500)
            {
                LOG(error) << "R3BLandDigitizer: Too many hits in one paddle: " << m << " hits"
                           << "entry " << l << " after sorting"
                           << "paddle " << paddle << "padlle mult " << PMmult[paddle] << "Sector " << scint << "eloss "
                           << eloss << "light yield " << light << "time " << time << "X In " << x << "Y In " << y
                           << "Z In " << z;
                return;
            }

            gGeoManager->FindNode(x, y, z);
            gGeoManager->CdUp();
            Double_t local_point[] = { 0., 0., 0. };
            Double_t global_point[3];
            gGeoManager->LocalToMaster(local_point, global_point);
            xpaddle[paddle] = global_point[0];
            ypaddle[paddle] = global_point[1];
            zpaddle[paddle] = global_point[2];

            if (fLandDigiPar->GetGeometryFileName().Contains("proto"))
            {
                // vertical paddles
                PM_res[paddle][m].Ltime = time + (plength - y) / cMedia;
                PM_res[paddle][m].LlightCFD = light * exp(-att * (plength - y));
                PM_res[paddle][m].LlightQDC = light * exp(-att * (plength - y));

                PM_res[paddle][m].Rtime = time + (plength + y) / cMedia;
                PM_res[paddle][m].RlightCFD = light * exp(-att * (plength + y));
                PM_res[paddle][m].RlightQDC = light * exp(-att * (plength + y));
            }
            else
            {
                if ((int)(((paddle - 1) / paddle_per_plane)) / 2. ==
                    (int)((int)(((paddle - 1) / paddle_per_plane)) / 2.))
                {
                    // horizontal paddles
                    PM_res[paddle][m].Ltime = time + (plength - x) / cMedia;
                    PM_res[paddle][m].LlightCFD = light * exp(-att * (plength - x));
                    PM_res[paddle][m].LlightQDC = light * exp(-att * (plength - x));

                    PM_res[paddle][m].Rtime = time + (plength + x) / cMedia;
                    PM_res[paddle][m].RlightCFD = light * exp(-att * (plength + x));
                    PM_res[paddle][m].RlightQDC = light * exp(-att * (plength + x));
                }
                else
                {
                    // vertical paddles
                    PM_res[paddle][m].Ltime = time + (plength - y) / cMedia;
                    PM_res[paddle][m].LlightCFD = light * exp(-att * (plength - y));
                    PM_res[paddle][m].LlightQDC = light * exp(-att * (plength - y));

                    PM_res[paddle][m].Rtime = time + (plength + y) / cMedia;
                    PM_res[paddle][m].RlightCFD = light * exp(-att * (plength + y));
                    PM_res[paddle][m].RlightQDC = light * exp(-att * (plength + y));
                }
            }
        } //! eloss
    }     //! MC hits

    Double_t temp[1500][3];

    /*  for (Int_t i=0;i<npaddles;i++){
        //cout<<"loop over paddles "<<i<<endl;

        // sort hits according to time; left PM
        for (Int_t j=0;j<PMmult[i];j++){
          temp[j][0]=PM_res[i][j+1].Ltime;
          temp[j][1]=PM_res[i][j+1].LlightCFD;
          temp[j][2]=PM_res[i][j+1].LlightQDC;
        }

        for (Int_t j=0;j<PMmult[i];j++){
          Double_t min = 100000;
          Int_t index = 0;
          for (Int_t k=0;k<PMmult[i];k++){
            if(temp[k][0] < min){
              min = temp[k][0];
              index = k;
            }
          }
          PM_res[i][j+1].Ltime=temp[index][0];
          PM_res[i][j+1].LlightCFD=temp[index][1];
          PM_res[i][j+1].LlightQDC=temp[index][2];
          temp[index][0] = 100000.;
        }

        // sort hits according to time; right PM
        for (Int_t j=0;j<PMmult[i];j++){
          temp[j][0]=PM_res[i][j+1].Rtime;
          temp[j][1]=PM_res[i][j+1].RlightCFD;
          temp[j][2]=PM_res[i][j+1].RlightQDC;
        }

        for (Int_t j=0;j<PMmult[i];j++){
          Double_t min = 100000;
          Int_t index = 0;
          for (Int_t k=0;k<PMmult[i];k++){
            if(temp[k][0] < min){
              min = temp[k][0];
              index = k;
            }
          }
          PM_res[i][j+1].Rtime=temp[index][0];
          PM_res[i][j+1].RlightCFD=temp[index][1];
          PM_res[i][j+1].RlightQDC=temp[index][2];
          temp[index][0] = 100000.;
        }
      }
      */

    // Check for leading edge
    Double_t triggerTime = 1e100;
    //  Double_t triggerTime=40.;
    for (Int_t i = 0; i < npaddles; i++)
    {
        if (fLandDigiPar->GetGeometryFileName().Contains("proto"))
        {
            if (0 == i || 9 == i || 140 == i || 149 == i)
            {
                continue;
            }
        }

        for (Int_t j = 0; j < PMmult[i]; j++)
        {
            PM_res[i][j + 1].Lenergy = lambda * PM_res[i][j + 1].LlightCFD;
            PM_res[i][j + 1].Lenergy =
                PM_res[i][j + 1].Lenergy +
                PM_res[i][j].Lenergy * exp(-lambda * (PM_res[i][j + 1].Ltime - PM_res[i][j].Ltime));

            if (PM_res[i][j + 1].Lenergy > Thresh)
            {
                if (PM_res[i][j + 1].Ltime < triggerTime)
                {
                    triggerTime = PM_res[i][j + 1].Ltime;
                } // find minimum
            }     // if above threshold

            PM_res[i][j + 1].Renergy = lambda * PM_res[i][j + 1].RlightCFD;
            PM_res[i][j + 1].Renergy =
                PM_res[i][j + 1].Renergy +
                PM_res[i][j].Renergy * exp(-lambda * (PM_res[i][j + 1].Rtime - PM_res[i][j].Rtime));

            if (PM_res[i][j + 1].Renergy > Thresh)
            {
                if (PM_res[i][j + 1].Rtime < triggerTime)
                {
                    triggerTime = PM_res[i][j + 1].Rtime;
                } // find minimum
            }     // if above threshold
        }         // digis
    }             // paddles

    //  Double_t temp[1500][3];
    for (Int_t i = 0; i < npaddles; i++)
    {
        if (fLandDigiPar->GetGeometryFileName().Contains("proto"))
        {
            if (0 == i || 9 == i || 140 == i || 149 == i)
            {
                continue;
            }
        }

        if (!PMmult[i])
        {
            continue;
        }

        // sort hits according to time; left PM
        for (Int_t j = 0; j < PMmult[i]; j++)
        {
            temp[j][0] = PM_res[i][j + 1].Ltime;
            temp[j][1] = PM_res[i][j + 1].LlightCFD;
            temp[j][2] = PM_res[i][j + 1].LlightQDC;
        }

        for (Int_t j = 0; j < PMmult[i]; j++)
        {
            Double_t min = 100000;
            Int_t index = 0;
            for (Int_t k = 0; k < PMmult[i]; k++)
            {
                if (temp[k][0] < min)
                {
                    min = temp[k][0];
                    index = k;
                }
            }
            PM_res[i][j + 1].Ltime = temp[index][0];
            PM_res[i][j + 1].LlightCFD = temp[index][1];
            PM_res[i][j + 1].LlightQDC = temp[index][2];
            temp[index][0] = 100000.;
        }

        // sort hits according to time; right PM
        for (Int_t j = 0; j < PMmult[i]; j++)
        {
            temp[j][0] = PM_res[i][j + 1].Rtime;
            temp[j][1] = PM_res[i][j + 1].RlightCFD;
            temp[j][2] = PM_res[i][j + 1].RlightQDC;
        }

        for (Int_t j = 0; j < PMmult[i]; j++)
        {
            Double_t min = 100000;
            Int_t index = 0;
            for (Int_t k = 0; k < PMmult[i]; k++)
            {
                if (temp[k][0] < min)
                {
                    min = temp[k][0];
                    index = k;
                }
            }
            PM_res[i][j + 1].Rtime = temp[index][0];
            PM_res[i][j + 1].RlightCFD = temp[index][1];
            PM_res[i][j + 1].RlightQDC = temp[index][2];
            temp[index][0] = 100000.;
        }

        Double_t tofl = 0.;
        Double_t tofr = 0.;
        Double_t lightr = 0.;
        Double_t lightl = 0.;
        Int_t multl = 0;
        Int_t multr = 0;
        PM_res[i][0].Lenergy = 0.;
        PM_res[i][0].Renergy = 0.;
        PM_res[i][0].Ltime = PM_res[i][1].Ltime;
        PM_res[i][0].Rtime = PM_res[i][1].Rtime;

        for (Int_t j = 0; j < PMmult[i]; j++)
        {
            if (PM_res[i][j + 1].Lenergy > Thresh && 0 == multl)
            {
                // This is supposed to mimic a QDC and a TDC
                // check if light is larger than threshold and register time
                // Take also time resolution of sigma=timeRes ps into account
                Double_t num = fRnd->Gaus(0., timeRes);
                tofl = PM_res[i][j + 1].Ltime + num;
                multl = multl + 1;
                lightl = 0.;
                for (Int_t k = 1; k <= PMmult[i]; k++)
                {
                    // add all times inside +- tofRange (ns) for QDC
                    hRLTimeToTrig->Fill(PM_res[i][k].Ltime - triggerTime);
                    if (TMath::Abs(PM_res[i][k].Ltime - triggerTime + fTOFRange / 2.) < fTOFRange)
                    {
                        //          if(TMath::Abs(PM_res[i][k].Ltime-triggerTime) < 4.3) {
                        lightl = lightl + PM_res[i][k].LlightQDC;
                        hPMl->Fill(PM_res[i][k].Ltime, PM_res[i][k].LlightCFD);
                    }
                }
            }

            if (PM_res[i][j + 1].Renergy > Thresh && 0 == multr)
            {
                // This is supposed to mimic a QDC and a TDC
                // check if light is larger than threshold and register time
                // Take also time resolution of sigma=timeRes ps into account
                Double_t num = fRnd->Gaus(0., timeRes);
                tofr = PM_res[i][j + 1].Rtime + num;
                multr = multr + 1;
                lightr = 0.;
                for (Int_t k = 1; k <= PMmult[i]; k++)
                {
                    // add all times inside +- tofRange (ns) for QDC
                    hRLTimeToTrig->Fill(PM_res[i][k].Rtime - triggerTime);
                    if (TMath::Abs(PM_res[i][k].Rtime - triggerTime + fTOFRange / 2.) < fTOFRange)
                    {
                        //          if(TMath::Abs(PM_res[i][k].Rtime-triggerTime) < 4.3) {
                        lightr = lightr + PM_res[i][k].RlightQDC;
                        hPMr->Fill(PM_res[i][k].Rtime, PM_res[i][k].RlightCFD);
                    }
                }
            }
        }

        // Multiplicity if only one PM has fired
        if (multl > 0 || multr > 0)
        {
            mult2 = mult2 + 1;
        }

        if (multl > 0 && multr > 0)
        {

            lightl *= TMath::Exp((2. * (plength)) * att / 2.);
            lightr *= TMath::Exp((2. * (plength)) * att / 2.);

            // PMT Saturation
            lightl = lightl / (1. + fSaturationCoefficient * lightl);
            lightr = lightr / (1. + fSaturationCoefficient * lightr);

            lightl = fRnd->Gaus(lightl, 0.05 * lightl);
            lightr = fRnd->Gaus(lightr, 0.05 * lightr);

            //      if(fBeamEnergy > 1000.) {
            //        lightl = fRnd->Gaus(lightl, 0.03*lightl);
            //        lightr = fRnd->Gaus(lightr, 0.03*lightr);
            //        Double_t k = 0.016;
            //        Double_t k1 = -5.5e-5;
            //        Double_t l = 10.;
            //        Double_t loss = fRnd->Uniform(0.96, 0.97);
            //        lightl *= loss;
            //        loss = fRnd->Uniform(0.96, 0.97);
            //        lightr *= loss;
            //        if(lightl > l) {
            //          lightl = l + (lightl-l) / (1. + k*(lightl-l) + k1*TMath::Power(lightl-l, 2));
            //        }
            //        if(lightr > l) {
            //          lightr = l + (lightr-l) / (1. + k*(lightr-l) + k1*TMath::Power(lightr-l, 2));
            //        }
            //      } else {
            //        lightl = fRnd->Gaus(lightl, 0.03*lightl);
            //        lightr = fRnd->Gaus(lightr, 0.03*lightr);
            //        Double_t k = 0.015;
            //        Double_t k1 = -5.2e-5;
            //        Double_t l = 10.;
            //        Double_t loss = fRnd->Uniform(0.85, 0.87);
            //        lightl *= loss;
            //        loss = fRnd->Uniform(0.85, 0.87);
            //        lightr *= loss;
            //        if(lightl > l) {
            //          lightl = l + (lightl-l) / (1. + k*(lightl-l) + k1*TMath::Power(lightl-l, 2));
            //        }
            //        if(lightr > l) {
            //          lightr = l + (lightr-l) / (1. + k*(lightr-l) + k1*TMath::Power(lightr-l, 2));
            //        }
            //      }

            if (lightl < threshL[i] || lightr < threshR[i])
            {
                continue;
            }

            // multiplicity if 2 PM's have fired
            mult1 = mult1 + 1;

            QDC_temp[mult1] = TMath::Sqrt(lightl * lightr);
            TDC_temp[mult1] = (tofl + tofr) / 2. - plength / cMedia;

            if (fLandDigiPar->GetGeometryFileName().Contains("proto"))
            {
                // vertical paddles
                xpos_temp[mult1] = xpaddle[i];
                ypos_temp[mult1] = (tofr - tofl) / 2. * cMedia;
                zpos_temp[mult1] = zpaddle[i];
            }
            else
            {
                if ((int)(((i - 1) / paddle_per_plane)) / 2. == (int)((int)(((i - 1) / paddle_per_plane)) / 2.))
                {
                    // horizontal paddles
                    xpos_temp[mult1] = (tofr - tofl) / 2. * cMedia;
                    ypos_temp[mult1] = ypaddle[i];
                    zpos_temp[mult1] = zpaddle[i];
                    //	    cout << "delta tof x " << (tofl - tofr) << endl;
                }
                else
                {
                    // vertical paddles
                    xpos_temp[mult1] = xpaddle[i];
                    ypos_temp[mult1] = (tofr - tofl) / 2. * cMedia;
                    zpos_temp[mult1] = zpaddle[i];
                    //	    cout << "delta tof y " << (tofl - tofr) << endl;
                }
            }

            // Here is an example how to fill the R3BLandDigi structure
            Double_t tdcL = tofl;
            Double_t tdcR = tofr;
            Double_t qdcL = lightl;
            Double_t qdcR = lightr;
            Int_t paddleNr = i + 1;
            Double_t qdc = QDC_temp[mult1];
            Double_t tdc = TDC_temp[mult1];
            Double_t xx = xpos_temp[mult1];
            Double_t yy = ypos_temp[mult1];
            Double_t zz = zpos_temp[mult1];

            AddHit(paddleNr, tdcL, tdcR, tdc, qdcL, qdcR, qdc, xx, yy, zz);
        }
    } // loop over paddles

    // control histograms
    if (mult2 > 0)
    {
        hMult2->Fill(mult2, 1.);
    }
    hMult1->Fill(mult1, 1.);

    if (fVerbose)
    {
        LOG(info) << "R3BLandDigitizer: produced " << fLandDigi->GetEntries() << " digis";
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void R3BLandDigitizer::Reset()
{
    if (fLandDigi)
    {
        fLandDigi->Clear();
    }
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void R3BLandDigitizer::Finish()
{
    hPMl->Write();
    hPMr->Write();
    hMult1->Write();
    hMult2->Write();
    hRLTimeToTrig->Write();
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
R3BLandDigi* R3BLandDigitizer::AddHit(Int_t paddleNr,
                                      Double_t tdcL,
                                      Double_t tdcR,
                                      Double_t tdc,
                                      Double_t qdcL,
                                      Double_t qdcR,
                                      Double_t qdc,
                                      Double_t xx,
                                      Double_t yy,
                                      Double_t zz)
{
    R3BLandDigi* digi = new ((*fLandDigi)[fLandDigi->GetEntriesFast()])
        R3BLandDigi(paddleNr, tdcL, tdcR, tdc, qdcL, qdcR, qdc, xx, yy, zz);
    return digi;
}
// ----------------------------------------------------------------------------

ClassImp(R3BLandDigitizer)
