/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

#include "R3BFi23aDigitizerHit.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"

// includes for modeling
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoShapeAssembly.h"
#include "TParticle.h"
#include "TVirtualMC.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "R3BFibPoint.h"
#include "R3BMCTrack.h"

using std::cout;
using std::endl;

R3BFi23aDigitizerHit::R3BFi23aDigitizerHit()
    : FairTask("R3B Fi23a Digitization scheme ")
    , fFi23aPoints(NULL)
{

    esigma = 0.001;
    tsigma = 0.01;
    ysigma = 1;
}

R3BFi23aDigitizerHit::R3BFi23aDigitizerHit(Double_t e, Double_t t, Double_t x, Double_t y)
    : FairTask("R3B Fi23a Digitization scheme ")
    , fFi23aPoints(NULL)
{

    esigma = e;
    tsigma = t;
    xsigma = x;
    ysigma = y;
}

R3BFi23aDigitizerHit::~R3BFi23aDigitizerHit() {}

void R3BFi23aDigitizerHit::SetEnergyResolution(Double_t e) { esigma = e; }
void R3BFi23aDigitizerHit::SetTimeResolution(Double_t t) { tsigma = t; }
void R3BFi23aDigitizerHit::SetXPositionResolution(Double_t x) { xsigma = x; }
void R3BFi23aDigitizerHit::SetYPositionResolution(Double_t y) { ysigma = y; }

InitStatus R3BFi23aDigitizerHit::Init()
{

    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init:No FairRootManager";
    fFi23aPoints = (TClonesArray*)ioman->GetObject("Fi23aPoint");

    fMCTrack = (TClonesArray*)ioman->GetObject("MCTrack");

    // Register output array DchDigi
    fFi23aHits = new TClonesArray("R3BFiberMAPMTHitData", 1000);
    ioman->Register("Fi23aHit", "Digital response in Fi23a", fFi23aHits, kTRUE);

    // for sigmas
    prnd = new TRandom3();

    return kSUCCESS;
}

void R3BFi23aDigitizerHit::Exec(Option_t* opt)
{
    Reset();

    auto Digitize = [this](TClonesArray* Points, TClonesArray* Hits, Int_t NumOfFibers) {
        Int_t entryNum = Points->GetEntries();

        if (!entryNum)
            return;

        // creating the storage for energy and time for each PMT
        std::vector<Double_t>* energy = new std::vector<Double_t>[NumOfFibers];
        std::vector<Double_t>* time = new std::vector<Double_t>[NumOfFibers];
        std::vector<Double_t>* x = new std::vector<Double_t>[NumOfFibers];
        std::vector<Double_t>* y = new std::vector<Double_t>[NumOfFibers];

        for (Int_t i = 0; i < NumOfFibers; ++i)
        {
            time[i].push_back(99999);
            energy[i].push_back(0);
            x[i].push_back(0);
            y[i].push_back(0);
        }

        struct TempHit
        {
            Int_t fiberID;
            Double_t Energy;
            Double_t Time;
            Double_t X;
            Double_t Y;

            TempHit(Int_t id, Double_t e, Double_t t, Double_t x, Double_t y)
                : fiberID(id)
                , Energy(e)
                , Time(t)
                , X(x)
                , Y(y)
            {
            }
        };

        // ordering the hits in time

        std::vector<TempHit> TempHits;

        for (Int_t i = 0; i < entryNum; ++i)
        {
            R3BFibPoint* data_element = (R3BFibPoint*)Points->At(i);

            TempHits.push_back(TempHit(data_element->GetDetectorID(),
                                       data_element->GetEnergyLoss(),
                                       data_element->GetTime(),
                                       data_element->GetXIn(),
                                       data_element->GetYIn()));
        }

        std::sort(TempHits.begin(), TempHits.end(), [](const TempHit& lhs, const TempHit& rhs) {
            return lhs.Time < rhs.Time;
        });

        // registering the hits in the PMTs

        for (TempHit& Hit : TempHits)
        {
            if (Hit.Energy < 1.e-12)
            {
                continue;
            }

            Int_t fiberID = Hit.fiberID;

            LOG(debug) << "Hit Fi23a in: fiber: " << Hit.fiberID << " x: " << Hit.X << " y: " << Hit.Y
                       << " Eloss: " << Hit.Energy << " t: " << Hit.Time;

            if (Hit.Time - time[fiberID].back() < 30)
            {
                energy[fiberID].back() += Hit.Energy;
                x[fiberID].back() = (time[fiberID].back() > Hit.Time) ? Hit.X : x[fiberID].back();
                y[fiberID].back() = (time[fiberID].back() > Hit.Time) ? Hit.Y : y[fiberID].back();
                time[fiberID].back() = (time[fiberID].back() > Hit.Time) ? Hit.Time : time[fiberID].back();
            }
            else
            {
                energy[fiberID].push_back(Hit.Energy);
                time[fiberID].push_back(Hit.Time);
                x[fiberID].push_back(Hit.X);
                y[fiberID].push_back(Hit.Y);
            }
        }

        // creating the final hits

        for (Int_t i = 0; i < NumOfFibers; ++i)
        {
            for (Double_t& energyl : energy[i])
            {
                if (energyl > 1.e-12)
                {
                    Double_t fiber_id = i;
                    LOG(debug) << "Hit Fi23a out: fiber: " << i << " x: " << (x[i].at(&energyl - energy[i].data()))
                               << " y: " << (y[i].at(&energyl - energy[i].data())) << " Eloss: " << energyl
                               << " t: " << time[i].at(&energyl - energy[i].data());

                    Double_t xx = x[i].at(&energyl - energy[i].data());

                    Bool_t granularity = true;

                    if (granularity)
                    {
                        LOG(debug) << "x before granularity: " << xx;
                        xx = -detector_width / 2. + fiber_thickness * (1 + air_layer) / 2. +
                             double(fiber_id) * (1 + air_layer) * fiber_thickness;

                        // Double_t fiber_width = 0.1; // cm
                        // xx = (int)((xx + fiber_width / 2.) / fiber_width) * fiber_width;
                        LOG(debug) << "x after granularity: " << xx;
                    }

                    new ((*Hits)[Hits->GetEntries()])
                        R3BFiberMAPMTHitData(1,
                                             prnd->Gaus(xx, xsigma),
                                             prnd->Gaus((y[i].at(&energyl - energy[i].data())), ysigma),
                                             prnd->Gaus(energyl, esigma),
                                             prnd->Gaus(time[i].at(&energyl - energy[i].data()), tsigma),
                                             i + 1,
                                             0.,
                                             0.,
                                             0.,
                                             0.);
                }
            }
        }

        delete[] energy;
        delete[] time;
        delete[] x;
        delete[] y;
    };

    // running the digitizer for the Fi detectors

    if (fFi23aPoints)
    {
        Digitize(fFi23aPoints, fFi23aHits, fiber_nbr);
    }
}
// -------------------------------------------------------------------------

void R3BFi23aDigitizerHit::Reset()
{
    if (fFi23aHits)
        fFi23aHits->Clear();
}

void R3BFi23aDigitizerHit::Finish() {}

ClassImp(R3BFi23aDigitizerHit)
