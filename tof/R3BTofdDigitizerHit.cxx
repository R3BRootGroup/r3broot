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

// TofdDigitzer 07 aug 2017 created by Alicia Wongel

#include "R3BTofdDigitizerHit.h"
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
#include "TRandom3.h"
#include "TVector3.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "R3BMCTrack.h"
#include "R3BTofdPoint.h"

using namespace std;

R3BTofdDigitizerHit::R3BTofdDigitizerHit()
    : FairTask("R3B Tofd Digitization scheme ")
    , fTofdPoints(NULL)
    , fTofdHits(NULL)
{

    // set default values for smearing
    fsigma_y = 0.001;
    fsigma_t = 0.001;
    fsigma_ELoss = 0.001;
}

R3BTofdDigitizerHit::~R3BTofdDigitizerHit()
{
    if (fTofdPoints)
        delete fTofdPoints;
    if (fTofdHits)
        delete fTofdHits;
}

InitStatus R3BTofdDigitizerHit::Init()
{

    //  cout<<"Init "<<endl;
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fTofdPoints = (TClonesArray*)ioman->GetObject("TOFdPoint");
    fMCTrack = (TClonesArray*)ioman->GetObject("MCTrack");

    // Register output array fTofdHits
    fTofdHits = new TClonesArray("R3BTofdHitData", 1000);
    ioman->Register("TofdHit", "Digital response in Tofd", fTofdHits, kTRUE);

    // Get random number for smearing in y, t, ELoss

    fRnd = new TRandom3();

    fHist1 = new TH1F("fHist1", "Energy loss histogram of Monte Carlo Points", 5000, 0., 5.);
    fHist2 = new TH1F("fHist2", "Energy loss histogram without merging", 5000, 0., 5.);
    fHist3 = new TH1F("fHist3", "Energy loss histogram with merging", 5000, 0., 5.);

    return kSUCCESS;
}

void R3BTofdDigitizerHit::Exec(Option_t* opt)
{

    Reset();

    //
    vector<R3BTofdPoint*> vPoints[1000];
    Int_t n_entries = fTofdPoints->GetEntries();
    Bool_t ChannelFired[1000] = { kFALSE };

    Double_t X_Pos[10000] = { 0 };
    Double_t Y_Pos[10000] = { 0 };
    Double_t Z_Pos[10000] = { 0 };
    Double_t Time[10000] = { 0 };
    Double_t Energy_Loss[10000] = { 0 };

    Int_t number_paddles_hit = 0;

    for (Int_t entry = 0; entry < n_entries; entry++)
    {

        R3BTofdPoint* data_element = (R3BTofdPoint*)fTofdPoints->At(entry);

        Int_t DetectorID = data_element->GetDetectorID();
        Double_t energy_loss = data_element->GetEnergyLoss();

        // discard all hits with an energy loss < cut
        if (energy_loss < 0.0000001)
            continue;

        // add hits to vector
        vPoints[DetectorID].push_back(data_element);

        ChannelFired[DetectorID] = kTRUE;
    }

    R3BTofdHitData* hit = nullptr;
    map<Int_t, R3BTofdHitData*> MapOfHits;

    for (Int_t channel = 0; channel < 1000; channel++)
    {
        if (ChannelFired[channel] == kTRUE)
        {

            sort(vPoints[channel].begin(),
                 vPoints[channel].end(),
                 [](const R3BTofdPoint* lhs, const R3BTofdPoint* rhs) { return lhs->GetTime() < rhs->GetTime(); });

            for (Int_t point = 0; point < vPoints[channel].size(); point++)
            {

                // energy threshold
                if (vPoints[channel].at(point)->GetEnergyLoss() < 0.0000001)
                    continue;
                fHist1->Fill(vPoints[channel].at(point)->GetEnergyLoss());

                if (0 == point || (vPoints[channel].at(point)->GetTime() - MapOfHits[channel]->GetTime()) > 30)
                { // add new hits
                    int layer_label = 0;
                    int paddle_number = 0;
                    if (channel < 200)
                    {
                        layer_label = 0;
                        paddle_number = channel - 100 - 1;
                        // paddle_number = channel;
                    }
                    else if (channel < 300)
                    {
                        layer_label = 1;
                        paddle_number = channel - 200 - 1;
                        // paddle_number = channel - 100;
                    }

                    Double_t xpaddle = -detector_width / 2. + paddle_width / 2. * (1 + layer_label) +
                                       paddle_number * paddle_width + paddle_number * air_gap_paddles;

                    X_Pos[channel] =
                        (vPoints[channel].at(point)->GetXIn() + vPoints[channel].at(point)->GetXOut()) / 2.;
                    Y_Pos[channel] = (vPoints[channel].at(point)->GetYIn() + vPoints[channel].at(point)->GetYOut()) /
                                     2.; // get y-position //local
                    Z_Pos[channel] = (vPoints[channel].at(point)->GetZIn() + vPoints[channel].at(point)->GetZOut()) /
                                     2.; // get y-position //local
                                         // position for PAS simulations
                                         /*
                                                             Double_t PositionX = -160.;
                                                             Double_t PositionY = 0.;
                                                             Double_t PositionZ = 700.;
                                                             Double_t RotationY = -14.;
                                         */
                                         // position for s494 simulations
                    Double_t PositionX = -129.6732;
                    Double_t PositionY = 1.492008;
                    Double_t PositionZ = 686.1278;
                    Double_t RotationY = -18.;

                    TVector3 posGlobal;
                    posGlobal.SetX(X_Pos[channel]);
                    posGlobal.SetY(Y_Pos[channel]);
                    posGlobal.SetZ(Z_Pos[channel]);

                    TVector3 pos0;
                    pos0.SetX(0);
                    pos0.SetY(0);
                    pos0.SetZ(0);
                    pos0.RotateY(RotationY * TMath::DegToRad());
                    TVector3 trans(PositionX, PositionY, PositionZ);
                    pos0 += trans;

                    TVector3 local = posGlobal - pos0;
                    local.RotateY(-RotationY * TMath::DegToRad());
                    Double_t x_local = local.X();
                    Double_t y_local = local.Y();

                    Bool_t granularity = true;
                    if (granularity)
                    {
                        LOG(debug) << "x before granularity: " << x_local << ", layer = " << layer_label + 1 << endl;

                        Double_t xp = -1000.;
                        // calculate x-position
                        if (layer_label + 1 == 1 || layer_label + 1 == 3)
                        {
                            xp = -detector_width / 2 + (paddle_width + air_gap_paddles) / 2 +
                                 paddle_number * (paddle_width + air_gap_paddles) - 0.04;
                        }
                        if (layer_label + 1 == 2 || layer_label + 1 == 4)
                        {
                            xp = -detector_width / 2 + (paddle_width + air_gap_paddles) +
                                 paddle_number * (paddle_width + air_gap_paddles) - 0.04;
                        }

                        x_local = xp;

                        LOG(debug) << "x after granularity: " << x_local << ", layer = " << layer_label + 1 << endl;
                    }

                    Time[channel] = vPoints[channel].at(point)->GetTime();
                    Energy_Loss[channel] = vPoints[channel].at(point)->GetEnergyLoss();
                    /*
                                        if (Energy_Loss[channel] < 0.03)
                                        {
                                            Energy_Loss[channel] = 2.;
                                        }
                                        else if (Energy_Loss[channel] < 0.08)
                                        {
                                            Energy_Loss[channel] = 6.;
                                        }
                                        else
                                        {
                                            Energy_Loss[channel] = 8.;
                                        }
                    */
                    // add to HitData and introduce smearing of y-position, time and energy loss

                    LOG(debug) << "Hit Tofd: ch = " << channel << " paddle = " << paddle_number << " x = " << x_local
                               << " y = " << y_local << " t = " << Time[channel] << " eloss = " << Energy_Loss[channel];

                    //                    if(channel < 200)

                    Int_t qcharge = (int)(88.339847 * Energy_Loss[channel] + 1.5039877 + 0.5);

                    gRandom->SetSeed(0);

                    Double_t ysigma = 2.;
                    if (qcharge < 3)
                        ysigma = 4.;
                    y_local = gRandom->Gaus(y_local, ysigma);

                    //  if(layer_label == 0)
                    {
                        MapOfHits.insert(pair<Int_t, R3BTofdHitData*>(
                            channel,
                            new R3BTofdHitData(0.,
                                               // fRnd->Uniform((X_Pos[channel] - 1.35), (X_Pos[channel] + 1.35)),
                                               x_local,
                                               y_local,
                                               0.,
                                               Time[channel],
                                               qcharge, // Energy_Loss[channel],
                                               layer_label + 1)));
                    }
                }

                else
                { // update hit
                    MapOfHits[channel]->SetTime(vPoints[channel].at(point)->GetTime());
                    MapOfHits[channel]->SetEloss(vPoints[channel].at(point)->GetEnergyLoss() +
                                                 MapOfHits[channel]->GetEloss());
                }
            }
        }
    }

    // merge energy loss of adjacent paddles in different layers

    Int_t counter = 0;
    vector<R3BTofdHitData*> vHits;

    for (auto& i : MapOfHits)
    {

        fHist2->Fill(i.second->GetEloss());

        if (i.first >= number_paddles)
        {
            continue;
        }

        vHits.push_back(i.second);

        Int_t k1 = i.first + number_paddles;
        Int_t k2 = i.first + number_paddles - 1;

        // check which paddle(s) in 2nd layer were hit
        if (MapOfHits.count(k1))
        {
            vHits.push_back(MapOfHits[k1]);
        }

        if (i.first != 0)
        {
            if (MapOfHits.count(k2))
            {
                vHits.push_back(MapOfHits[k2]);
            }
        }

        // sum over energy loss
        Double_t sum_eloss = 0;
        for (Int_t j = 0; j < vHits.size(); j++)
        {
            sum_eloss += vHits[j]->GetEloss();
        }

        // average over energy loss
        Double_t merged_eloss = sum_eloss / (vHits.size());

        // set merged energy loss as new energy loss for hits
        for (Int_t j = 0; j < vHits.size(); j++)
        {
            vHits[j]->SetEloss(merged_eloss);
        }

        vHits.clear();
    }

    for (auto& i : MapOfHits)
    {
        new ((*fTofdHits)[counter]) R3BTofdHitData(*(i.second));
        fHist3->Fill(i.second->GetEloss());
        delete i.second;
        counter++;
    }
}
// -------------------------------------------------------------------------

void R3BTofdDigitizerHit::Reset()
{
    // Clear the structure
    //   cout << " -I- Digit Reset() called " << endl;

    if (fTofdHits)
        fTofdHits->Clear();
}

void R3BTofdDigitizerHit::Finish()
{
    fHist1->Write();
    fHist2->Write();
    fHist3->Write();
}

ClassImp(R3BTofdDigitizerHit)
