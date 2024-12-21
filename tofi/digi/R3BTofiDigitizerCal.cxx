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

#include "R3BTofiDigitizerCal.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "R3BMCTrack.h"
#include "R3BTofi.h"
#include "R3BTofiPoint.h"
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
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using std::cout;
using std::endl;

R3BTofiDigitizerCal::R3BTofiDigitizerCal()
    : FairTask("R3B Tofi Digitization scheme ")
    , fTofiPoints(NULL)
{
}

R3BTofiDigitizerCal::~R3BTofiDigitizerCal() {}

InitStatus R3BTofiDigitizerCal::Init()
{
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";

    fTofiPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("TOFiPoint"));

    // fMCTrack = (TClonesArray*)ioman->GetObject("MCTrack");
    // ioman->Register("MCTrack", "Monte Carlo data", fMCTrack, kTRUE);

    // Register output array fTofiCals
    fTofiCals = new TClonesArray("R3BTofiCalData", 1000);
    ioman->Register("TofiCal", "Digital response in Tofi", fTofiCals, kTRUE);
    maxevent = ioman->CheckMaxEventNo();

    fCalTriggerItems = new TClonesArray("R3BTofiCalData", 1000);
    ioman->Register("TofiTriggerCal", "Land", fCalTriggerItems, kTRUE);

    // Get random number for smearing in y, t, ELoss
    prnd = new TRandom3();

    fhMultTofi = new TH1F("multTofi", "Multiplicity ToFI", 20, 0, 20);
    fhMultTofi->GetXaxis()->SetTitle("Multiplicity");

    return kSUCCESS;
}

void R3BTofiDigitizerCal::Exec(Option_t* opt)
{
    //	cout<<"R3BTofiDigitizerCal Exec Entry"<<endl;
    if (counter / 10000. == (int)counter / 10000)
        std::cout << "\rEvents: " << counter << " / " << maxevent << " (" << (int)(counter * 100. / maxevent) << " %) "
                  << std::flush;
    counter += 1;

    Reset();

    //    cout<<"R3BTofiDigitizerCal Exec Before Digitize"<<endl;

    auto Digitize = [this](TClonesArray* Points, TClonesArray* Hits, TClonesArray* Trigger, Int_t NumOfChannels)
    {
        Int_t entryNum = Points->GetEntries();

        if (!entryNum)
            return;

        // creating the storage for energy and time for each paddle
        std::vector<Double_t>* energy = new std::vector<Double_t>[NumOfChannels];
        std::vector<Double_t>* time = new std::vector<Double_t>[NumOfChannels];
        std::vector<Double_t>* y = new std::vector<Double_t>[NumOfChannels];
        // std::vector<Double_t>* energy_temp = new std::vector<Double_t>[NumOfChannels];
        // std::vector<Double_t>* time_temp = new std::vector<Double_t>[NumOfChannels];
        // std::vector<Double_t>* y_temp = new std::vector<Double_t>[NumOfChannels];

        for (Int_t i = 0; i < NumOfChannels; ++i)
        {
            time[i].push_back(-99999);
            energy[i].push_back(0);
            y[i].push_back(0);
            // time_temp[i].push_back(-99999);
            // energy_temp[i].push_back(0);
            // y_temp[i].push_back(0);
        }
        // Sim data:
        struct TempHit
        {
            Int_t ChannelID;
            Double_t Energy;
            Double_t Time;
            Double_t Y;

            TempHit(Int_t id, Double_t e, Double_t t, Double_t y)
                : ChannelID(id)
                , Energy(e)
                , Time(t)
                , Y(y)
            {
            }
        };

        // ordering the hits in time

        std::vector<TempHit> TempHits;

        for (Int_t i = 0; i < entryNum; ++i)
        {
            R3BTofiPoint* data_element = dynamic_cast<R3BTofiPoint*>(Points->At(i));

            TempHits.push_back(TempHit(data_element->GetDetectorID(), // channel nummer
                                       data_element->GetEnergyLoss(),
                                       data_element->GetTime(),
                                       data_element->GetYIn()));
        }

        std::sort(TempHits.begin(),
                  TempHits.end(),
                  [](const TempHit& lhs, const TempHit& rhs) { return lhs.Time < rhs.Time; });

        Int_t number_paddles_hit = 0;
        for (TempHit& Hit : TempHits) // loop over full range of TempHits; Hit is reference to each TempHit
        {

            // discard all hits with an energy loss < cut
            if (Hit.Energy < 0.000001 || Hit.Time > 1.e9)
                continue;

            Int_t ChannelID = Hit.ChannelID;
            // pile-up simulation:
            if (Hit.Time - time[ChannelID].back() > 200.)
            {
                y[ChannelID].push_back(Hit.Y); // get y-position //local
                time[ChannelID].push_back(Hit.Time);
                energy[ChannelID].push_back(Hit.Energy);
                //     cout<<"t>200: "<<Hit.Time <<"; "<<Hit.Energy*1000<<endl;
            }

            else
            {
                energy[ChannelID].back() += Hit.Energy;
                y[ChannelID].back() = (time[ChannelID].back() > Hit.Time) ? Hit.Y : y[ChannelID].back();
                time[ChannelID].back() = (time[ChannelID].back() > Hit.Time) ? Hit.Time : time[ChannelID].back();
                //       cout<<"t<200: "<<Hit.Time <<"; "<<time[ChannelID].back()<<",
                //       "<<energy[ChannelID].back()*1000.<<
                //       "; "<<Hit.Energy*1000.<<endl;
            }
        }

        // creating the final hits

        int layer_label = -100;
        int paddle_number = -100;
        Double_t yrnd, yns, ToT_up, ToT_down, ernd, ens, timernd, timeL_up = -1., timeL_down = -1., timeT_up = -1.,
                                                                  timeT_down = -1.;
        Int_t mult = 0;

        for (Int_t i = 0; i < NumOfChannels; ++i)
        {
            for (Double_t& energyl : energy[i])
            {
                if (energyl > 0.0001)
                {
                    mult++;
                    if (i < 200) // ch: 101-144, layer: 0, paddle_number: 0-43
                    {
                        layer_label = 0;
                        paddle_number = i - 101;
                    }
                    else if (i < 300 && i > 200) // ch: 201-244, layer: 1, paddle_number: 0-43
                    {
                        layer_label = 1;
                        paddle_number = i - 201;
                    }
                    else if (i < 400 && i > 300) // ch: 301-344, layer: 2, paddle_number: 0-43
                    {
                        layer_label = 2;
                        paddle_number = i - 301;
                    }
                    else if (i < 500 && i > 400) // ch: 401-444, layer: 3, paddle_number: 0-43
                    {
                        layer_label = 3;
                        paddle_number = i - 401;
                    }

                    /* From y-position get ToT_up and ToT_down only considering absorption;
                     * ToT_up = energyl*exp(-l/labs), labs assume 100cm
                     */
                    yrnd = prnd->Gaus((y[i].at(&energyl - energy[i].data())), ysigma);
                    ernd = prnd->Gaus(energyl, esigma) * 1000.; // GeV->MeV
                    ToT_up =
                        ernd * exp(-(50. - yrnd) / 100.); // *10 just to come to the tot region we measure in experiment
                    ToT_down = ernd * exp(-(50 + yrnd) / 100.);

                    /* Now, make from ToT and time leading and trailing times:
                     * time=(timeL_up+timeL_down)/2, and timeL_up-timeL_down=y => timeL_up, timeL_down
                     * time_trail = time_lead+ToT
                     */
                    timernd = prnd->Gaus(time[i].at(&energyl - energy[i].data()),
                                         tsigma); //*3. just to get some larger values
                    Double_t veff = 12.;
                    timeL_up = timernd - 50. / veff + (50 - yrnd) / veff;
                    timeL_down = timernd - 50. / veff + (50 + yrnd) / veff;

                    timeT_up = timeL_up + ToT_up;
                    timeT_down = timeL_down + ToT_down;

                    /*cout<<"layer= "<<layer_label<<" paddle= "<<paddle_number<<" elos= "<<ernd<<"  time= "
                        <<timernd<< " tL_up= "<<timeL_up<<" tT_up= "<<timeT_up<<" tL_down= "<<timeL_down
                        <<" tT_down= "<<timeT_down<<endl;
                    */
                    /*   CalData format:      --> in Tree:
                        fPlane(detector) 1-4  --> TofiCal.fDetector
                      , fBar(channel)    1-44 --> TofiCal.fBar
                      , fSide(side)      1-2  --> TofiCal.fSide
                      , fLeading_ns(lead)     --> TofiCal.fLeading_ns
                      , fTrailing_ns(trail)   --> TofiCal.fTrailing_ns
                    */
                    new ((*Hits)[Hits->GetEntries()])
                        R3BTofiCalData(layer_label + 1, paddle_number + 1, 1, timeL_up, timeT_up);

                    new ((*Hits)[Hits->GetEntries()])
                        R3BTofiCalData(layer_label + 1, paddle_number + 1, 2, timeL_down, timeT_down);

                    // Int_t card = (int)paddle_number/8.+layer_label*6;
                    for (Int_t j = 0; j < 12; j++)
                    {
                        new ((*Trigger)[Trigger->GetEntries()]) R3BTofiCalData(5, j + 1, 1, 0., 0.);
                    }
                }
            }
        }
        fhMultTofi->Fill(mult);

        delete[] energy;
        delete[] time;
        delete[] y;
    };

    // running the digitizer for the Tofi

    if (fTofiPoints)
    {
        Digitize(fTofiPoints, fTofiCals, fCalTriggerItems, number_channels);
    }
}

// -------------------------------------------------------------------------

void R3BTofiDigitizerCal::Reset()
{
    // Clear the structure
    //   cout << " -I- Digit Reset() called " << endl;

    if (fTofiCals)
        fTofiCals->Clear();
    if (fCalTriggerItems)
        fCalTriggerItems->Clear();
}

void R3BTofiDigitizerCal::Finish() { fhMultTofi->Write(); }

ClassImp(R3BTofiDigitizerCal)
