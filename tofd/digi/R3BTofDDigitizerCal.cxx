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

#include "R3BTofDDigitizerCal.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "R3BMCTrack.h"
#include "R3BTofdCalData.h"
#include "R3BTofdPoint.h"
#include "TClonesArray.h"

#include "TMath.h"
#include "TRandom3.h"
#include <vector>

R3BTofDDigitizerCal::R3BTofDDigitizerCal()
    : FairTask("R3BTofDDigitizerCal")
    , fTofdPoints(NULL)
    , fTofdCals(NULL)
    , fCalTriggerItems(NULL)
    , fMCTrack(NULL)
{
    // detector parameter
    number_layers = 4;
    number_paddles = 44; // per layer
    number_channels = number_layers * number_paddles;
}

R3BTofDDigitizerCal::~R3BTofDDigitizerCal()
{
    if (fTofdPoints)
        delete fTofdPoints;
    if (fTofdCals)
        delete fTofdCals;
    if (fCalTriggerItems)
        delete fCalTriggerItems;
}

InitStatus R3BTofDDigitizerCal::Init()
{
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    LOG_IF(fatal, !ioman) << "Init: No FairRootManager";

    fTofdPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("TOFdPoint"));

    // fMCTrack = (TClonesArray*)ioman->GetObject("MCTrack");
    // ioman->Register("MCTrack", "Monte Carlo data", fMCTrack, kTRUE);

    // Register output array fTofdCals
    fTofdCals = new TClonesArray("R3BTofdCalData");
    ioman->Register("TofdCal", "Digital response in TofD", fTofdCals, kTRUE);
    maxevent = ioman->CheckMaxEventNo();

    fCalTriggerItems = new TClonesArray("R3BTofdCalData");
    ioman->Register("TofdTriggerCal", "TofdTriggerCal data", fCalTriggerItems, kTRUE);

    // Get random number for smearing in y, t, ELoss
    prnd = new TRandom3();

    return kSUCCESS;
}

void R3BTofDDigitizerCal::Exec(Option_t* opt)
{
    Reset();

    auto Digitize = [this](TClonesArray* Points, TClonesArray* Hits, TClonesArray* Trigger, Int_t NumOfChannels)
    {
        Int_t entryNum = Points->GetEntries();

        if (!entryNum)
            return;

        // creating the storage for energy and time for each paddle
        std::vector<Double_t>* energy = new std::vector<Double_t>[NumOfChannels];
        std::vector<Double_t>* time = new std::vector<Double_t>[NumOfChannels];
        std::vector<Double_t>* y = new std::vector<Double_t>[NumOfChannels];

        for (Int_t i = 0; i < NumOfChannels; ++i)
        {
            time[i].push_back(-99999);
            energy[i].push_back(0);
            y[i].push_back(0);
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
            R3BTofdPoint* data_element = dynamic_cast<R3BTofdPoint*>(Points->At(i));
            // Plane and paddle data members in 1-base
            auto channel = (data_element->GetPlane() - 1) * number_paddles + data_element->GetPaddle() - 1;
            TempHits.push_back(TempHit(channel, // channel nummer
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
            // discard all hits with an energy loss < cut -> 1 keV
            if (Hit.Energy < 0.000001)
                continue;

            Int_t ChannelID = Hit.ChannelID;
            // pile-up simulation:
            if (Hit.Time - time[ChannelID].back() > 200.)
            {
                y[ChannelID].push_back(Hit.Y); // get y-position //local
                time[ChannelID].push_back(Hit.Time);
                energy[ChannelID].push_back(Hit.Energy);
            }

            else
            {
                energy[ChannelID].back() += Hit.Energy;
                y[ChannelID].back() = (time[ChannelID].back() > Hit.Time) ? Hit.Y : y[ChannelID].back();
                time[ChannelID].back() = (time[ChannelID].back() > Hit.Time) ? Hit.Time : time[ChannelID].back();
            }
        }

        // creating the final hits

        int layer_label = 0;
        int paddle_number = 0;
        Double_t yrnd, yns, ToT_up, ToT_down, ernd, ens, timernd, timeL_up = -1., timeL_down = -1., timeT_up = -1.,
                                                                  timeT_down = -1.;

        for (Int_t i = 0; i < NumOfChannels; ++i)
        {
            for (Double_t& energyl : energy[i])
            {
                if (energyl > 0.0001)
                {
                    if (i < number_paddles) // plane 1, paddle_number: 0-43
                    {
                        layer_label = 0;
                        paddle_number = i;
                    }
                    else if (i < number_paddles * 2) // plane 2, paddle_number: 0-43
                    {
                        layer_label = 1;
                        paddle_number = i - number_paddles;
                    }
                    else if (i < number_paddles * 3) // plane 3, paddle_number: 0-43
                    {
                        layer_label = 2;
                        paddle_number = i - number_paddles * 2;
                    }
                    else if (i < number_paddles * 4) // plane 4, paddle_number: 0-43
                    {
                        layer_label = 3;
                        paddle_number = i - number_paddles * 3;
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

                    new ((*Hits)[Hits->GetEntries()])
                        R3BTofdCalData(layer_label + 1, paddle_number + 1, 1, timeL_up, timeT_up);

                    new ((*Hits)[Hits->GetEntries()])
                        R3BTofdCalData(layer_label + 1, paddle_number + 1, 2, timeL_down, timeT_down);

                    for (Int_t j = 0; j < 12; j++)
                    {
                        new ((*Trigger)[Trigger->GetEntries()]) R3BTofdCalData(5, j + 1, 1, 0., 0.);
                    }
                }
            }
        }

        delete[] energy;
        delete[] time;
        delete[] y;
    };

    // running the digitizer for the tofd
    if (fTofdPoints)
    {
        Digitize(fTofdPoints, fTofdCals, fCalTriggerItems, number_channels);
    }
}

// -------------------------------------------------------------------------

void R3BTofDDigitizerCal::Reset()
{
    // Clear the structure
    if (fTofdCals)
        fTofdCals->Clear();
    if (fCalTriggerItems)
        fCalTriggerItems->Clear();
}

ClassImp(R3BTofDDigitizerCal);
