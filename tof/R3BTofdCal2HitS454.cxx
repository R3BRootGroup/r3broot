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
// -----                 R3BTofdCal2HitS454                   -----
// -----            Created May 2016 by M.Heil            -----
// -----           Modified Mar 2020 by L.Bott            -----
// ------------------------------------------------------------

#include "R3BTofdCal2HitS454.h"
#include "R3BEventHeader.h"
#include "R3BTCalEngine.h"
#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include "R3BTofdHitModulePar.h"
#include "R3BTofdHitPar.h"
#include <FairRootManager.h>

#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"

#include "TClonesArray.h"
#include "TMath.h"

#include <iostream>
#include <sstream>
#include <string>

// This header is included in several places, and would re-define the mapping
// table every time -> linker nor happy.
// So do some forward declaration for now, and keep the include in just one
// place, for now R3BTofdCal2Histo.cxx.
// #include "mapping_tofd_trig.hh"
extern unsigned g_tofd_trig_map[4][2][48];
void tofd_trig_map_setup();

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

#define N_TOFD_HIT_PLANE_MAX 4
#define N_TOFD_HIT_PADDLE_MAX 44

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
} // namespace

R3BTofdCal2HitS454::R3BTofdCal2HitS454()
    : R3BTofdCal2HitS454("TofdCal2Hit", 1)
{
}

R3BTofdCal2HitS454::R3BTofdCal2HitS454(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems(NULL)
    , fCalTriggerItems(NULL)
    , fHitItems(new TClonesArray("R3BTofdHitData"))
    , fNofHitItems(0)
    , fNofHitPars(0)
    , fHitPar(NULL)
    , fTrigger(-1)
    , fTpat(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fTofdQ(1)
    , fTofdHisto(true)
    , fTofdTotPos(true)
    , fnEvents(0)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , maxevent(0)
    , wrongtrigger(0)
    , wrongtpat(0)
    , headertpat(0)
    , events_in_cal_level(0)
    , inbarcoincidence(0)
    , eventstore(0)
    , singlehit(0)
    , multihit(0)
    , bars_with_multihit(0)
    , events_wo_tofd_hits(0)
{
    if (fTofdHisto)
    {
        for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
        {
            fhQ[i] = NULL;
            fhxy[i] = NULL;
            fhQvsEvent[i] = NULL;
            fhTdiff[i] = NULL;
            fhTsync[i] = NULL;
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {
                fhQvsPos[i][j] = NULL;
                // fhQvsTHit[i][j] = NULL;
                // fhTvsTHit[i][j] = NULL;
            }
        }
    }
}

R3BTofdCal2HitS454::~R3BTofdCal2HitS454()
{
    if (fTofdHisto)
    {
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            if (fhQ[i])
                delete fhQ[i];
            if (fhxy[i])
                delete fhxy[i];
            if (fhQvsEvent[i])
                delete fhQvsEvent[i];
            if (fhTdiff[i])
                delete fhTdiff[i];
            if (fhTsync[i])
                delete fhTsync[i];
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {
                if (fhQvsPos[i][j])
                    delete fhQvsPos[i][j];
            }
        }
    }
    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

InitStatus R3BTofdCal2HitS454::Init()
{
    fHitPar = dynamic_cast<R3BTofdHitPar*>(FairRuntimeDb::instance()->getContainer("TofdHitPar"));
    if (!fHitPar)
    {
        LOG(error) << "Could not get access to TofdHitPar-Container.";
        fNofHitPars = 0;
        return kFATAL;
    }
    fNofHitPars = fHitPar->GetNumModulePar();
    if (fNofHitPars == 0)
    {
        LOG(error) << "There are no Hit parameters in container TofdHitPar";
        return kFATAL;
    }

    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdCal"));
    if (NULL == fCalItems)
        LOG(fatal) << "Branch TofdCal not found";
    fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdTriggerCal"));
    if (NULL == fCalTriggerItems)
        LOG(fatal) << "Branch TofdTriggerCal not found";
    maxevent = mgr->CheckMaxEventNo();

    // request storage of Hit data in output tree
    mgr->Register("TofdHit", "Land", fHitItems, kTRUE);

    tofd_trig_map_setup();

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BTofdCal2HitS454::SetParContainers()
{
    fHitPar = dynamic_cast<R3BTofdHitPar*>(FairRuntimeDb::instance()->getContainer("TofdHitPar"));
    if (!fHitPar)
    {
        LOG(error) << "Could not get access to TofdHitPar-Container.";
        fNofHitPars = 0;
        return;
    }
}

InitStatus R3BTofdCal2HitS454::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

namespace
{
    uint64_t n1, n2;
};

void R3BTofdCal2HitS454::Exec(Option_t* option)
{
    static uint32_t counter = 0;
    if (0 == counter % 10000)
        std::cout << "\rEvents: " << counter << " / " << maxevent << " (" << (int)(counter * 100. / maxevent) << " %) "
                  << std::flush;
    ++counter;

    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
    {
        wrongtrigger++;
        return;
    }
    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    if (fTpat_bit >= 0)
    {
        Int_t itpat = header->GetTpat();
        Int_t tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if ((header) && (tpatvalue == 0))
        {
            wrongtpat++;
            return;
        }
    }
    headertpat++;
    Double_t timeP0 = 0.;
    Double_t randx;
    std::vector<std::vector<std::vector<Double_t>>> q;
    std::vector<std::vector<std::vector<Double_t>>> thit;
    std::vector<std::vector<std::vector<Double_t>>> x;
    std::vector<std::vector<std::vector<Double_t>>> y;
    std::vector<std::vector<std::vector<Double_t>>> yToT;
    UInt_t vmultihits[N_PLANE_MAX + 1][N_TOFD_HIT_PADDLE_MAX + 1];
    for (Int_t i = 0; i <= fNofPlanes; i++)
    {
        q.push_back(std::vector<std::vector<Double_t>>());
        thit.push_back(std::vector<std::vector<Double_t>>());
        x.push_back(std::vector<std::vector<Double_t>>());
        y.push_back(std::vector<std::vector<Double_t>>());
        yToT.push_back(std::vector<std::vector<Double_t>>());
        for (Int_t j = 0; j <= N_TOFD_HIT_PADDLE_MAX; j++)
        {
            vmultihits[i][j] = 0;
            q[i].push_back(std::vector<Double_t>());
            thit[i].push_back(std::vector<Double_t>());
            x[i].push_back(std::vector<Double_t>());
            y[i].push_back(std::vector<Double_t>());
            yToT[i].push_back(std::vector<Double_t>());
        }
    }

    //    std::cout<<"new event!*************************************\n";

    Int_t nHits = fCalItems->GetEntries();
    if (nHits == 0)
        events_wo_tofd_hits++;

    Int_t nHitsEvent = 0;
    // Organize cals into bars.
    struct Entry
    {
        std::vector<R3BTofdCalData*> top;
        std::vector<R3BTofdCalData*> bot;
    };
    std::map<size_t, Entry> bar_map;
    // puts("Event");
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto* hit = dynamic_cast<R3BTofdCalData*>(fCalItems->At(ihit));
        size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();

        /*std::cout << "Hits: " << hit->GetDetectorId() << ' ' << hit->GetBarId() << ' ' << hit->GetSideId() << '  '
                  << hit->GetTimeLeading_ns() << ' ' << hit->GetTimeTrailing_ns() << ' '
                  << hit->GetTimeTrailing_ns() - hit->GetTimeLeading_ns() << '\n';
        */
        auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
        auto& vec = 1 == hit->GetSideId() ? ret.first->second.top : ret.first->second.bot;
        vec.push_back(hit);
        events_in_cal_level++;
    }

    static bool s_was_trig_missing = false;
    auto trig_num = fCalTriggerItems->GetEntries();
    // Find coincident PMT hits.
    // std::cout << "Print:\n";
    for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
    {
        //    reset:
        // for (auto it2 = it->second.top.begin(); it->second.top.end() != it2; ++it2) {
        // std::cout << "Top: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        // }
        // for (auto it2 = it->second.bot.begin(); it->second.bot.end() != it2; ++it2) {
        // std::cout << "Bot: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        // }
        auto const& top_vec = it->second.top;
        auto const& bot_vec = it->second.bot;
        size_t top_i = 0;
        size_t bot_i = 0;
        for (; top_i < top_vec.size() && bot_i < bot_vec.size();)
        {
            auto top = top_vec.at(top_i);
            auto bot = bot_vec.at(bot_i);
            auto top_trig_i = g_tofd_trig_map[top->GetDetectorId() - 1][top->GetSideId() - 1][top->GetBarId() - 1];
            auto bot_trig_i = g_tofd_trig_map[bot->GetDetectorId() - 1][bot->GetSideId() - 1][bot->GetBarId() - 1];
            Double_t top_trig_ns = 0, bot_trig_ns = 0;
            if (top_trig_i < trig_num && bot_trig_i < trig_num)
            {
                auto top_trig = dynamic_cast<R3BTofdCalData const*>(fCalTriggerItems->At(top_trig_i));
                auto bot_trig = dynamic_cast<R3BTofdCalData const*>(fCalTriggerItems->At(bot_trig_i));
                top_trig_ns = top_trig->GetTimeLeading_ns();
                bot_trig_ns = bot_trig->GetTimeLeading_ns();
                /*
                                std::cout << "Top: " << top->GetDetectorId() << ' ' << top->GetSideId() << ' ' <<
                   top->GetBarId() << ' '
                                << top_trig_i << ' ' << top_trig->GetTimeLeading_ns() << std::endl;
                                std::cout << "Bot: " <<
                                bot->GetDetectorId() << ' ' << bot->GetSideId() << ' ' << bot->GetBarId() << ' ' <<
                   bot_trig_i << ' '
                                << bot_trig->GetTimeLeading_ns() << std::endl;
                */
                ++n1;
            }
            else
            {
                if (!s_was_trig_missing)
                {
                    LOG(error) << "R3BTofdCal2HitS454Par::Exec() : Missing trigger information!";
                    s_was_trig_missing = true;
                }
                ++n2;
            }

            // Shift the cyclic difference window by half a window-length and move it back,
            // this way the trigger time will be at 0.
            auto top_ns =
                fmod(top->GetTimeLeading_ns() - top_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) - c_range_ns / 2;
            auto bot_ns =
                fmod(bot->GetTimeLeading_ns() - bot_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) - c_range_ns / 2;
            /*
                        if(top_ns>2000 || bot_ns>2000){
                            std::cout << top->GetTimeLeading_ns() << ' ' << top_trig_ns << ' ' << top_ns << std::endl;
                            std::cout << bot->GetTimeLeading_ns() << ' ' << bot_trig_ns << ' ' << bot_ns << std::endl;
                        }
            */
            auto dt = top_ns - bot_ns;
            // Handle wrap-around.
            auto dt_mod = fmod(dt + c_range_ns, c_range_ns);
            if (dt < 0)
            {
                // We're only interested in the short time-differences, so we
                // want to move the upper part of the coarse counter range close
                // to the lower range, i.e. we cut the middle of the range and
                // glue zero and the largest values together.
                dt_mod -= c_range_ns;
            }
            // std::cout << top_i << ' ' << bot_i << ": " << top_ns << ' ' << bot_ns << " = " << dt << ' ' <<
            // std::abs(dt_mod) << '\n';
            if (std::abs(dt_mod) < c_bar_coincidence_ns)
            {
                inbarcoincidence++;
                // Hit!
                // std::cout << "Hit!\n";
                Int_t iPlane = top->GetDetectorId(); // 1..n
                Int_t iBar = top->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    LOG(error) << "R3BTofdCal2HitS454Par::Exec() : more detectors than expected! Det: " << iPlane
                               << " allowed are 1.." << fNofPlanes;
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    LOG(error) << "R3BTofdCal2HitS454Par::Exec() : more bars then expected! Det: " << iBar
                               << " allowed are 1.." << fPaddlesPerPlane;
                    continue;
                }

                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                // std::cout<<"ToT: "<<top_tot << " "<<bot_tot<<"\n";

                // register multi hits
                vmultihits[iPlane][iBar] += 1;

                nHitsEvent += 1;
                R3BTofdHitModulePar* par = fHitPar->GetModuleParAt(iPlane, iBar);
                if (!par)
                {
                    LOG(info) << "R3BTofdCal2HitS454::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                              << ", Bar: " << top->GetBarId();
                    continue;
                }
                // walk corrections
                if (par->GetPar1Walk() == 0. || par->GetPar2Walk() == 0. || par->GetPar3Walk() == 0. ||
                    par->GetPar4Walk() == 0. || par->GetPar5Walk() == 0.)
                    LOG(info) << "Walk correction not found!";
                auto bot_ns_walk = bot_ns - walk(bot_tot,
                                                 par->GetPar1Walk(),
                                                 par->GetPar2Walk(),
                                                 par->GetPar3Walk(),
                                                 par->GetPar4Walk(),
                                                 par->GetPar5Walk());
                auto top_ns_walk = top_ns - walk(top_tot,
                                                 par->GetPar1Walk(),
                                                 par->GetPar2Walk(),
                                                 par->GetPar3Walk(),
                                                 par->GetPar4Walk(),
                                                 par->GetPar5Walk());
                // calculate tdiff
                auto tdiff = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2()));

                // calculate time of hit
                Double_t THit = (bot_ns + top_ns) / 2. - par->GetSync();
                if (std::isnan(THit))
                {
                    LOG(fatal) << "ToFD THit not found";
                }
                if (timeP0 == 0.)
                    timeP0 = THit;

                thit[iPlane][iBar].push_back(THit);

                // calculate y-position
                auto pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();
                y[iPlane][iBar].push_back(pos);

                // calculate y-position from ToT
                auto posToT =
                    par->GetLambda() * log((top_tot * par->GetToTOffset2()) / (bot_tot * par->GetToTOffset1()));
                yToT[iPlane][iBar].push_back(posToT);

                if (fTofdTotPos)
                    pos = posToT;

                Float_t paddle_width = 2.70000;
                Float_t paddle_thickness = 0.50000;
                Float_t air_gap_paddles = 0.04;
                Float_t air_gap_layer = 5.;
                // define number of layers and paddles with sizes of the detector
                Int_t number_layers = 2;   // number of layers
                Int_t number_paddles = 44; // number of paddles per layer
                Float_t detector_width =
                    number_paddles * paddle_width + (number_paddles - 1) * air_gap_paddles + paddle_width;
                Float_t detector_thickness = (number_layers - 1) * air_gap_layer + number_layers * paddle_thickness;

                // calculate x-position
                if (iPlane == 1 || iPlane == 3)
                {
                    // x[iPlane][iBar].push_back(iBar * 2.8 - 23. * 2.8);
                    x[iPlane][iBar].push_back(-detector_width / 2 + (paddle_width + air_gap_paddles) / 2 +
                                              (iBar - 1) * (paddle_width + air_gap_paddles) - 0.04);
                    //					cout << "Test: " << iBar * 2.8 - 23. * 2.8 << "  " <<
                    //					-detector_width/2 + (paddle_width+air_gap_paddles)/2 +
                    //(iBar-1)*(paddle_width+air_gap_paddles) - 0.04<<endl;
                }
                if (iPlane == 2 || iPlane == 4)
                {
                    // x[iPlane][iBar].push_back(iBar * 2.8 - 23. * 2.8 + 1.4);
                    x[iPlane][iBar].push_back(-detector_width / 2 + (paddle_width + air_gap_paddles) +
                                              (iBar - 1) * (paddle_width + air_gap_paddles) - 0.04);
                }

                Double_t para[4];
                Double_t qb = 0.;
                if (fTofdTotPos)
                {
                    // via pol3
                    para[0] = par->GetPola();
                    para[1] = par->GetPolb();
                    para[2] = par->GetPolc();
                    para[3] = par->GetPold();
                    qb = TMath::Sqrt(top_tot * bot_tot) /
                         (para[0] + para[1] * pos + para[2] * pow(pos, 2) + para[3] * pow(pos, 3));
                    qb = qb * fTofdQ;
                }
                else
                {
                    // via double exponential:
                    para[0] = par->GetPar1a();
                    para[1] = par->GetPar1b();
                    para[2] = par->GetPar1c();
                    para[3] = par->GetPar1d();
                    auto q1 =
                        bot_tot / (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                    para[0] = par->GetPar2a();
                    para[1] = par->GetPar2b();
                    para[2] = par->GetPar2c();
                    para[3] = par->GetPar2d();
                    auto q2 =
                        top_tot / (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                    q1 = q1 * fTofdQ;
                    q2 = q2 * fTofdQ;
                    qb = (q1 + q2) / 2.;
                }

                Double_t parz[3];
                parz[0] = par->GetPar1za();
                parz[1] = par->GetPar1zb();
                parz[2] = par->GetPar1zc();

                if (parz[0] > 0 && parz[2] > 0)
                    LOG(debug) << "Charges in this event " << parz[0] * TMath::Power(qb, parz[2]) + parz[1] << " plane "
                               << iPlane << " ibar " << iBar;
                else
                    LOG(debug) << "Charges in this event " << qb << " plane " << iPlane << " ibar " << iBar;
                LOG(debug) << "Times in this event " << THit << " plane " << iPlane << " ibar " << iBar;
                if (iPlane == 1 || iPlane == 3)
                    LOG(debug) << "x in this event "
                               << -detector_width / 2 + (paddle_width + air_gap_paddles) / 2 +
                                      (iBar - 1) * (paddle_width + air_gap_paddles) - 0.04
                               << " ibar " << iBar;
                if (iPlane == 2 || iPlane == 4)
                    LOG(debug) << "x in this event "
                               << -detector_width / 2 + (paddle_width + air_gap_paddles) +
                                      (iBar - 1) * (paddle_width + air_gap_paddles) - 0.04
                               << " plane " << iPlane << " ibar " << iBar;
                LOG(debug) << "y in this event " << pos << " plane " << iPlane << " ibar " << iBar << "\n";

                if (parz[0] > 0 && parz[2] > 0)
                {
                    q[iPlane][iBar].push_back(parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                }
                else
                {
                    q[iPlane][iBar].push_back(qb);
                    parz[0] = 1.;
                    parz[1] = 0.;
                    parz[2] = 1.;
                }

                if (fTofdHisto)
                {
                    // fill control histograms
                    CreateHistograms(iPlane, iBar);
                    fhTsync[iPlane - 1]->Fill(iBar, THit);
                    fhTdiff[iPlane - 1]->Fill(iBar, tdiff);
                    fhQvsPos[iPlane - 1][iBar - 1]->Fill(pos, parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                    // fhQvsTHit[iPlane - 1][iBar - 1]->Fill(qb, THit);
                    // fhTvsTHit[iPlane - 1][iBar - 1]->Fill(dt_mod, THit);
                }

                ++top_i;
                ++bot_i;
            }
            else if (dt < 0 && dt > -c_range_ns / 2)
            {
                ++top_i;
                LOG(debug) << "Not in bar coincidence increase top counter";
            }
            else
            {
                ++bot_i;
                LOG(debug) << "Not in bar coincidence increase bot counter";
            }
        }
    }

    // Now all hits in this event are analyzed

    LOG(debug) << "Hits in this event: " << nHitsEvent;
    //    if (nHitsEvent == 0)
    //        events_wo_tofd_hits++;

    // init arrays to store hits
    Double_t tArrQ[nHitsEvent + 1];
    Double_t tArrT[nHitsEvent + 1];
    Double_t tArrX[nHitsEvent + 1];
    Double_t tArrY[nHitsEvent + 1];
    Double_t tArrYT[nHitsEvent + 1];
    Double_t tArrP[nHitsEvent + 1];
    Double_t tArrB[nHitsEvent + 1];
    Bool_t tArrU[nHitsEvent + 1];
    for (int i = 0; i < (nHitsEvent + 1); i++)
    {
        tArrQ[i] = -1.;
        tArrT[i] = -1.;
        tArrX[i] = -1.;
        tArrY[i] = -1.;
        tArrYT[i] = -1.;
        tArrP[i] = -1.;
        tArrB[i] = -1.;
        tArrU[i] = kFALSE;
    }

    for (Int_t i = 1; i <= fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane + 1; j++)
        {
            if (vmultihits[i][j] > 1)
            {
                bars_with_multihit++;
                multihit += vmultihits[i][j] - 1;
            }
        }
    }

    // sort hits for time
    for (Int_t i = 1; i <= fNofPlanes; i++)
    { // loop over planes i
        for (Int_t j = 1; j < fPaddlesPerPlane + 1; j++)
        { // loop over virtual paddles j
            if (thit[i][j].empty() == false)
            { // check paddle for entries
                for (Int_t m = 0; m < thit[i][j].size(); m++)
                { // loop over multihits m
                    Int_t p = 0;
                    if (tArrT[0] == -1.)
                    { // first entry
                        LOG(debug) << "First entry plane/bar " << i << "/" << j;
                        tArrQ[0] = q[i][j].at(m);
                        tArrT[0] = thit[i][j].at(m);
                        tArrX[0] = x[i][j].at(m);
                        tArrY[0] = y[i][j].at(m);
                        tArrYT[0] = yToT[i][j].at(m);
                        tArrP[0] = i;
                        tArrB[0] = j;
                    }
                    else
                    {
                        if (thit[i][j].at(m) < tArrT[0])
                        { // new first entry with smaller time
                            LOG(debug) << "Insert new first " << i << "/" << j;
                            insertX(nHitsEvent, tArrQ, q[i][j].at(m), 1);
                            insertX(nHitsEvent, tArrT, thit[i][j].at(m), 1);
                            insertX(nHitsEvent, tArrX, x[i][j].at(m), 1);
                            insertX(nHitsEvent, tArrY, y[i][j].at(m), 1);
                            insertX(nHitsEvent, tArrYT, yToT[i][j].at(m), 1);
                            insertX(nHitsEvent, tArrP, i, 1);
                            insertX(nHitsEvent, tArrB, j, 1);
                        }
                        else
                        {
                            while (thit[i][j].at(m) > tArrT[p] && tArrT[p] != -1.)
                            {
                                p++; // find insert position
                                if (p > nHitsEvent + 1)
                                    LOG(fatal) << "Insert position oor"; // should not happen
                            }

                            LOG(debug) << "Will insert at " << p;
                            if (p > 0 && thit[i][j].at(m) > tArrT[p - 1] && thit[i][j].at(m) != tArrT[p])
                            { // insert at right position
                                LOG(debug) << "Insert at " << p << " " << i << "/" << j;
                                insertX(nHitsEvent, tArrQ, q[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrT, thit[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrX, x[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrY, y[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrYT, yToT[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrP, i, p + 1);
                                insertX(nHitsEvent, tArrB, j, p + 1);
                            }
                            else
                            {

                                LOG(error) << "Insert event with exact same time " << fnEvents;
                                p = p + 1;
                                LOG(debug) << "Insert at " << p << " " << i << "/" << j;
                                insertX(nHitsEvent, tArrQ, q[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrT, thit[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrX, x[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrY, y[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrYT, yToT[i][j].at(m), p + 1);
                                insertX(nHitsEvent, tArrP, i, p + 1);
                                insertX(nHitsEvent, tArrB, j, p + 1);
                            }
                        }
                    }
                }
            }
        }
    }

    // Now we have all hits in this event time sorted

    // print events
    std::stringstream ss;
    ss << "Time sorted hits:\n [CHARGE] ";
    if (tArrT[0] != -1.)
    {
        for (Int_t a = 0; a < nHitsEvent; a++)
            ss << tArrQ[a] << " ";
        ss << "\n [TIME] ";
        for (Int_t a = 0; a < nHitsEvent; a++)
            ss << tArrT[a] << " ";
        ss << "\n [XPOS] ";
        for (Int_t a = 0; a < nHitsEvent; a++)
            ss << tArrX[a] << " ";
        ss << "\n [YPOS] ";
        for (Int_t a = 0; a < nHitsEvent; a++)
        {
            if (fTofdTotPos)
                ss << tArrYT[a] << " ";
            else
                ss << tArrY[a] << " ";
        }
        ss << "\n [PLANE] ";
        for (Int_t a = 0; a < nHitsEvent; a++)
            ss << tArrP[a] << " ";
        ss << "\n [BAR] ";
        for (Int_t a = 0; a < nHitsEvent; a++)
            ss << tArrB[a] << " ";
    }
    LOG(debug) << ss.str() << "\n";

    // Now we can analyze the hits in this event

    if (fTofdHisto)
    {
        for (Int_t a = 0; a < nHitsEvent; a++)
        { // loop over all hits
            eventstore++;
            fhQ[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrQ[a]);        // charge per plane
            fhQvsEvent[((Int_t)tArrP[a]) - 1]->Fill(fnEvents, tArrQ[a]); // charge vs event #
            if (fTofdTotPos)
            {
                fhxy[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrYT[a]); // xy of plane
                // std::cout<<"Use ToT y Position\n";
            }
            else
            {
                fhxy[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrY[a]); // xy of plane
                // std::cout<<"Use DEXP y Position\n";
            }
        }
    }

    for (Int_t hit = 0; hit < nHitsEvent; hit++)
    { // loop over hits
        if (tArrU[hit] == false)
        {
            // cout << "Single Hit for Plane " << tArrP[hit] << " " << tArrB[hit] << endl;
            ;
            tArrU[hit] = true;
            // store single hits
            singlehit++;
            if (fTofdTotPos)
            {
                new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(
                    tArrT[hit], tArrX[hit], tArrYT[hit], tArrQ[hit], -1., tArrQ[hit], tArrP[hit], tArrB[hit]);
            }
            else
            {
                new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(
                    tArrT[hit], tArrX[hit], tArrY[hit], tArrQ[hit], -1., tArrQ[hit], tArrP[hit], tArrB[hit]);
            }
        }
    }

    LOG(debug) << "Used up hits in this event:";
    for (Int_t a = 0; a < nHitsEvent; a++)
    {
        LOG(debug) << "Event " << a << " " << tArrU[a] << " ";
        if (tArrU[a] != true)
            LOG(fatal) << "Not all events analyzed!";
    }

    LOG(debug) << "------------------------------------------------------\n";
    fnEvents++;
}

void R3BTofdCal2HitS454::CreateHistograms(Int_t iPlane, Int_t iBar)
{
    Double_t max_charge = 80.;
    // create histograms if not already existing

    if (NULL == fhTsync[iPlane - 1])
    {
        char strName[255];
        sprintf(strName, "Time_Sync_Plane_%d", iPlane);
        fhTsync[iPlane - 1] = new TH2F(strName, "", 50, 0, 50, 10000, -2000., 2000.);
        fhTsync[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTsync[iPlane - 1]->GetYaxis()->SetTitle("THit in ns");
    }

    if (NULL == fhTdiff[iPlane - 1])
    {
        char strName1[255];
        char strName2[255];
        sprintf(strName1, "Time_Diff_Plane_%d", iPlane);
        sprintf(strName2, "Time Diff Plane %d", iPlane);
        fhTdiff[iPlane - 1] = new TH2F(strName1, strName2, 50, 0, 50, 400, -8., 8.);
        fhTdiff[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTdiff[iPlane - 1]->GetYaxis()->SetTitle("Time difference (PM1 - PM2) in ns");
    }

    if (NULL == fhQvsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 400, -100, 100, max_charge * 10, 0., max_charge);
        fhQvsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Position in cm");
    }
    /*
    if (NULL == fhQvsTHit[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_THit_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsTHit[iPlane - 1][iBar - 1] = new TH2F(strName, "", 1000, 0., max_charge, 1000, -10, 40);
        fhQvsTHit[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("THit in ns");
        fhQvsTHit[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Charge");
    }
    */
    /*
    if (NULL == fhTvsTHit[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "T_vs_ToF_Plane_%d_Bar_%d", iPlane, iBar);
        fhTvsTHit[iPlane - 1][iBar - 1] = new TH2F(strName, "", 625, -25, 25, 1000, -10, 40);
        fhTvsTHit[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("THit in ns");
        fhTvsTHit[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("T1-T2 in ns");
    }
    */
    if (NULL == fhQ[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "Q_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Q Plane %d ", iPlane);
        fhQ[iPlane - 1] = new TH2F(strName1, strName2, 44, 0, 44, max_charge * 10, 0., max_charge);
        fhQ[iPlane - 1]->GetYaxis()->SetTitle("Charge");
        fhQ[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
    }
    /*
    if (NULL == fhxy[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "xy_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "xy of Plane %d ", iPlane);
        fhxy[iPlane - 1] = new TH2F(strName1, strName2, 320, -80, 80, 400, -100., 100.);
        fhxy[iPlane - 1]->GetYaxis()->SetTitle("y-position in cm");
        fhxy[iPlane - 1]->GetXaxis()->SetTitle("x-position in cm");
    }
    */
    if (NULL == fhxy[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "xy_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "xy of Plane %d ", iPlane);
        fhxy[iPlane - 1] = new TH2F(strName1, strName2, 46, 0, 46, 400, -100., 100.);
        fhxy[iPlane - 1]->GetYaxis()->SetTitle("y-position in cm");
        fhxy[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
    }

    if (NULL == fhQvsEvent[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "QvsEvent_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Charge vs Event # Plane %d ", iPlane);
        fhQvsEvent[iPlane - 1] = new TH2F(strName1, strName2, 2e5, 0, 2e9, max_charge * 10, 0., max_charge);
        fhQvsEvent[iPlane - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsEvent[iPlane - 1]->GetXaxis()->SetTitle("Event #");
    }
}
void R3BTofdCal2HitS454::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
        fNofHitItems = 0;
    }
}

void R3BTofdCal2HitS454::FinishTask()
{
    if (fTofdHisto)
    {
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            if (fhQ[i])
                fhQ[i]->Write();
            if (fhxy[i])
                fhxy[i]->Write();
            if (fhQvsEvent[i])
                fhQvsEvent[i]->Write();
            if (fhTdiff[i])
                fhTdiff[i]->Write();
            if (fhTsync[i])
                fhTsync[i]->Write();
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {

                // control histogram time particles
                if (fhQvsPos[i][j])
                    fhQvsPos[i][j]->Write();
                /*
                if (fhQvsTHit[i][j])
                    fhQvsTHit[i][j]->Write();
                */
                /*
                if (fhTvsTHit[i][j])
                    fhTvsTHit[i][j]->Write();
                */
            }
        }
    }
    std::cout << "\n\nSome statistics:\n"
              << "Total number of events in tree   " << maxevent << "\n"
              << "Max Event analyzed               " << fnEvents + wrongtrigger + wrongtpat << "\n"
              << "Events skipped due to trigger    " << wrongtrigger << "\n"
              << "Events skipped due to tpat       " << wrongtpat << "\n"
              << "Events with correct header&tpat  " << headertpat << "\n"
              << "Events without ToFd hits         " << events_wo_tofd_hits << "\n"
              << "Events in cal level              " << events_in_cal_level << "\n"
              << "Hits in bar coincidence          " << inbarcoincidence << "\n"
              << "Bars with multihits              " << bars_with_multihit << "\n"
              << "Multihits                        " << multihit << "\n"
              << "Events stored                    " << eventstore << " <-> " << inbarcoincidence
              << " (Events in bar coincidence)\n"
              << "Events in single planes          " << singlehit << "\n"
              << "Good events in total             " << eventstore << " <-> " << singlehit << " = singlehit \n";

    std::cout << "n1=" << n1 << " n2=" << n2 << std::endl;
}

Double_t R3BTofdCal2HitS454::betaCorr(Double_t delta)
{
    //    Double_t corr=-3.*delta;  //correction for Ag

    Double_t corr = -1. * delta; // correction for 12C
    corr = 0.;
    return corr;
}
/* old method
Double_t R3BTofdCal2HitS454::walk(Double_t q)
{
    Double_t y;
    //
    //   Double_t p0 = 18.;
    //    Double_t p1 = -0.5;
    //
    //    y = p0 * TMath::Power(q,p1);

    Double_t par1, par2, par3, par4, par5;
    Int_t voltage = 444;

    if (voltage == 444)
    { // voltage set in s444
        par1 = 2.178871e+01;
        par2 = -3.565959e-03;
        par3 = 5.713045e+01;
        par4 = 4.007571e-02;
        par5 = -9.537515e-05;
    }
    if (voltage == 500)
    {
        par1 = 1.64344e+01;
        par2 = 2.84000e-01;
        par3 = 3.47659e+02;
        par4 = -2.70050e-01;
        par5 = 3.61515e-04;
    }
    if (voltage == 600)
    {
        par1 = 1.22606e+01;
        par2 = 3.12697e-01;
        par3 = 4.40109e+02;
        par4 = -1.86328e-01;
        par5 = 1.49519e-04;
    }
    y = -30.2 + par1 * TMath::Power(q, par2) + par3 / q + par4 * q + par5 * q * q;
    return y;
}
*/
Double_t R3BTofdCal2HitS454::walk(Double_t Q,
                                  Double_t par1,
                                  Double_t par2,
                                  Double_t par3,
                                  Double_t par4,
                                  Double_t par5) // new method
{
    Double_t y = 0;
    y = -30.2 + par1 * TMath::Power(Q, par2) + par3 / Q + par4 * Q + par5 * Q * Q;
    return y;
}
Double_t R3BTofdCal2HitS454::saturation(Double_t x)
{
    Double_t kor;
    Int_t voltage = 600;
    if (voltage == 600)
    {
        if (x < 173)
        {
            kor = 0.;
        }
        else if (x > 208)
        {
            kor = -1.73665e+03 + 2.82009e+01 * 208. - 1.53846e-01 * (208. * 208.) + 2.82425e-04 * (208. * 208. * 208.);
        }
        else
        {
            kor = -1.73665e+03 + 2.82009e+01 * x - 1.53846e-01 * (x * x) + 2.82425e-04 * (x * x * x);
        }
    }
    if (voltage == 500)
    {
        if (x < 95.5)
        {
            kor = 0.;
        }
        else if (x > 124)
        {
            kor = 1.08 * x - 112.44;
        }
        else
        {
            kor = 643.257 - 16.7823 * x + 0.139822 * (x * x) - 0.000362154 * (x * x * x);
        }
    }
    if (voltage == 700)
    {
        if (x < 198)
        {
            kor = 0.;
        }
        else if (x > 298)
        {
            kor = 0.21 * x - 45.54;
        }
        else
        {
            kor = -19067 + 383.93 * x - 3.05794 * (x * x) + 0.0120429 * (x * x * x) - 2.34619e-05 * (x * x * x * x) +
                  1.81076e-08 * (x * x * x * x * x);
        }
    }
    return kor;
}

Double_t* R3BTofdCal2HitS454::insertX(Int_t n, Double_t arr[], Double_t x, Int_t pos)
{
    Int_t i;

    // increase the size by 1
    n++;

    // shift elements forward
    for (i = n; i >= pos; i--)
        arr[i] = arr[i - 1];

    // insert x at pos
    arr[pos - 1] = x;

    return arr;
}

ClassImp(R3BTofdCal2HitS454)
