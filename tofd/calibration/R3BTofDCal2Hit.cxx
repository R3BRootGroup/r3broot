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
// -----                 R3BTofDCal2Hit                   -----
// -----    Created 18/03/22 by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------

#include "R3BTofDCal2Hit.h"
#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTCalEngine.h"
#include "R3BTofDHitModulePar.h"
#include "R3BTofDHitPar.h"
#include "R3BTofDMappingPar.h"
#include "R3BTofdCalData.h"
#include "R3BTofdHitData.h"
#include <FairRootManager.h>

#include "FairLogger.h"
#include "FairRuntimeDb.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"

#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom3.h"

#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
    uint64_t n1, n2;
} // namespace

R3BTofDCal2Hit::R3BTofDCal2Hit()
    : R3BTofDCal2Hit("TofdCal2Hit", 1)
{
}

R3BTofDCal2Hit::R3BTofDCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems(nullptr)
    , fCalTriggerItems(nullptr)
    , fHitItems(new TClonesArray("R3BTofdHitData"))
    , fNofHitPars(0)
    , fHitPar(nullptr)
    , fMapPar(nullptr)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(44)
    , fTofdQ(0)
    , fTofdHisto(false)
    , fTofdTotPos(false)
    , fnEvents(0)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , maxevent(0)
    , wrongtrigger(0)
    , wrongtpat(0)
    , notpat(0)
    , headertpat(0)
    , events_in_cal_level(0)
    , inbarcoincidence(0)
    , eventstore(0)
    , singlehit(0)
    , multihit(0)
    , bars_with_multihit(0)
    , events_wo_tofd_hits(0)
    , lasttpatevent(0)
    , goodpair(0)
    , goodpair4(0)
    , goodpair3(0)
    , goodpair1(0)
    , goodpair2(0)
    , goodpair5(0)
    , goodpair6(0)
    , goodpair7(0)
{
    fhNoTpat = nullptr;
    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        fhQ[i] = nullptr;
        fhxy[i] = nullptr;
        fhQvsEvent[i] = nullptr;
        fhTdiff[i] = nullptr;
        fhTsync[i] = nullptr;
        fhQ0Qt[i] = nullptr;
        fhTvsQ[i] = nullptr;
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            fhQvsPos[i][j] = nullptr;
        }
    }
}

R3BTofDCal2Hit::~R3BTofDCal2Hit()
{
    if (fTofdHisto)
    {
        if (fhNoTpat)
            delete fhNoTpat;
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
            if (fhQ0Qt[i])
                delete fhQ0Qt[i];
            if (fhTvsQ[i])
                delete fhTvsQ[i];
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
        fHitItems = nullptr;
    }
}

void R3BTofDCal2Hit::SetParContainers()
{
    fMapPar = dynamic_cast<R3BTofDMappingPar*>(FairRuntimeDb::instance()->getContainer("tofdMappingPar"));
    R3BLOG_IF(warn, !fMapPar, "Could not get access to tofdMappingPar container");

    fHitPar = dynamic_cast<R3BTofDHitPar*>(FairRuntimeDb::instance()->getContainer("tofdHitPar"));
    if (!fHitPar)
    {
        R3BLOG(error, "Could not get access to tofdHitPar container");
        fNofHitPars = 0;
        return;
    }
}

void R3BTofDCal2Hit::SetParameter()
{
    R3BLOG_IF(info, fMapPar, "Nb of planes " << fMapPar->GetNbPlanes() << " and paddles " << fMapPar->GetNbPaddles());

    if (fHitPar)
        fNofHitPars = fHitPar->GetNumModulePar();
    else
        fNofHitPars = 0;

    R3BLOG(info, "Parameters in the tofdHitPar container: " << fNofHitPars);

    return;
}

InitStatus R3BTofDCal2Hit::Init()
{
    R3BLOG(info, "");
    auto* mgr = FairRootManager::Instance();
    if (!mgr)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    R3BLOG_IF(fatal, nullptr == header, "EventHeader. not found");

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdCal"));
    R3BLOG_IF(fatal, nullptr == fCalItems, "TofdCal not found");

    fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdTriggerCal"));
    R3BLOG_IF(fatal, nullptr == fCalTriggerItems, "TofdTriggerCal not found");

    maxevent = mgr->CheckMaxEventNo();

    // request storage of Hit data in output tree
    mgr->Register("TofdHit", "TofdHit data", fHitItems, !fOnline);

    SetParameter();

    if (fTofdHisto)
    {
        for (Int_t i = 1; i <= fNofPlanes; i++)
            for (Int_t j = 1; j <= N_TOFD_HIT_PADDLE_MAX; j++)
                CreateHistograms(i, j);
    }

    // Definition of a time stich object to correlate times coming from different systems
    fTimeStitch = std::make_unique<R3BCoarseTimeStitch>();

    return kSUCCESS;
}

InitStatus R3BTofDCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BTofDCal2Hit::Exec(Option_t* option)
{
    static uint32_t counter = 0;
    ++counter;

    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
    {
        wrongtrigger++;
        return;
    }
    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit1 = fTpat1 - 1;
    Int_t fTpat_bit2 = fTpat2 - 1;
    Int_t tpatbin;
    Int_t tpatsum = 0;

    if (header && fTpat1 >= 0 && fTpat2 >= 0)
    {
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            // LOG(debug)<<"tpatbin "tpatbin;
            if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2))
            {
                wrongtpat++;
                return;
            }
            if (tpatbin != 0)
            {
                //  fhTpat->Fill(i+1);
                LOG(debug) << "Accepted Tpat: " << i + 1;
            }
            tpatsum += tpatbin;
        }
        if (tpatsum < 1)
        {
            if (fhNoTpat)
                fhNoTpat->Fill(counter - lasttpatevent);

            lasttpatevent = counter;
            notpat++;
            // LOG(fatal)<<"No Tpat info";
            return;
        }
    }

    headertpat++;
    Double_t timeP0 = 0.;
    Double_t randx;

    UInt_t vmultihits[N_PLANE_MAX + 1][N_TOFD_HIT_PADDLE_MAX + 1];
    for (Int_t i = 0; i <= fNofPlanes; i++)
    {
        for (Int_t j = 0; j <= N_TOFD_HIT_PADDLE_MAX; j++)
        {
            vmultihits[i][j] = 0;
        }
    }
    struct hit
    {
        Double_t charge;
        Double_t time;
        Double_t xpos;
        Double_t ypos;
        Int_t plane;
        Int_t bar;
        Double_t time_raw;
        Double_t tof;
    };

    //    std::cout<<"new event!*************************************\n";

    std::vector<hit> event;

    Int_t nHits = fCalItems->GetEntriesFast();
    LOG(debug) << "Leading and trailing edges in this event: " << nHits;
    if (nHits == 0)
        events_wo_tofd_hits++;

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
        size_t idx = (hit->GetDetectorId() - 1) * fPaddlesPerPlane + (hit->GetBarId() - 1);
        auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
        auto& vec = 1 == hit->GetSideId() ? ret.first->second.bot : ret.first->second.top;
        vec.push_back(hit);
        events_in_cal_level++;
    }

    // Build trigger map.
    std::vector<R3BTofdCalData const*> trig_map;
    for (int i = 0; i < fCalTriggerItems->GetEntriesFast(); ++i)
    {
        auto trig = dynamic_cast<R3BTofdCalData const*>(fCalTriggerItems->At(i));
        if (trig_map.size() < trig->GetBarId())
        {
            trig_map.resize(trig->GetBarId());
        }
        trig_map.at(trig->GetBarId() - 1) = trig;
    }

    bool s_was_trig_missing = false;
    // Find coincident PMT hits.
    // std::cout << "Print:\n";
    for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
    {
        auto const& top_vec = it->second.top;
        auto const& bot_vec = it->second.bot;
        size_t top_i = 0;
        size_t bot_i = 0;
        for (; top_i < top_vec.size() && bot_i < bot_vec.size();)
        {
            auto top = top_vec.at(top_i);
            auto bot = bot_vec.at(bot_i);

            Int_t top_trig_i = 0;
            Int_t bot_trig_i = 0;
            if (fMapPar)
            {
                top_trig_i = fMapPar->GetTrigMap(top->GetDetectorId(), top->GetBarId(), top->GetSideId());
                bot_trig_i = fMapPar->GetTrigMap(bot->GetDetectorId(), bot->GetBarId(), bot->GetSideId());
            }

            Double_t top_trig_ns = 0, bot_trig_ns = 0;
            if (top_trig_i < trig_map.size() && trig_map.at(top_trig_i) && bot_trig_i < trig_map.size() &&
                trig_map.at(bot_trig_i))
            {
                auto top_trig = trig_map.at(top_trig_i);
                auto bot_trig = trig_map.at(bot_trig_i);
                top_trig_ns = top_trig->GetTimeLeading_ns();
                bot_trig_ns = bot_trig->GetTimeLeading_ns();

                ++n1;
            }
            else
            {
                if (!s_was_trig_missing)
                {
                    R3BLOG(error, "Missing trigger information!");
                    R3BLOG(error, "Top: " << top->GetDetectorId() << ' ' << top->GetSideId() << ' ' << top->GetBarId());
                    R3BLOG(error, "Bot: " << bot->GetDetectorId() << ' ' << bot->GetSideId() << ' ' << bot->GetBarId());
                    s_was_trig_missing = true;
                }
                ++n2;
            }

            // Shift the cyclic difference window by half a window-length and move it back,
            // this way the trigger time will be at 0.
            auto top_ns = fTimeStitch->GetTime(top->GetTimeLeading_ns() - top_trig_ns);
            auto bot_ns = fTimeStitch->GetTime(bot->GetTimeLeading_ns() - bot_trig_ns);

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
                Int_t iPlane = top->GetDetectorId(); // 1..n
                Int_t iBar = top->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    R3BLOG(error, "More detectors than expected! Det: " << iPlane << " allowed are 1.." << fNofPlanes);
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    R3BLOG(error, "More bars then expected! Det: " << iBar << " allowed are 1.." << fPaddlesPerPlane);
                    continue;
                }

                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                auto THit_raw = (bot->GetTimeLeading_ns() + top->GetTimeLeading_ns()) / 2.; // needed for TOF for ROLUs

                // std::cout<<"ToT: "<<top_tot << " "<<bot_tot<<"\n";

                // register multi hits
                vmultihits[iPlane][iBar] += 1;

                auto par = fHitPar->GetModuleParAt(iPlane, iBar);
                if (!par)
                {
                    R3BLOG(error, "Hit par not found, Plane: " << top->GetDetectorId() << ", Bar: " << top->GetBarId());
                    continue;
                }

                // walk corrections
                if (par->GetPar1Walk() == 0. || par->GetPar2Walk() == 0. || par->GetPar3Walk() == 0. ||
                    par->GetPar4Walk() == 0. || par->GetPar5Walk() == 0.)
                {
                    R3BLOG(debug, "TofD walk correction not found");
                }
                else
                {
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
                }

                // calculate tdiff
                auto tdiff = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2()));

                // calculate time of hit
                Double_t THit = (bot_ns + top_ns) / 2. - par->GetSync();
                if (std::isnan(THit))
                {
                    R3BLOG(fatal, "TofD THit not found");
                }
                if (timeP0 == 0.)
                    timeP0 = THit;

                // calculate y-position
                auto pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();

                // calculate y-position from ToT
                auto posToT =
                    par->GetLambda() * log((top_tot * par->GetToTOffset2()) / (bot_tot * par->GetToTOffset1()));

                if (fTofdTotPos)
                {
                    pos = posToT;
                }

                Float_t paddle_width = 2.700;
                Float_t air_gap_paddles = 0.04;
                Float_t air_gap_layer = 5.;
                Float_t detector_width =
                    fPaddlesPerPlane * paddle_width + (fPaddlesPerPlane - 1) * air_gap_paddles + paddle_width;
                Double_t xp = -1000.;
                // calculate x-position
                if (iPlane == 1 || iPlane == 3)
                {
                    xp = -detector_width / 2 + (paddle_width + air_gap_paddles) / 2 +
                         (iBar - 1) * (paddle_width + air_gap_paddles) +
                         gRandom->Uniform(-paddle_width / 2., paddle_width / 2.);
                }
                if (iPlane == 2 || iPlane == 4)
                {
                    xp = -detector_width / 2 + (paddle_width + air_gap_paddles) +
                         (iBar - 1) * (paddle_width + air_gap_paddles) +
                         gRandom->Uniform(-paddle_width / 2., paddle_width / 2.);
                }

                Double_t para[4];
                Double_t qb = 0.;
                if (fTofdQ > 0)
                {
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
                        auto q1 = bot_tot /
                                  (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                        para[0] = par->GetPar2a();
                        para[1] = par->GetPar2b();
                        para[2] = par->GetPar2c();
                        para[3] = par->GetPar2d();
                        auto q2 = top_tot /
                                  (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                        q1 = q1 * fTofdQ;
                        q2 = q2 * fTofdQ;
                        qb = (q1 + q2) / 2.;
                    }
                }
                else
                {
                    qb = TMath::Sqrt(top_tot * bot_tot);
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
                               << " plane " << iPlane << " ibar " << iBar;
                if (iPlane == 2 || iPlane == 4)
                    LOG(debug) << "x in this event "
                               << -detector_width / 2 + (paddle_width + air_gap_paddles) +
                                      (iBar - 1) * (paddle_width + air_gap_paddles) - 0.04
                               << " plane " << iPlane << " ibar " << iBar;
                LOG(debug) << "y in this event " << pos << " plane " << iPlane << " ibar " << iBar << "\n";

                // Tof with respect LOS detector
                auto tof = fTimeStitch->GetTime((bot_ns + top_ns) / 2. - header->GetTStart(), "tamex", "vftx");
                auto tof_corr = tof - par->GetTofSyncOffset();

                event.push_back(
                    { parz[0] + parz[1] * qb + parz[2] * qb * qb, THit, xp, pos, iPlane, iBar, THit_raw, tof_corr });

                if (fTofdHisto)
                {
                    // fill control histograms
                    fhTsync[iPlane - 1]->Fill(iBar, THit);
                    fhTdiff[iPlane - 1]->Fill(iBar, tdiff);
                    fhQvsPos[iPlane - 1][iBar - 1]->Fill(pos, parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                    // fhQvsTHit[iPlane - 1][iBar - 1]->Fill(qb, THit);
                    // fhTvsTHit[iPlane - 1][iBar - 1]->Fill(dt_mod, THit);
                }

                for (Int_t e = 0; e < event.size(); e++)
                {
                    LOG(debug) << event[e].charge << " " << event[e].time << " " << event[e].xpos << " "
                               << event[e].ypos << " " << event[e].plane << " " << event[e].bar;
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

    LOG(debug) << "Hits in this event: " << event.size();
    Bool_t tArrU[event.size() + 1];
    for (int i = 0; i < (event.size() + 1); i++)
        tArrU[i] = kFALSE;

    for (Int_t i = 1; i <= fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane + 1; j++)
        {
            if (vmultihits[i][j] > 1)
            {
                bars_with_multihit++;
                multihit += vmultihits[i][j] - 1;
                // if(vmultihits[i][j]>3)cout<<vmultihits[i][j] - 1<<"\n";
            }
        }
    }

    std::sort(event.begin(), event.end(), [](hit const& a, hit const& b) { return a.time < b.time; });
    // Now we have all hits in this event time sorted

    if (fTofdHisto)
    {
        LOG(debug) << "Charge Time xpos ypos plane bar";
        for (Int_t hit = 0; hit < event.size(); hit++)
        {
            LOG(debug) << event[hit].charge << " " << event[hit].time << " " << event[hit].xpos << " "
                       << event[hit].ypos << " " << event[hit].plane << " " << event[hit].bar;
            // if (event[hit].plane == 2 && (event[hit].bar < 21 || event[hit].bar > 24)) fhTvsQ[event[hit].plane -
            // 1]->Fill(event[hit].time-event[0].time,event[hit].charge);
            if (event[hit].plane == 2 && (event[hit].bar == 18))
                fhTvsQ[event[hit].plane - 1]->Fill(event[hit].time - event[0].time, event[hit].charge);
        }
    }

    // Now we can analyze the hits in this event

    // select events with feasible times
    Double_t hit_coinc = 20.; // coincidence window for hits in one event in ns. physics says max 250 ps
    Double_t time0;
    for (Int_t ihit = 0; ihit < event.size();)
    { // loop over all hits in this event
        LOG(debug) << "Set new coincidence window: " << event[ihit].plane << " " << event[ihit].bar << " "
                   << event[ihit].time << " " << event[ihit].charge;
        time0 = event[ihit].time;              // time of first hit in coincidence window
        Double_t charge0 = event[ihit].charge; // charge of first hit in coincidence window
        Int_t plane0 = event[ihit].plane;      // plane of first hit in coincidence window
        std::vector<Double_t> goodcharge;
        std::vector<Double_t> goodplane;
        std::vector<Double_t> goodbar;
        struct goodhit
        {
            Double_t goodq;
            Double_t goodp;
            Double_t goodb;
        };
        struct by_charge
        {
            bool operator()(goodhit const& a, goodhit const& b) const noexcept { return a.goodq < b.goodq; }
        };
        std::vector<goodhit> goodevents;

        while (event[ihit].time < time0 + hit_coinc)
        { // check if in coincidence window
            if (fTofdHisto)
            {
                if (event[ihit].plane == plane0 && charge0 != event[ihit].charge)
                {
                    fhQ0Qt[event[ihit].plane - 1]->Fill(charge0, event[ihit].charge);
                }
            }

            if ((event[ihit].charge > 5.5 && event[ihit].charge < 6.5) ||
                (event[ihit].charge > 1.5 && event[ihit].charge < 2.5))
            {
                goodcharge.push_back(event[ihit].charge);
                goodplane.push_back(event[ihit].plane);
                goodbar.push_back(event[ihit].bar);
            }

            LOG(debug) << "Hit in coincidence window: " << event[ihit].plane << " " << event[ihit].bar << " "
                       << event[ihit].time << " " << event[ihit].charge;

            ihit++;
            if (ihit >= event.size())
                break;
        }
        if (goodcharge.size() > 3)
        {
            if (goodcharge.size() == 4)
            {
                if (std::accumulate(goodcharge.begin(), goodcharge.end(), 0.0) > 14. &&
                    std::accumulate(goodcharge.begin(), goodcharge.end(), 0.0) < 18.)
                {
                    if (std::accumulate(goodplane.begin(), goodplane.end(), 0) == 6.)
                    {
                        LOG(debug) << "Found good pair 2 times in all planes";
                        for (Int_t g = 0; g < goodcharge.size(); g++)
                        {
                            LOG(debug) << goodcharge.at(g);
                            LOG(debug) << goodplane.at(g);
                            LOG(debug) << goodbar.at(g);
                            goodevents.push_back({ goodcharge.at(g), goodplane.at(g), goodbar.at(g) });
                        }
                        goodpair++;
                        goodpair4++;
                    }
                    else
                    {
                        goodpair3++;
                    }
                }
            }
            else if (std::accumulate(goodcharge.begin(), goodcharge.end(), 0.0) > 14.)
                goodpair6++;
        }
        if (goodcharge.size() == 3)
        {
            if (std::accumulate(goodcharge.begin(), goodcharge.end(), 0.0) > 8.5 &&
                std::accumulate(goodcharge.begin(), goodcharge.end(), 0.0) < 15.5)
            {
                if (std::accumulate(goodplane.begin(), goodplane.end(), 0) == 4 ||
                    std::accumulate(goodplane.begin(), goodplane.end(), 0) == 5)
                {
                    LOG(debug) << "Found good pair at least once in all planes";
                    for (Int_t g = 0; g < goodcharge.size(); g++)
                    {
                        LOG(debug) << goodcharge.at(g);
                    }
                    goodpair++;
                    goodpair5++;
                }
                else
                {
                    goodpair7++;
                }
            }
        }
        if (goodcharge.size() == 2)
        {
            if (std::accumulate(goodcharge.begin(), goodcharge.end(), 0.0) > 7. &&
                std::accumulate(goodcharge.begin(), goodcharge.end(), 0.0) < 9.)
            {
                if (std::accumulate(goodplane.begin(), goodplane.end(), 0) == 2. ||
                    std::accumulate(goodplane.begin(), goodplane.end(), 0) == 4.)
                {
                    LOG(debug) << "Found good pair in one plane";
                    for (Int_t g = 0; g < goodcharge.size(); g++)
                    {
                        LOG(debug) << goodcharge.at(g);
                    }
                    goodpair++;
                    goodpair1++;
                }
                if (std::accumulate(goodplane.begin(), goodplane.end(), 0) == 3.)
                {
                    LOG(debug) << "Found good pair in different planes";
                    for (Int_t g = 0; g < goodcharge.size(); g++)
                    {
                        LOG(debug) << goodcharge.at(g);
                    }
                    goodpair++;
                    goodpair2++;
                }
            }
        }
        std::sort(
            goodevents.begin(), goodevents.end(), [](goodhit const& a, goodhit const& b) { return a.goodq < b.goodq; });
        for (Int_t g = 0; g < goodevents.size(); g++)
        {
            LOG(debug) << goodevents[g].goodq;
            LOG(debug) << goodevents[g].goodp;
            LOG(debug) << goodevents[g].goodb;
        }
    }

    if (fTofdHisto)
    {
        size_t ihit = 0;
        for (; ihit < event.size();)
        {                                                                          // loop over all hits
            fhQ[event[ihit].plane - 1]->Fill(event[ihit].bar, event[ihit].charge); // charge per plane
            fhQvsEvent[event[ihit].plane - 1]->Fill(fnEvents, event[ihit].charge); // charge vs event #
            fhxy[event[ihit].plane - 1]->Fill(event[ihit].bar, event[ihit].ypos);  // xy of plane
            ihit++;
        }
    }

    // store events
    for (Int_t hit = 0; hit < event.size(); hit++)
    { // loop over hits
        if (tArrU[hit] == false)
        {
            eventstore++;
            tArrU[hit] = true;
            // store single hits
            singlehit++;
            new ((*fHitItems)[fHitItems->GetEntriesFast()]) R3BTofdHitData(event[hit].time,
                                                                           event[hit].xpos,
                                                                           event[hit].ypos,
                                                                           event[hit].charge,
                                                                           -1.,
                                                                           event[hit].charge,
                                                                           event[hit].plane,
                                                                           event[hit].bar,
                                                                           event[hit].time_raw,
                                                                           event[hit].tof);
        }
    }

    LOG(debug) << "Used up hits in this event:";
    for (Int_t a = 0; a < event.size(); a++)
    {
        LOG(debug) << "Event " << a << " " << tArrU[a] << " ";
        if (tArrU[a] != true)
            LOG(fatal) << "Not all events analyzed!";
    }

    LOG(debug) << "------------------------------------------------------\n";
    fnEvents++;
}

void R3BTofDCal2Hit::CreateHistograms(Int_t iPlane, Int_t iBar)
{
    Double_t max_charge = 80.;
    // create histograms if not already existing
    if (NULL == fhNoTpat)
    {
        fhNoTpat = new TH1F("NoTpat", "NoTpat", 200, 0, 200);
        fhNoTpat->GetXaxis()->SetTitle("No Tpat event dist");
    }

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

    if (NULL == fhTvsQ[iPlane - 1])
    {
        char strName[255];
        char strName2[255];
        sprintf(strName, "TvsQ_Plane_%d", iPlane);
        sprintf(strName2, "Time Hit vs Q per event Plane %d", iPlane);
        fhTvsQ[iPlane - 1] = new TH2F(strName, "", 10000, 0, 1000, 1000, 0., max_charge);
        fhTvsQ[iPlane - 1]->GetXaxis()->SetTitle("Time between hits in ns");
        fhTvsQ[iPlane - 1]->GetYaxis()->SetTitle("Charge");
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
    if (NULL == fhQ0Qt[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "QvsQt0_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Q vs Q_time0 Plane %d ", iPlane);
        fhQ0Qt[iPlane - 1] =
            new TH2F(strName1, strName2, max_charge * 10, 0., max_charge, max_charge * 10, 0., max_charge);
        fhQ0Qt[iPlane - 1]->GetYaxis()->SetTitle("Charge particle i");
        fhQ0Qt[iPlane - 1]->GetXaxis()->SetTitle("Charge first particle");
    }
}
void R3BTofDCal2Hit::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
    }
}

void R3BTofDCal2Hit::FinishTask()
{
    if (fTofdHisto)
    {
        //   if (fhTpat)
        //     fhTpat->Write();
        if (fhNoTpat)
            fhNoTpat->Write();
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
            if (fhQ0Qt[i])
                fhQ0Qt[i]->Write();
            if (fhTvsQ[i])
                fhTvsQ[i]->Write();
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

    std::stringstream sprint;
    sprint << "\n\nSome statistics:\n";
    sprint << "Total number of events in tree    " << maxevent << "\n";
    sprint << "Max Event analyzed                " << fnEvents + wrongtrigger + wrongtpat + notpat << "\n";
    sprint << "Events skipped due to trigger     " << wrongtrigger << "\n";
    sprint << "Events skipped due to tpat        " << wrongtpat << "\n";
    sprint << "Events with no tpat               " << notpat << "\n";
    sprint << "Events with correct header&tpat   " << headertpat << "\n";
    sprint << "Events without ToFd hits          " << events_wo_tofd_hits << "\n";
    sprint << "Events in cal level               " << events_in_cal_level << " leading and trailing edges\n";
    sprint << "Hits in bar coincidence           " << inbarcoincidence;
    sprint << " leading and trailing edges glued together\n";
    sprint << "Could not match                   " << events_in_cal_level - inbarcoincidence << " edges\n";
    sprint << "Bars with multihits               " << bars_with_multihit << "\n";
    sprint << "# Multihits in planes and bars    " << multihit << "\n";
    sprint << "Events stored                     " << eventstore << " <-> " << inbarcoincidence;
    sprint << " (Events in bar coincidence)\n";
    sprint << "Events in single planes           " << singlehit;
    sprint << "\n";
    //<< "Good events in total            " << eventstore << " <-> " << singlehit << " = singlehit \n";
    sprint << "Really good events                " << goodpair4 << " 2 particles 2 times in 2 planes \n";
    sprint << "Good events                       " << goodpair1 << " 2 particles in 1 plane \n";
    sprint << "Good events                       " << goodpair2 << " 2 particles in different planes \n";
    sprint << "Good events                       " << goodpair5 << " 2 particles at least once in 2 planes\n";
    sprint << "Good events                       " << goodpair7 << " 3 in 2 planes\n";
    sprint << "Good events                       " << goodpair3 << " 2 particles in 2 planes odd\n";
    sprint << "Good events                       " << goodpair6 << " more than 2 particles in 2 planes\n";
    sprint << "Good events in total              " << goodpair << "\n";

    sprint << "n1=" << n1 << " n2=" << n2;

    R3BLOG(info, sprint.str());
}

Double_t R3BTofDCal2Hit::walk(Double_t Q,
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

ClassImp(R3BTofDCal2Hit)
