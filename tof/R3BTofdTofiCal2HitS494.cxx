/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum f�r Schwerionenforschung GmbH    *
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
// ------------------------------------------------------------
// -----                 R3BTofdTofiCal2HitS494                   -----
// -----            Created May 2016 by M.Heil            -----
// -----           Modified Mar 2020 by L.Bott            -----
// ------------------------------------------------------------

#include "R3BTofdTofiCal2HitS494.h"

#include "R3BEventHeader.h"
#include "R3BTofdCalData.h"

#include "R3BTofdHitData.h"
#include "R3BTofdHitModulePar.h"
#include "R3BTofdHitPar.h"

#include "R3BTofiCalData.h"
#include "R3BTofiHitData.h"
#include "R3BTofiHitModulePar.h"
#include "R3BTofiHitPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRuntimeDb.h"
#include "R3BLogger.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"

#include "TClonesArray.h"
#include "TMath.h"

#include <iostream>
#include <numeric>
#include <sstream>
#include <string>

// This header is included in several places, and would re-define the mapping
// table every time -> linker nor happy.
// So do some forward declaration for now, and keep the include in just one
// place, for now R3BTofdCal2Histo.cxx.
//#include "mapping_tofd_trig.hh"
extern unsigned g_tofd_trig_map[4][2][44];
extern void tofd_trig_map_setup();
extern unsigned g_tofi_trig_map[2][24];
void tofi_trig_map_setup();

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
} // namespace

R3BTofdTofiCal2HitS494::R3BTofdTofiCal2HitS494()
    : FairTask("TofdTofiCal2Hit", 1)
    , fHitItems(new TClonesArray("R3BTofdHitData"))
    , fNofHitItems(0)
    , fNofHitPars(0)
    , fHitPar(NULL)
    , fNofHitParsTofi(0)
    , fHitParTofi(NULL)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fNofPlanes(2)
    , fPaddlesPerPlane(44)
    , fNofPlanesTofi(2)
    , fPaddlesPerPlaneTofi(24)
    , fTofdQ(1)
    , fTofiQ(1)
    , fTofiGap(0)
    , fTofdHisto(true)
    , fTofdTotPos(true)
    , fTofiTotPos(true)
    , fShowProgress(true)
    , fYOffset(true)
    , fOffsetFile("")
    , fSimu(false)
{
    if (fTofdHisto)
    {
        // fhTpat = NULL;
        fhNoTpat = NULL;
        for (Int_t i = 0; i < fNofPlanes + fNofPlanesTofi; i++)
        {
            fhQ[i] = NULL;
            fhxy[i] = NULL;
            fhQvsEvent[i] = NULL;
            fhTdiff[i] = NULL;
            fhTsync[i] = NULL;
            fhQ0Qt[i] = NULL;
            fhTvsQ[i] = NULL;
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                fhQvsPos[i][j] = NULL;
            }
        }
    }
}

R3BTofdTofiCal2HitS494::R3BTofdTofiCal2HitS494(const char* name, Int_t iVerbose)
    : FairTask("TofdTofiCal2Hit", 1)
    , fHitItems(new TClonesArray("R3BTofdHitData"))
    , fNofHitItems(0)
    , fNofHitPars(0)
    , fHitPar(NULL)
    , fNofHitParsTofi(0)
    , fHitParTofi(NULL)
    , fTrigger(-1)
    , fTpat1(-1)
    , fTpat2(-1)
    , fNofPlanes(2)
    , fPaddlesPerPlane(44)
    , fNofPlanesTofi(2)
    , fPaddlesPerPlaneTofi(24)
    , fTofdQ(1)
    , fTofiQ(1)
    , fTofiGap(0)
    , fTofdHisto(true)
    , fTofdTotPos(true)
    , fTofiTotPos(true)
    , fShowProgress(true)
    , fYOffset(true)
    , fOffsetFile("")
    , fSimu(false)
{
    if (fTofdHisto)
    {
        //  fhTpat = NULL;
        fhNoTpat = NULL;
        for (Int_t i = 0; i < fNofPlanes + fNofPlanesTofi; i++)
        {
            fhQ[i] = NULL;
            fhxy[i] = NULL;
            fhQvsEvent[i] = NULL;
            fhTdiff[i] = NULL;
            fhTsync[i] = NULL;
            fhQ0Qt[i] = NULL;
            fhTvsQ[i] = NULL;
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                fhQvsPos[i][j] = NULL;
            }
        }
    }
}

R3BTofdTofiCal2HitS494::~R3BTofdTofiCal2HitS494()
{
    if (fTofdHisto)
    {
        // if (fhTpat)
        //   delete fhTpat;
        if (fhNoTpat)
            delete fhNoTpat;
        for (Int_t i = 0; i < fNofPlanes + fNofPlanesTofi; i++)
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
                /*
                if (fhQvsTHit[i][j])
                    delete fhQvsTHit[i][j];
                */
                /*
                if (fhTvsTHit[i][j])
                    delete fhTvsTHit[i][j];
                */
            }
        }
    }
    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

InitStatus R3BTofdTofiCal2HitS494::Init()
{
    fHitPar = (R3BTofdHitPar*)FairRuntimeDb::instance()->getContainer("TofdHitPar");
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

    // if(!fSimu)
    {
        fHitParTofi = (R3BTofiHitPar*)FairRuntimeDb::instance()->getContainer("TofiHitPar");
        if (!fHitParTofi)
        {
            LOG(error) << "Could not get access to TofiHitPar-Container.";
            fNofHitParsTofi = 0;
            return kFATAL;
        }
        fNofHitParsTofi = fHitParTofi->GetNumModulePar();
        if (fNofHitParsTofi == 0)
        {
            LOG(error) << "There are no Hit parameters in container TofiHitPar";
            return kFATAL;
        }
    }

    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";
    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (header)
        R3BLOG(info, "EventHeader. was found");
    else
        R3BLOG(info, "EventHeader. was not found");

    for (int det = 0; det < DET_MAX; det++)
    {
        cout << "Reading detector " << det << ", " << fDetectorNames[det] << endl;
        if (det == 0)
            maxevent = mgr->CheckMaxEventNo();

        fCalItems.push_back((TClonesArray*)mgr->GetObject(Form("%sCal", fDetectorNames[det])));
        if (NULL == fCalItems.at(det))
        {
            printf("Could not find Cal data for '%s'.\n", fDetectorNames[det]);
        }
    }

    // request storage of Hit data in output tree
    mgr->Register("TofdHit", "Land", fHitItems, kTRUE);

    // Getting y offsets *****
    std::string foff_param_file = fOffsetFile;
    ifstream infile(foff_param_file.c_str());
    Int_t iplane, ibar;
    if (!(infile.is_open()) && fYOffset)
        cout << "** ERROR ** Y-OFFSET FILE NOT FOUND! Yoffsets SET to ZERO " << foff_param_file.c_str() << endl;
    if (infile.is_open())
    {
        LOG(info) << "TOFD/TOFI y-OFFSET PARAMS READ FROM FILE " << foff_param_file.c_str() << endl;
        for (Int_t ivec = 0; ivec < 132; ivec++)
        {
            infile >> iplane >> ibar;
            infile >> ytofd_offsetZ2[iplane - 1][ibar - 1] >> ytofd_offsetZ3[iplane - 1][ibar - 1] >>
                ytofd_offsetZ4[iplane - 1][ibar - 1] >> ytofd_offsetZ5[iplane - 1][ibar - 1] >>
                ytofd_offsetZ6[iplane - 1][ibar - 1] >> ytofd_offsetZ7[iplane - 1][ibar - 1] >>
                ytofd_offsetZ8[iplane - 1][ibar - 1];

            /*  cout<<iplane <<", "<< ibar <<", "<<
                    ytofd_offsetZ2[iplane-1][ibar-1] <<", "<<
                    ytofd_offsetZ6[iplane-1][ibar-1] <<", "<<
                    ytofd_offsetZ8[iplane-1][ibar-1]<<endl;*/
        }
    }
    // *******
    // fill control histograms
    if (fTofdHisto)
    {
        for (Int_t i = 1; i < fNofPlanes + fNofPlanesTofi + 1; i++)
        {
            for (Int_t k = 1; k < N_TOFD_HIT_PADDLE_MAX + 1; k++)
            {
                CreateHistograms(i, k);
            }
        }
    }

    tofd_trig_map_setup();
    tofi_trig_map_setup();

    cout << "" << endl;
    cout << fNofPlanes << " planes and " << N_TOFD_HIT_PADDLE_MAX << " bars from ToFD will be treated" << endl;
    cout << fNofPlanesTofi << " planes and " << N_TOFI_HIT_PADDLE_MAX << " bars from ToFI will be treated" << endl;

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BTofdTofiCal2HitS494::SetParContainers()
{
    fHitPar = (R3BTofdHitPar*)FairRuntimeDb::instance()->getContainer("TofdHitPar");
    if (!fHitPar)
    {
        LOG(error) << "Could not get access to TofdHitPar-Container.";
        fNofHitPars = 0;
        return;
    }
    // if(!fSimu)
    {
        fHitParTofi = (R3BTofiHitPar*)FairRuntimeDb::instance()->getContainer("TofiHitPar");
        if (!fHitParTofi)
        {
            LOG(error) << "Could not get access to TofiHitPar-Container.";
            fNofHitParsTofi = 0;
            return;
        }
    }
}

InitStatus R3BTofdTofiCal2HitS494::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

namespace
{
    uint64_t n1, n2;
};

void R3BTofdTofiCal2HitS494::Exec(Option_t* option)
{
    static uint32_t counter = 0;
    if (0 == counter % 10000)
        std::cout << "\rEvents: " << counter << " / " << maxevent << " (" << (int)(counter * 100. / maxevent) << " %)"
                  << std::flush;
    ++counter;
    // cout<<"trigger soll: "<<fTrigger <<", triggert ist: "<<header->GetTrigger()<<endl;
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
    if (header)
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
            if (NULL == fhNoTpat)
            {
                fhNoTpat = new TH1F("NoTpat", "NoTpat", 200, 0, 200);
                fhNoTpat->GetXaxis()->SetTitle("No Tpat event dist");
            }
            fhNoTpat->Fill(counter - lasttpatevent);
            lasttpatevent = counter;
            notpat++;
            // LOG(fatal)<<"No Tpat info";
            return;
        }
    }
    /*
    Int_t fTpat_bit = fTpat - 1;
    if (fTpat_bit >= 0)
    {
        Int_t itpat = header->GetTpat();

        Int_t tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        Int_t t1 = 1 << fTpat_bit;
        Int_t t2 = itpat & (1 << fTpat_bit);
        Int_t t3 = tpatvalue;
        std::cout<<itpat<<" "<<t1<<" "<<t2<<" "<<t3<<"\n";
        if ((header) && (tpatvalue == 0))
        {
            wrongtpat++;
            return;
        }
    }
    */
    headertpat++;
    struct finhit
    {
        Double_t charge;
        Double_t time;
        Double_t xpos;
        Double_t ypos;
        Int_t plane;
        Int_t bar;
        Double_t time_raw;
    };

    //    std::cout<<"new event!*************************************\n";

    std::vector<finhit> finevent;
    // beginning of Tofd
    if (1 == 1)
    {
        Double_t timeP0 = 0.;
        Double_t randx;

        UInt_t vmultihits[fNofPlanes + 1][N_TOFD_HIT_PADDLE_MAX + 1];
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
        };

        //    std::cout<<"new event!*************************************\n";

        std::vector<hit> event;

        auto detTofd = fCalItems.at(DET_TOFD);
        Int_t nHits = detTofd->GetEntriesFast();

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
            auto* hit = (R3BTofdCalData*)detTofd->At(ihit);
            size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();

            auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
            auto& vec = 1 == hit->GetSideId() ? ret.first->second.bot : ret.first->second.top;
            vec.push_back(hit);
            tofd_events_in_cal_level++;
        }

        // Build trigger map.
        std::vector<R3BTofdCalData const*> trig_map;
        auto detTofdTrigg = fCalItems.at(DET_TOFD_TRIGG);
        for (int i = 0; i < detTofdTrigg->GetEntries(); ++i)
        {
            auto trig = (R3BTofdCalData const*)detTofdTrigg->At(i);
            if (trig_map.size() < trig->GetBarId())
            {
                trig_map.resize(trig->GetBarId());
            }
            trig_map.at(trig->GetBarId() - 1) = trig;
        }

        bool s_was_trig_missing = false;
        // Find coincident PMT hits.
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
                        LOG(error) << "R3BTofdCal2HitS494Par::Exec() : Missing trigger information!";
                        LOG(error) << "Top: " << top->GetDetectorId() << ' ' << top->GetSideId() << ' '
                                   << top->GetBarId();
                        LOG(error) << "Bot: " << bot->GetDetectorId() << ' ' << bot->GetSideId() << ' '
                                   << bot->GetBarId();
                        s_was_trig_missing = true;
                    }
                    ++n2;
                }

                // Shift the cyclic difference window by half a window-length and move it back,
                // this way the trigger time will be at 0.
                auto top_ns = fmod(top->GetTimeLeading_ns() - top_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                              c_range_ns / 2;
                auto bot_ns = fmod(bot->GetTimeLeading_ns() - bot_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                              c_range_ns / 2;

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
                if (std::abs(dt_mod) < c_bar_coincidence_ns)
                {
                    tofd_inbarcoincidence++;
                    // Hit!
                    // std::cout << "Hit!\n";
                    Int_t iPlane = top->GetDetectorId(); // 1..n
                    Int_t iBar = top->GetBarId();        // 1..n
                    if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                    {
                        //   LOG(error) << "R3BTofdCal2HitS494Par::Exec() : more planes than expected! Det: " << iPlane
                        //            << " allowed are 1.." << fNofPlanes;
                        continue;
                    }
                    if (iBar > fPaddlesPerPlane) // same here
                    {
                        LOG(error) << "R3BTofdCal2HitS494Par::Exec() : more bars then expected! Det: " << iBar
                                   << " allowed are 1.." << fPaddlesPerPlane;
                        continue;
                    }

                    auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                    auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                    auto THit_raw =
                        (bot->GetTimeLeading_ns() + top->GetTimeLeading_ns()) / 2.; // needed for TOF for ROLUs

                    // std::cout<<"ToT: "<<top_tot << " "<<bot_tot<<"\n";

                    // register multi hits
                    vmultihits[iPlane][iBar] += 1;

                    R3BTofdHitModulePar* par = fHitPar->GetModuleParAt(iPlane, iBar);
                    if (!par)
                    {
                        LOG(info) << "R3BTofdCal2HitS494::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                                  << ", Bar: " << top->GetBarId();
                        continue;
                    }

                    // walk corrections
                    if (par->GetPar1Walk() == 0. || par->GetPar2Walk() == 0. || par->GetPar3Walk() == 0. ||
                        par->GetPar4Walk() == 0. || par->GetPar5Walk() == 0.)
                        LOG(debug) << "Walk correction not found!";
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

                    // calculate y-position
                    auto pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();

                    //  cout << "Test1: " << pos << " bot: " << bot_ns << " top: " << top_ns << " offsets: " <<
                    //  par->GetOffset1() << "  " << par->GetOffset2() << " veff: " << par->GetVeff() <<  endl;

                    // calculate y-position from ToT
                    auto posToT =
                        par->GetLambda() * log((top_tot * par->GetToTOffset2()) / (bot_tot * par->GetToTOffset1()));

                    //    cout<<"TOFD POS: "<<posToT<<"; "<<par->GetLambda() <<", "<<par->GetToTOffset2()<<";
                    //    "<<par->GetToTOffset1()<<", "<<top_tot<<", "<<bot_tot<<endl;

                    //    cout << "Test2: " << posToT << " bot: " << bot_tot << " top: " << top_tot << " offsets: " <<
                    //   par->GetToTOffset1() << "  " << par->GetToTOffset2() << " lambda: " << par->GetLambda() <<
                    //   endl;

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
                    Double_t xp = -1000.;
                    // calculate x-position
                    if (iPlane == 1 || iPlane == 3)
                    {
                        xp = -detector_width / 2 + (paddle_width + air_gap_paddles) / 2 +
                             (iBar - 1) * (paddle_width + air_gap_paddles) - 0.04;
                    }
                    if (iPlane == 2 || iPlane == 4)
                    {
                        xp = -detector_width / 2 + (paddle_width + air_gap_paddles) +
                             (iBar - 1) * (paddle_width + air_gap_paddles) - 0.04;
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

                            //   cout<<"TOFD Q: "<<para[0]<<"; "<<para[1]<<"; "<<para[2]<<", "<<para[3]<<", "<<qb<<endl;
                        }
                        else
                        {
                            // via double exponential:
                            para[0] = par->GetPar1a();
                            para[1] = par->GetPar1b();
                            para[2] = par->GetPar1c();
                            para[3] = par->GetPar1d();
                            auto q1 =
                                bot_tot /
                                (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                            para[0] = par->GetPar2a();
                            para[1] = par->GetPar2b();
                            para[2] = par->GetPar2c();
                            para[3] = par->GetPar2d();
                            auto q2 =
                                top_tot /
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

                    /*   if(fSimu)
                       {
                           pos = (bot_ns - top_ns) * par->GetVeff();
                           if(par->GetVeff() == 0.) cout<< "tofd veff != 6 "<<par->GetVeff()<<"; "<<iPlane<<";
                       "<<iBar<<endl;
                       }*/

                    Double_t parz[3];
                    parz[0] = par->GetPar1za();
                    parz[1] = par->GetPar1zb();
                    parz[2] = par->GetPar1zc();

                    // cout<<"parz: "<<parz[0]<<", "<<parz[1]<<", "<<parz[2]<<endl;

                    if (parz[0] > 0 && parz[2] > 0)
                        LOG(debug) << "Charges in this event " << parz[0] * TMath::Power(qb, parz[2]) + parz[1]
                                   << " plane " << iPlane << " ibar " << iBar;
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

                    if (parz[0] > 0 && parz[2] > 0)
                    {
                        event.push_back(
                            { parz[0] * TMath::Power(qb, parz[2]) + parz[1], THit, xp, pos, iPlane, iBar, THit_raw });
                    }
                    else
                    {
                        parz[0] = 1.;
                        parz[1] = 0.;
                        parz[2] = 1.;

                        event.push_back({ qb, THit, xp, pos, iPlane, iBar, THit_raw });
                    }

                    if (fTofdHisto)
                    {
                        // fill control histograms
                        fhTsync[iPlane - 1]->Fill(iBar, THit);
                        fhTdiff[iPlane - 1]->Fill(iBar, tdiff);
                        // fhQvsPos[iPlane - 1][iBar - 1]->Fill(pos, parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
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
                    tofd_bars_with_multihit++;
                    tofd_multihit += vmultihits[i][j] - 1;
                    // if(vmultihits[i][j]>3)cout<<vmultihits[i][j] - 1<<"\n";
                }
            }
        }

        std::sort(event.begin(), event.end(), [](hit const& a, hit const& b) { return a.time < b.time; });
        // Now we have all hits in this event time sorted

        LOG(debug) << "Tofd: Charge Time xpos ypos plane bar";
        for (Int_t hit = 0; hit < event.size(); hit++)
        {
            LOG(debug) << event[hit].charge << " " << event[hit].time << " " << event[hit].xpos << " "
                       << event[hit].ypos << " " << event[hit].plane << " " << event[hit].bar;
            // if (event[hit].plane == 2 && (event[hit].bar < 21 || event[hit].bar > 24)) fhTvsQ[event[hit].plane -
            // 1]->Fill(event[hit].time-event[0].time,event[hit].charge);
            if (event[hit].plane == 2 && (event[hit].bar == 18))
                fhTvsQ[event[hit].plane - 1]->Fill(event[hit].time - event[0].time, event[hit].charge);
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
                    if (accumulate(goodcharge.begin(), goodcharge.end(), 0.0) > 14. &&
                        accumulate(goodcharge.begin(), goodcharge.end(), 0.0) < 18.)
                    {
                        if (accumulate(goodplane.begin(), goodplane.end(), 0) == 6.)
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
                else if (accumulate(goodcharge.begin(), goodcharge.end(), 0.0) > 14.)
                    goodpair6++;
            }
            if (goodcharge.size() == 3)
            {
                if (accumulate(goodcharge.begin(), goodcharge.end(), 0.0) > 8.5 &&
                    accumulate(goodcharge.begin(), goodcharge.end(), 0.0) < 15.5)
                {
                    if (accumulate(goodplane.begin(), goodplane.end(), 0) == 4 ||
                        accumulate(goodplane.begin(), goodplane.end(), 0) == 5)
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
                if (accumulate(goodcharge.begin(), goodcharge.end(), 0.0) > 7. &&
                    accumulate(goodcharge.begin(), goodcharge.end(), 0.0) < 9.)
                {
                    if (accumulate(goodplane.begin(), goodplane.end(), 0) == 2. ||
                        accumulate(goodplane.begin(), goodplane.end(), 0) == 4.)
                    {
                        LOG(debug) << "Found good pair in one plane";
                        for (Int_t g = 0; g < goodcharge.size(); g++)
                        {
                            LOG(debug) << goodcharge.at(g);
                        }
                        goodpair++;
                        goodpair1++;
                    }
                    if (accumulate(goodplane.begin(), goodplane.end(), 0) == 3.)
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
            std::sort(goodevents.begin(), goodevents.end(), [](goodhit const& a, goodhit const& b) {
                return a.goodq < b.goodq;
            });
            for (Int_t g = 0; g < goodevents.size(); g++)
            {
                LOG(debug) << goodevents[g].goodq;
                LOG(debug) << goodevents[g].goodp;
                LOG(debug) << goodevents[g].goodb;
            }
        }

        // store events
        for (Int_t hit = 0; hit < event.size(); hit++)
        { // loop over hits
            if (tArrU[hit] == false)
            {
                tofd_eventstore++;
                tArrU[hit] = true;
                // store single hits
                tofd_singlehit++;
                if (fTofdTotPos && !fSimu)
                {
                    Int_t iplane = event[hit].plane;
                    Int_t ibar = event[hit].bar;

                    event[hit].ypos = event[hit].ypos * sqrt(event[hit].charge) / sqrt(8.);

                    if (fYOffset)
                    {
                        if (event[hit].charge > 1.4 && event[hit].charge < 2.6)
                        {
                            // cout<<"Zfrag=2 orig: "<<iplane<<", "<<ibar<<", "<<event[hit].ypos<<",
                            // "<<ytofd_offsetZ2[iplane-1][ibar-1]<<endl;
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ2[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=2 corr: "<<iplane<<", "<<ibar<<", "<<event[hit].ypos<<",
                            // "<<ytofd_offsetZ2[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 2.6 && event[hit].charge < 3.5)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ3[iplane - 1][ibar - 1];
                            //  cout<<"Zfrag=3: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ3[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 3.5 && event[hit].charge < 4.5)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ4[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=4: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ4[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 4.5 && event[hit].charge < 5.5)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ5[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=5: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ5[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 5.5 && event[hit].charge < 6.6)
                        {
                            //   cout<<"Zfrag=6 orig: "<<iplane<<", "<<ibar<<", "<<event[hit].ypos<<",
                            //   "<<ytofd_offsetZ6[iplane-1][ibar-1]<<endl;
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ6[iplane - 1][ibar - 1];
                            //  cout<<"Zfrag=6 corr: "<<iplane<<", "<<ibar<<", "<<event[hit].ypos<<",
                            //  "<<ytofd_offsetZ6[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 6.6 && event[hit].charge <= 7.2)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ7[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=7: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ7[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge > 7.2)
                        { //&& event[hit].charge < 8.8) {
                          //  cout<<"Zfrag=8 orig: "<<iplane<<", "<<ibar<<", "<<event[hit].ypos<<",
                          //  "<<ytofd_offsetZ8[iplane-1][ibar-1]<<endl;
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ8[iplane - 1][ibar - 1];
                            //   cout<<"Zfrag=8 corr: "<<iplane<<", "<<ibar<<", "<<event[hit].ypos<<",
                            //   "<<ytofd_offsetZ8[iplane-1][ibar-1]<<endl;
                        }
                    }
                }
                if (fTofdHisto)
                {
                    fhQ[event[hit].plane - 1]->Fill(event[hit].bar, event[hit].charge);  // charge per plane
                    fhQvsEvent[event[hit].plane - 1]->Fill(fnEvents, event[hit].charge); // charge vs event #
                    fhxy[event[hit].plane - 1]->Fill(event[hit].bar, event[hit].ypos);   // xy of plane
                    fhQvsPos[event[hit].plane - 1][event[hit].bar - 1]->Fill(event[hit].ypos, event[hit].charge);
                }
                /*AKH     new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(event[hit].time,
                                                                       event[hit].xpos,
                                                                       event[hit].ypos,
                                                                       event[hit].charge,
                                                                       -1.,
                                                                       event[hit].charge,
                                                                       event[hit].plane,
                                                                       event[hit].bar,
                                                                       event[hit].time_raw); AKH*/
                finevent.push_back({ event[hit].charge,
                                     event[hit].time,
                                     event[hit].xpos,
                                     event[hit].ypos,
                                     event[hit].plane,
                                     event[hit].bar });
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
    } // End of ToFD

    // ***** TOFI *****************
    if (1 == 1) //&& !fSimu)
    {
        Double_t timeP0 = 0.;
        Double_t randx;

        UInt_t vmultihits[fNofPlanesTofi + 1][N_TOFI_HIT_PADDLE_MAX + 1];
        for (Int_t i = 0; i <= fNofPlanesTofi; i++)
        {
            for (Int_t j = 0; j <= N_TOFI_HIT_PADDLE_MAX; j++)
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
        };

        //    std::cout<<"new event!*************************************\n";
        std::vector<hit> event;

        auto detTofi = fCalItems.at(DET_TOFI);
        Int_t nHits = detTofi->GetEntries();
        // Organize cals into bars.
        struct Entry
        {
            std::vector<R3BTofiCalData*> top;
            std::vector<R3BTofiCalData*> bot;
        };
        std::map<size_t, Entry> bar_map;
        // puts("Event");
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto* hit = (R3BTofiCalData*)detTofi->At(ihit);
            size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();
            /*
                    std::cout << "Hits: " << hit->GetDetectorId() << ' ' << hit->GetBarId() << ' ' << hit->GetSideId()
               << '
               '
                              << hit->GetTimeLeading_ns() << ' ' << hit->GetTimeTrailing_ns()
                              << ' ' << hit->GetTimeTrailing_ns() - hit->GetTimeLeading_ns() << '\n';
            */
            auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
            // auto& vec = 1 == hit->GetSideId() ? ret.first->second.bot : ret.first->second.top;
            auto& vec = 1 == hit->GetSideId() ? ret.first->second.top : ret.first->second.bot;
            vec.push_back(hit);
            tofi_events_in_cal_level++;
        }

        // Build trigger map.
        std::vector<R3BTofiCalData const*> trig_map_tofi;
        auto detTofiTrigger = fCalItems.at(DET_TOFI_TRIGG);
        for (int i = 0; i < detTofiTrigger->GetEntries(); ++i)
        {
            auto trig_tofi = (R3BTofiCalData const*)detTofiTrigger->At(i);
            if (trig_map_tofi.size() < trig_tofi->GetBarId())
            {
                trig_map_tofi.resize(trig_tofi->GetBarId());
            }
            trig_map_tofi.at(trig_tofi->GetBarId() - 1) = trig_tofi;
        }

        bool s_was_trig_missing = false;
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
                auto top_trig_i = g_tofi_trig_map[top->GetSideId() - 1][top->GetBarId() - 1];
                auto bot_trig_i = g_tofi_trig_map[bot->GetSideId() - 1][bot->GetBarId() - 1];
                Double_t top_trig_ns = 0, bot_trig_ns = 0;
                if (top_trig_i < trig_map_tofi.size() && trig_map_tofi.at(top_trig_i) &&
                    bot_trig_i < trig_map_tofi.size() && trig_map_tofi.at(bot_trig_i))
                {
                    auto top_trig = trig_map_tofi.at(top_trig_i);
                    auto bot_trig = trig_map_tofi.at(bot_trig_i);
                    top_trig_ns = top_trig->GetTimeLeading_ns();
                    bot_trig_ns = bot_trig->GetTimeLeading_ns();

                    ++n1;
                }
                else
                {
                    if (!s_was_trig_missing)
                    {
                        LOG(error) << "R3BTofiCal2HitS494Par::Exec() : Missing trigger information!";
                        s_was_trig_missing = true;
                    }
                    ++n2;
                }

                // Shift the cyclic difference window by half a window-length and move it back,
                // this way the trigger time will be at 0.
                auto top_ns = fmod(top->GetTimeLeading_ns() - top_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                              c_range_ns / 2;
                auto bot_ns = fmod(bot->GetTimeLeading_ns() - bot_trig_ns + c_range_ns + c_range_ns / 2, c_range_ns) -
                              c_range_ns / 2;

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

                if (std::abs(dt_mod) < c_bar_coincidence_ns)
                {
                    tofi_inbarcoincidence++;
                    // Hit!
                    // std::cout << "Hit!\n";
                    Int_t iPlane = top->GetDetectorId(); // 1..n
                    Int_t iBar = top->GetBarId();        // 1..n
                    if (iPlane > fNofPlanesTofi)         // this also errors for iDetector==0
                    {

                        ++top_i;
                        ++bot_i;
                        continue;
                    }
                    if (iBar > fPaddlesPerPlaneTofi) // same here
                    {

                        ++top_i;
                        ++bot_i;
                        continue;
                    }

                    auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                    auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                    // register multi hits
                    vmultihits[iPlane][iBar] += 1;

                    R3BTofiHitModulePar* parTofi = fHitParTofi->GetModuleParAt(iPlane, iBar);
                    if (!parTofi)
                    {
                        LOG(info) << "R3BTofiCal2HitS494::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                                  << ", Bar: " << top->GetBarId();
                        ++top_i;
                        ++bot_i;
                        continue;
                    }
                    // walk corrections
                    if (parTofi->GetPar1Walk() == 0. || parTofi->GetPar2Walk() == 0. || parTofi->GetPar3Walk() == 0. ||
                        parTofi->GetPar4Walk() == 0. || parTofi->GetPar5Walk() == 0.)
                        LOG(debug) << "Walk correction not found!";
                    auto bot_ns_walk = bot_ns - walk(bot_tot,
                                                     parTofi->GetPar1Walk(),
                                                     parTofi->GetPar2Walk(),
                                                     parTofi->GetPar3Walk(),
                                                     parTofi->GetPar4Walk(),
                                                     parTofi->GetPar5Walk());
                    auto top_ns_walk = top_ns - walk(top_tot,
                                                     parTofi->GetPar1Walk(),
                                                     parTofi->GetPar2Walk(),
                                                     parTofi->GetPar3Walk(),
                                                     parTofi->GetPar4Walk(),
                                                     parTofi->GetPar5Walk());
                    // calculate tdiff
                    auto tdiff = ((bot_ns + parTofi->GetOffset1()) - (top_ns + parTofi->GetOffset2()));

                    // calculate time of hit
                    Double_t THit = (bot_ns + top_ns) / 2. - parTofi->GetSync();
                    if (std::isnan(THit))
                    {
                        LOG(fatal) << "ToFi THit not found";
                    }
                    if (timeP0 == 0.)
                        timeP0 = THit;

                    // calculate y-position
                    auto pos =
                        ((bot_ns + parTofi->GetOffset1()) - (top_ns + parTofi->GetOffset2())) * parTofi->GetVeff();

                    // calculate y-position from ToT
                    auto posToT = parTofi->GetLambda() *
                                  log((top_tot * parTofi->GetToTOffset2()) / (bot_tot * parTofi->GetToTOffset1()));

                    //  cout<<"TOFI POS: "<<posToT<<", "<<parTofi->GetLambda()<<", "<<parTofi->GetToTOffset2()<<",
                    //  "<<parTofi->GetToTOffset1()<<"; "<<top_tot<<
                    //", "<<bot_tot<<endl;

                    if (fTofiTotPos)
                        pos = posToT;

                    /*	if(fSimu)
                        {
                             pos = (bot_ns - top_ns) * parTofi->GetVeff();
                             if(parTofi->GetVeff() == 0) cout<<"tofi veff=0 "<<parTofi->GetVeff()<<", "<<iPlane<<";
                       "<<iBar<<endl;
                         }*/

                    Float_t paddle_width = 0.50000;
                    Float_t paddle_thickness = 0.50000;
                    Float_t air_gap_paddles = 0.04;
                    Float_t air_gap_layer = 5.;
                    Float_t gap_center_layer = 6.;
                    if (fTofiGap > 0.)
                        gap_center_layer = fTofiGap;
                    // define number of layers and paddles with sizes of the detector
                    Int_t number_layers = fNofPlanesTofi;        // Sabina 2;   // number of layers
                    Int_t number_paddles = fPaddlesPerPlaneTofi; // number of paddles per layer
                    Float_t detector_width =
                        number_paddles * paddle_width + (number_paddles - 1) * air_gap_paddles + gap_center_layer;
                    Float_t detector_thickness = (number_layers - 1) * air_gap_layer + number_layers * paddle_thickness;

                    // calculate x-position
                    randx = (std::rand() / (float)RAND_MAX);
                    Double_t xp = -1000.;
                    if (iPlane == 1 || iPlane == 3)
                    {
                        //	-detector_width / 2 + (paddle_width + air_gap_paddles) / 2 +
                        //                            (iBar - 1) * (paddle_width + air_gap_paddles) - 0.04);

                        if (iBar <= number_paddles / 2)
                            xp = -detector_width / 2 + paddle_width / 2 + (iBar - 1) * (paddle_width + air_gap_paddles);
                        else
                            xp = -detector_width / 2 + paddle_width / 2 - air_gap_paddles +
                                 (iBar - 1) * (paddle_width + air_gap_paddles) + gap_center_layer;
                    }
                    if (iPlane == 2 || iPlane == 4)
                    {
                        if (iBar <= number_paddles / 2)
                            xp = -detector_width / 2 + paddle_width + (iBar - 1) * (paddle_width + air_gap_paddles);
                        else
                            xp = -detector_width / 2 + paddle_width - air_gap_paddles +
                                 (iBar - 1) * (paddle_width + air_gap_paddles) + gap_center_layer;
                    }

                    Double_t para[4];
                    Double_t qb = 0.;
                    if (fTofiTotPos)
                    {
                        // via pol3
                        para[0] = parTofi->GetPola();
                        para[1] = parTofi->GetPolb();
                        para[2] = parTofi->GetPolc();
                        para[3] = parTofi->GetPold();

                        qb = TMath::Sqrt(top_tot * bot_tot) /
                             (para[0] + para[1] * pos + para[2] * pow(pos, 2) + para[3] * pow(pos, 3));
                        qb = qb * fTofiQ;

                        // cout<<"TOFI Q: "<<para[0]<<"; "<<para[1]<<"; "<<para[2]<<", "<<para[3]<<", "<<qb<<endl;
                    }
                    else
                    {
                        // via double exponential:
                        para[0] = parTofi->GetPar1a();
                        para[1] = parTofi->GetPar1b();
                        para[2] = parTofi->GetPar1c();
                        para[3] = parTofi->GetPar1d();
                        auto q1 = bot_tot /
                                  (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                        para[0] = parTofi->GetPar2a();
                        para[1] = parTofi->GetPar2b();
                        para[2] = parTofi->GetPar2c();
                        para[3] = parTofi->GetPar2d();
                        auto q2 = top_tot /
                                  (para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3]);
                        q1 = q1 * fTofiQ;
                        q2 = q2 * fTofiQ;
                        qb = (q1 + q2) / 2.;
                    }

                    Double_t parz[3];
                    parz[0] = parTofi->GetPar1za();
                    parz[1] = parTofi->GetPar1zb();
                    parz[2] = parTofi->GetPar1zc();

                    if (parz[0] > 0 && parz[2] > 0)
                        LOG(debug) << "Charges in this event " << parz[0] * TMath::Power(qb, parz[2]) + parz[1]
                                   << " plane " << iPlane << " ibar " << iBar;
                    else
                        LOG(debug) << "Charges in this event " << qb << " plane " << iPlane << " ibar " << iBar;
                    LOG(debug) << "Times in this event " << THit << " plane " << iPlane << " ibar " << iBar;
                    if (iPlane == 1 || iPlane == 3)
                    {
                        if (iBar <= number_paddles / 2)
                            LOG(debug) << "x in this event "
                                       << -detector_width / 2 + paddle_width / 2 +
                                              (iBar - 1) * (paddle_width + air_gap_paddles)
                                       << " plane " << iPlane << " ibar " << iBar;
                        else
                            LOG(debug) << "x in this event "
                                       << -detector_width / 2 + paddle_width / 2 - air_gap_paddles +
                                              (iBar - 1) * (paddle_width + air_gap_paddles) + gap_center_layer
                                       << " plane " << iPlane << " ibar " << iBar;
                    }
                    if (iPlane == 2 || iPlane == 4)
                    {
                        if (iBar <= number_paddles / 2)
                            LOG(debug) << "x in this event "
                                       << -detector_width / 2 + paddle_width +
                                              (iBar - 1) * (paddle_width + air_gap_paddles)
                                       << " plane " << iPlane << " ibar " << iBar;
                        else
                            LOG(debug) << "x in this event "
                                       << -detector_width / 2 + paddle_width - air_gap_paddles +
                                              (iBar - 1) * (paddle_width + air_gap_paddles) + gap_center_layer
                                       << " plane " << iPlane << " ibar " << iBar;
                    }
                    LOG(debug) << "y in this event " << pos << " plane " << iPlane << " ibar " << iBar << "\n";

                    // cout<<"parz: "<<parz[0]<<", "<<parz[1]<<", "<<parz[2]<<endl;
                    if (parz[0] > 0 && parz[2] > 0)
                    {
                        event.push_back({ parz[0] * TMath::Power(qb, parz[2]) + parz[1], THit, xp, pos, iPlane, iBar });
                    }
                    else
                    {
                        parz[0] = 1.;
                        parz[1] = 0.;
                        parz[2] = 1.;
                        event.push_back({ qb, THit, xp, pos, iPlane, iBar });
                    }

                    if (fTofdHisto)
                    {
                        // fill control histograms
                        // a  fhTsync[iPlane+fNofPlanes - 1]->Fill(iBar, THit);
                        fhTdiff[iPlane + fNofPlanes - 1]->Fill(iBar, tdiff);
                        // a fhQvsPos[iPlane+fNofPlanes - 1][iBar - 1]->Fill(pos, parz[0] * TMath::Power(qb, parz[2]) +
                        // parz[1]); fhQvsTHit[iPlane+fNofPlanes - 1][iBar - 1]->Fill(qb, THit);
                        // fhTvsTHit[iPlane+fNofPlanes - 1][iBar - 1]->Fill(dt_mod, THit);
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
        if (event.empty())
            events_wo_tofi_hits++;

        // init arrays to store hits
        Bool_t tArrU[event.size() + 1];
        for (int i = 0; i < (event.size() + 1); i++)
            tArrU[i] = kFALSE;

        for (Int_t i = 1; i <= fNofPlanesTofi; i++)
        {
            for (Int_t j = 0; j < fPaddlesPerPlaneTofi + 1; j++)
            {
                if (vmultihits[i][j] > 1)
                {
                    tofi_bars_with_multihit++;
                    tofi_multihit += vmultihits[i][j] - 1;
                }
            }
        }

        std::sort(event.begin(), event.end(), [](hit const& a, hit const& b) { return a.time < b.time; });
        // Now we have all hits in this event time sorted

        LOG(debug) << "Tofi: Charge Time xpos ypos plane bar";
        for (Int_t hit = 0; hit < event.size(); hit++)
        {
            LOG(debug) << event[hit].charge << " " << event[hit].time << " " << event[hit].xpos << " "
                       << event[hit].ypos << " " << event[hit].plane << " " << event[hit].bar;
        }

        // Now we can analyze the hits in this event

        for (Int_t hit = 0; hit < event.size(); hit++)
        { // loop over hits
            if (tArrU[hit] == false)
            {
                tArrU[hit] = true;
                // store single hits
                if (fTofiTotPos && !fSimu)
                {
                    Int_t iplane = event[hit].plane + fNofPlanes;
                    Int_t ibar = event[hit].bar;

                    event[hit].ypos = event[hit].ypos * sqrt(event[hit].charge) / sqrt(8.);

                    if (fYOffset)
                    {
                        // cout<<"Tofi y in: "<<event[hit].charge<<", "<<event[hit].ypos<<endl;
                        if (event[hit].charge > 1.4 && event[hit].charge < 2.6)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ2[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=2: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ2[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 2.6 && event[hit].charge < 3.5)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ3[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=3: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ3[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 3.5 && event[hit].charge < 4.5)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ4[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=4: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ4[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 4.5 && event[hit].charge < 5.5)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ5[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=5: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ5[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 5.5 && event[hit].charge < 6.6)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ6[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=6: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ6[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge >= 6.6 && event[hit].charge <= 7.2)
                        {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ7[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=7: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ7[iplane-1][ibar-1]<<endl;
                        }
                        if (event[hit].charge > 7.2)
                        { //&& event[hit].charge < 8.8) {
                            event[hit].ypos = event[hit].ypos - ytofd_offsetZ8[iplane - 1][ibar - 1];
                            // cout<<"Zfrag=8: "<<iplane<<", "<<ibar<<", "<<ytofd_offsetZ8[iplane-1][ibar-1]<<endl;
                        }
                        // cout<<"Tofi y first out: "<<event[hit].charge<<", "<<event[hit].ypos<<endl;
                    }
                }
                if (fTofdHisto)
                {
                    tofi_eventstore++;

                    // if (abs(event[ihit].time) < 50)
                    {
                        fhQ[event[hit].plane + fNofPlanes - 1]->Fill(event[hit].bar,
                                                                     event[hit].charge); // charge per plane
                        fhQvsEvent[event[hit].plane + fNofPlanes - 1]->Fill(fnEvents,
                                                                            event[hit].charge); // charge vs event #
                        fhTsync[event[hit].plane + fNofPlanes - 1]->Fill(event[hit].bar, event[hit].time);
                        fhQvsPos[event[hit].plane + fNofPlanes - 1][event[hit].bar - 1]->Fill(event[hit].ypos,
                                                                                              event[hit].charge);
                        fhxy[event[hit].plane + fNofPlanes - 1]->Fill(event[hit].bar, event[hit].ypos); // xy of plane
                    }
                }

                tofi_singlehit++;

                /*AKH     new ((*fHitItems)[fNofHitItems++])
                         R3BTofdHitData(event[hit].time,
                                        event[hit].xpos,
                                        event[hit].ypos,
                                        event[hit].charge,
                                        -1.,
                                        event[hit].charge,
                                        event[hit].plane + fNofPlanes, // +2 in order to add Tofi as 3rd Tofd plane
                                        event[hit].bar); AKH*/

                finevent.push_back({ event[hit].charge,
                                     event[hit].time,
                                     event[hit].xpos,
                                     event[hit].ypos,
                                     event[hit].plane + 2,
                                     event[hit].bar });
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
    } // end Tofi

    LOG(debug) << "Tofi+Tofd: Charge Time xpos ypos plane bar";
    for (Int_t hit = 0; hit < finevent.size(); hit++)
    {
        LOG(debug) << finevent[hit].charge << " " << finevent[hit].time << " " << finevent[hit].xpos << " "
                   << finevent[hit].ypos << " " << finevent[hit].plane << " " << finevent[hit].bar;

        new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(finevent[hit].time,
                                                          finevent[hit].xpos,
                                                          finevent[hit].ypos,
                                                          finevent[hit].charge,
                                                          -1.,
                                                          finevent[hit].charge,
                                                          finevent[hit].plane,
                                                          finevent[hit].bar);
    }

    fnEvents++;
}

void R3BTofdTofiCal2HitS494::CreateHistograms(Int_t iPlane, Int_t iBar)
{
    Double_t max_charge = 20.;
    // create histograms if not already existing

    /*    if (NULL == fhTpat)
        {
            fhTpat = new TH1F("Tpat", "Tpat", 20, 0, 20);
            fhTpat->GetXaxis()->SetTitle("Tpat value");
        }
    */
    //  cout<<"Creating histogramm: "<<iPlane<<", "<<iBar<<endl;

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
        fhTvsQ[iPlane - 1] = new TH2F(strName, "", 10000, 0, 1000, max_charge * 10, 0., max_charge);
        fhTvsQ[iPlane - 1]->GetXaxis()->SetTitle("Time between hits in ns");
        fhTvsQ[iPlane - 1]->GetYaxis()->SetTitle("Charge");
    }

    if (NULL == fhQvsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 400, -100, 100, max_charge * 10, 0., max_charge);
        fhQvsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("y Position in cm");
    }

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
void R3BTofdTofiCal2HitS494::FinishEvent()
{
    for (Int_t det = 0; det < DET_MAX; det++)
    {
        if (fCalItems.at(det))
        {
            fCalItems.at(det)->Clear();
        }
    }
    if (fHitItems)
    {
        fHitItems->Clear();
        fNofHitItems = 0;
    }
}

void R3BTofdTofiCal2HitS494::FinishTask()
{
    if (fTofdHisto)
    {
        //   if (fhTpat)
        //     fhTpat->Write();
        if (fhNoTpat)
            fhNoTpat->Write();
        for (Int_t i = 0; i < fNofPlanes + fNofPlanesTofi; i++)
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
    std::cout << "\n\nSome statistics:\n"
              << "Total number of events in tree    " << maxevent << "\n"
              << "Max Event analyzed                " << fnEvents + wrongtrigger + wrongtpat + notpat << "\n"
              << "Events skipped due to trigger     " << wrongtrigger << "\n"
              << "Events skipped due to tpat        " << wrongtpat << "\n"
              << "Events with no tpat               " << notpat << "\n"
              << "Events with correct header&tpat   " << headertpat << "\n"
              << "Events without ToFd hits          " << events_wo_tofd_hits << "\n"
              << "tofd Events in cal level          " << tofd_events_in_cal_level << " leading and trailing edges\n"
              << "tofi Events in cal level          " << tofi_events_in_cal_level << " leading and trailing edges\n"
              << "tofd Hits in bar coincidence      " << tofd_inbarcoincidence << " \n"
              << "tofi Hits in bar coincidence      " << tofi_inbarcoincidence << " \n"
              << "tofd leading and trailing edges glued together\n"
              << "Could not match                   " << tofd_events_in_cal_level - tofd_inbarcoincidence << " edges\n"
              << "tofi leading and trailing edges glued together\n"
              << "Could not match                   " << tofi_events_in_cal_level - tofi_inbarcoincidence << " edges\n"
              << "Tofd Bars with multihits          " << tofd_bars_with_multihit << "\n"
              << "Tofi Bars with multihits          " << tofi_bars_with_multihit << "\n"
              << "tofd # Multihits in planes and bars    " << tofd_multihit << "\n"
              << "tofi # Multihits in planes and bars    " << tofi_multihit << "\n"
              << "tofd Events stored                     " << tofd_eventstore << " <-> " << tofd_inbarcoincidence
              << " (Events in bar coincidence)\n"
              << "tofi Events stored                     " << tofi_eventstore << " <-> " << tofi_inbarcoincidence
              << " (Events in bar coincidence)\n"
              << "tofd Events in single planes           " << tofd_singlehit << "\n"
              << "tofi Events in single planes           " << tofi_singlehit
              << "\n"
              //<< "Good events in total            " << eventstore << " <-> " << singlehit << " = singlehit \n";
              << "Really good events                " << goodpair4 << " 2 particles 2 times in 2 planes \n"
              << "Good events                       " << goodpair1 << " 2 particles in 1 plane \n"
              << "Good events                       " << goodpair2 << " 2 particles in different planes \n"
              << "Good events                       " << goodpair5 << " 2 particles at least once in 2 planes\n"
              << "Good events                       " << goodpair7 << " 3 in 2 planes\n"
              << "Good events                       " << goodpair3 << " 2 particles in 2 planes odd\n"
              << "Good events                       " << goodpair6 << " more than 2 particles in 2 planes\n"
              << "Good events in total              " << goodpair << "\n";

    std::cout << "n1=" << n1 << " n2=" << n2 << std::endl;
}

Double_t R3BTofdTofiCal2HitS494::betaCorr(Double_t delta)
{
    //    Double_t corr=-3.*delta;  //correction for Ag

    Double_t corr = -1. * delta; // correction for 12C
    corr = 0.;
    return corr;
}
/* old method
Double_t R3BTofdTofiCal2HitS494::walk(Double_t q)
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
Double_t R3BTofdTofiCal2HitS494::walk(Double_t Q,
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
Double_t R3BTofdTofiCal2HitS494::saturation(Double_t x)
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

ClassImp(R3BTofdTofiCal2HitS494)
