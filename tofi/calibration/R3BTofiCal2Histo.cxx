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
// ----------------------------------------------------------
// ----- Create histograms for parameters for Tofi      -----
// -----     Created July 2019 by L.Bott                -----
// ----------------------------------------------------------
#include "R3BTofiCal2Histo.h"
#include "R3BEventHeader.h"
#include "R3BTofiCalData.h"
#include "R3BTofiHitModulePar.h"
#include "R3BTofiHitPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <stdlib.h>

#include "mapping_tofi_trig.hh"
extern unsigned g_tofi_trig_map[2][24];
void tofi_trig_map_setup();

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

namespace
{
    double c_range_ns = 2048. * 5.;
    double c_bar_coincidence_ns = 20.; // nanoseconds.
} // namespace

R3BTofiCal2Histo::R3BTofiCal2Histo()
    : FairTask("R3BTofiCal2Histo", 1)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fTpat(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fNEvents(0)
    , fCal_Par(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fTofiY(0.)
    , fTofiQ(0.)
    , fTofiSmiley(true)
    , fTofiZ(false)
    , fParaFile("")
{
    for (Int_t i = 0; i < N_TOFI_HIT_PLANE_MAX; i++)
    {
        fh_Tofi_TotPm[i] = NULL;
        fhTdiff[i] = NULL;
        fhTsync[i] = NULL;
        for (Int_t j = 0; j < N_TOFI_HIT_PADDLE_MAX; j++)
        {
            fhLogTot1vsLogTot2[i][j] = NULL;
            fhSqrtQvsPosToT[i][j] = NULL;
            fhQvsPos[i][j] = NULL;
            fhTot1vsPos[i][j] = NULL;
            fhTot2vsPos[i][j] = NULL;
            // fhTot1vsTot2[i][j] = NULL;
        }
    }
}

R3BTofiCal2Histo::R3BTofiCal2Histo(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fUpdateRate(1000000)
    , fMinStats(100000)
    , fTrigger(-1)
    , fTpat(-1)
    , fNofPlanes(5)
    , fPaddlesPerPlane(6)
    , fNEvents(0)
    , fCal_Par(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fTofiY(0.)
    , fTofiQ(0.)
    , fTofiSmiley(true)
    , fTofiZ(false)
    , fParaFile("")
    , maxevent(0)
{
    for (Int_t i = 0; i < N_TOFI_HIT_PLANE_MAX; i++)
    {
        fh_Tofi_TotPm[i] = NULL;
        fhTdiff[i] = NULL;
        fhTsync[i] = NULL;
        for (Int_t j = 0; j < N_TOFI_HIT_PADDLE_MAX; j++)
        {
            fhLogTot1vsLogTot2[i][j] = NULL;
            fhSqrtQvsPosToT[i][j] = NULL;
            fhQvsPos[i][j] = NULL;
            fhTot1vsPos[i][j] = NULL;
            fhTot2vsPos[i][j] = NULL;
            // fhTot1vsTot2[i][j] = NULL;
        }
    }
}

R3BTofiCal2Histo::~R3BTofiCal2Histo()
{
    for (Int_t i = 0; i < N_TOFI_HIT_PLANE_MAX; i++)
    {
        if (fh_Tofi_TotPm[i])
            delete fh_Tofi_TotPm[i];
        if (fhTdiff[i])
            delete fhTdiff[i];
        if (fhTsync[i])
            delete fhTsync[i];
        for (Int_t j = 0; j < N_TOFI_HIT_PADDLE_MAX; j++)
        {
            if (fhLogTot1vsLogTot2[i][j])
                delete fhLogTot1vsLogTot2[i][j];
            if (fhSqrtQvsPosToT[i][j])
                delete fhSqrtQvsPosToT[i][j];
            if (fhQvsPos[i][j])
                delete fhQvsPos[i][j];
            if (fhTot1vsPos[i][j])
                delete fhTot1vsPos[i][j];
            if (fhTot2vsPos[i][j])
                delete fhTot2vsPos[i][j];
            /*
            if (fhTot1vsTot2[i][j])
                delete fhTot1vsTot2[i][j];
            */
        }
    }
    if (fCal_Par)
    {
        delete fCal_Par;
    }
}

InitStatus R3BTofiCal2Histo::Init()
{
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        return kFATAL;
    }

    header = dynamic_cast<R3BEventHeader*>(rm->GetObject("R3BEventHeader"));
    // may be = NULL!

    fCalData = dynamic_cast<TClonesArray*>(rm->GetObject("TofiCal"));
    if (NULL == fCalData)
    {
        LOG(fatal) << "Branch TofiCal not found";
    }
    fCalTriggerItems = dynamic_cast<TClonesArray*>(rm->GetObject("TofiTriggerCal"));
    if (NULL == fCalTriggerItems)
        LOG(fatal) << "Branch TofiTriggerCal not found";

    maxevent = rm->CheckMaxEventNo();

    if (!fNofModules)
    {
        LOG(error) << "R3BTofiCal2Histo::Init() Number of modules not set. ";
        return kFATAL;
    }

    tofi_trig_map_setup();

    return kSUCCESS;
}

void R3BTofiCal2Histo::SetParContainers()
{
    // container needs to be created in tcal/R3BTCalContFact.cxx AND R3BTCal needs
    // to be set as dependency in CMakelists.txt (in this case in the tof directory)
    fCal_Par = dynamic_cast<R3BTofiHitPar*>(FairRuntimeDb::instance()->getContainer("TofiHitPar"));
    if (!fCal_Par)
    {
        LOG(error) << "R3BTofiCal2Histo::Init() Couldn't get handle on TofiHitPar. ";
    }
    //    fCal_Par->setChanged();
}

namespace
{
    uint64_t n1, n2;
};

void R3BTofiCal2Histo::Exec(Option_t* option)
{

    if (fNEvents / 10000. == fNEvents / 10000)
        std::cout << "Events: " << fNEvents << " / " << maxevent << " (" << (Int_t)(fNEvents * 100. / maxevent)
                  << " %)             \r" << std::flush;

    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    if (fTpat_bit >= 0)
    {
        Int_t itpat = header->GetTpat();
        Int_t tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if ((header) && (tpatvalue == 0))
            return;
    }

    // Tofi detector
    std::vector<std::vector<UInt_t>> multihits(N_TOFI_HIT_PLANE_MAX, std::vector<UInt_t>(N_TOFI_HIT_PADDLE_MAX));
    UInt_t nHits = fCalData->GetEntries();
    // Organize cals into bars.
    struct Entry
    {
        std::vector<R3BTofiCalData*> top;
        std::vector<R3BTofiCalData*> bot;
    };
    std::map<size_t, Entry> bar_map;
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto* hit = dynamic_cast<R3BTofiCalData*>(fCalData->At(ihit));
        size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();
        // std::cout << "Hits: " << hit->GetDetectorId() << ' ' << hit->GetBarId() << ' ' << hit->GetSideId() << ' '
        //          << hit->GetTimeLeading_ns() << ' ' << hit->GetTimeTrailing_ns() << '\n';
        auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
        auto& vec = 1 == hit->GetSideId() ? ret.first->second.top : ret.first->second.bot;
        vec.push_back(hit);
    }

    static bool s_was_trig_missing = false;
    UInt_t trig_num = fCalTriggerItems->GetEntriesFast();
    // if (trig_num>0) cout << "Trigger information: " << trig_num << endl;
    // Find coincident PMT hits.
    // std::cout << "Print:\n";
    for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
    {
        // for (auto it2 = it->second.top.begin(); it->second.top.end() != it2; ++it2) {
        // std::cout << "Top: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        //}
        // for (auto it2 = it->second.bot.begin(); it->second.bot.end() != it2; ++it2) {
        // std::cout << "Bot: " << (*it2)->GetDetectorId() << ' ' << (*it2)->GetBarId() << ' ' <<
        // (*it2)->GetTimeLeading_ns() << '\n';
        //}
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
            // cout << "Triger mapping: " << top_trig_i << "  "  << bot_trig_i << endl;
            // cout << "Trigger num: " << trig_num << endl;
            Double_t top_trig_ns = 0, bot_trig_ns = 0;
            if (top_trig_i < trig_num && bot_trig_i < trig_num)
            {
                auto top_trig = dynamic_cast<R3BTofiCalData const*>(fCalTriggerItems->At(top_trig_i));
                auto bot_trig = dynamic_cast<R3BTofiCalData const*>(fCalTriggerItems->At(bot_trig_i));
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
                    LOG(error) << "R3BTofiCal2Histo::Exec() : Missing trigger information!";
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
                        if(top_ns>0 || bot_ns>0){
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
                // Hit!
                // std::cout << "Hit!\n";
                Int_t iPlane = top->GetDetectorId(); // 1..n
                Int_t iBar = top->GetBarId();        // 1..n
                if (iPlane > fNofPlanes)             // this also errors for iDetector==0
                {
                    // LOG(error) << "R3BTofiCal2HitPar::Exec() : more detectors than expected! Det: " << iPlane
                    //           << " allowed are 1.." << fNofPlanes;
                    ++top_i;
                    ++bot_i;
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    // LOG(error) << "R3BTofiCal2HitPar::Exec() : more bars then expected! Det: " << iBar
                    //           << " allowed are 1.." << fPaddlesPerPlane;
                    ++top_i;
                    ++bot_i;
                    continue;
                }
                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                // register multi hits
                multihits[iPlane - 1][iBar - 1] += 1;

                // prepare offset and sync calculation
                if (fTofiQ < 0.1)
                {
                    LOG(debug) << "Fill histo for offset and sync calculation plane " << iPlane << " bar " << iBar;
                    // calculate tdiff
                    auto tdiff = bot_ns - top_ns;

                    // create histograms
                    CreateHistograms(iPlane, iBar);

                    // fill control histograms
                    // fhTot1vsTot2[iPlane - 1][iBar - 1]->Fill(top_tot, bot_tot);
                    fhLogTot1vsLogTot2[iPlane - 1][iBar - 1]->Fill(TMath::Log(top_tot), TMath::Log(bot_tot));
                    fh_Tofi_TotPm[iPlane - 1]->Fill(iBar, top_tot);
                    fh_Tofi_TotPm[iPlane - 1]->Fill(-iBar - 1, bot_tot);

                    // Time differences of one paddle; offset  histo
                    fhTdiff[iPlane - 1]->Fill(iBar, tdiff);

                    // offset histo via ToT
                    auto posToT = 0.;
                    if (fTofiY == 0.)
                        posToT = TMath::Log(top_tot / bot_tot);
                    else
                    {
                        R3BTofiHitModulePar* par = fCal_Par->GetModuleParAt(iPlane, iBar);
                        if (!par)
                        {
                            LOG(info) << "R3BTofiCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                                      << ", Bar: " << top->GetBarId();
                            ++top_i;
                            ++bot_i;
                            continue;
                        }
                        posToT = log((top_tot + par->GetToTOffset2()) / (bot_tot + par->GetToTOffset1()));
                    }
                    fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->Fill(posToT, sqrt(top_tot * bot_tot));
                    // fhTot1vsPos[iPlane - 1][iBar - 1]->Fill(posToT, bot_tot);
                    // fhTot2vsPos[iPlane - 1][iBar - 1]->Fill(posToT, top_tot);

                    // Time of hit
                    auto THit = (top_ns + bot_ns) / 2.;
                    fhTsync[iPlane - 1]->Fill(iBar, THit);
                }
                else
                {
                    // get parameter
                    R3BTofiHitModulePar* para = fCal_Par->GetModuleParAt(iPlane, iBar);
                    if (!para)
                    {
                        LOG(info) << "R3BTofiCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                                  << ", Bar: " << top->GetBarId();
                        ++top_i;
                        ++bot_i;
                        continue;
                    }

                    // calculate tdiff with offest
                    auto tdiff = (bot_ns + para->GetOffset1()) - (top_ns + para->GetOffset2());

                    // walk corrections
                    if (para->GetPar1Walk() == 0. || para->GetPar2Walk() == 0. || para->GetPar3Walk() == 0. ||
                        para->GetPar4Walk() == 0. || para->GetPar5Walk() == 0.)
                        LOG(info) << "Walk correction not found!";
                    bot_ns = bot_ns - walk(bot_tot,
                                           para->GetPar1Walk(),
                                           para->GetPar2Walk(),
                                           para->GetPar3Walk(),
                                           para->GetPar4Walk(),
                                           para->GetPar5Walk());
                    top_ns = top_ns - walk(top_tot,
                                           para->GetPar1Walk(),
                                           para->GetPar2Walk(),
                                           para->GetPar3Walk(),
                                           para->GetPar4Walk(),
                                           para->GetPar5Walk());

                    auto posToT =
                        para->GetLambda() * log((top_tot + para->GetToTOffset2()) / (bot_tot + para->GetToTOffset1()));

                    // create histograms
                    CreateHistograms(iPlane, iBar);

                    // fill control histograms
                    // fhTot1vsTot2[iPlane - 1][iBar - 1]->Fill(top_tot, bot_tot);
                    fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->Fill(posToT, sqrt(top_tot * bot_tot));

                    // Time differences of one paddle
                    fhTdiff[iPlane - 1]->Fill(iBar, tdiff);

                    // Time of hit
                    auto THit = (bot_ns + top_ns) / 2. - para->GetSync();

                    // Sync of one plane
                    fhTsync[iPlane - 1]->Fill(iBar, THit);
                }

                // prepare double exponential fit
                if (!fTofiSmiley)
                {

                    LOG(debug) << "Prepare histo for double exponential fit";
                    R3BTofiHitModulePar* par = fCal_Par->GetModuleParAt(iPlane, iBar);
                    if (!par)
                    {
                        LOG(info) << "R3BTofiCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                                  << ", Bar: " << top->GetBarId();
                        ++top_i;
                        ++bot_i;
                        continue;
                    }

                    // calculate y position
                    auto pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();

                    // fill fitting histograms and smiley histogram
                    fhTot1vsPos[iPlane - 1][iBar - 1]->Fill(pos, bot_tot);
                    fhTot2vsPos[iPlane - 1][iBar - 1]->Fill(pos, top_tot);
                }
                // prepare charge fit / quench correction
                if (fTofiZ == true)
                {
                    LOG(debug) << "Prepare histo for quenching correction";
                    // get parameter
                    R3BTofiHitModulePar* par = fCal_Par->GetModuleParAt(iPlane, iBar);
                    if (!par)
                    {
                        LOG(info) << "R3BTofiCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
                                  << ", Bar: " << top->GetBarId();
                        ++top_i;
                        ++bot_i;
                        continue;
                    }

                    // calculate y position
                    auto posToT =
                        par->GetLambda() * log((top_tot * par->GetToTOffset2()) / (bot_tot * par->GetToTOffset1()));

                    Double_t parq[4];
                    parq[0] = par->GetPar1a();
                    parq[1] = par->GetPar1b();
                    parq[2] = par->GetPar1c();
                    parq[3] = par->GetPar1d();

                    // calculate charge
                    Double_t qb = 0.;
                    if (fTofiSmiley)
                    {
                        qb = TMath::Sqrt(top_tot * bot_tot) /
                             (parq[0] + parq[1] * posToT + parq[2] * pow(posToT, 2) + parq[3] * pow(posToT, 3));
                        qb = qb *
                             fTofiQ; // theory says: dE ~ Z^2 but we see quenching -> just use linear and fit the rest
                    }
                    else
                    {
                        // double exponential
                        auto pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();
                        auto q1 = bot_tot /
                                  (parq[0] * (exp(-parq[1] * (pos + 100.)) + exp(-parq[2] * (pos + 100.))) + parq[3]);
                        parq[0] = par->GetPar2a();
                        parq[1] = par->GetPar2b();
                        parq[2] = par->GetPar2c();
                        parq[3] = par->GetPar2d();
                        auto q2 = top_tot /
                                  (parq[0] * (exp(-parq[1] * (pos + 100.)) + exp(-parq[2] * (pos + 100.))) + parq[3]);
                        q1 = q1 *
                             fTofiQ; // theory says: dE ~ Z^2 but we see quenching -> just use linear and fit the rest
                        q2 = q2 * fTofiQ;
                        qb = (q1 + q2) / 2.;
                    }

                    // fill control histograms and Q vs Pos without multihits
                    if (multihits[iPlane - 1][iBar - 1] < 2 && (qb > 0.))
                    {
                        fhQvsPos[iPlane - 1][iBar - 1]->Fill(posToT, qb);
                    }
                }

                ++top_i;
                ++bot_i;

                // Increment events
                fNEvents += 1;
            }
            else if (dt < 0 && dt > -c_range_ns / 2)
            {
                ++top_i;
            }
            else
            {
                ++bot_i;
            }
        }
    }
}

void R3BTofiCal2Histo::CreateHistograms(Int_t iPlane, Int_t iBar)
{
    Double_t max_charge = 60.;
    if (NULL == fhTot1vsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Tot1_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        if (iPlane < 3)
            fhTot1vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 200, -100, 100, 400, 0., 200.);
        if (iPlane > 2)
            fhTot1vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 200, -100, 100, 400, 0., 200.);
        fhTot1vsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Pos in cm");
        fhTot1vsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToT of PM1 in ns");
    }
    if (NULL == fhTot2vsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Tot2_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        if (iPlane < 3)
            fhTot2vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 200, -100, 100, 400, 0., 200.);
        if (iPlane > 2)
            fhTot2vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 200, -100, 100, 400, 0., 200.);
        fhTot2vsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Pos in cm");
        fhTot2vsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToT of PM2 in ns");
    }
    if (NULL == fhTdiff[iPlane - 1])
    {
        char strName1[255];
        char strName2[255];
        sprintf(strName1, "Time_Diff_Plane_%d", iPlane);
        sprintf(strName2, "Time Diff Plane %d", iPlane);
        fhTdiff[iPlane - 1] = new TH2F(strName1, strName2, 50, 0, 50, 4000, -20., 20.);
        fhTdiff[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTdiff[iPlane - 1]->GetYaxis()->SetTitle("Time difference (PM1 - PM2) in ns");
    }
    if (NULL == fhTsync[iPlane - 1])
    {
        char strName[255];
        char strName2[255];
        sprintf(strName, "Time_Sync_Plane_%d", iPlane);
        sprintf(strName2, "Time Sync Plane %d", iPlane);
        fhTsync[iPlane - 1] = new TH2F(strName, strName2, 50, 0, 50, 4000, -2000, 2000.);
        fhTsync[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTsync[iPlane - 1]->GetYaxis()->SetTitle("THit in ns");
    }
    if (NULL == fh_Tofi_TotPm[iPlane - 1])
    {
        char strName[255];
        sprintf(strName, "Tofi_ToT_plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Tofi ToT plane %d", iPlane);
        fh_Tofi_TotPm[iPlane - 1] = new TH2F(strName, strName2, 90, -45, 45, 300, 0., 300.);
        fh_Tofi_TotPm[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fh_Tofi_TotPm[iPlane - 1]->GetYaxis()->SetTitle("ToT / ns");
    }
    if (NULL == fhLogTot1vsLogTot2[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Plane_%d_Bar_%d_LogToT1vsLogToT2", iPlane, iBar);
        fhLogTot1vsLogTot2[iPlane - 1][iBar - 1] = new TH2F(strName, "", 400, 2., 6., 400, 2., 6.);
        fhLogTot1vsLogTot2[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Log(ToT) of PM2");
        fhLogTot1vsLogTot2[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Log(ToT) of PM1");
    }
    if (NULL == fhSqrtQvsPosToT[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "SqrtQ_vs_PosToT_Plane_%d_Bar_%d", iPlane, iBar);
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1] =
            new TH2F(strName, "", 20000, -100, 100, max_charge * 4, 0., max_charge * 4);
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("sqrt(PM1*PM2)");
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Position from ToT in cm");
    }
    if (NULL == fhQvsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 20000, -100, 100, max_charge * 10, 0., max_charge);
        fhQvsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Position in cm");
    }
    /*
    if (NULL == fhTot1vsTot2[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Plane_%d_Bar_%d_ToT1vsToT2", iPlane, iBar);
        fhTot1vsTot2[iPlane - 1][iBar - 1] = new TH2F(strName, "", 300, 0., 300., 300, 0., 300.);
        fhTot1vsTot2[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT of PM2 in ns");
        fhTot1vsTot2[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToT of PM1 in ns");
    }
    */
}

void R3BTofiCal2Histo::FinishEvent() {}

void R3BTofiCal2Histo::FinishTask()
{
    for (Int_t i = 0; i < N_TOFI_HIT_PLANE_MAX; i++)
    {
        if (fh_Tofi_TotPm[i])
            fh_Tofi_TotPm[i]->Write(); // control histogram for ToT
        if (fhTsync[i])
            fhTsync[i]->Write(); // histogram for sync calculation
        if (fhTdiff[i])
            fhTdiff[i]->Write(); // histogram for offset and veff calculation
        for (Int_t j = 0; j < N_TOFI_HIT_PADDLE_MAX; j++)
        {
            if (fhLogTot1vsLogTot2[i][j])
                fhLogTot1vsLogTot2[i][j]->Write(); // control histogram Log(ToT) Pm1 vs Log(ToT) Pm2
            if (fhSqrtQvsPosToT[i][j])
                fhSqrtQvsPosToT[i][j]->Write(); // histogram for ToT offset calculation
            if (fhQvsPos[i][j])
                fhQvsPos[i][j]->Write(); // histogram for charge fit
            if (fhTot1vsPos[i][j])
                fhTot1vsPos[i][j]->Write(); // histogram for position dependence of charge 1
            if (fhTot2vsPos[i][j])
                fhTot2vsPos[i][j]->Write(); // histogram for position dependence of charge 2
            /*
            if (fhTot1vsTot2[i][j])
                fhTot1vsTot2[i][j]->Write(); // control histogram ToT Pm1 vs ToT Pm2
            */
        }
    }
}

Double_t R3BTofiCal2Histo::walk(Double_t Q, Double_t par1, Double_t par2, Double_t par3, Double_t par4, Double_t par5)
{
    Double_t y = 0;
    y = -30.2 + par1 * TMath::Power(Q, par2) + par3 / Q + par4 * Q + par5 * Q * Q;
    return y;
}

Double_t R3BTofiCal2Histo::saturation(Double_t x)
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

ClassImp(R3BTofiCal2Histo)
