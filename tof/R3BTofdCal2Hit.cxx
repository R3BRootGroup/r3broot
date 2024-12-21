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
// -----                 R3BTofdCal2Hit                   -----
// -----            Created May 2016 by M.Heil            -----
// -----           Modified Dec 2019 by L.Bott            -----
// ------------------------------------------------------------

#include "R3BTofdCal2Hit.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
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

using namespace std;
#define IS_NAN(x) TMath::IsNaN(x)

#define N_TOFD_HIT_PLANE_MAX 4
#define N_TOFD_HIT_PADDLE_MAX 44

namespace
{
    double c_range_ns = 2048 * 5;
    double c_bar_coincidence_ns = 20; // nanoseconds.
} // namespace

R3BTofdCal2Hit::R3BTofdCal2Hit()
    : R3BTofdCal2Hit("R3BTofdCal2Hit", 1)
{
}

R3BTofdCal2Hit::R3BTofdCal2Hit(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItems(NULL)
    , fCalItemsLos(NULL)
    , fHitItemsLos(NULL)
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
    , countloshit(0)
    , wrongtrigger(0)
    , wrongtpat(0)
    , headertpat(0)
    , events_in_cal_level(0)
    , inbarcoincidence(0)
    , countreset(0)
    , hitsbeforereset(0)
    , eventstore(0)
    , incoincidence(0)
    , inaverage12(0)
    , inaverage34(0)
    , singlehit(0)
    , multihit(0)
    , bars_with_multihit(0)
    , events_wo_tofd_hits(0)
{
    fhLosXYP = NULL;
    fhChargeLosTofD = NULL;
    fh_los_pos = NULL;
    if (fTofdHisto)
    {
        fhxy12 = NULL;
        fhxy12tot = NULL;
        fhxy34 = NULL;
        fhxy34tot = NULL;
        fhCharge = NULL;
        fhAverageCharge = NULL;
        for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
        {
            fhQ[i] = NULL;
            fhxy[i] = NULL;
            fhQvsEvent[i] = NULL;
            fhQM[i] = NULL;
            fhMvsQ[i] = NULL;
            fhTdiff[i] = NULL;
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {
                fhQvsPos[i][j] = NULL;
                fhTdiffvsQ[i][2 * j] = NULL;
                fhTdiffvsQ[i][2 * j + 1] = NULL;
                fhQvsQ[i][2 * j] = NULL;
                fhQvsQ[i][2 * j + 1] = NULL;
                fhQvsTof[i][j] = NULL;
                fhTvsTof[i][j] = NULL;
                fhToTvsTofw[i][j] = NULL;
            }
        }
    }
}

R3BTofdCal2Hit::~R3BTofdCal2Hit()
{
    if (fhLosXYP)
        delete fhLosXYP;
    if (fhChargeLosTofD)
        delete fhChargeLosTofD;
    if (fh_los_pos)
        delete fh_los_pos;
    if (fTofdHisto)
    {
        if (fhxy12)
            delete fhxy12;
        if (fhxy12tot)
            delete fhxy12tot;
        if (fhxy34)
            delete fhxy34;
        if (fhxy34tot)
            delete fhxy34tot;
        if (fhAverageCharge)
            delete fhAverageCharge;
        if (fhCharge)
            delete fhCharge;
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            if (fhQ[i])
                delete fhQ[i];
            if (fhxy[i])
                delete fhxy[i];
            if (fhQvsEvent[i])
                delete fhQvsEvent[i];
            if (fhQM[i])
                delete fhQM[i];
            if (fhMvsQ[i])
                delete fhMvsQ[i];
            if (fhTdiff[i])
                delete fhTdiff[i];
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {
                if (fhQvsPos[i][j])
                    delete fhQvsPos[i][j];
                if (fhTdiffvsQ[i][2 * j])
                    delete fhTdiffvsQ[i][2 * j];
                if (fhTdiffvsQ[i][2 * j + 1])
                    delete fhTdiffvsQ[i][2 * j + 1];
                if (fhQvsQ[i][2 * j])
                    delete fhQvsQ[i][2 * j];
                if (fhQvsQ[i][2 * j + 1])
                    delete fhQvsQ[i][2 * j + 1];
                if (fhQvsTof[i][j])
                    delete fhQvsTof[i][j];
                if (fhTvsTof[i][j])
                    delete fhTvsTof[i][j];
                if (fhToTvsTofw[i][j])
                    delete fhQvsTof[i][j];
            }
        }
    }
    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

InitStatus R3BTofdCal2Hit::Init()
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
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));
    if (header == nullptr)
    {
        header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));
        R3BLOG(warn, "R3BEventHeader was found instead of EventHeader.");
    }

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("TofdCal"));
    if (NULL == fCalItems)
        LOG(fatal) << "Branch TofdCal not found";
    maxevent = mgr->CheckMaxEventNo();
    fCalItemsLos = dynamic_cast<TClonesArray*>(mgr->GetObject("LosCal"));
    if (NULL == fCalItemsLos)
        LOG(warn) << "Branch LosCal not found";
    fHitItemsLos = dynamic_cast<TClonesArray*>(mgr->GetObject("LosHit"));
    if (NULL == fHitItemsLos)
        LOG(warn) << "Branch LosHit not found";

    // request storage of Hit data in output tree
    mgr->Register("TofdHit", "Tofd hit data", fHitItems, kTRUE);

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BTofdCal2Hit::SetParContainers()
{
    fHitPar = dynamic_cast<R3BTofdHitPar*>(FairRuntimeDb::instance()->getContainer("TofdHitPar"));
    if (!fHitPar)
    {
        LOG(error) << "Could not get access to TofdHitPar-Container.";
        fNofHitPars = 0;
        return;
    }
}

InitStatus R3BTofdCal2Hit::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BTofdCal2Hit::Exec(Option_t* option)
{
    if (fnEvents / 100000. == (int)fnEvents / 100000)
        std::cout << "\rEvents: " << fnEvents << " / " << maxevent << " (" << (int)(fnEvents * 100. / maxevent)
                  << " %) " << std::flush;

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
    Double_t timeRef = 0.;
    Double_t timeLos = 0;
    Double_t timeP0 = 0.;
    Double_t LosTresM = 0;
    Double_t LosQ = 0;
    Double_t xLosP = 1000;
    Double_t yLosP = 1000;
    Double_t randx;
    std::vector<std::vector<std::vector<Double_t>>> q;
    std::vector<std::vector<std::vector<Double_t>>> tof;
    std::vector<std::vector<std::vector<Double_t>>> x;
    std::vector<std::vector<std::vector<Double_t>>> y;
    std::vector<std::vector<std::vector<Double_t>>> yToT;
    UInt_t vmultihits[N_PLANE_MAX + 1][N_TOFD_HIT_PADDLE_MAX * 2 + 1];
    for (Int_t i = 0; i <= fNofPlanes; i++)
    {
        q.push_back(std::vector<std::vector<Double_t>>());
        tof.push_back(std::vector<std::vector<Double_t>>());
        x.push_back(std::vector<std::vector<Double_t>>());
        y.push_back(std::vector<std::vector<Double_t>>());
        yToT.push_back(std::vector<std::vector<Double_t>>());
        for (Int_t j = 0; j <= 2 * N_TOFD_HIT_PADDLE_MAX; j++)
        {
            vmultihits[i][j] = 0;
            q[i].push_back(std::vector<Double_t>());
            tof[i].push_back(std::vector<Double_t>());
            x[i].push_back(std::vector<Double_t>());
            y[i].push_back(std::vector<Double_t>());
            yToT[i].push_back(std::vector<Double_t>());
        }
    }
    if (fHitItemsLos)
    {
        Int_t nHits = fHitItemsLos->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            countloshit++;
            LOG(warn) << "LOS Ihit  " << ihit << " " << nHits;
            R3BLosHitData* hitData = dynamic_cast<R3BLosHitData*>(fHitItemsLos->At(ihit));
            if (ihit == 0)
                timeLos = hitData->fTime_ns;
            LOG(warn) << "LOS Time " << timeLos;
            if (std::isnan(timeLos))
                return; /// in s444 run2930 Event 22333208
            if (NULL == fh_los_pos)
            {
                char strName[255];
                sprintf(strName, "LOS_X_vs_Y_MCFD");
                //     fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);
                fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);
            }
            fh_los_pos->Fill(hitData->fX_cm, hitData->fY_cm);

            if (ihit == 0)
                LosQ = hitData->fZ;
            if (NULL == fhChargeLosTofD)
            {
                char strName[255];
                sprintf(strName, "LosQ_vs_TofdQ");
                fhChargeLosTofD = new TH2F(strName, "", 400, 0., 100., 800, 0., 200.);
                fhChargeLosTofD->GetYaxis()->SetTitle("Charge LOS");
                fhChargeLosTofD->GetXaxis()->SetTitle("Charge ToFD");
            }
        }
    }
    /*
        if(fCalItemsLos)
        {
            Int_t nHits = fCalItemsLos->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                LOG(warn) << "LOS Ihit  "<< ihit<<" "<<nHits<<FairLogger::endl;
                R3BLosCalData *calData = (R3BLosCalData*)fCalItemsLos->At(ihit);
                timeLos=(calData->fTimeV_r_ns+calData->fTimeV_l_ns+calData->fTimeV_t_ns+calData->fTimeV_b_ns)/4.;
                LosTresM=(calData->fTimeV_r_ns+calData->fTimeV_l_ns)/2.-(calData->fTimeV_t_ns+calData->fTimeV_b_ns)/2.;
                LOG(warn) << "LOS MCFD  "<< LosTresM<<" "<<timeLos<<FairLogger::endl;
                if (NULL == fh_los_pos)
                {
                    char strName[255];
                    sprintf(strName, "LOS_X_vs_Y_MCFD");
                    //     fh_los_pos = new TH2F(strName, "", 2000, -10., 10., 2000, -10., 10.);
                    fh_los_pos = new TH2F(strName, "", 2000, -100., 100., 2000, -100., 100.);
                }
                fh_los_pos->Fill((calData->fTimeV_r_ns-calData->fTimeV_l_ns+2.15)*5.3,
                                 (calData->fTimeV_b_ns-calData->fTimeV_t_ns+1.28)*5.2);   // Xe: (2.15,5.3), (1.28,5.2)
                if (NULL == fhLosTimeM )
                {
                    char strName[255];
                    sprintf(strName, "LOS_Time_Resolution_MCFD");
                    //     fhLosTimeM = new TH1F(strName, "", 4000, -2., 2.);
                    fhLosTimeM = new TH1F(strName, "", 4000, -2000., 2000.);
                }
                fhLosTimeM->Fill(LosTresM);
            }
        }
    */

    // std::cout<<"new event!*************************************\n";

    Int_t nHits = fCalItems->GetEntries();
    Int_t nHitsEvent = 0;
    // Organize cals into bars.
    struct Entry
    {
        std::vector<R3BTofdCalData*> top;
        std::vector<R3BTofdCalData*> bot;
    };
    std::map<size_t, Entry> bar_map;
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto* hit = dynamic_cast<R3BTofdCalData*>(fCalItems->At(ihit));
        size_t idx = hit->GetDetectorId() * fPaddlesPerPlane * hit->GetBarId();
        // std::cout << "Hits: " << hit->GetDetectorId() << ' ' << hit->GetBarId() << ' ' << hit->GetSideId() << ' '
        //          << hit->GetTimeLeading_ns() << ' ' << hit->GetTimeTrailing_ns() << '\n';
        auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
        auto& vec = 1 == hit->GetSideId() ? ret.first->second.top : ret.first->second.bot;
        vec.push_back(hit);
        events_in_cal_level++;
    }

    // Find coincident PMT hits.
    // std::cout << "Print:\n";
    for (auto it = bar_map.begin(); bar_map.end() != it; ++it)
    {
    reset:
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
            auto top_ns = top->GetTimeLeading_ns();
            auto bot_ns = bot->GetTimeLeading_ns();
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
                    LOG(error) << "R3BTofdCal2HitPar::Exec() : more detectors than expected! Det: " << iPlane
                               << " allowed are 1.." << fNofPlanes;
                    continue;
                }
                if (iBar > fPaddlesPerPlane) // same here
                {
                    LOG(error) << "R3BTofdCal2HitPar::Exec() : more bars then expected! Det: " << iBar
                               << " allowed are 1.." << fPaddlesPerPlane;
                    continue;
                }

                auto top_tot = fmod(top->GetTimeTrailing_ns() - top->GetTimeLeading_ns() + c_range_ns, c_range_ns);
                auto bot_tot = fmod(bot->GetTimeTrailing_ns() - bot->GetTimeLeading_ns() + c_range_ns, c_range_ns);

                // register multi hits
                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
                if (iPlane == 1 || iPlane == 3)
                    vmultihits[iPlane][iBar * 2 - 2] += 1;
                if (iPlane == 2 || iPlane == 4)
                    vmultihits[iPlane][iBar * 2] += 1;
                vmultihits[iPlane][iBar * 2 - 1] += 1;

                nHitsEvent += 1;
                R3BTofdHitModulePar* par = fHitPar->GetModuleParAt(iPlane, iBar);
                if (!par)
                {
                    LOG(info) << "R3BTofdCal2Hit::Exec : Hit par not found, Plane: " << top->GetDetectorId()
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

                // calculate time-of-flight
                if (timeLos == 0)
                    LOG(warn) << "Los Time is zero! ";
                Double_t ToF = (bot_ns + top_ns) / 2. - timeLos - par->GetSync();
                if (std::isnan(ToF))
                {
                    LOG(fatal) << "ToFD ToF not found";
                }
                if (timeP0 == 0.)
                    timeP0 = ToF;

                if (timeLos == 0)
                { // no LOS in s454
                    /// What to do here:
                    /// check if all ToF in one event are in a range of 3000ns (readout window) and shift times
                    /// according to that
                    ///
                    ///       this is the first hit
                    ///       I
                    /// e.g. 171; 9439; 179; 1117; 175 -->> 171+c_range_ns; 9439; 179+c_range_ns; 1117+c_range_ns;
                    /// 175+c_range_ns
                    ///             I
                    ///             this should be the first hit -> counter resets -> other hits follow
                    if (ToF - timeP0 < -3000.)
                    {
                        ToF += c_range_ns;
                    }
                    if (ToF - timeP0 > 3000.)
                    {
                        timeP0 = ToF;
                        it = bar_map.begin();
                        countreset++;
                        hitsbeforereset += nHitsEvent;
                        for (Int_t i = 0; i <= fNofPlanes; i++)
                        {
                            for (Int_t j = 0; j <= 2 * N_TOFD_HIT_PADDLE_MAX; j++)
                            {
                                tof[i][j].clear();
                                x[i][j].clear();
                                y[i][j].clear();
                                yToT[i][j].clear();
                                q[i][j].clear();
                                vmultihits[i][j] = 0;
                                nHitsEvent = 0;
                            }
                        }
                        LOG(warn) << "Found new first hit -> will reset";
                        goto reset; /// TODO: how to do without goto?
                    }
                }

                if (timeLos != 0)
                {
                    LOG(debug) << "Found LOS detector";
                    LOG(debug) << "check for coincidence range: TOF: " << ToF << " c_range: " << c_range_ns << "\n";
                    while (ToF < -c_range_ns / 2)
                    {
                        ToF += c_range_ns;
                        LOG(debug) << "Shift up\n";
                    }
                    while (ToF > c_range_ns / 2)
                    {
                        ToF -= c_range_ns;
                        LOG(debug) << "Shift down\n";
                    }
                }

                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
                if (iPlane == 1 || iPlane == 3)
                    tof[iPlane][iBar * 2 - 2].push_back(ToF);
                if (iPlane == 2 || iPlane == 4)
                    tof[iPlane][iBar * 2].push_back(ToF);
                tof[iPlane][iBar * 2 - 1].push_back(ToF);

                // calculate y-position
                auto pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();

                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
                if (iPlane == 1 || iPlane == 3)
                    y[iPlane][iBar * 2 - 2].push_back(pos);
                if (iPlane == 2 || iPlane == 4)
                    y[iPlane][iBar * 2].push_back(pos);
                y[iPlane][iBar * 2 - 1].push_back(pos);

                // calculate y-position from ToT
                auto posToT =
                    par->GetLambda() * log((top_tot * par->GetToTOffset2()) / (bot_tot * par->GetToTOffset1()));
                if (iPlane == 1 || iPlane == 3)
                    yToT[iPlane][iBar * 2 - 2].push_back(posToT);
                if (iPlane == 2 || iPlane == 4)
                    yToT[iPlane][iBar * 2].push_back(posToT);
                yToT[iPlane][iBar * 2 - 1].push_back(posToT);

                if (fTofdTotPos)
                    pos = posToT;

                // calculate x-position
                randx = (std::rand() / (float)RAND_MAX);
                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
                if (iPlane == 1 || iPlane == 3)
                    x[iPlane][iBar * 2 - 2].push_back(iBar * 2.8 - 21. * 2.8 - 1.4 - 1.4 * randx);
                if (iPlane == 2 || iPlane == 4)
                    x[iPlane][iBar * 2].push_back(iBar * 2.8 - 21. * 2.8 + 1.4 - 1.4 * randx);
                x[iPlane][iBar * 2 - 1].push_back(iBar * 2.8 - 21. * 2.8 - 1.4 * randx);
                // correct for position dependence and calculate nuclear charge Z
                Double_t para[4];
                para[0] = par->GetPar1a();
                para[1] = par->GetPar1b();
                para[2] = par->GetPar1c();
                para[3] = par->GetPar1d();

                Double_t qb = 0.;
                if (fTofdTotPos)
                {
                    // via pol3
                    qb = TMath::Sqrt(top_tot * bot_tot) /
                         (para[0] + para[1] * pos + para[2] * pow(pos, 2) + para[3] * pow(pos, 3));
                    qb = qb * fTofdQ;
                }
                else
                {
                    // via double exponential:
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

                /*
                if(parz[0]>0 && parz[2]>0) std::cout<<"Charges in this event
                "<<parz[0]*TMath::Power((q1+q2)/2.,parz[2])+parz[1]<<" plane "<<iPlane<<" ibar "<<iBar<<"\n"; else
                std::cout<<"Charges in this event "<<(q1+q2)/2.<<" plane "<<iPlane<<" ibar "<<iBar<<"\n";
                std::cout<<"Times in this event "<<ToF<<" plane "<<iPlane<<" ibar "<<iBar<<endl;
                if(iPlane==1 || iPlane==3)std::cout<<"x in this event "<<iBar*2.8-21.*2.8-1.4-1.4*randx<<" plane
                "<<iPlane<<" ibar "<<iBar<<" vbars "<<iBar*2-2<<" "<<iBar*2-1<<"\n"; if(iPlane==2 ||
                iPlane==4)std::cout<<"x in this event "<<iBar*2.8-21.*2.8+1.4-1.4*randx<<" plane "<<iPlane<<" ibar
                "<<iBar<<" vbars "<<iBar*2-1<<" "<<iBar*2<<"\n"; std::cout<<"y in this event "<<pos<<" plane
                "<<iPlane<<" ibar "<<iBar<<endl<<endl;
                */

                // we increase the number of bars by a factor of 2 in order to compare planes with half bar width
                // overlap
                if (parz[0] > 0 && parz[2] > 0)
                {
                    if (iPlane == 1 || iPlane == 3)
                        q[iPlane][iBar * 2 - 2].push_back(parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                    if (iPlane == 2 || iPlane == 4)
                        q[iPlane][iBar * 2].push_back(parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                    q[iPlane][iBar * 2 - 1].push_back(parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                }
                else
                {
                    if (iPlane == 1 || iPlane == 3)
                        q[iPlane][iBar * 2 - 2].push_back(qb);
                    if (iPlane == 2 || iPlane == 4)
                        q[iPlane][iBar * 2].push_back(qb);
                    q[iPlane][iBar * 2 - 1].push_back(qb);
                    parz[0] = 1.;
                    parz[1] = 0.;
                    parz[2] = 1.;
                }

                if (fTofdHisto)
                {
                    // fill control histograms
                    CreateHistograms(iPlane, iBar);
                    // fhTof[iPlane-1]->Fill(iBar,ToF);
                    // fhTsync[iPlane-1]->Fill(iBar,ToF);
                    fhTdiff[iPlane - 1]->Fill(iBar, tdiff);
                    fhQvsPos[iPlane - 1][iBar - 1]->Fill(pos, parz[0] * TMath::Power(qb, parz[2]) + parz[1]);
                    fhQvsTof[iPlane - 1][iBar - 1]->Fill(qb, ToF);
                    fhTvsTof[iPlane - 1][iBar - 1]->Fill(dt_mod, ToF);
                    fhToTvsTofw[iPlane - 1][iBar - 1]->Fill((bot_tot + top_tot) / 2.,
                                                            ToF); // needed to get TOF w/o walk correction
                }

                // Time reference in case on has the master signal in one of the TDC channels.
                // Not used at the moment.
                timeRef = 0;

                ++top_i;
                ++bot_i;
            }
            else if (dt < 0 && dt > -c_range_ns / 2)
            {
                ++top_i;
                LOG(warn) << "Not in bar coincidence increase top counter";
            }
            else
            {
                ++bot_i;
                LOG(warn) << "Not in bar coincidence increase bot counter";
            }
        }
    }

    // Now all hits in this event are analyzed

    Double_t hit_coinc = 5000000.; // coincidence window for hits in one event in ns. physics says max 250 ps
    Double_t maxChargeDiff = 1.;   // maximum charge difference between two planes for averaged hits

    LOG(debug) << "Hits in this event: " << nHitsEvent;
    if (nHitsEvent == 0)
        events_wo_tofd_hits++;

    // init arrays to store hits
    Double_t tArrQ[2 * nHitsEvent + 1];
    Double_t tArrT[2 * nHitsEvent + 1];
    Double_t tArrX[2 * nHitsEvent + 1];
    Double_t tArrY[2 * nHitsEvent + 1];
    Double_t tArrYT[2 * nHitsEvent + 1];
    Double_t tArrP[2 * nHitsEvent + 1];
    Double_t tArrB[2 * nHitsEvent + 1];
    Bool_t tArrU[2 * nHitsEvent + 1];
    for (int i = 0; i < (2 * nHitsEvent + 1); i++)
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
        for (Int_t j = 0; j < fPaddlesPerPlane * 2 + 1; j++)
        {
            if (vmultihits[i][j] > 1)
            {
                bars_with_multihit++;
                multihit += vmultihits[i][j] - 1;
            }
        }
    }

    // order events for time
    for (Int_t i = 1; i <= fNofPlanes; i++)
    { // loop over planes i
        for (Int_t j = 0; j < fPaddlesPerPlane * 2 + 1; j++)
        { // loop over virtual paddles j
            if (tof[i][j].empty() == false)
            { // check paddle for entries
                for (Int_t m = 0; m < tof[i][j].size(); m++)
                { // loop over multihits m
                    Int_t p = 0;
                    if (tArrT[0] == -1.)
                    { // first entry
                        LOG(debug) << "First entry plane/bar " << i << "/" << j;
                        tArrQ[0] = q[i][j].at(m);
                        tArrT[0] = tof[i][j].at(m);
                        tArrX[0] = x[i][j].at(m);
                        tArrY[0] = y[i][j].at(m);
                        tArrYT[0] = yToT[i][j].at(m);
                        tArrP[0] = i;
                        tArrB[0] = j;
                    }
                    else
                    {
                        if (tof[i][j].at(m) < tArrT[0])
                        { // new first entry with smaller time
                            LOG(debug) << "Insert new first " << i << " " << j;
                            insertX(2 * nHitsEvent, tArrQ, q[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrT, tof[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrX, x[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrY, y[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrYT, yToT[i][j].at(m), 1);
                            insertX(2 * nHitsEvent, tArrP, i, 1);
                            insertX(2 * nHitsEvent, tArrB, j, 1);
                        }
                        else
                        {
                            while (tof[i][j].at(m) > tArrT[p] && tArrT[p] != -1.)
                            {
                                p++; // find insert position
                                if (p > 2 * nHitsEvent + 1)
                                    LOG(fatal) << "Insert position oor"; // should not happen
                            }

                            LOG(debug) << "Will insert at " << p;
                            if (p > 0 && tof[i][j].at(m) > tArrT[p - 1] && tof[i][j].at(m) != tArrT[p])
                            { // insert at right position
                                LOG(debug) << "Insert at " << p << " " << i << " " << j;
                                insertX(2 * nHitsEvent, tArrQ, q[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrT, tof[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrX, x[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrY, y[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrYT, yToT[i][j].at(m), p + 1);
                                insertX(2 * nHitsEvent, tArrP, i, p + 1);
                                insertX(2 * nHitsEvent, tArrB, j, p + 1);
                            }
                            else
                            {
                                if (tof[i][j].at(m) == tArrT[p])
                                { // handle virtual bars
                                    LOG(debug) << "Insert virtual bar " << i << " " << j;
                                    insertX(2 * nHitsEvent, tArrQ, q[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrT, tof[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrX, x[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrY, y[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrYT, yToT[i][j].at(m), p + 2);
                                    insertX(2 * nHitsEvent, tArrP, i, p + 2);
                                    insertX(2 * nHitsEvent, tArrB, j, p + 2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // print time sorted events
    /*
    if(tArrT[0]!=-1.){

            for (Int_t a = 0; a < 2*nHitsEvent; a++)
                std::cout << tArrQ[a] << " ";
            std::cout << "\n";
            for (Int_t a = 0; a < 2*nHitsEvent; a++)
                std::cout << tArrT[a] << " ";
            std::cout << "\n";
            for (Int_t a = 0; a < 2*nHitsEvent; a++)
                std::cout << tArrX[a] << " ";
            std::cout << "\n";
            for (Int_t a = 0; a < 2*nHitsEvent; a++)
                std::cout << tArrY[a] << " ";
            std::cout << "\n";
            for (Int_t a = 0; a < 2*nHitsEvent; a++)
                std::cout << tArrP[a] << " ";
            std::cout << "\n";
            for (Int_t a = 0; a < 2*nHitsEvent; a++)
                std::cout << tArrB[a] << " ";
            std::cout << "\n";
    }
    */

    // Now we can analyze the hits in this event

    if (fTofdHisto)
    {
        for (Int_t a = 0; a < 2 * nHitsEvent; a++)
        { // loop over all hits
            eventstore++;
            fhQ[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrQ[a]);        // charge per plane
            fhQvsEvent[((Int_t)tArrP[a]) - 1]->Fill(fnEvents, tArrQ[a]); // charge vs event #
            if (fTofdTotPos)
            {
                fhxy[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrYT[a]); // xy of plane
            }
            else
            {
                fhxy[((Int_t)tArrP[a]) - 1]->Fill(tArrB[a], tArrY[a]); // xy of plane
            }
        }
    }

    // select events with feasible times
    Double_t time0;
    for (Int_t ihit = 0; ihit < 2 * nHitsEvent;)
    { // loop over all hits in this event
        LOG(warn) << "\nSet new coincidence window: " << tArrP[ihit] << " " << tArrB[ihit] << " " << tArrT[ihit] << " "
                  << tArrQ[ihit];
        time0 = tArrT[ihit];            // time of first hit in coincidence window
        Double_t charge0 = tArrQ[ihit]; // charge of first hit in coincidence window
        Double_t plane0 = tArrP[ihit];  // plane of first hit in coincidence window
        Int_t hitscoinc = 0;
        Int_t nAverage = 0;
        Int_t nAverage12 = 0;
        Int_t nAverage34 = 0;
        Int_t nAverage14[100] = { 0 };
        Double_t sumQ12 = 0;
        Double_t sumQ34 = 0;
        Double_t sumQ14[100] = { 0. };

        while (tArrT[ihit] < time0 + hit_coinc)
        { // check if in coincidence window
            incoincidence++;
            /*
            std::cout<<"Used up hits in this coincidence window:\n";
            for(Int_t a=0; a<2*nHitsEvent; a++)
                std::cout << tArrP[a] << " ";
            std::cout << "\n";
            for(Int_t a=0; a<2*nHitsEvent; a++)
                std::cout << tArrB[a] << " ";
            std::cout << "\n";
            for(Int_t a=0; a<2*nHitsEvent; a++)
                std::cout << tArrU[a] << " ";
            std::cout << "\n";
            */

            hitscoinc++; // number of hits in coincidence

            if (fTofdHisto)
            {
                if (tArrP[ihit] == plane0 && charge0 != tArrQ[ihit])
                {
                    fhQM[(Int_t)tArrP[ihit] - 1]->Fill(charge0, tArrQ[ihit]);
                }
            }

            LOG(debug) << "Hit in coincidence window: " << tArrP[ihit] << " " << tArrB[ihit] << " " << tArrT[ihit]
                       << " " << tArrQ[ihit];

            // try to average plane 1&2
            for (Int_t i = 1; i < hitscoinc; i++)
            { // loop over hits in coincidence
                if (tArrP[ihit] <= 2 && tArrP[ihit - i] <= 2)
                {
                    // cout << i << "  " << tArrP[ihit] << "  " << tArrB[ihit] << " "
                    // << tArrP[ihit-i] << "  " << tArrB[ihit-i] << endl;
                    if (hitscoinc > 0 && (Int_t)(tArrP[ihit] - tArrP[ihit - i]) != 0)
                    { // check if planes differ
                        /// find overlapping virtualbars && similar charge in both planes? && bar wasn't used for other
                        /// average?
                        if (tArrB[ihit - i] == tArrB[ihit] && abs(tArrQ[ihit - i] - tArrQ[ihit]) < maxChargeDiff &&
                            tArrU[ihit] == false && tArrU[ihit - i] == false)
                        {
                            inaverage12++;
                            // cout << "Try to average " << tArrQ[ihit] << " plane: " << tArrP[ihit] << " Bar: " <<
                            // tArrB[ihit]
                            //             << " and " << tArrQ[ihit - i] << " plane: " << tArrP[ihit - i] << " Bar: "
                            //             << tArrB[ihit - i] << endl;

                            nAverage12++; // number of averaged hits in this coincidence window
                            sumQ12 += (tArrQ[ihit] + tArrQ[ihit - i]) / 2.; // average charges and add to sum

                            tArrU[ihit] = tArrU[ihit - i] = true; // set involved bars as used

                            if ((ihit - i) % 2 != 0)
                                tArrU[ihit - (i + 1)] =
                                    true; // set the associated virtual bars of the used bars as used
                            else
                                tArrU[ihit - (i - 1)] = true;
                            if ((ihit - i) % 2 != 0)
                                tArrU[ihit + 1] = true;
                            else
                                tArrU[ihit - 1] = true;

                            if (fTofdHisto)
                            {
                                fhCharge->Fill((tArrQ[ihit] + tArrQ[ihit - i]) / 2.); // Fill charge histogram

                                if (((Int_t)tArrP[ihit]) == 2)
                                {
                                    fhQvsQ[((Int_t)tArrP[ihit]) - 1][((Int_t)tArrB[ihit]) - 1]->Fill(
                                        tArrQ[ihit - i], tArrQ[ihit]); // Fill charge vs charge histogram
                                    fhTdiffvsQ[((Int_t)tArrP[ihit]) - 2][((Int_t)tArrB[ihit]) - 1]->Fill(
                                        tArrT[ihit] - tArrT[ihit - i],
                                        (tArrQ[ihit] + tArrQ[ihit - i]) / 2.); // Fill tdiff planes histogram
                                }
                                else
                                {
                                    fhQvsQ[((Int_t)tArrP[ihit]) - 1][((Int_t)tArrB[ihit]) - 1]->Fill(tArrQ[ihit],
                                                                                                     tArrQ[ihit - i]);
                                    if ((Int_t)tArrP[ihit] == 1)
                                        fhTdiffvsQ[((Int_t)tArrP[ihit]) - 1][((Int_t)tArrB[ihit]) - 1]->Fill(
                                            -(tArrT[ihit] - tArrT[ihit - i]),
                                            (tArrQ[ihit] + tArrQ[ihit - i]) / 2.); // Fill tdiff planes histogram
                                }

                                if ((tArrQ[ihit] || tArrQ[ihit - i] > 7.5) && (tArrQ[ihit] || tArrQ[ihit - i] < 8.5))
                                    fhxy12->Fill((tArrX[ihit] + tArrX[ihit - i]) / 2.,
                                                 (tArrY[ihit] + tArrY[ihit - i]) / 2.); // Fill average xy histogram
                                fhxy12tot->Fill((tArrX[ihit] + tArrX[ihit - i]) / 2.,
                                                (tArrYT[ihit] + tArrYT[ihit - i]) / 2.); // Fill average xy histogram
                            }

                            // store average
                            if (fTofdTotPos)
                            {
                                new ((*fHitItems)[fNofHitItems++])
                                    R3BTofdHitData((tArrT[ihit] + tArrT[ihit - i]) / 2.,
                                                   (tArrX[ihit] + tArrX[ihit - i]) / 2.,
                                                   (tArrYT[ihit] + tArrYT[ihit - i]) / 2.,
                                                   (tArrQ[ihit] + tArrQ[ihit - i]) / 2.,
                                                   abs(tArrT[ihit] - tArrT[ihit - i]),
                                                   (tArrQ[ihit] + tArrQ[ihit - i]) / 2.,
                                                   12);
                            }
                            else
                            {
                                new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData((tArrT[ihit] + tArrT[ihit - i]) / 2.,
                                                                                  (tArrX[ihit] + tArrX[ihit - i]) / 2.,
                                                                                  (tArrY[ihit] + tArrY[ihit - i]) / 2.,
                                                                                  (tArrQ[ihit] + tArrQ[ihit - i]) / 2.,
                                                                                  abs(tArrT[ihit] - tArrT[ihit - i]),
                                                                                  (tArrQ[ihit] + tArrQ[ihit - i]) / 2.,
                                                                                  12);
                            }
                        }

                        // std::cout<<"Used up averaged hits in this coincidence window:\n";
                        // for(Int_t a=0; a<2*nHitsEvent; a++)
                        //    std::cout << tArrU[a] << " ";
                        // std::cout << "\n";
                    }
                }
            }

            // try to average plane 3&4
            for (Int_t i = 1; i < hitscoinc; i++)
            { // loop over hits in coincidence
                if (tArrP[ihit] >= 3 && tArrP[ihit - i] >= 3)
                {
                    // std::cout<<i<<" "<<tArrP[ihit]<<" "<<tArrB[ihit]<<" "<<tArrP[ihit-i]<<" "<<tArrB[ihit-i]<<"\n";
                    if (hitscoinc > 0 && (Int_t)(tArrP[ihit] - tArrP[ihit - i]) != 0)
                    { // check if planes differ
                        /// find overlapping virtualbars         similar charge in both planes?               bar wasn't
                        /// used for other average?
                        if (tArrB[ihit - i] == tArrB[ihit] && abs(tArrQ[ihit - i] - tArrQ[ihit]) < maxChargeDiff &&
                            tArrU[ihit] == false && tArrU[ihit - i] == false)
                        {
                            inaverage34++;
                            LOG(warn) << "Try to average " << tArrQ[ihit] << " " << tArrP[ihit] << " " << tArrB[ihit]
                                      << " and " << tArrQ[ihit - i] << " " << tArrP[ihit - i] << " " << tArrB[ihit - i];

                            nAverage34++; // number of averaged hits in this coincidence window
                            sumQ34 += (tArrQ[ihit] + tArrQ[ihit - i]) / 2.; // average charges and add to sum

                            tArrU[ihit] = tArrU[ihit - i] = true; // set involved bars as used

                            if ((ihit - i) % 2 != 0)
                                tArrU[ihit - (i + 1)] =
                                    true; // set the associated virtual bars of the used bars as used
                            else
                                tArrU[ihit - (i - 1)] = true;
                            if ((ihit - i) % 2 != 0)
                                tArrU[ihit + 1] = true;
                            else
                                tArrU[ihit - 1] = true;

                            if (fTofdHisto)
                            {
                                fhCharge->Fill((tArrQ[ihit] + tArrQ[ihit - i]) / 2.); // Fill charge histogram

                                if (((Int_t)tArrP[ihit]) == fNofPlanes)
                                { /// TODO: maybe get first plane here?
                                    fhQvsQ[((Int_t)tArrP[ihit]) - 1][((Int_t)tArrB[ihit]) - 1]->Fill(
                                        tArrQ[ihit - i], tArrQ[ihit]); // Fill charge vs charge histogram
                                    fhTdiffvsQ[((Int_t)tArrP[ihit]) - 2][((Int_t)tArrB[ihit]) - 1]->Fill(
                                        tArrT[ihit] - tArrT[ihit - i],
                                        (tArrQ[ihit] + tArrQ[ihit - i]) / 2.); // Fill tdiff planes histogram
                                }
                                else
                                {
                                    fhQvsQ[((Int_t)tArrP[ihit]) - 1][((Int_t)tArrB[ihit]) - 1]->Fill(tArrQ[ihit],
                                                                                                     tArrQ[ihit - i]);
                                    if ((Int_t)tArrP[ihit] == 3)
                                        fhTdiffvsQ[((Int_t)tArrP[ihit]) - 1][((Int_t)tArrB[ihit]) - 1]->Fill(
                                            -(tArrT[ihit] - tArrT[ihit - i]),
                                            (tArrQ[ihit] + tArrQ[ihit - i]) / 2.); // Fill tdiff planes histogram
                                }

                                fhxy34->Fill((tArrX[ihit] + tArrX[ihit - i]) / 2.,
                                             (tArrY[ihit] + tArrY[ihit - i]) / 2.); // Fill average xy histogram
                                fhxy34tot->Fill((tArrX[ihit] + tArrX[ihit - i]) / 2.,
                                                (tArrYT[ihit] + tArrYT[ihit - i]) / 2.); // Fill average xy histogram
                            }

                            // store average
                            if (fTofdTotPos)
                            {
                                new ((*fHitItems)[fNofHitItems++])
                                    R3BTofdHitData((tArrT[ihit] + tArrT[ihit - i]) / 2.,
                                                   (tArrX[ihit] + tArrX[ihit - i]) / 2.,
                                                   (tArrYT[ihit] + tArrYT[ihit - i]) / 2.,
                                                   (tArrQ[ihit] + tArrQ[ihit - i]) / 2.,
                                                   abs(tArrT[ihit] - tArrT[ihit - i]),
                                                   (tArrQ[ihit] + tArrQ[ihit - i]) / 2.,
                                                   34);
                            }
                            else
                            {
                                new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData((tArrT[ihit] + tArrT[ihit - i]) / 2.,
                                                                                  (tArrX[ihit] + tArrX[ihit - i]) / 2.,
                                                                                  (tArrY[ihit] + tArrY[ihit - i]) / 2.,
                                                                                  (tArrQ[ihit] + tArrQ[ihit - i]) / 2.,
                                                                                  abs(tArrT[ihit] - tArrT[ihit - i]),
                                                                                  (tArrQ[ihit] + tArrQ[ihit - i]) / 2.,
                                                                                  34);
                            }
                        }

                        // std::cout<<"Used up averaged hits in this coincidence window:\n";
                        // for(Int_t a=0; a<2*nHitsEvent; a++)
                        //    std::cout << tArrU[a] << " ";
                        // std::cout << "\n";
                    }
                }
            }

            // try to average all 4 planes
            if (hitscoinc > 0)
            { // loop over hits in coincidence

                Int_t bar = tArrB[ihit];
                if (bar > 27 && bar < 41)
                {
                    // cout << "ihit: " << ihit << " plane: " << tArrP[ihit] << " bar: " << tArrB[ihit]
                    // 	<< "  Q: " << tArrQ[ihit] << endl;

                    // cout << "so far: " << nAverage14[bar] << " Q: " << sumQ14[bar] / nAverage14[bar] << endl;
                    if (abs(tArrQ[ihit] - sumQ14[bar] / nAverage14[bar]) < maxChargeDiff || nAverage14[bar] == 0)
                    {

                        // cout << "Average " << nAverage14[bar] << " Try to average "  << sumQ14[bar] / nAverage14[bar]
                        // << " with " << tArrQ[ihit] << " plane: " << tArrP[ihit] << " Bar: " << tArrB[ihit] << endl;

                        nAverage14[bar]++;          // number of averaged hits in this coincidence window
                        sumQ14[bar] += tArrQ[ihit]; // average charges and add to sum
                    }
                }
                if (fTofdHisto && nAverage14[bar] == 4)
                {
                    fhAverageCharge->Fill(bar, sumQ14[bar] / nAverage14[bar]);
                }
            }

            ihit++;
            if (ihit >= 2 * nHitsEvent)
                break;
        }
        nAverage = nAverage12 + nAverage34;
        if (fTofdHisto)
        {
            if (nAverage > 0)
            {
                // std::cout<<nAverage12<<"/"<<nAverage34<<" Events in coincidence window averaged\nCombined Charge"
                // <<sumQ12<<"/"<<sumQ34<<"\n";
                if (nAverage12 > 0)
                    fhMvsQ[0]->Fill(sumQ12, nAverage12); // Fill histogram number of averaged hits vs summed up charge
                if (nAverage34 > 0)
                    fhMvsQ[2]->Fill(sumQ34, nAverage34); // Fill histogram number of averaged hits vs summed up charge
            }
        }
    }
    for (Int_t hit = 0; hit < 2 * nHitsEvent; hit++)
    { // loop over not averaged hits
        if (tArrU[hit] == false)
        {
            LOG(debug) << "Single Hit for Plane " << tArrP[hit] << " " << tArrB[hit];
            tArrU[hit] = tArrU[hit + 1] = true;
            // store single hits only seen in planes
            singlehit++;
            if (fTofdTotPos)
            {
                new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(tArrT[hit],
                                                                  (tArrX[hit] + tArrX[hit + 1]) / 2.,
                                                                  tArrYT[hit],
                                                                  tArrQ[hit],
                                                                  -1.,
                                                                  tArrQ[hit],
                                                                  tArrP[hit]);
            }
            else
            {
                new ((*fHitItems)[fNofHitItems++]) R3BTofdHitData(tArrT[hit],
                                                                  (tArrX[hit] + tArrX[hit + 1]) / 2.,
                                                                  tArrY[hit],
                                                                  tArrQ[hit],
                                                                  -1.,
                                                                  tArrQ[hit],
                                                                  tArrP[hit]);
            }
            hit++;
        }
    }

    // std::cout<<"Used up hits in this event:\n";
    for (Int_t a = 0; a < 2 * nHitsEvent; a++)
    {
        // std::cout << tArrU[a] << " ";
        if (tArrU[a] != true)
            LOG(fatal);
    }
    // std::cout << "\n";

    LOG(debug) << "------------------------------------------------------\n";

    fnEvents++;
}

void R3BTofdCal2Hit::CreateHistograms(Int_t iPlane, Int_t iBar)
{
    Double_t max_charge = 80.;
    // create histograms if not already existing

    /*
    if (NULL == fhTsync[iPlane - 1])
    {
        char strName[255];
        sprintf(strName, "Time_Sync_Plane_%d", iPlane);
        fhTsync[iPlane - 1] = new TH2F(strName, "", 50,0,50,10000, -10, 90.);
        fhTsync[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTsync[iPlane - 1]->GetYaxis()->SetTitle("ToF in ns");
    }
    */
    /*
    if (NULL == fhTof[iPlane - 1])
    {
        char strName[255];
        sprintf(strName, "ToF_Plane_%d", iPlane);
        fhTof[iPlane - 1] = new TH2F(strName, "", 50,0,50,40000, -200, 200.);
        fhTof[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
        fhTof[iPlane - 1]->GetYaxis()->SetTitle("ToF in ns");
    }
    */
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

    if (NULL == fhQvsTof[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_ToF_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsTof[iPlane - 1][iBar - 1] = new TH2F(strName, "", 1000, 0., max_charge, 1000, -10, 40);
        fhQvsTof[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToF in ns");
        fhQvsTof[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Charge");
    }

    if (NULL == fhTvsTof[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "T_vs_ToF_Plane_%d_Bar_%d", iPlane, iBar);
        fhTvsTof[iPlane - 1][iBar - 1] = new TH2F(strName, "", 625, -25, 25, 1000, -10, 40);
        fhTvsTof[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToF in ns");
        fhTvsTof[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("T1-T2 in ns");
    }

    if (NULL == fhToTvsTofw[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "ToT_vs_ToF_Plane_%d_Bar_%d_w", iPlane, iBar);
        fhToTvsTofw[iPlane - 1][iBar - 1] = new TH2F(strName, "", 1000, 0., 200, 1000, -10, 40);
        fhToTvsTofw[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToT in ns");
        fhToTvsTofw[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToF in ns");
    }

    if (NULL == fhQvsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 200, -100, 100, 500, 0., max_charge);
        fhQvsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Position in cm");
    }

    if (NULL == fhQ[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "Q_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Q Plane %d ", iPlane);
        fhQ[iPlane - 1] = new TH2F(strName1, strName2, 90, 0, 90, max_charge * 10, 0., max_charge);
        fhQ[iPlane - 1]->GetYaxis()->SetTitle("Charge");
        fhQ[iPlane - 1]->GetXaxis()->SetTitle("Paddle number");
    }
    /*
    if (NULL == fhxy[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "xy_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "xy of Plane %d ", iPlane);
        fhxy[iPlane - 1] = new TH2F(strName1, strName2, 160, -80, 80, 400, -100., 100.);
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
        fhxy[iPlane - 1] = new TH2F(strName1, strName2, 90, 0, 90, 400, -100., 100.);
        fhxy[iPlane - 1]->GetYaxis()->SetTitle("y-position in cm");
        fhxy[iPlane - 1]->GetXaxis()->SetTitle("Bar #");
    }
    if (NULL == fhQvsEvent[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "QvsEvent_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Charge vs Event # Plane %d ", iPlane);
        fhQvsEvent[iPlane - 1] = new TH2F(strName1, strName2, 2e2, 0, 2e9, max_charge * 10, 0., max_charge); // 2e4 2e9
        fhQvsEvent[iPlane - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsEvent[iPlane - 1]->GetXaxis()->SetTitle("Event #");
    }
    // Multiplicity
    if (NULL == fhQM[iPlane - 1])
    {
        char strName1[255];
        sprintf(strName1, "QvsQt0_Plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Q vs Q_time0 Plane %d ", iPlane);
        fhQM[iPlane - 1] =
            new TH2F(strName1, strName2, max_charge * 10, 0., max_charge, max_charge * 10, 0., max_charge);
        fhQM[iPlane - 1]->GetYaxis()->SetTitle("Charge particle i");
        fhQM[iPlane - 1]->GetXaxis()->SetTitle("Charge first particle");
    }

    if (iPlane == 1 || iPlane == 3)
    {
        if (NULL == fhMvsQ[iPlane - 1])
        {
            char strName1[255];
            sprintf(strName1, "QvsHits_Plane_%d", iPlane);
            char strName2[255];
            sprintf(strName2, "Q vs Hit # Plane %d ", iPlane);
            char strName3[255];
            sprintf(strName3, "Hits in planes %d %d in coincidence window", iPlane, iPlane + 1);
            fhMvsQ[iPlane - 1] = new TH2F(strName1, strName2, max_charge * 10, 0., max_charge, 20, 0., 20);
            fhMvsQ[iPlane - 1]->GetYaxis()->SetTitle(strName3);
            fhMvsQ[iPlane - 1]->GetXaxis()->SetTitle("#sum Charge");
        }
    }

    if (iPlane == 1 || iPlane == 3)
    {
        if (NULL == fhTdiffvsQ[iPlane - 1][2 * iBar - 2])
        {
            char strName[255];
            sprintf(strName, "Tdiff_Plane_%dand%d_Bar_%dvsQ", iPlane, iPlane + 1, iBar * 2 - 1);
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 2] = new TH2F(strName, "", 1000, -10, 10, 1200, 0., 60.);
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 2]->GetYaxis()->SetTitle("charge");
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 2]->GetXaxis()->SetTitle("dt in ns");
        }
        if (NULL == fhTdiffvsQ[iPlane - 1][2 * iBar - 3])
        {
            char strName[255];
            sprintf(strName, "Tdiff_Plane_%dand%d_Bar_%dvsQ", iPlane, iPlane + 1, iBar * 2 - 2);
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 3] = new TH2F(strName, "", 1000, -10, 10, 1200, 0., 60.);
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 3]->GetYaxis()->SetTitle("charge");
            fhTdiffvsQ[iPlane - 1][iBar * 2 - 3]->GetXaxis()->SetTitle("dt in ns");
        }
    }

    if (NULL == fhQvsQ[iPlane - 1][iBar * 2 - 2])
    {
        char strName[255];
        sprintf(strName, "Q_vs_Q_Plane_%d_Bar_%d", iPlane, iBar * 2 - 1);
        fhQvsQ[iPlane - 1][iBar * 2 - 2] = new TH2F(strName, "", 1000, 0, max_charge, 1000, 0., max_charge);
        char strNamex[255];
        char strNamey[255];
        if (iPlane % 2 == 0)
        {
            sprintf(strNamex, "Charge plane %d", iPlane);
            sprintf(strNamey, "Charge plane %d", iPlane + 1);
        }
        else
        {
            sprintf(strNamex, "Charge plane %d", iPlane + 1);
            sprintf(strNamey, "Charge plane %d", iPlane);
        }
        fhQvsQ[iPlane - 1][iBar * 2 - 2]->GetYaxis()->SetTitle(strNamey);
        fhQvsQ[iPlane - 1][iBar * 2 - 2]->GetXaxis()->SetTitle(strNamex);
    }

    if (iPlane == 1 || iPlane == 3)
    {
        if (NULL == fhQvsQ[iPlane - 1][iBar * 2 - 3])
        {
            char strName[255];
            sprintf(strName, "Q_vs_Q_Plane_%d_Bar_%d", iPlane, iBar * 2 - 2);
            fhQvsQ[iPlane - 1][iBar * 2 - 3] = new TH2F(strName, "", 1000, 0, max_charge, 1000, 0., max_charge);
            char strNamex[255];
            char strNamey[255];
            if (iPlane % 2 == 0)
            {
                sprintf(strNamex, "Charge plane %d", iPlane);
                sprintf(strNamey, "Charge plane %d", iPlane + 1);
            }
            else
            {
                sprintf(strNamex, "Charge plane %d", iPlane + 1);
                sprintf(strNamey, "Charge plane %d", iPlane);
            }
            fhQvsQ[iPlane - 1][iBar * 2 - 3]->GetYaxis()->SetTitle(strNamey);
            fhQvsQ[iPlane - 1][iBar * 2 - 3]->GetXaxis()->SetTitle(strNamex);
        }
    }

    if (iPlane == 2 || iPlane == 4)
    {
        if (NULL == fhQvsQ[iPlane - 1][iBar * 2 - 1])
        {
            char strName[255];
            sprintf(strName, "Q_vs_Q_Plane_%d_Bar_%d", iPlane, iBar * 2);
            fhQvsQ[iPlane - 1][iBar * 2 - 1] = new TH2F(strName, "", 1000, 0, max_charge, 1000, 0., max_charge);
            char strNamex[255];
            char strNamey[255];
            if (iPlane % 2 == 0)
            {
                sprintf(strNamex, "Charge plane %d", iPlane);
                sprintf(strNamey, "Charge plane %d", iPlane + 1);
            }
            else
            {
                sprintf(strNamex, "Charge plane %d", iPlane + 1);
                sprintf(strNamey, "Charge plane %d", iPlane);
            }
            fhQvsQ[iPlane - 1][iBar * 2 - 1]->GetYaxis()->SetTitle(strNamey);
            fhQvsQ[iPlane - 1][iBar * 2 - 1]->GetXaxis()->SetTitle(strNamex);
        }
    }

    if (NULL == fhxy12)
    {
        char strName[255];
        sprintf(strName, "xy_of_ToFD_plane_12");
        fhxy12 = new TH2F(strName, "", 200, -100, 100, 200, -100., 100.);
        fhxy12->GetYaxis()->SetTitle("y-position in cm");
        fhxy12->GetXaxis()->SetTitle("x-position in cm");
    }

    if (NULL == fhxy12tot)
    {
        char strName[255];
        sprintf(strName, "xyToT_of_ToFD_plane_12");
        fhxy12tot = new TH2F(strName, "", 200, -100, 100, 200, -100., 100.);
        fhxy12tot->GetYaxis()->SetTitle("y-position in cm");
        fhxy12tot->GetXaxis()->SetTitle("x-position in cm");
    }

    if (NULL == fhxy34)
    {
        char strName[255];
        sprintf(strName, "xy_of_ToFD_plane_34");
        fhxy34 = new TH2F(strName, "", 200, -100, 100, 200, -100., 100.);
        fhxy34->GetYaxis()->SetTitle("y-position in cm");
        fhxy34->GetXaxis()->SetTitle("x-position in cm");
    }

    if (NULL == fhxy34tot)
    {
        char strName[255];
        sprintf(strName, "xyToT_of_ToFD_plane_34");
        fhxy34tot = new TH2F(strName, "", 200, -100, 100, 200, -100., 100.);
        fhxy34tot->GetYaxis()->SetTitle("y-position in cm");
        fhxy34tot->GetXaxis()->SetTitle("x-position in cm");
    }

    if (NULL == fhCharge)
    {
        char strName[255];
        sprintf(strName, "Charge_of_TofD");
        fhCharge = new TH1F(strName, "Charge of ToFD", 1000, 0., max_charge);
        fhCharge->GetYaxis()->SetTitle("Counts");
        fhCharge->GetXaxis()->SetTitle("Charge");
    }

    if (NULL == fhAverageCharge)
    {
        char strName[255];
        sprintf(strName, "Average_Charge");
        fhAverageCharge = new TH2F(strName, "", 90, 0, 90, 1200, 0., 60.);
        fhAverageCharge->GetYaxis()->SetTitle("Q");
        fhAverageCharge->GetXaxis()->SetTitle("Bar #");
    }

    /*
    if (NULL == fhChargevsTof)
    {
        char strName[255];
        sprintf(strName, "Charge_vs_Tof_of_TofD");
        fhChargevsTof = new TH2F(strName, "",  2000, -10., 40.,1000, 0, 100);
        fhChargevsTof->GetXaxis()->SetTitle("ToF in ns");
        fhChargevsTof->GetYaxis()->SetTitle("Charge");
    }
    */
    /*
    if (NULL == fhChargevsPos)
    {
        char strName[255];
        sprintf(strName, "Charge_vs_Pos_of_TofD");
        fhChargevsPos = new TH2F(strName, "", 100, 0, 100, 1000, 0., 100.);
        fhChargevsPos->GetYaxis()->SetTitle("Charge");
        fhChargevsPos->GetXaxis()->SetTitle("Bar number");
    }
    */
    /*
    if (NULL == fhQp12)
    {
        char strName[255];
        sprintf(strName, "Charge_vs_Pos_p12");
        fhQp12 = new TH2F(strName, "", 100, 0, 100, 5000, 0., max_charge);
        fhQp12->GetYaxis()->SetTitle("Average charge of plane 1 and 2");
        fhQp12->GetXaxis()->SetTitle("Bar number");
    }
    */
    /*
    if (NULL == fhQp34)
    {
        char strName[255];
        sprintf(strName, "Charge_vs_Pos_p34");
        fhQp34 = new TH2F(strName, "", 100, 0, 100, 1000, 0., max_charge);
        fhQp34->GetYaxis()->SetTitle("Average charge of plane 3 and 4");
        fhQp34->GetXaxis()->SetTitle("Bar number");
    }
    */
}
void R3BTofdCal2Hit::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
        fNofHitItems = 0;
    }
    if (fCalItemsLos)
    {
        fCalItemsLos->Clear();
    }
    if (fHitItemsLos)
    {
        fHitItemsLos->Clear();
    }
}

void R3BTofdCal2Hit::FinishTask()
{
    if (fhLosXYP)
        fhLosXYP->Write();
    // if (fhChargeLosTofD) fhChargeLosTofD->Write();
    if (fh_los_pos)
        fh_los_pos->Write();
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
            if (fhQM[i])
                fhQM[i]->Write();
            if (fhMvsQ[i])
                fhMvsQ[i]->Write();
            // if (fhTof[i]) fhTof[i]->Write();
            if (fhTdiff[i])
                fhTdiff[i]->Write();
            // if (fhTsync[i]) fhTsync[i]->Write();
            for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
            {

                // control histogram time particles
                if (fhQvsPos[i][j])
                    fhQvsPos[i][j]->Write();
                if (fhTdiffvsQ[i][2 * j])
                    fhTdiffvsQ[i][2 * j]->Write();
                if (fhTdiffvsQ[i][2 * j + 1])
                    fhTdiffvsQ[i][2 * j + 1]->Write();
                if (fhQvsQ[i][2 * j])
                    fhQvsQ[i][2 * j]->Write();
                if (fhQvsQ[i][2 * j + 1])
                    fhQvsQ[i][2 * j + 1]->Write();
                if (fhQvsTof[i][j])
                    fhQvsTof[i][j]->Write();
                if (fhTvsTof[i][j])
                    fhTvsTof[i][j]->Write();
                if (fhToTvsTofw[i][j])
                    fhToTvsTofw[i][j]->Write();
            }
        }
        if (fhxy12)
            fhxy12->Write();
        if (fhxy12tot)
            fhxy12tot->Write();
        if (fhxy34)
            fhxy34->Write();
        if (fhxy34tot)
            fhxy34tot->Write();
        if (fhCharge)
            fhCharge->Write();
        if (fhAverageCharge)
            fhAverageCharge->Write();
        // if (fhChargevsTof) fhChargevsTof->Write();
        // if (fhChargevsPos) fhChargevsPos->Write();
        // if (fhQp12) fhQp12->Write();
        // if (fhQp34) fhQp34->Write();
    }
    std::cout << "\n\nSome statistics:\n"
              << "Total number of events in tree  " << maxevent << "\n"
              << "Max Event analyzed              " << fnEvents + wrongtrigger + wrongtpat << "\n"
              << "Events in LOS                   " << countloshit << "\n"
              << "Events skipped due to trigger   " << wrongtrigger << "\n"
              << "Events skipped due to tpat      " << wrongtpat << "\n"
              << "Events with correct header&tpat " << headertpat << "\n"
              << "Events without ToFd hits        " << events_wo_tofd_hits << "\n"
              << "Events in cal level             " << events_in_cal_level << "\n"
              << "Hits in bar coincidence         " << inbarcoincidence << "\n"
              << "Number of resets                " << countreset << "\n"
              << "Hits before reset               " << hitsbeforereset << "\n"
              << "Bars with multihits             " << bars_with_multihit << "\n"
              << "Multihits                       " << multihit << "\n"
              << "Events stored                   " << eventstore / 2 << " <-> " << inbarcoincidence - hitsbeforereset
              << " = " << inbarcoincidence << " - " << hitsbeforereset
              << " (Events in bar coincidence - Hits before reset)\n"
              << "Events in coincidence window    " << incoincidence / 2 << "\n"
              << "Events in average p1&2          " << inaverage12 * 2 << "\n"
              << "Events in average p3&4          " << inaverage34 * 2 << "\n"
              << "Events in single planes         " << singlehit << "\n"
              << "Good events in total            " << eventstore / 2 << " <-> "
              << inaverage12 * 2 + inaverage34 * 2 + singlehit << " = " << inaverage12 * 2 << " + " << inaverage34 * 2
              << " + " << singlehit << "\n";
}

Double_t R3BTofdCal2Hit::betaCorr(Double_t delta)
{
    //    Double_t corr=-3.*delta;  //correction for Ag

    Double_t corr = -1. * delta; // correction for 12C
    corr = 0.;
    return corr;
}
/* old method
Double_t R3BTofdCal2Hit::walk(Double_t q)
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
Double_t R3BTofdCal2Hit::walk(Double_t Q,
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
Double_t R3BTofdCal2Hit::saturation(Double_t x)
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

Double_t* R3BTofdCal2Hit::insertX(Int_t n, Double_t arr[], Double_t x, Int_t pos)
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

ClassImp(R3BTofdCal2Hit)
