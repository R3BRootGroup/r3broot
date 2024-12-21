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

#include "R3BTofDCal2HitPar.h"
#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BTofDHitModulePar.h"
#include "R3BTofDHitPar.h"
#include "R3BTofDMappingPar.h"
#include "R3BTofdCalData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRtdbRun.h"
#include "FairRunIdGenerator.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TProfile.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <stdlib.h>

#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

namespace
{
    double c_range_ns = 2048. * 5.;
    double c_bar_coincidence_ns = 20.; // nanoseconds.
} // namespace

R3BTofDCal2HitPar::R3BTofDCal2HitPar()
    : R3BTofDCal2HitPar("R3BTofDCal2HitPar", 1)
{
}

R3BTofDCal2HitPar::R3BTofDCal2HitPar(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        fh_tofd_TotPm[i] = nullptr;
        fhTdiff[i] = nullptr;
        fhTsync[i] = nullptr;
        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
        {
            fhLogTot1vsLogTot2[i][j] = nullptr;
            fhSqrtQvsPosToT[i][j] = nullptr;
            fhQvsPos[i][j] = nullptr;
            fhTot1vsPos[i][j] = nullptr;
            fhTot2vsPos[i][j] = nullptr;
            fh1_tofsync[i][j] = nullptr;
        }
    }
}

R3BTofDCal2HitPar::~R3BTofDCal2HitPar()
{
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        if (fh_tofd_TotPm[i])
            delete fh_tofd_TotPm[i];
        if (fhTdiff[i])
            delete fhTdiff[i];
        if (fhTsync[i])
            delete fhTsync[i];
        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
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
            if (fh1_tofsync[i][j])
                delete fh1_tofsync[i][j];
        }
    }
    if (fHitPar)
    {
        delete fHitPar;
    }
}

InitStatus R3BTofDCal2HitPar::Init()
{
    R3BLOG(info, "");
    FairRootManager* rm = FairRootManager::Instance();
    if (!rm)
    {
        R3BLOG(fatal, "FairRootManager not found");
        return kFATAL;
    }

    fHeader = dynamic_cast<R3BEventHeader*>(rm->GetObject("EventHeader."));
    R3BLOG_IF(fatal, nullptr == fHeader, "EventHeader. not found");

    fCalData = dynamic_cast<TClonesArray*>(rm->GetObject("TofdCal"));
    R3BLOG_IF(fatal, nullptr == fCalData, "TofdCal not found");

    fCalTriggerItems = dynamic_cast<TClonesArray*>(rm->GetObject("TofdTriggerCal"));
    R3BLOG_IF(fatal, nullptr == fCalTriggerItems, "TofdTriggerCal not found");

    fHitPar = dynamic_cast<R3BTofDHitPar*>(FairRuntimeDb::instance()->getContainer("tofdHitPar"));
    if (!fHitPar)
    {
        R3BLOG(error, "Could not get access to tofdHitPar container");
        return kFATAL;
    }

    for (Int_t i = 1; i <= fNofPlanes; i++)
        for (Int_t j = 1; j <= fPaddlesPerPlane; j++)
            CreateHistograms(i, j);

    // Definition of a time stich object to correlate times coming from different systems
    fTimeStitch = std::make_unique<R3BCoarseTimeStitch>();

    return kSUCCESS;
}

void R3BTofDCal2HitPar::SetParContainers()
{
    fMapPar = dynamic_cast<R3BTofDMappingPar*>(FairRuntimeDb::instance()->getContainer("tofdMappingPar"));
    R3BLOG_IF(warn, !fMapPar, "Could not get access to tofdMappingPar container");
}

namespace
{
    uint64_t n1, n2;
};

void R3BTofDCal2HitPar::Exec(Option_t* /*option*/)
{
    // test for requested trigger (if possible)
    if ((fTrigger >= 0) && (fHeader) && (fHeader->GetTrigger() != fTrigger))
        return;

    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit = fTpat - 1;
    if (fTpat_bit >= 0)
    {
        Int_t itpat = fHeader->GetTpat();
        Int_t tpatvalue = (itpat & (1 << fTpat_bit)) >> fTpat_bit;
        if ((fHeader) && (tpatvalue == 0))
            return;
    }

    // ToFD detector
    std::vector<std::vector<UInt_t>> multihits(fNofPlanes, std::vector<UInt_t>(fPaddlesPerPlane));
    for (Int_t i = 0; i < fNofPlanes; ++i)
        for (Int_t j = 0; j < fPaddlesPerPlane; ++j)
            multihits[i][j] = 0;

    UInt_t nHits = fCalData->GetEntriesFast();
    // Organize cals into bars.
    struct Entry
    {
        std::vector<R3BTofdCalData*> top;
        std::vector<R3BTofdCalData*> bot;
    };
    std::map<size_t, Entry> bar_map;
    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto* hit = dynamic_cast<R3BTofdCalData*>(fCalData->At(ihit));
        size_t idx = (hit->GetDetectorId() - 1) * fPaddlesPerPlane + (hit->GetBarId() - 1);
        auto ret = bar_map.insert(std::pair<size_t, Entry>(idx, Entry()));
        auto& vec = 1 == hit->GetSideId() ? ret.first->second.bot : ret.first->second.top;
        vec.push_back(hit);
    }

    static bool s_was_trig_missing = false;
    auto trig_num = fCalTriggerItems->GetEntriesFast();
    // Find coincident PMT hits.
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
            if (top_trig_i < trig_num && bot_trig_i < trig_num)
            {
                auto top_trig = dynamic_cast<R3BTofdCalData const*>(fCalTriggerItems->At(top_trig_i));
                auto bot_trig = dynamic_cast<R3BTofdCalData const*>(fCalTriggerItems->At(bot_trig_i));
                top_trig_ns = top_trig->GetTimeLeading_ns();
                bot_trig_ns = bot_trig->GetTimeLeading_ns();
                ++n1;
            }
            else
            {
                if (!s_was_trig_missing)
                {
                    R3BLOG(error, "Missing trigger information!");
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
            if (std::abs(dt_mod) < c_bar_coincidence_ns)
            {
                // Hit
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
                // register multi hits
                multihits[iPlane - 1][iBar - 1] += 1;

                if (fParameter == 1)
                {
                    // calculate tdiff
                    auto tdiff = (bot_ns) - (top_ns);

                    // fill control histograms
                    fhLogTot1vsLogTot2[iPlane - 1][iBar - 1]->Fill(TMath::Log(top_tot), TMath::Log(bot_tot));
                    fh_tofd_TotPm[iPlane - 1]->Fill(iBar, top_tot);
                    fh_tofd_TotPm[iPlane - 1]->Fill(-iBar - 1, bot_tot);

                    // Time differences of one paddle; offset  histo
                    fhTdiff[iPlane - 1]->Fill(iBar, tdiff);

                    auto posToT = TMath::Log(top_tot / bot_tot);
                    fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->Fill(posToT, sqrt(top_tot * bot_tot));

                    // Time of hit
                    auto THit = (top_ns + bot_ns) / 2.;
                    fhTsync[iPlane - 1]->Fill(iBar, THit);

                    // Tof with respect LOS detector
                    auto tof = fTimeStitch->GetTime((top_ns + bot_ns) / 2. - fHeader->GetTStart(), "tamex", "vftx");
                    fh1_tofsync[iPlane - 1][iBar - 1]->Fill(tof);
                    // std::cout << "top" << top_ns << " bot"<<bot_ns << " start" << header->GetTStart() << std::endl;
                }

                // prepare offset and sync calculation
                if (fTofdQ == 0.0 && fParameter > 1)
                {
                    LOG(debug) << "Fill histo for offset and sync calculation plane " << iPlane << " bar " << iBar;
                    auto par = fHitPar->GetModuleParAt(iPlane, iBar);
                    if (!par)
                    {
                        R3BLOG(error,
                               "Hit par not found, Plane: " << top->GetDetectorId() << ", Bar: " << top->GetBarId());
                        continue;
                    }
                    // calculate tdiff
                    auto tdiff = (bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2());

                    // fill control histograms
                    // fhTot1vsTot2[iPlane - 1][iBar - 1]->Fill(top_tot, bot_tot);
                    fhLogTot1vsLogTot2[iPlane - 1][iBar - 1]->Fill(TMath::Log(top_tot), TMath::Log(bot_tot));
                    fh_tofd_TotPm[iPlane - 1]->Fill(iBar, top_tot);
                    fh_tofd_TotPm[iPlane - 1]->Fill(-iBar - 1, bot_tot);

                    // Time differences of one paddle; offset  histo
                    fhTdiff[iPlane - 1]->Fill(iBar, tdiff);

                    // offset histo via ToT
                    /*if (fTofdY == 0.)
                    {
                        R3BLOG(debug, "Will prepare for offset and sync calculation");
                        posToT = TMath::Log(top_tot / bot_tot);
                    }
                    else
                    {*/
                    auto posToT = TMath::Log((top_tot + par->GetToTOffset2()) / (bot_tot + par->GetToTOffset1()));
                    //}
                    fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->Fill(posToT, sqrt(top_tot * bot_tot));
                    // fhTot1vsPos[iPlane - 1][iBar - 1]->Fill(posToT, bot_tot);
                    // fhTot2vsPos[iPlane - 1][iBar - 1]->Fill(posToT, top_tot);

                    // Time of hit
                    auto THit = (top_ns + bot_ns) / 2. - par->GetSync();
                    fhTsync[iPlane - 1]->Fill(iBar, THit);

                    // Tof with respect LOS detector
                    auto tof = fTimeStitch->GetTime((top_ns + bot_ns) / 2. - fHeader->GetTStart(), "tamex", "vftx");
                    fh1_tofsync[iPlane - 1][iBar - 1]->Fill(tof - par->GetTofSyncOffset());
                }
                else if (fTofdQ > 0 && fParameter > 1)
                {
                    // get parameter
                    auto para = fHitPar->GetModuleParAt(iPlane, iBar);
                    if (!para)
                    {
                        R3BLOG(error,
                               "Hit par not found, Plane: " << top->GetDetectorId() << ", Bar: " << top->GetBarId());
                        continue;
                    }

                    // calculate tdiff with offest
                    auto tdiff = (bot_ns + para->GetOffset1()) - (top_ns + para->GetOffset2());
                    /*
                                        // walk corrections
                                        if (para->GetPar1Walk() == 0. || para->GetPar2Walk() == 0. ||
                       para->GetPar3Walk() == 0. || para->GetPar4Walk() == 0. || para->GetPar5Walk() == 0.){
                                            R3BLOG(warn, "Walk correction not found!");
                                        }else{

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
                                        }*/

                    auto posToT =
                        para->GetLambda() * log((top_tot + para->GetToTOffset2()) / (bot_tot + para->GetToTOffset1()));

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
                if (!fTofdSmiley && fTofdQ > 0.1 && fParameter > 1)
                {
                    LOG(debug) << "Prepare histo for double exponential fit";
                    auto par = fHitPar->GetModuleParAt(iPlane, iBar);
                    if (!par)
                    {
                        R3BLOG(error,
                               "Hit par not found, Plane: " << top->GetDetectorId() << ", Bar: " << top->GetBarId());
                        continue;
                    }

                    // calculate y position
                    auto pos = ((bot_ns + par->GetOffset1()) - (top_ns + par->GetOffset2())) * par->GetVeff();

                    // fill fitting histograms and smiley histogram
                    fhTot1vsPos[iPlane - 1][iBar - 1]->Fill(pos, bot_tot);
                    fhTot2vsPos[iPlane - 1][iBar - 1]->Fill(pos, top_tot);
                }

                // prepare charge fit / quench correction
                if (fTofdZ == true && fParameter > 1)
                {
                    LOG(debug) << "Prepare histo for quenching correction";
                    // get parameter
                    auto par = fHitPar->GetModuleParAt(iPlane, iBar);
                    if (!par)
                    {
                        R3BLOG(error,
                               "Hit par not found, Plane: " << top->GetDetectorId() << ", Bar: " << top->GetBarId());
                        continue;
                    }

                    // Calculate y position
                    auto posToT =
                        par->GetLambda() * log((top_tot + par->GetToTOffset2()) / (bot_tot + par->GetToTOffset1()));

                    Double_t parq[4];
                    parq[0] = par->GetPar1a();
                    parq[1] = par->GetPar1b();
                    parq[2] = par->GetPar1c();
                    parq[3] = par->GetPar1d();

                    // Calculate charge
                    Double_t qb = 0.;
                    if (fTofdSmiley)
                    {
                        qb = TMath::Sqrt(top_tot * bot_tot); //
                                                             //(parq[0] + parq[1] * posToT + parq[2] * pow(posToT, 2) +
                                                             // parq[3] * pow(posToT, 3));
                        // qb = qb * fTofdQ; // theory says: dE ~ Z^2 but we see quenching -> just use linear and fit
                        // the rest
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
                             fTofdQ; // theory says: dE ~ Z^2 but we see quenching -> just use linear and fit the rest
                        q2 = q2 * fTofdQ;
                        qb = (q1 + q2) / 2.;
                    }

                    // Fill control histograms and Q vs Pos without multihits
                    if (multihits[iPlane - 1][iBar - 1] < 2 && (qb > 0.))
                    {
                        // std::cout<< " h: " << posToT << " , " << qb << std::endl;
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

void R3BTofDCal2HitPar::CreateHistograms(Int_t iPlane, Int_t iBar)
{
    Double_t max_charge = fMaxQ;
    if (NULL == fhTdiff[iPlane - 1])
    {
        char strName1[255];
        char strName2[255];
        sprintf(strName1, "Time_Diff_Plane_%d", iPlane);
        sprintf(strName2, "Time Diff Plane %d", iPlane);
        fhTdiff[iPlane - 1] = new TH2F(strName1, strName2, 50, 0, 50, 4000, -20., 20.);
        fhTdiff[iPlane - 1]->GetXaxis()->SetTitle("Bar ");
        fhTdiff[iPlane - 1]->GetYaxis()->SetTitle("Time difference (PM1 - PM2) in ns");
    }
    if (NULL == fhTsync[iPlane - 1])
    {
        char strName[255];
        sprintf(strName, "Time_Sync_Plane_%d", iPlane);
        fhTsync[iPlane - 1] = new TH2F(strName, strName, 45, 0, 45, 10000, -3500, 1000.);
        fhTsync[iPlane - 1]->GetXaxis()->SetTitle("Bar ");
        fhTsync[iPlane - 1]->GetYaxis()->SetTitle("THit in ns");
    }
    if (NULL == fh_tofd_TotPm[iPlane - 1])
    {
        char strName[255];
        sprintf(strName, "Tofd_ToT_plane_%d", iPlane);
        char strName2[255];
        sprintf(strName2, "Tofd ToT plane %d", iPlane);
        fh_tofd_TotPm[iPlane - 1] = new TH2F(strName, strName2, 90, -45, 45, 1500, 0., 500.);
        fh_tofd_TotPm[iPlane - 1]->GetXaxis()->SetTitle("Bar ");
        fh_tofd_TotPm[iPlane - 1]->GetYaxis()->SetTitle("ToT / ns");
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
            new TH2F(strName, "", 2000, -100, 100, max_charge * 4, 0., max_charge * 4);
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("sqrt(PM1*PM2)");
        fhSqrtQvsPosToT[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Position from ToT in cm");
    }
    if (NULL == fhQvsPos[iPlane - 1][iBar - 1])
    {
        char strName[255];
        sprintf(strName, "Q_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        fhQvsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 2000, -100, 100, max_charge / 2, 0., max_charge);
        fhQvsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("Charge");
        fhQvsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Position in cm");
    }

    if (NULL == fhTot1vsPos[iPlane - 1][iBar - 1] && !fTofdSmiley)
    {
        char strName[255];
        sprintf(strName, "Tot1_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        fhTot1vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 200, -100, 100, 400, 0., 200.);
        fhTot1vsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Pos in cm");
        fhTot1vsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToT of PM1 in ns");
    }
    if (NULL == fhTot2vsPos[iPlane - 1][iBar - 1] && !fTofdSmiley)
    {
        char strName[255];
        sprintf(strName, "Tot2_vs_Pos_Plane_%d_Bar_%d", iPlane, iBar);
        fhTot2vsPos[iPlane - 1][iBar - 1] = new TH2F(strName, "", 200, -100, 100, 400, 0., 200.);
        fhTot2vsPos[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("Pos in cm");
        fhTot2vsPos[iPlane - 1][iBar - 1]->GetYaxis()->SetTitle("ToT of PM2 in ns");
    }
    if (fh1_tofsync[iPlane - 1][iBar - 1] == NULL)
    {
        char strName[255];
        sprintf(strName, "tofdiff_plane_%d_bar_%d", iPlane, iBar);
        fh1_tofsync[iPlane - 1][iBar - 1] = new TH1F(strName, strName, 5000, -250, -100);
        fh1_tofsync[iPlane - 1][iBar - 1]->GetXaxis()->SetTitle("ToF [ns]");
    }
}

void R3BTofDCal2HitPar::calcOffset()
{
    TCanvas* cOffset = new TCanvas("cOffset", "cOffset", 10, 10, 1000, 900);
    cOffset->Divide(2, 2);
    R3BTofDHitModulePar* mpar;
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        if (fhTdiff[i])
        {
            LOG(warn) << "Found histo Time_Diff_Plane_" << i + 1;
            // auto* h = (TH2F*)fhTdiff[i]->Clone();
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                mpar = new R3BTofDHitModulePar();
                Double_t offset = 0.;
                cOffset->cd(i + 1);
                fhTdiff[i]->Draw("colz");
                auto* histo_py = fhTdiff[i]->ProjectionY("histo_py", j + 2, j + 2, "");
                histo_py->Rebin(4);
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t Max = histo_py->GetXaxis()->GetBinCenter(binmax);
                TF1* fgaus = new TF1("fgaus", "gaus(0)", Max - 0.3, Max + 0.3);
                histo_py->Fit("fgaus", "QR0");
                offset = fgaus->GetParameter(1); // histo_py->GetXaxis()->GetBinCenter(binmax);
                LOG(warn) << " Plane  " << i + 1 << " Bar " << j + 1 << " Offset  " << offset;
                mpar->SetPlane(i + 1);
                mpar->SetPaddle(j + 1);
                mpar->SetOffset1(-offset / 2.);
                mpar->SetOffset2(offset / 2.);
                fHitPar->AddModulePar(mpar);
            }
        }
    }
    fHitPar->setChanged();
}
void R3BTofDCal2HitPar::calcToTOffset(Double_t totLow, Double_t totHigh)
{
    TCanvas* cToTOffset = new TCanvas("cToTOffset", "cToTOffset", 10, 10, 1000, 900);
    cToTOffset->Divide(1, 2);
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
        {
            Double_t offset = 0.;
            R3BTofDHitModulePar* par = fHitPar->GetModuleParAt(i + 1, j + 1);
            if (fhSqrtQvsPosToT[i][j])
            {
                LOG(info) << "Found histo SqrtQ_vs_PosToT_Plane_" << i + 1 << "_Bar_" << j + 1;
                // auto* h = fhSqrtQvsPosToT[i][j]->Clone();
                cToTOffset->cd(1);
                fhSqrtQvsPosToT[i][j]->Draw("colz");
                auto histo_py = fhSqrtQvsPosToT[i][j]->ProjectionX("histo_py", totLow, totHigh, "");
                cToTOffset->cd(2);
                histo_py->Rebin(2);
                histo_py->Draw();
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t Max = histo_py->GetXaxis()->GetBinCenter(binmax);
                TF1* fgaus = new TF1(
                    "fgaus", "gaus(0)", Max - 0.2, Max + 0.2); // new TF1("fgaus", "gaus(0)", Max - 0.06, Max + 0.06);
                histo_py->Fit("fgaus", "QR0");
                offset = fgaus->GetParameter(1);
                fgaus->Draw("SAME");
                histo_py->SetAxisRange(Max - .5, Max + .5, "X");
                fhSqrtQvsPosToT[i][j]->SetAxisRange(Max - .5, Max + .5, "X");
                fhSqrtQvsPosToT[i][j]->SetAxisRange(totLow, totHigh, "Y");
                cToTOffset->Update();
                delete fgaus;
            }
            LOG(warn) << " Plane  " << i + 1 << " Bar " << j + 1 << " ToT Offset  " << offset << "\n";
            par->SetToTOffset1(sqrt(exp(offset)));
            par->SetToTOffset2(1. / sqrt(exp(offset)));
        }
    }
    fHitPar->setChanged();
}

void R3BTofDCal2HitPar::calcSync()
{
    TCanvas* cSync = new TCanvas("cSync", "cSync", 10, 10, 1000, 900);
    cSync->Divide(2, 2);
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        if (fhTsync[i])
        {
            LOG(info) << "Found histo Time_Sync_Plane_" << i + 1;
            // auto h = fhTsync[i]->Clone();
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                cSync->cd(i + 1);
                fhTsync[i]->Draw("colz");
                auto* histo_py = fhTsync[i]->ProjectionY("histo_py", j + 2, j + 2, "");
                R3BTofDHitModulePar* par = fHitPar->GetModuleParAt(i + 1, j + 1);
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t Max = histo_py->GetXaxis()->GetBinCenter(binmax);
                Double_t MaxEntry = histo_py->GetBinContent(binmax);
                TF1* fgaus = new TF1("fgaus", "gaus(0)", Max - 10., Max + 10.);
                fgaus->SetParameters(MaxEntry, Max, 20);
                histo_py->Fit("fgaus", "QR0");
                Double_t sync = fgaus->GetParameter(1); // histo_py->GetXaxis()->GetBinCenter(binmax);
                par->SetSync(sync);
                LOG(info) << " Plane  " << i + 1 << " Bar " << j + 1 << " Sync  " << sync;
            }
        }
    }
    fHitPar->setChanged();
}

void R3BTofDCal2HitPar::zcorr(TH2F* histo, Int_t min, Int_t max, Double_t* pars, Int_t index)
{
    if (histo->GetEntries() < 100)
    {
        R3BLOG(warn, "Nb of events below 100 for histo with index" << index);
        return;
    }

    Double_t par[3000] = { 0 };
    Int_t nPeaks = fNbZPeaks;
    Double_t x[3000] = { 0 };
    char strName[255];
    sprintf(strName, "canvas_%d", index);
    TCanvas* c1 = new TCanvas(strName, "", 100, 100, 800, 800);
    c1->Divide(1, 3);
    c1->cd(1);
    auto* h = dynamic_cast<TH2F*>(histo->Clone());
    h->Draw("colz");
    h->SetAxisRange(min, max, "Y");
    // Projection of charge axis
    auto* h1 = h->ProjectionY("p_y");
    c1->cd(2);
    h1->Draw();
    // Use TSpectrum to find the peak candidates
    TSpectrum* s = new TSpectrum(nPeaks);
    Int_t nfound = s->Search(h1, 10, "", 0.001);
    std::cout << "Found " << nfound << " candidate peaks to fit\n";

    if (nfound == 0)
    {
        delete s;
        delete c1;
        delete h;
        delete h1;
        return;
    }

    c1->Update();
    // Eliminate background peaks
    nPeaks = 0;
    Double_t* xpeaks = s->GetPositionX();
    for (Int_t p = 0; p <= nfound; p++)
    {
        Float_t xp = xpeaks[p];
        Int_t bin = h1->GetXaxis()->FindBin(xp);
        Float_t yp = h1->GetBinContent(bin);
        if (yp - TMath::Sqrt(yp) < 1.)
            continue;
        par[2 * nPeaks] = yp;
        par[2 * nPeaks + 1] = xp;
        nPeaks++;
    }

    Double_t peaks[nPeaks];
    for (Int_t i = 0; i < nPeaks; i++)
    {
        // printf("Found peak @ %f\n",xpeaks[i]);
        peaks[i] = par[2 * i + 1];
    }
    c1->Update();

    if (nPeaks < 2)
    {
        pars[0] = 0.;
        pars[1] = fTofdQ / peaks[0];
        pars[2] = 0.;

        delete s;
        delete c1;
        delete h;
        delete h1;
        return;
    }

    // select useful peaks
    sort(peaks, peaks + nPeaks);
    Double_t zpeaks[3000] = { 0 };
    string doagain = "y";
    Int_t nfp;
doagainfit:
    do
    {
        nfp = 0;
        for (Int_t i = 0; i < nPeaks; i++)
        {
            std::cout << "Peak @ " << peaks[i];
            Int_t z = 0;
            while ((std::cout << " corresponds to Z=") && !(std::cin >> z))
            {
                std::cout << "That's not a number;";
                std::cin.clear();
                std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            }
            if (z == 0)
                continue;
            x[nfp] = (Double_t)z;
            zpeaks[nfp] = peaks[i];
            // std::cout<<"z real " << x[nfp] << " z found " << zpeaks[nfp] <<"\n";
            nfp++;
        }
        std::cout << "Do again? (y/n) ";
        std::cin >> doagain;
    } while (doagain != "n");

    // fit charge axis
    std::cout << "Selected " << nfp << " useful peaks to fit\nStart fitting...\n";
    c1->cd(3)->Clear();
    c1->Update();
    c1->cd(3);
    auto gr1 = new TGraph(nfp, zpeaks, x);
    gr1->SetMarkerColor(4);
    gr1->SetMarkerSize(1.5);
    gr1->SetMarkerStyle(20);
    gr1->Draw("AP");
    // TF1* fitz = new TF1("fitz", "[0]*TMath::Power(x,[2])+[1]", min, max);
    if (fZfitType != "pol1" && fZfitType != "pol2")
    {
        R3BLOG(error, "Fit " << fZfitType << " is not allowed, use pol1 or pol2 ");
        return;
    }
    auto fitz = new TF1("fitz", fZfitType, min, max);
    fitz->SetLineColor(2);
    fitz->SetLineWidth(2);
    fitz->SetLineStyle(1);
    // fitz->SetParameters(1.5, 2., .1);
    gr1->Fit("fitz", "Q");
    fitz->Draw("lsame");
    c1->Update();
    // write parameters
    std::cout << "Is OK? (y/n) ";
    std::cin >> doagain;
    if (doagain == "n")
        goto doagainfit;

    int nbpars = 2;
    if (fZfitType == "pol2")
        nbpars = 3;

    for (Int_t j = 0; j < nbpars; j++)
    {
        pars[j] = fitz->GetParameter(j);
        // std::cout<<Form("par%i= ",j)<<pars[j]<<"\n";
    }
    delete s;
    delete h;
    delete h1;
    delete gr1;
    delete c1;
    delete fitz;
}

void R3BTofDCal2HitPar::calcVeff()
{
    TCanvas* cVeff = new TCanvas("cVeff", "cVeff", 10, 10, 1000, 900);
    cVeff->Divide(2, 2);
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
        {
            Double_t max = 0.;
            Double_t veff = 7.;
            if (fhTdiff[i])
            {
                LOG(info) << "Found histo Time_Diff_Plane_" << i + 1;
                auto par = fHitPar->GetModuleParAt(i + 1, j + 1);
                if (!par)
                {
                    LOG(warn) << "Hit par not found, Plane: " << i + 1 << ", Bar: " << j + 1;
                    continue;
                }
                // auto* h = (TH2F*)histofilename->Get(Form("Time_Diff_Plane_%i", i + 1))->Clone();
                cVeff->cd(i + 1);
                // h->Draw("colz");
                auto* histo_py = fhTdiff[i]->ProjectionY("histo_py", j + 2, j + 2, "");
                Int_t binmax = histo_py->GetMaximumBin();
                max = histo_py->GetXaxis()->GetBinCenter(binmax);
                Double_t maxEntry = histo_py->GetBinContent(binmax);
                auto* fgaus = new TF1("fgaus", "gaus(0)", max - 0.3, max + 0.3); /// TODO: find best range
                fgaus->SetParameters(maxEntry, max, 20);
                histo_py->Fit("fgaus", "QR0");
                Double_t offset1 = par->GetOffset1();
                Double_t offset2 = par->GetOffset2();
                max = fgaus->GetParameter(1) + offset1 - offset2; /// TODO: needs to be tested
                // max = max+offset1-offset2;
                veff = fTofdY / max; // effective speed of light in [cm/s]
                LOG(info) << "Plane  " << i + 1 << " Bar " << j + 1 << " offset  " << par->GetOffset1();
                LOG(info) << "Plane  " << i + 1 << " Bar " << j + 1 << " max  " << max;
                LOG(info) << "Plane  " << i + 1 << " Bar " << j + 1 << " veff  " << veff;
                par->SetVeff(veff);
            }
        }
    }
    fHitPar->setChanged();
}

void R3BTofDCal2HitPar::calcLambda(Double_t totLow, Double_t totHigh)
{
    TCanvas* cToTOffset = new TCanvas("cLambda", "cLambda", 10, 10, 1000, 900);
    cToTOffset->Divide(1, 2);
    for (Int_t i = 0; i < fNofPlanes; i++)
    {
        for (Int_t j = 0; j < fPaddlesPerPlane; j++)
        {
            Double_t offset = 0.;
            auto par = fHitPar->GetModuleParAt(i + 1, j + 1);
            if (fhSqrtQvsPosToT[i][j])
            {
                LOG(info) << "Found histo SqrtQ_vs_PosToT_Plane_" << i + 1 << "_Bar_" << j + 1;
                // auto* h = (TH2F*)histofilename->Get(Form("SqrtQ_vs_PosToT_Plane_%i_Bar_%i", i + 1, j + 1))->Clone();
                cToTOffset->cd(1);
                fhSqrtQvsPosToT[i][j]->Draw("colz");
                auto* histo_py = fhSqrtQvsPosToT[i][j]->ProjectionX("histo_py", totLow, totHigh, "");
                cToTOffset->cd(2);
                histo_py->Draw();
                Int_t binmax = histo_py->GetMaximumBin();
                Double_t Max = histo_py->GetXaxis()->GetBinCenter(binmax);
                TF1* fgaus = new TF1("fgaus", "gaus(0)", Max - 0.06, Max + 0.06);
                histo_py->Fit("fgaus", "QR0");
                offset = fgaus->GetParameter(1);
                fgaus->Draw("SAME");
                histo_py->SetAxisRange(Max - .5, Max + .5, "X");
                fhSqrtQvsPosToT[i][j]->SetAxisRange(Max - .5, Max + .5, "X");
                fhSqrtQvsPosToT[i][j]->SetAxisRange(totLow, totHigh, "Y");
                cToTOffset->Update();
                delete fgaus;
                delete histo_py;
            }
            else
                LOG(error) << "Missing histo plane " << i + 1 << " bar " << j + 1;
            Double_t lambda = fTofdY / offset;
            LOG(info) << " Plane  " << i + 1 << " Bar " << j + 1 << " ToT Offset  " << offset << " Lambda " << lambda
                      << "\n";
            par->SetLambda(lambda);
        }
    }
    fHitPar->setChanged();
}

void R3BTofDCal2HitPar::smiley(TH2F* histo, Double_t min, Double_t max, Double_t* para)
{
    // This fits the smiley: Sqrt(q1*q2) returns position dependent charge, we fit that via pol3 and try to correct
    Double_t y[1000], x[1000];
    Int_t n = 0;
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = 0;
    }
    TGraph* gr1 = new TGraph();
    TGraph* gr2 = new TGraph();
    TCanvas* cfit_smiley = new TCanvas("cfit_smiley", "fit smiley", 100, 100, 800, 800);
    cfit_smiley->Clear();
    cfit_smiley->Divide(1, 4);
    cfit_smiley->cd(1);
    TH2F* histo1 = dynamic_cast<TH2F*>(histo->Clone());
    histo1->Draw("colz");
    TH2F* histo2 = dynamic_cast<TH2F*>(histo->Clone());
    histo2->RebinX(50);
    histo2->GetYaxis()->SetRangeUser(fTofdTotLow, fTofdTotHigh);
    // histo2->SetAxisRange(fTofdTotLow,fTofdTotHigh,"Y");
    cfit_smiley->cd(2);
    histo2->Draw("colz");
    std::cout << "Searching for points to fit...\n";
    for (Int_t i = 1; i < histo2->GetNbinsX(); i++)
    {
        // std::cout<<"Bin "<<i<<" of "<<histo2->GetNbinsX()<<" with cut: "<<fTofdTotLow<<" < sqrt(q1*q2) <
        // "<<fTofdTotHigh<<"\n";
        cfit_smiley->cd(2);
        TLine* l = new TLine(
            histo2->GetXaxis()->GetBinCenter(i), fTofdTotLow, histo2->GetXaxis()->GetBinCenter(i), fTofdTotHigh);
        l->SetLineColor(kRed);
        l->SetLineWidth(2.);
        l->Draw();
        cfit_smiley->cd(3);
        auto* histo_py = histo2->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        // cfit_smiley->Update();
        x[n] = histo2->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);

        if ((x[n] < min || x[n] > max) || (y[n] < fTofdTotLow || y[n] > fTofdTotHigh))
        {
            delete histo_py;
            continue;
        }
        if (histo_py->GetMaximum() > 5)
        {
            n++;
            delete l;
        }
        delete histo_py;
    }
    gr1 = new TGraph(n, x, y);
    gr1->SetTitle("Points found for fitting; x position in cm; sqrt(tot1*tot2)");
    gr1->Draw("A*");
    std::cout << "Start fitting\n";
    TF1* f1 = new TF1("f1", "pol3", min, max);
    f1->SetLineColor(2);
    gr1->Fit("f1", "Q", "", min, max);
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f1->GetParameter(j);
        std::cout << "Parameter: " << para[j] << "\n";
    }
    // fit again but with more information and better cuts
    std::cout << "Fit again with more information\n";
    n = 0;
    cfit_smiley->cd(4);
    for (Int_t i = 1; i < histo2->GetNbinsX(); i++)
    {
        Double_t pos = histo2->GetXaxis()->GetBinCenter(i);
        Double_t ymean = f1->Eval(pos);
        histo2->SetAxisRange(ymean - 5., ymean + 5., "Y");
        histo2->Draw("colz");
        auto* histo_py = histo2->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        x[n] = histo2->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);
        if (histo_py->GetMaximum() > 2)
            n++;
        delete histo_py;
    }
    gr2 = new TGraph(n, x, y);
    gr2->SetTitle("More information;x position in cm;sqrt(q1*q2)");
    gr2->Draw("A*");
    f1->DrawCopy("SAME");
    TF1* f2 = new TF1("f2", "pol3", min, max);
    f2->SetParameters(para[0], para[1], para[2], para[3]);
    f2->SetLineColor(3);
    gr2->Fit("f2", "0Q", "", min, max);
    f2->Draw("SAME");
    std::cout << "Will write:\n";
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f2->GetParameter(j);
        std::cout << "Parameter: " << para[j] << "\n";
    }
    histo2->GetYaxis()->SetRangeUser(fTofdTotLow, fTofdTotHigh);
    auto legend = new TLegend(.9, 0.7, .99, 0.9);
    legend->AddEntry("f1", "First Fit", "l");
    legend->AddEntry("f2", "Second Fit", "l");
    legend->Draw();
    cfit_smiley->Update();
    // gPad->WaitPrimitive();
    // gSystem->Sleep(3000);
    delete histo1;
    delete histo2;
    delete gr1;
    delete gr2;
    delete f1;
    delete f2;
    delete cfit_smiley;
}

void R3BTofDCal2HitPar::doubleExp(TH2F* histo, Double_t min, Double_t max, Double_t* para)
{
    // This fits the exponential decay of the light in a paddle. The 2 PMTs are fit with the same function but one
    // side will deliver negative attenuation parameters and the other positive.
    Double_t y[1000], x[1000];
    Int_t n = 0;
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = 0;
    }
    TGraph* gr1 = new TGraph();
    TGraph* gr2 = new TGraph();
    TCanvas* cfit_exp = new TCanvas("cfit_exp", "fit exponential", 100, 100, 800, 800);
    cfit_exp->Clear();
    cfit_exp->Divide(1, 3);
    cfit_exp->cd(1);
    TH2F* histo1 = dynamic_cast<TH2F*>(histo->Clone());
    TH2F* histo2 = dynamic_cast<TH2F*>(histo->Clone());
    histo1->Draw("colz");
    cfit_exp->cd(2);
    for (Int_t i = 1; i < histo1->GetNbinsX() - 1; i++)
    {
        auto* histo_py = histo1->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        x[n] = histo1->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);
        if ((x[n] < -40. || x[n] > 40.) || y[n] < 50.)
        {
            delete histo_py;
            continue;
        }
        if (histo_py->GetMaximum() > 5)
            n++;
        delete histo_py;
    }
    gr1 = new TGraph(n, x, y);
    gr1->Draw("A*");
    TF1* f1 = new TF1("f1", "[0]*(exp(-[1]*(x+100.))+exp(-[2]*(x+100.)))+[3]", min, max);
    f1->SetParameters(520., 0.001, 17234, -485.);
    f1->SetLineColor(2);
    gr1->Fit("f1", "", "", min, max);
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f1->GetParameter(j);
        std::cout << "Parameter: " << para[j] << "\n";
    }
    // fit again but with more information and better cuts
    n = 0;
    cfit_exp->cd(3);
    for (Int_t i = 1; i < histo2->GetNbinsX(); i++)
    {
        Double_t pos = histo2->GetXaxis()->GetBinCenter(i);
        Double_t ymean = para[0] * (exp(-para[1] * (pos + 100.)) + exp(-para[2] * (pos + 100.))) + para[3];
        histo2->SetAxisRange(ymean - 5., ymean + 5., "Y");
        histo2->Draw("colz");
        auto* histo_py = histo2->ProjectionY("histo_py", i, i, "");
        histo_py->Draw();
        x[n] = histo2->GetXaxis()->GetBinCenter(i);
        Int_t binmax = histo_py->GetMaximumBin();
        y[n] = histo_py->GetXaxis()->GetBinCenter(binmax);
        if (histo_py->GetMaximum() > 2)
            n++;
        delete histo_py;
    }
    gr2 = new TGraph(n, x, y);
    gr2->Draw("A*");
    TF1* f2 = new TF1("f2", "[0]*(exp(-[1]*(x+100.))+exp(-[2]*(x+100.)))+[3]", min, max);
    f2->SetParameters(para[0], para[1], para[2], para[3]);
    f2->SetLineColor(2);
    gr2->Fit("f2", "", "", min, max);
    for (Int_t j = 0; j <= 3; j++)
    {
        para[j] = f2->GetParameter(j);
        std::cout << "Parameter: " << para[j] << "\n";
    }
    cfit_exp->Update();
    // gPad->WaitPrimitive();
    // gSystem->Sleep(3000);
    delete gr1;
    delete gr2;
    delete f1;
    delete f2;
}

void R3BTofDCal2HitPar::FinishTask()
{
    if (fParameter == 1)
    {
        // Determine time offset of the 2 PMTs of one paddle. This procedure
        // assumes a sweep run in the middle of the ToF wall horizontally.
        // Since all paddles are mounted vertically one can determine the offset.
        // Half of the offset is added to PM1 and half to PM2.
        LOG(info) << "Calling function calcOffset";
        calcOffset();
        // Determine ToT offset between top and bottom PMT
        LOG(info) << "Calling function calcToTOffset";
        calcToTOffset(fTofdTotLow, fTofdTotHigh);
        // Determine sync offset between paddles
        LOG(info) << "Calling function calcSync";
        calcSync();
        LOG(error) << "Call walk correction before next step!";

        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                auto par = fHitPar->GetModuleParAt(i + 1, j + 1);
                Int_t binmax = fh1_tofsync[i][j]->GetMaximumBin();
                auto tofsync = fh1_tofsync[i][j]->GetXaxis()->GetBinCenter(binmax);

                TF1* fgauss = new TF1("fgaus", "gaus(0)", tofsync - 0.25, tofsync + 0.25);
                fh1_tofsync[i][j]->Fit("fgaus", "QR");
                auto tof_offset = fgauss->GetParameter(1);

                par->SetTofSyncOffset(tof_offset - fMeanTof);
                LOG(info) << " Plane  " << i + 1 << " Bar " << j + 1 << " Tof-Sync  " << tof_offset;
            }
        }
    }
    else if (fParameter == 2)
    {
        // Determine effective speed of light in [cm/s] for each paddle
        LOG(info) << "Calling function";
        calcVeff();
        // Determine light attenuation lambda for each paddle
        LOG(info) << "Calling function calcLambda";
        calcLambda(fTofdTotLow, fTofdTotHigh);
    }
    else if (fParameter == 3)
    {
        // calculation of position dependend charge
        if (fTofdSmiley)
        {
            LOG(info) << "Calling function smiley";
            Double_t para2[4];
            for (Int_t i = 0; i < 4; i++)
                para2[i] = 0.;
            Double_t min2 = -40.; // -40 effective bar length
            Double_t max2 = 40.;  // 40 effective bar length = 80 cm
                                  // we will use 50 here for some fit safety margin
            for (Int_t i = 0; i < fNofPlanes; i++)
            {
                for (Int_t j = 0; j < fPaddlesPerPlane; j++)
                {
                    if (fhSqrtQvsPosToT[i][j])
                    {
                        LOG(info) << "Calling Plane " << i + 1 << " Bar " << j + 1;
                        auto par = fHitPar->GetModuleParAt(i + 1, j + 1);
                        smiley(fhSqrtQvsPosToT[i][j], min2, max2, para2);
                        par->SetPola(para2[0]);
                        par->SetPolb(para2[1]);
                        par->SetPolc(para2[2]);
                        par->SetPold(para2[3]);
                    }
                }
            }
            fHitPar->setChanged();
        }
        else
        {
            LOG(info) << "Calling function doubleExp";
            Double_t para[4];
            for (Int_t i = 0; i < 4; i++)
                para[i] = 0.;
            Double_t min = -40.; // effective bar length
            Double_t max = 40.;  // effective bar length = 80 cm

            for (Int_t i = 0; i < fNofPlanes; i++)
            {
                for (Int_t j = 0; j < fPaddlesPerPlane; j++)
                {
                    if (fhTot1vsPos[i][j])
                    {
                        auto par = fHitPar->GetModuleParAt(i + 1, j + 1);
                        doubleExp(fhTot1vsPos[i][j], min, max, para);
                        Double_t offset1 = par->GetOffset1();
                        Double_t offset2 = par->GetOffset2();
                        Double_t veff = par->GetVeff();
                        Double_t sync = par->GetSync();
                        par->SetPar1a(para[0]);
                        par->SetPar1b(para[1]);
                        par->SetPar1c(para[2]);
                        par->SetPar1d(para[3]);
                    }
                    if (fhTot2vsPos[i][j])
                    {
                        auto par = fHitPar->GetModuleParAt(i + 1, j + 1);
                        doubleExp(fhTot2vsPos[i][j], min, max, para);
                        Double_t offset1 = par->GetOffset1();
                        Double_t offset2 = par->GetOffset2();
                        Double_t veff = par->GetVeff();
                        Double_t sync = par->GetSync();
                        par->SetPar2a(para[0]);
                        par->SetPar2b(para[1]);
                        par->SetPar2c(para[2]);
                        par->SetPar2d(para[3]);
                    }
                }
            }
            fHitPar->setChanged();
        }
    }

    if (fParameter == 4)
    {
        // Z correction for each plane
        LOG(warn) << "Calling function zcorr";
        Double_t para[8];
        Double_t pars[3];
        pars[0] = 0.;
        pars[1] = 0.;
        pars[2] = 0.;
        Int_t min = 0, max = fMaxQ; // select range for peak search
        for (Int_t i = 0; i < fNofPlanes; i++)
        {
            for (Int_t j = 0; j < fPaddlesPerPlane; j++)
            {
                if (fhQvsPos[i][j])
                {
                    auto par = fHitPar->GetModuleParAt(i + 1, j + 1);
                    std::cout << "Calling Plane: " << i + 1 << " Bar " << j + 1 << "\n";
                    Int_t index = i * fPaddlesPerPlane + j;
                    zcorr(fhQvsPos[i][j], min, max, pars, index);
                    std::cout << "Write parameter: " << pars[0] << " " << pars[1] << " " << pars[2] << "\n";
                    par->SetPar1za(pars[0]);
                    par->SetPar1zb(pars[1]);
                    par->SetPar1zc(pars[2]);
                }
            }
        }
        fHitPar->setChanged();
    }

    for (Int_t i = 0; i < N_TOFD_HIT_PLANE_MAX; i++)
    {
        if (fh_tofd_TotPm[i])
            fh_tofd_TotPm[i]->Write(); // control histogram for ToT
        if (fhTsync[i])
            fhTsync[i]->Write(); // histogram for sync calculation
        if (fhTdiff[i])
            fhTdiff[i]->Write(); // histogram for offset and veff calculation
        for (Int_t j = 0; j < N_TOFD_HIT_PADDLE_MAX; j++)
        {
            if (fh1_tofsync[i][j])
                fh1_tofsync[i][j]->Write(); // histogram for ToF sync calculation
            if (fhLogTot1vsLogTot2[i][j])
                fhLogTot1vsLogTot2[i][j]->Write(); // control histogram Log(ToT) Pm1 vs Log(ToT) Pm2
            if (fhSqrtQvsPosToT[i][j])
                fhSqrtQvsPosToT[i][j]->Write(); // histogram for ToT offset calculation
            if (fhQvsPos[i][j])
                fhQvsPos[i][j]->Write(); // histogram for charge fit
            if (!fTofdSmiley)
            {
                if (fhTot1vsPos[i][j])
                    fhTot1vsPos[i][j]->Write(); // histogram for position dependence of charge 1
                if (fhTot2vsPos[i][j])
                    fhTot2vsPos[i][j]->Write(); // histogram for position dependence of charge 2
            }
            /*
            if (fhTot1vsTot2[i][j])
                fhTot1vsTot2[i][j]->Write(); // control histogram ToT Pm1 vs ToT Pm2
            */
        }
    }
}

Double_t R3BTofDCal2HitPar::walk(Double_t Q, Double_t par1, Double_t par2, Double_t par3, Double_t par4, Double_t par5)
{
    Double_t y = 0;
    y = -30.2 + par1 * TMath::Power(Q, par2) + par3 / Q + par4 * Q + par5 * Q * Q;
    return y;
}

Double_t R3BTofDCal2HitPar::saturation(Double_t x)
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

ClassImp(R3BTofDCal2HitPar)
