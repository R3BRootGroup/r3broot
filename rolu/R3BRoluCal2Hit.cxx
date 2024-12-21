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

// ----------------------------------------------------------------------------
// -----                      R3BRoluCal2Hit                             ------
// -----            Created July 2019 by A. Kelic-Heil                   ------
// -----    Modified July 2021 by A. Kelic-Heil to include 2 ROLUs       ------
// ----------------------------------------------------------------------------

#include "R3BRoluCal2Hit.h"
#include "FairLogger.h"
#include "R3BRoluCalData.h"
#include "R3BRoluHitData.h"
#include "R3BRoluMapped2Cal.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"
#include <FairRootManager.h>

#include "TClonesArray.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include <fstream>
#include <iomanip>
#include <iostream>

#include <boost/multi_array.hpp>
#include <cassert>
#include <cmath> // NAN

using namespace std;

template <class T>
static void init_array(T& array, double init)
{
    // boost::multi_array prefers cascaded for loops over index sets
    // (e.g. python's dict().keys())
    // thus we can not really use a range based for loop
    // and instead wrangle pointers like it is 1980 again.
    for (double* eee = array.origin(); eee < array.origin() + array.num_elements(); ++eee)
        *eee = init;
}

namespace
{
    double c_period = 2048 * 5;
} // namespace
#define IS_NAN(x) TMath::IsNaN(x)

R3BRoluCal2Hit::R3BRoluCal2Hit()
    : R3BRoluCal2Hit("RoluCal2Hit", 1)
{
}

R3BRoluCal2Hit::R3BRoluCal2Hit(const char* name, int iVerbose)
    : FairTask(name, iVerbose)
    , fHitItems(new TClonesArray("R3BRoluHitData"))
{
}

R3BRoluCal2Hit::~R3BRoluCal2Hit()
{
    for (int iDet = 0; iDet < fNofDetectors; iDet++)
    {
        if (fhQ_L[iDet])
            delete (fhQ_L[iDet]);
        if (fhQ_R[iDet])
            delete (fhQ_R[iDet]);
        if (fhQ_O[iDet])
            delete (fhQ_O[iDet]);
        if (fhQ_U[iDet])
            delete (fhQ_U[iDet]);
    }

    if (fHitItems)
    {
        delete fHitItems;
        fHitItems = NULL;
    }
}

InitStatus R3BRoluCal2Hit::Init()
{
    // get access to Cal data
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "R3BRoluCal2Hit::Init() FairRootManager not found";

    fCalItems = dynamic_cast<TClonesArray*>(mgr->GetObject("RoluCal"));
    if (NULL == fCalItems)
        LOG(fatal) << "R3BRoluCal2Hit::Init() Branch RoluCal not found";

    if (fSkipTrigger == false)
    {
        fCalTriggerItems = dynamic_cast<TClonesArray*>(mgr->GetObject("RoluTriggerCal"));
        if (NULL == fCalTriggerItems)
            LOG(warn) << "R3BRoluCal2Hit::Init() Branch RoluTriggerCal not found";
    }

    // request storage of Hit data in output tree
    mgr->Register("RoluHit", "RoluHitData", fHitItems, !fOnline);
    maxevent = mgr->CheckMaxEventNo();

    return kSUCCESS;
}

InitStatus R3BRoluCal2Hit::ReInit() { return kSUCCESS; }

/*
 * Calculate a single hit time for Rolu detector
 */
void R3BRoluCal2Hit::Exec(Option_t* option)
{
    if (fnEvents / 100000. == (int)fnEvents / 100000)
        std::cout << "\rEvents: " << fnEvents << " / " << maxevent << " (" << (int)(fnEvents * 100. / maxevent)
                  << " %) " << std::flush;

    // min,max,Nbins for ToT spectra
    double fhQmin = 0.;
    double fhQmax = 200; // 300.; //150
    int fhQbin = 4000;   // 0; //3000; //1500
    for (int iDet = 0; iDet < fNofDetectors; iDet++)
    {
        if (NULL == fhQ_R[iDet])
        {
            fhQ_R[iDet] = new TH1F(Form("TOT_R%d", iDet + 1), Form("TOTr%d", iDet + 1), fhQbin, fhQmin, fhQmax);
        }

        if (NULL == fhQ_O[iDet])
        {
            fhQ_O[iDet] = new TH1F(Form("TOT_O%d", iDet + 1), Form("TOTo%d", iDet + 1), fhQbin, fhQmin, fhQmax);
        }

        if (NULL == fhQ_L[iDet])
        {
            fhQ_L[iDet] = new TH1F(Form("TOT_L%d", iDet + 1), Form("TOTl%d", iDet + 1), fhQbin, fhQmin, fhQmax);
        }

        if (NULL == fhQ_U[iDet])
        {
            fhQ_U[iDet] = new TH1F(Form("TOT_U%d", iDet + 1), Form("TOTu%d", iDet + 1), fhQbin, fhQmin, fhQmax);
        }
    }

    int nTrig = 0;
    double lead_trig_ns = 0. / 0.;
    // Calibrate the trigger
    if (fSkipTrigger == false)
    {
        nTrig = fCalTriggerItems->GetEntriesFast();
        for (unsigned int iTrig = 0; iTrig < nTrig; ++iTrig)
        {
            auto cur_cal = dynamic_cast<R3BRoluCalData*>(fCalTriggerItems->At(iTrig));
            lead_trig_ns = cur_cal->GetTimeL_ns(0);
            // cout<<"Trigger: "<<lead_trig_ns<<endl;
        }
    }

    int nParts = fCalItems->GetEntriesFast();

    if (nParts < 1 || nTrig < 1)
        return;

    int iDet = 0;
    using A = boost::multi_array<double, 3>;
    auto dims = boost::extents[nParts][fNofDetectors][4];
    A timeRolu_L(dims);
    init_array(timeRolu_L, NAN);
    A timeRolu_T(dims);
    init_array(timeRolu_T, NAN);
    A totRolu(dims);
    init_array(totRolu, NAN);

    for (int iPart = 0; iPart < nParts; iPart++)
    {
        // Parts is the number of particle passing through detector in one event
        R3BRoluCalData* calItem = dynamic_cast<R3BRoluCalData*>(fCalItems->At(iPart));
        iDet = calItem->GetDetector();

        for (int iCha = 0; iCha < 4; iCha++)
        {

            if (!(IS_NAN(calItem->GetTimeL_ns(iCha))))
            { // TAMEX leading
                timeRolu_L[iPart][iDet - 1][iCha] = calItem->GetTimeL_ns(iCha);
            }
            if (!(IS_NAN(calItem->GetTimeT_ns(iCha))))
            { // TAMEX trailing
                timeRolu_T[iPart][iDet - 1][iCha] = calItem->GetTimeT_ns(iCha);
            }

            if (timeRolu_T[iPart][iDet - 1][iCha] > 0. && timeRolu_L[iPart][iDet - 1][iCha] > 0. &&
                !(IS_NAN(timeRolu_T[iPart][iDet - 1][iCha])) && !(IS_NAN(timeRolu_L[iPart][iDet - 1][iCha])))
            {
                while (timeRolu_T[iPart][iDet - 1][iCha] - timeRolu_L[iPart][iDet - 1][iCha] <= 0.)
                {
                    timeRolu_T[iPart][iDet - 1][iCha] = timeRolu_T[iPart][iDet - 1][iCha] + 2048. * fClockFreq;
                }

                totRolu[iPart][iDet - 1][iCha] = timeRolu_T[iPart][iDet - 1][iCha] - timeRolu_L[iPart][iDet - 1][iCha];
            }

            double time_to_trig =
                fmod(timeRolu_L[iPart][iDet - 1][iCha] - lead_trig_ns + c_period + c_period / 2, c_period) -
                c_period / 2;

            // cout<<"ROLU cal2hit: "<<iDet<<", "<<iCha+1<<"; "<<timeRolu_L[iPart][iDet - 1][iCha]<<",
            // "<<totRolu[iPart][iDet - 1][iCha]<<endl;
            new ((*fHitItems)[fHitItems->GetEntriesFast()])
                R3BRoluHitData(iDet, iCha + 1, time_to_trig, totRolu[iPart][iDet - 1][iCha]);
        }

        fhQ_R[iDet - 1]->Fill(totRolu[iPart][iDet - 1][0]);
        fhQ_O[iDet - 1]->Fill(totRolu[iPart][iDet - 1][1]);
        fhQ_L[iDet - 1]->Fill(totRolu[iPart][iDet - 1][2]);
        fhQ_U[iDet - 1]->Fill(totRolu[iPart][iDet - 1][3]);
    }

    fnEvents++;
}

void R3BRoluCal2Hit::FinishEvent()
{
    if (fHitItems)
    {
        fHitItems->Clear();
    }
}

void R3BRoluCal2Hit::FinishTask()
{
    for (int iDet = 0; iDet < fNofDetectors; iDet++)
    {
        if (fhQ_L[iDet])
            fhQ_L[iDet]->Write();
        if (fhQ_R[iDet])
            fhQ_R[iDet]->Write();
        if (fhQ_O[iDet])
            fhQ_O[iDet]->Write();
        if (fhQ_U[iDet])
            fhQ_U[iDet]->Write();
    }
}

ClassImp(R3BRoluCal2Hit);
ClassImp(R3BRoluCal2Hit);
