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
// -----                  R3BPtofMapped2TCal                -----
// -----          Created Feb 9th 2017 by R.Plag          -----
// ------------------------------------------------------------

/* March 2017
 * Reconstructing data into PaddleItems. Each paddle holds two times,
 * (for PM1 and PM2). This allows to plot
 * the time differences e.g. via evt->Draw(...) interactively (aka without
 * looping over all channels) which is crucial for a quick check of the
 * detector status during the experiment.
 *
 * This code was copied from Ptof. Differences:
 * - how many planes?
 *
 */

#include "R3BPtofMapped2Cal.h"

#include "R3BEventHeader.h"
#include "R3BPaddleTamexMappedData.h"
#include "R3BTCalEngine.h"
#include "R3BTCalPar.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include "TClonesArray.h"
#include "TMath.h"

#define IS_NAN(x) TMath::IsNaN(x)

R3BPtofMapped2Cal::R3BPtofMapped2Cal()
    : FairTask("PtofTcal", 1)
    , fMappedItems(NULL)
    , fCalItems(new TClonesArray("R3BPaddleCalData"))
    , fNofCalItems(0)
    , fNofTcalPars(0)
    , fTcalPar(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{
}

R3BPtofMapped2Cal::R3BPtofMapped2Cal(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fMappedItems(NULL)
    , fCalItems(new TClonesArray("R3BPaddleCalData"))
    , fNofCalItems(0)
    , fNofTcalPars(0)
    , fTcalPar(NULL)
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
{
}

R3BPtofMapped2Cal::~R3BPtofMapped2Cal()
{
    if (fCalItems)
    {
        delete fCalItems;
        fCalItems = NULL;
        fNofCalItems = 0;
    }
}

InitStatus R3BPtofMapped2Cal::Init()
{
    if (fTcalPar == NULL)
    {
        LOG(error) << "Have no TCal parameter container";
        return kFATAL;
    }

    fNofTcalPars = fTcalPar->GetNumModulePar();
    if (fNofTcalPars == 0)
    {
        LOG(error) << "No TCal parameters in container PtofTCalPar";
        return kFATAL;
    }

    LOG(info) << "R3BPtofMapped2Cal::Init : read " << fNofTcalPars << " modules";

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    // get access to Mapped data
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("PtofMapped"));
    if (NULL == fMappedItems)
        LOG(fatal) << "Branch PtofMapped not found";

    // request storage of Cal data in output tree
    mgr->Register("PtofCal", "Land", fCalItems, kTRUE);

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BPtofMapped2Cal::SetParContainers()
{
    fTcalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("PtofTCalPar"));
    if (!fTcalPar)
        LOG(error) << "Could not get access to PtofTCalPar-Container.";
}

InitStatus R3BPtofMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BPtofMapped2Cal::Exec(Option_t* option)
{
    Int_t nHits = fMappedItems->GetEntriesFast();

    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        R3BPaddleTamexMappedData* mapped = dynamic_cast<R3BPaddleTamexMappedData*>(fMappedItems->At(ihit));

        Int_t iPlane = mapped->GetPlaneId(); // 1..n; no need to check range
        Int_t iBar = mapped->GetBarId();     // 1..n
        //		LOG(error) << "R3BPtofMapped2Cal::Exec : Plane: " << 	iPlane << ", Bar: " << iBar ;

        R3BPaddleCalData* cal = new ((*fCalItems)[fNofCalItems++]) R3BPaddleCalData(iPlane, iBar);

        // convert times to ns
        for (int tube = 0; tube < 2; tube++)     // PM1 and PM2
            for (int edge = 0; edge < 2; edge++) // loop over leading and trailing edge
            {
                // check if there is indeed data for this tube and edge:

                // fetch calib params:
                R3BTCalModulePar* par = fTcalPar->GetModuleParAt(iPlane, iBar, tube * 2 + edge + 1); // 1..4
                if (!par)
                {
                    LOG(info) << "R3BPtofMapped2Cal::Exec : Tcal par not found, Plane: " << iPlane << ", Bar: " << iBar
                              << ", Tube: " << (tube + 1);
                    continue;
                }

                // Convert TDC to [ns] ...
                if (mapped->GetFineTime(tube, edge) == -1)
                {
                    continue;
                }
                Double_t time_ns = par->GetTimeVFTX(mapped->GetFineTime(tube, edge));
                if (time_ns < 0. || time_ns > fClockFreq)
                {
                    LOG(error) << "R3BPtofMapped2Cal::Exec : bad time in ns: Plane= " << iPlane << ", bar= " << iBar
                               << ",tube= " << (tube + 1) << ",edge= " << edge
                               << ", time in channels = " << mapped->GetFineTime(tube, edge)
                               << ", time in ns = " << time_ns;
                    continue;
                }

                // ... and add clock time
                time_ns = fClockFreq - time_ns + mapped->GetCoarseTime(tube, edge) * fClockFreq;
                if (iBar == 40)
                {
                    LOG(error) << "Test: Plane: " << iPlane << "  Bar: " << iBar << "  tube= " << tube
                               << "  edge= " << edge << "  time in ns = " << time_ns;
                }
                cal->SetTime(tube, edge, time_ns);
            }
    }
}

void R3BPtofMapped2Cal::FinishEvent()
{
    if (fCalItems)
    {
        fCalItems->Clear();
        fNofCalItems = 0;
    }
}

void R3BPtofMapped2Cal::FinishTask() {}

ClassImp(R3BPtofMapped2Cal)
