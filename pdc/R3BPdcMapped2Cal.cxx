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

#include <FairRootManager.h>
#include <assert.h>

#include "R3BPdcMapped2Cal.h"

#include "TClonesArray.h"
#include "TMath.h"

#include "FairLogger.h"
#include "FairRuntimeDb.h"

#include "R3BPdcCalData.h"
#include "R3BPdcMappedData.h"
#include "R3BTCalEngine.h"

#include "R3BEventHeader.h"

extern "C"
{
// #include "ext_data_client.h"
#include "ext_h101_pdc.h"
}
#define LENGTH(x) (sizeof x / sizeof *x)
#define IS_NAN(x) TMath::IsNaN(x)

R3BPdcMapped2Cal::R3BPdcMapped2Cal()
    : FairTask("R3BPdcMapped2Cal", 1)
    , fMappedItems(nullptr)
    , fCalItems(new TClonesArray("R3BPdcCalData"))
    , fCalTriggerItems(new TClonesArray("R3BPdcCalData"))
    , fTcalPar(0)
    , fNofTcalPars(0)
    , fTrigger(-1)
    , fClockFreq(1000. / CTDC_16_CLOCK_MHZ)
{
}

R3BPdcMapped2Cal::R3BPdcMapped2Cal(Int_t iVerbose)
    : FairTask("R3BPdcMapped2Cal", iVerbose)
    , fMappedItems(nullptr)
    , fCalItems(new TClonesArray("R3BPdcCalData"))
    , fCalTriggerItems(new TClonesArray("R3BPdcCalData"))
    , fTcalPar(0)
    , fNofTcalPars(0)
    , fTrigger(-1)
    , fClockFreq(1000. / CTDC_16_CLOCK_MHZ)
{
}

R3BPdcMapped2Cal::~R3BPdcMapped2Cal()
{
    delete fCalItems;
    delete fCalTriggerItems;
}

InitStatus R3BPdcMapped2Cal::Init()
{
    fNofTcalPars = fTcalPar->GetNumModulePar();
    if (fNofTcalPars == 0)
    {
        LOG(error) << "There are no TCal parameters in container PdcTCalPar";
        return kFATAL;
    }

    LOG(info) << "R3BPdcMapped2Cal::Init : read " << fNofTcalPars << " modules";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.
    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    // get access to Mapped data
    fMappedItems = dynamic_cast<TClonesArray*>(mgr->GetObject("PdcMapped"));
    if (NULL == fMappedItems)
        LOG(fatal) << "Branch PdcMapped not found";

    // request storage of Cal data in output tree
    mgr->Register("PdcCal", "Land", fCalItems, kTRUE);
    mgr->Register("PdcTriggerCal", "Land", fCalTriggerItems, kTRUE);

    return kSUCCESS;
}

// Note that the container may still be empty at this point.
void R3BPdcMapped2Cal::SetParContainers()
{
    fTcalPar = dynamic_cast<R3BTCalPar*>(FairRuntimeDb::instance()->getContainer("PdcTCalPar"));
    if (!fTcalPar)
    {
        LOG(error) << "Could not get access to PdcTCalPar-Container.";
        fNofTcalPars = 0;
    }
}

InitStatus R3BPdcMapped2Cal::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BPdcMapped2Cal::Exec(Option_t* option)
{
    auto mapped_num = fMappedItems->GetEntriesFast();
    LOG(debug) << "R3BPdcMapped2Cal::Exec:fMappedItems=" << fMappedItems->GetName() << '.';
    for (auto i = 0; i < mapped_num; i++)
    {
        auto mapped = dynamic_cast<R3BPdcMappedData*>(fMappedItems->At(i));
        assert(mapped);

        auto wire = mapped->GetWireId();
        LOG(debug) << " R3BPdcMapped2Cal::Exec:wire=" << wire
                   << ":Edge=" << (mapped->GetEdgeId() == 1 ? "Leading" : "Trailing") << '.';

        // Fetch tcal parameters.
        R3BTCalModulePar* par;
        par = fTcalPar->GetModuleParAt(mapped->GetPlaneId(), mapped->GetWireId(), mapped->GetEdgeId());
        if (!par)
        {
            LOG(warn) << "R3BPdcMapped2Cal::Exec (" << fName << "): Wire=" << wire << ": TCal par not found.";
            continue;
        }
        Double_t time_ns = -1;

        // Calibrate fine time.
        auto fine_raw = mapped->GetTimeFine();
        if (-1 == fine_raw)
        {
            // TODO: Is this really ok?
            continue;
        }
        auto fine_ns = par->GetTimeClockTDC(fine_raw);
        LOG(debug) << " R3BPdcMapped2Cal::Exec: Fine raw=" << fine_raw << " -> ns=" << fine_ns << '.';

        if (fine_ns < 0. || fine_ns >= fClockFreq)
        {
            LOG(error) << "R3BPdcMapped2Cal::Exec (" << fName << "): Wire=" << wire
                       << ": Bad CTDC fine time (raw=" << fine_raw << ",ns=" << fine_ns << ").";
            continue;
        }

        time_ns = mapped->GetTimeCoarse() * fClockFreq - fine_ns;

        // LOG(debug) << " R3BPdcMapped2Cal::Exec (" << fName << "): wire=" << wire
        //               << ": Time=" << time_ns << "ns.";

        // cout << "mapped2cal plane: " << mapped->GetPlaneId() << " Wire: " << mapped->GetWireId()
        //     << " Edge: " << mapped->GetEdgeId() << " Time: " << time_ns << endl;

        if (LENGTH(EXT_STR_h101_PDC_onion::PDC_P) + 1 == mapped->GetPlaneId())
        {
            new ((*fCalTriggerItems)[fCalTriggerItems->GetEntriesFast()])
                R3BPdcCalData(mapped->GetPlaneId(), mapped->GetWireId(), mapped->GetEdgeId(), time_ns);
        }
        else
        {
            new ((*fCalItems)[fCalItems->GetEntriesFast()])
                R3BPdcCalData(mapped->GetPlaneId(), mapped->GetWireId(), mapped->GetEdgeId(), time_ns);
        }
    }
    fnEvents++;
}

void R3BPdcMapped2Cal::FinishEvent()
{
    if (fCalItems)
    {
        fCalItems->Clear();
    }
    if (fCalTriggerItems)
    {
        fCalTriggerItems->Clear();
    }
}

void R3BPdcMapped2Cal::FinishTask() {}

ClassImp(R3BPdcMapped2Cal)
