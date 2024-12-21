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

#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"

#include <iomanip>

#include "R3BRpcMapped2PreCal.h"

#include "R3BRpcMappedData.h"

#include "R3BRpcPreCalData.h"
double counter_bad = 0;
double counter_good = 0;

// R3BRpcMapped2PreCal: Constructor
R3BRpcMapped2PreCal::R3BRpcMapped2PreCal()
    : FairTask("R3B RPC Calibrator")
    , fTCalPar(NULL)
    , fOnline(kFALSE)
    , fFpgaCorrelationFile("")
    , fClockFreq(1. / VFTX_CLOCK_MHZ * 1000.)
    , fMappedDataCA(NULL)
    , fRpcPreCalDataCA(NULL)
{
}

R3BRpcMapped2PreCal::~R3BRpcMapped2PreCal()
{
    LOG(info) << "R3BRpcMapped2PreCal: Delete instance";
    delete fRpcPreCalDataCA;
}

InitStatus R3BRpcMapped2PreCal::Init()
{
    LOG(info) << "R3BRpcMapped2PreCal::Init()";

    // INPUT DATA
    FairRootManager* rootManager = FairRootManager::Instance();
    if (!rootManager)
    {
        LOG(fatal) << "R3BRpcMapped2PreCal::FairRootManager not found";
        return kFATAL;
    }

    fMappedDataCA = dynamic_cast<TClonesArray*>(rootManager->GetObject("R3BRpcMappedData"));
    if (!fMappedDataCA)
    {
        LOG(error) << "R3BRpcMapped2PreCalPar::Init() fMappedDataCA not found";
        return kFATAL;
    }

    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    if (!rtdb)
    {
        LOG(error) << "R3BRpcMapped2PreCal:: FairRuntimeDb not opened";
    }

    fTCalPar = dynamic_cast<R3BTCalPar*>(rtdb->getContainer("RpcTCalPar"));
    if (!fTCalPar)
    {
        LOG(error) << "R3BRpcMapped2PreCalPar::Init() Couldn't get handle on RpcTCalPar container";
    }
    else
    {
        LOG(info) << "R3BRpcMapped2PreCal:: RPCTCalPar container open";
    }

    // FPGA mapping
    int lines = 0;

    ifstream in;

    string chn_id;
    int side_lut;

    in.open(fFpgaCorrelationFile, ios::in);

    if (fFpgaCorrelationFile == "")
    {
        LOG(error) << "R3BRpcMapped2PreCalPar::Init() Couldn't get handle on fFpgaCorrelationFile";
    }

    while (in >> chn_id >> side_lut)
    {

        lines++;
        if (chn_id.find("NC") != string::npos || chn_id.find("NB") != string::npos)
        {
            continue;
        }
        int chn = stoi(chn_id);
        int fpga = (lines - 1) / 32.;
        lut[chn - 1][side_lut].push_back(fpga);
    }
    in.close();

    // OUTPUT DATA
    // Calibrated data
    fRpcPreCalDataCA = new TClonesArray("R3BRpcPreCalData", 50);

    rootManager->Register("R3BRpcPreCalData", "RPC Pre Cal", fRpcPreCalDataCA, !fOnline);

    return kSUCCESS;
}

InitStatus R3BRpcMapped2PreCal::ReInit() { return kSUCCESS; }

void R3BRpcMapped2PreCal::Exec(Option_t* option)
{
    // Reset entries in output arrays, local arrays
    Reset();

    if (!fTCalPar)
    {
        LOG(warn) << "R3BRpcMapped2PreCal::Parameter container not found";
    }

    Entry Ref_vec[10];

    Double_t times_Ref = 0. / 0.;
    Double_t times_Pmt = 0. / 0.;
    Double_t times_Strip = 0. / 0.;

    // loop over the 3 mapped structures

    // loop over ref data
    Int_t nHits = fMappedDataCA->GetEntries();
    for (Int_t i = 0; i < nHits; i++)
    {
        auto map1 = dynamic_cast<R3BRpcMappedData*>(fMappedDataCA->At(i));

        UInt_t iDetector = map1->GetDetId();
        if (iDetector != 2)
        {
            continue;
        }

        UInt_t iFpga;
        if (map1->GetChannelId() == 7)
        {
            iFpga = 5;
        }
        else
        {
            iFpga = map1->GetChannelId(); // now 1..4
        }

        UInt_t iEdge_Side = 0; // 0
        R3BTCalModulePar* par_Refs = fTCalPar->GetModuleParAt(iDetector + 1, iFpga, iEdge_Side + 1);

        if (NULL == par_Refs)
        {
            continue;
        }

        // Convert TDC to [ns] ...

        times_Ref = par_Refs->GetTimeVFTX(map1->GetFineTime());

        if (times_Ref < 0. || times_Ref > fClockFreq) //|| IS_NAN(times_Ref))
        {

            LOG(info) << "R3BRpcMapped2PreCal::Exec : Bad time in ns: det= " << iDetector << ", FPGA= " << iFpga
                      << ", time in channels = " << map1->GetFineTime() << ", time in ns = " << times_Ref;
            continue;
        }

        auto& entry = Ref_vec[map1->GetChannelId() - 1];
        entry.time = 5. * map1->GetCoarseTime() - times_Ref;
        entry.Mapped = map1;
    }

    std::vector<Entry> strip_vec[41][2][2];

    std::vector<Entry> bar_vec[5][2][2];
    double leading_strip = 0;
    double leading_pmt = 0;

    double c = 2048. * 5.;
    nHits = fMappedDataCA->GetEntries();
    for (Int_t i = 0; i < nHits; i++)
    {
        auto map1 = dynamic_cast<R3BRpcMappedData*>(fMappedDataCA->At(i));
        UInt_t iDetector = map1->GetDetId();
        UInt_t iFpga;
        if (map1->GetChannelId() == 7)
        {
            iFpga = 5;
        }
        else
        {
            iFpga = map1->GetChannelId(); // now 1..4
        }
        if (iDetector == 2)
        {
            TClonesArray& clref = *fRpcPreCalDataCA;
            Int_t size = clref.GetEntriesFast();
            new (clref[size]) R3BRpcPreCalData(2, iFpga, Ref_vec[iFpga - 1].time, 0, 0);
        }
        // loop over strip data
        if (iDetector == 0)
        {
            UInt_t iStrip = map1->GetChannelId();                      // now 1..41
            UInt_t iEdge_Side = map1->GetEdge() * 2 + map1->GetSide(); // 0,3
            R3BTCalModulePar* par_Strips = fTCalPar->GetModuleParAt(iDetector + 1, iStrip, iEdge_Side + 1);

            if (NULL == par_Strips)
            {
                continue;
            }

            // Convert TDC to [ns] ...

            times_Strip = par_Strips->GetTimeVFTX(map1->GetFineTime());

            if (times_Strip < 0. || times_Strip > fClockFreq) //|| IS_NAN(times_Strip))
            {

                LOG(info) << "R3BRpcMapped2PreCal::Exec : Bad time in ns: det= " << iDetector << ", strip= " << iStrip
                          << ", side= " << map1->GetSide() << ", edge= " << map1->GetEdge()
                          << ", time in channels = " << map1->GetFineTime() << ", time in ns = " << times_Strip;
                continue;
            }
            Entry entry;
            entry.time =
                fmod((5. * map1->GetCoarseTime() - times_Strip - Ref_vec[lut[iStrip - 1][map1->GetSide()][0]].time -
                      (Ref_vec[8].time - Ref_vec[6].time) + 4 * c + c / 2),
                     c) -
                c / 2;
            entry.Mapped = map1;
            strip_vec[map1->GetChannelId() - 1][map1->GetSide()][map1->GetEdge()].push_back(entry);
        }

        // loop over pmt data

        if (iDetector == 1)
        {
            UInt_t iPmt = map1->GetChannelId();                        // now 1..41
            UInt_t iEdge_Side = map1->GetEdge() * 2 + map1->GetSide(); // 0,3
            R3BTCalModulePar* par_Pmts = fTCalPar->GetModuleParAt(iDetector + 1, iPmt, iEdge_Side + 1);

            if (NULL == par_Pmts)
            {
                continue;
            }

            // Convert TDC to [ns] ...

            times_Pmt = par_Pmts->GetTimeVFTX(map1->GetFineTime());

            if (times_Pmt < 0. || times_Pmt > fClockFreq) //|| IS_NAN(times_Pmt))
            {

                LOG(info) << "R3BRpcMapped2PreCal::Exec : Bad time in ns: det= " << iDetector << ", pmt= " << iPmt
                          << ", side= " << map1->GetSide() << ", edge= " << map1->GetEdge()
                          << ", time in channels = " << map1->GetFineTime() << ", time in ns = " << times_Pmt << " "
                          << fClockFreq;
                continue;
            }

            Entry entry;
            entry.time = fmod((5. * map1->GetCoarseTime() - times_Pmt - Ref_vec[1].time -
                               (Ref_vec[8].time - Ref_vec[6].time) + 4 * c + c / 2),
                              c) -
                         c / 2;
            entry.Mapped = map1;
            bar_vec[map1->GetChannelId() - 1][map1->GetSide()][map1->GetEdge()].push_back(entry);
        }
    }

    for (Int_t strip = 0; strip < 41; strip++)
    {
        for (Int_t side = 0; side < 2; side++)
        {

            size_t lead_i = 0, trail_i = 0;
            for (;;)
            {
                if (lead_i >= strip_vec[strip][side][1].size() || trail_i >= strip_vec[strip][side][0].size())
                {
                    break;
                }

                auto const& entry_lead = strip_vec[strip][side][1].at(lead_i);
                auto const& entry_trail = strip_vec[strip][side][0].at(trail_i);

                double tot = fmod((entry_trail.time - entry_lead.time + c + c / 2), c) - c / 2;
                if (tot < 0)
                {
                    ++trail_i;
                    continue;
                }
                // It fills the R3BRpcStripCalData
                TClonesArray& clref = *fRpcPreCalDataCA;
                Int_t size = clref.GetEntriesFast();
                new (clref[size]) R3BRpcPreCalData(0, strip + 1, entry_lead.time, tot, side);
                ++lead_i;
                ++trail_i;

                leading_strip = entry_lead.time;
            }
        }
    }

    for (Int_t bar = 0; bar < 5; bar++)
    {
        for (Int_t side = 0; side < 2; side++)
        {

            size_t lead_i = 0, trail_i = 0;
            for (;;)
            {
                if (lead_i >= bar_vec[bar][side][1].size() || trail_i >= bar_vec[bar][side][0].size())
                {
                    break;
                }

                auto const& entry_lead = bar_vec[bar][side][1].at(lead_i);
                auto const& entry_trail = bar_vec[bar][side][0].at(trail_i);

                double tot = fmod((entry_trail.time - entry_lead.time + c + c / 2), c) - c / 2;
                if (tot < 0)
                {
                    ++trail_i;
                    continue;
                }
                // It fills the R3BRpcPmtCalData
                TClonesArray& clref = *fRpcPreCalDataCA;
                Int_t size = clref.GetEntriesFast();
                new (clref[size]) R3BRpcPreCalData(1, bar + 1, entry_lead.time, tot, side);
                ++lead_i;
                ++trail_i;

                leading_pmt = entry_lead.time;
            }
        }
    }
    return;
}

void R3BRpcMapped2PreCal::Finish() {}

void R3BRpcMapped2PreCal::Reset()
{
    LOG(debug) << "Clearing RPCCalData Structure";
    if (fRpcPreCalDataCA)
    {
        fRpcPreCalDataCA->Clear();
    }
}

ClassImp(R3BRpcMapped2PreCal)
