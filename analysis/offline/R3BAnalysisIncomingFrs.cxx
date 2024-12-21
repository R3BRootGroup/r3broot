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

// ---------------------------------------------------------------------
// -----           R3BAnalysisIncomingFrs source file              -----
// -----        Created 03/05/21  by J.L. Rodriguez-Sanchez        -----
// ---------------------------------------------------------------------

/*
 * This task should make the analysis of the incoming projectiles from FRS
 *
 */

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "R3BAnalysisIncomingFrs.h"
#include "R3BEventHeader.h"
#include "R3BLosCalData.h"
#include "R3BLosMappedData.h"
#include "R3BMusicHitData.h"
#include "R3BSamplerMappedData.h"
#include "R3BSci2TcalData.h"
#include "R3BTCalEngine.h"

#include "TClonesArray.h"
#include "TMath.h"
#include <TRandom3.h>
#include <TRandomGen.h>
#include <algorithm>
#include <array>
#include <vector>
#define IS_NAN(x) TMath::IsNaN(x)
using namespace std;

#define SPEED_OF_LIGHT_MNS 0.299792458

R3BAnalysisIncomingFrs::R3BAnalysisIncomingFrs()
    : FairTask("AnalysisIncomingFrs", 1)
    , fTrigger(-1)
    , fTpat(-1)
    , fClockFreq(5.)
    , fNEvents(0)
    , fTcalSci2(NULL)
    , fHitItemsMus(NULL)
    , fFrsDataCA(NULL)
    , fToFoffset(0)
    , fToFmin(-5000)
    , fToFmax(5000)
    , fTof2InvV_p0(-7.8)
    , fTof2InvV_p1(0.0073)
    , fL2(137)
    , fPos_p0(-11)
    , fPos_p1(54.7)
    , fDispersionS2(7000)
    , fBrho0_S2toCC(9.458)
    , fP0(-2.12371e7)
    , fP1(4.9473e7)
    , fP2(-2.87635e7)
    , fZprimary(50.)
    , fZoffset(-1.3)
    , fPos2min(-200.)
    , fPos2max(200.)
    , ZMUSIC_cut(50)
    , ZMUSIC_wcut(0.5)
    , fOnline(kFALSE)
{
}

R3BAnalysisIncomingFrs::R3BAnalysisIncomingFrs(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrigger(-1)
    , fTpat(-1)
    , fClockFreq(5.)
    , fNEvents(0)
    , fTcalSci2(NULL)
    , fHitItemsMus(NULL)
    , fFrsDataCA(NULL)
    , fToFoffset(0)
    , fToFmin(-5000)
    , fToFmax(5000)
    , fTof2InvV_p0(-7.8)
    , fTof2InvV_p1(0.0073)
    , fL2(137)
    , fPos_p0(-11)
    , fPos_p1(54.7)
    , fDispersionS2(7000.)
    , fBrho0_S2toCC(9.458)
    , fP0(-2.12371e7)
    , fP1(4.9473e7)
    , fP2(-2.87635e7)
    , fZprimary(50.)
    , fZoffset(-1.3)
    , ZMUSIC_cut(50)
    , ZMUSIC_wcut(0.5)
    , fOnline(kFALSE)
{
}

R3BAnalysisIncomingFrs::~R3BAnalysisIncomingFrs()
{
    if (fTcalSci2)
        delete fTcalSci2;
    if (fHitItemsMus)
        delete fHitItemsMus;
    if (fFrsDataCA)
        delete fFrsDataCA;
}

InitStatus R3BAnalysisIncomingFrs::Init()
{
    // Initialize random number:
    std::srand(std::time(0)); // use current time as seed for random generator

    LOG(info) << "R3BAnalysisIncomingFrs::Init ";

    // try to get a handle on the EventHeader. EventHeader may not be
    // present though and hence may be null. Take care when using.

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
        LOG(fatal) << "FairRootManager not found";

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("R3BEventHeader"));

    // Get objects for detectors on all levels
    assert(DET_MAX + 1 == sizeof(fDetectorNames) / sizeof(fDetectorNames[0]));
    for (int det = 0; det < DET_MAX; det++)
    {
        fMappedItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sMapped", fDetectorNames[det]))));
        if (NULL == fMappedItems.at(det))
        {
            printf("Could not find mapped data for '%s'.\n", fDetectorNames[det]);
        }
        fCalItems.push_back(dynamic_cast<TClonesArray*>(mgr->GetObject(Form("%sCal", fDetectorNames[det]))));
    }

    // --- Get access to Sci2 data at Tcal level --- //
    fTcalSci2 = dynamic_cast<TClonesArray*>(mgr->GetObject("Sci2Tcal"));
    if (!fTcalSci2)
    {
        LOG(info) << "R3BAnalysisIncomingFrs::Init()  Could not find ci2Tcal";
    }

    // get access to hit data of the MUSIC
    fHitItemsMus = dynamic_cast<TClonesArray*>(mgr->GetObject("MusicHitData"));
    if (!fHitItemsMus)
        LOG(warn) << "R3BAnalysisIncomingFrs: MusicHitData not found";

    // OUTPUT DATA
    fFrsDataCA = new TClonesArray("R3BFrsData", 5);
    if (!fOnline)
    {
        mgr->Register("FrsData", "Analysis FRS", fFrsDataCA, kTRUE);
    }
    else
    {
        mgr->Register("FrsData", "Analysis FRS", fFrsDataCA, kFALSE);
    }

    return kSUCCESS;
}

void R3BAnalysisIncomingFrs::Exec(Option_t* option)
{
    //  cout << "fNEvents " << fNEvents << endl;

    FairRootManager* mgr = FairRootManager::Instance();
    if (NULL == mgr)
    {
        // FairLogger::GetLogger()->Fatal(MESSAGE_ORIGIN, "FairRootManager not found");
        LOG(error) << "FairRootManager not found";
        return;
    }

    double Zmusic = 0., Music_ang = 0.;
    if (fHitItemsMus && fHitItemsMus->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitItemsMus->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BMusicHitData* hit = dynamic_cast<R3BMusicHitData*>(fHitItemsMus->At(ihit));
            if (!hit)
                continue;
            Zmusic = hit->GetZcharge();
            Music_ang = hit->GetTheta() * 1000.; // mrad
        }
    }

    // --- --------------- --- //
    // --- local variables --- //
    // --- --------------- --- //
    Int_t multTcal[3];
    Double_t iRawTimeNs[3][64];
    UInt_t nHits, iCh, cpt = 0;
    Double_t ToFraw = 0., PosRaw = -10000., ToFcal = 0., PosCal = -10000.;
    Double_t Velo = 0., Beta = 0., Gamma = 0., Brho = 0., AoQ = 0.;
    Double_t ToFraw_m1 = 0., PosRaw_m1 = -10000., ToFcal_m1 = 0., PosCal_m1 = -10000.;
    Double_t Velo_m1 = 0., Beta_m1 = 0., Gamma_m1 = 0., Brho_m1 = 0., AoQ_m1 = 0.;
    for (UShort_t pmt = 0; pmt < 3; pmt++)
    {
        multTcal[pmt] = 0;
        for (UShort_t m = 0; m < 64; m++)
            iRawTimeNs[pmt][m] = 0.;
    }

    // --- -------------- --- //
    // --- read tcal data --- //
    // --- -------------- --- //
    if (fTcalSci2 && fTcalSci2->GetEntriesFast())
    {
        nHits = fTcalSci2->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSci2TcalData* hittcal = dynamic_cast<R3BSci2TcalData*>(fTcalSci2->At(ihit));
            if (!hittcal)
                continue;
            iCh = hittcal->GetChannel() - 1;
            iRawTimeNs[iCh][multTcal[iCh]] = hittcal->GetRawTimeNs();
            multTcal[iCh]++;
        } // --- end of loop over Tcal data --- //
    }

    time = header->GetTimeStamp(); // / 1.6; // divided by 1.6 for stand alone daq with Vulom generating time stamps

    // time = 0;

    if (time_start == 0 && time > 0)
    {
        time_start = time;
        fNEvents_start = fNEvents;
    }

    if (header->GetTrigger() == 12)
        time_spill_start = time; // header->GetTimeStamp();    // spill start in nsec
    if (header->GetTrigger() == 13)
        time_spill_end = time; // header->GetTimeStamp();    // spill end  in nsec

    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;

    if (fMappedItems.at(DET_SAMPLER))
    {

        auto det = fMappedItems.at(DET_SAMPLER);
        Int_t nHitsSamp = det->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHitsSamp; ihit++)
        {
            auto hit = dynamic_cast<R3BSamplerMappedData*>(det->At(ihit));
            // time is in steps of 10 ns
            // is is a 34 bit number, so max 1073741823
            samplerCurr = hit->GetTime();
            Double_t dt = ((double)(samplerCurr - samplerPrev)) / 100.;

            samplerPrev = samplerCurr;
        }
    }

    //----------------------------------------------------------------------
    // LOS detector
    //----------------------------------------------------------------------
    // Time nomenclature for LOS:
    // lt=0, l=1,lb=2,b=3,rb=4,r=5,rt=6,t=7
    //   t / \lt
    //   rt| |l
    //   r | |lb
    //   rb\ /b

    Double_t timeTofd = 0;

    Double_t timeLosV[fNofLosDetectors][32]; // mean time from VFTX
    Double_t LosTresV[fNofLosDetectors][32];
    Double_t timeLosT[fNofLosDetectors][32];
    Double_t LosTresT[fNofLosDetectors][32];
    Double_t timeLos[fNofLosDetectors][32];
    Double_t totsum[fNofLosDetectors][32];
    Double_t xT_cm[fNofLosDetectors][32];
    Double_t yT_cm[fNofLosDetectors][32];
    Double_t xToT_cm[fNofLosDetectors][32];
    Double_t yToT_cm[fNofLosDetectors][32];
    Double_t xV_cm[fNofLosDetectors][32];
    Double_t yV_cm[fNofLosDetectors][32];

    Double_t time_V[fNofLosDetectors][32][8]; // [det][multihit][pm]
    Double_t time_L[fNofLosDetectors][32][8];
    Double_t time_T[fNofLosDetectors][32][8];
    Double_t tot[fNofLosDetectors][32][8];
    Double_t time_MTDC[32][8] = { { 0. } };
    Double_t LosTresMTDC[32];

    for (Int_t idet = 0; idet < fNofLosDetectors; idet++)
    {
        for (Int_t imult = 0; imult < 32; imult++)
        {
            timeLosV[idet][imult] = 0.0;
            LosTresV[idet][imult] = 0.0 / 0.0;
            timeLosT[idet][imult] = 0.0;
            LosTresT[idet][imult] = 0.0 / 0.0;
            timeLos[idet][imult] = 0.0;
            totsum[idet][imult] = 0.0;
            xT_cm[idet][imult] = 0.0 / 0.0;
            yT_cm[idet][imult] = 0.0 / 0.0;
            xToT_cm[idet][imult] = -100000.;
            yToT_cm[idet][imult] = -100000.;
            xV_cm[idet][imult] = 0.0 / 0.0;
            yV_cm[idet][imult] = 0.0 / 0.0;
            for (Int_t icha = 0; icha < 8; icha++)
            {
                time_V[idet][imult][icha] = 0.0 / 0.0; // [det][multihit][pm]
                time_L[idet][imult][icha] = 0.0 / 0.0;
                time_T[idet][imult][icha] = 0.0 / 0.0;
                tot[idet][imult][icha] = 0.0 / 0.0;
            }
        }
    }
    Int_t nPartLOS = 0;
    Int_t nPartc[fNofLosDetectors];
    for (Int_t d = 0; d < fNofLosDetectors; d++)
        nPartc[d] = 0;

    Bool_t iLOSType[fNofLosDetectors][32];
    Bool_t iLOSPileUp[fNofLosDetectors][32];
    for (Int_t idet = 0; idet < fNofLosDetectors; idet++)
    {
        for (Int_t imult = 0; imult < 32; imult++)
        {
            iLOSType[idet][imult] = false;
            iLOSPileUp[idet][imult] = false;
        }
    }

    Int_t Multip;

    if (fMappedItems.at(DET_LOS))
    {
        auto det = fMappedItems.at(DET_LOS);
        nHits = det->GetEntriesFast();

        //       cout<<"nHits in LOS: "<<nHits<<endl;
        if (nHits > 0)
            nLosEvents += 1;
        Multip = nHits;

        //    cout<<"LOS: nHits"<<nHits<<endl;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BLosMappedData* hit = dynamic_cast<R3BLosMappedData*>(det->At(ihit));
            if (!hit)
                continue;

            // channel numbers are stored 1-based (1..n)
            Int_t iDet = hit->GetDetector(); // 1..
            Int_t iCha = hit->GetChannel();  // 1..
            Int_t iTyp = hit->GetType();     // 0,1,2,3
        }
    }

    if (fCalItems.at(DET_LOS))
    {
        auto det = fCalItems.at(DET_LOS);
        nPartLOS = det->GetEntriesFast();

        Int_t iDet = 0;
        Double_t time_V_LOS1[32][8] = { { 0. } };
        Double_t time_V_LOS2[32][8] = { { 0. } };

        for (Int_t iPart = 0; iPart < nPartLOS; iPart++)
        {
            /*
             * nPart is the number of particle passing through LOS detector in one event
             */
            R3BLosCalData* calData = dynamic_cast<R3BLosCalData*>(det->At(iPart));
            iDet = calData->GetDetector();

            Double_t sumvtemp = 0, sumltemp = 0, sumttemp = 0;
            for (Int_t iCha = 0; iCha < 8; iCha++)
            {
                sumvtemp += calData->GetTimeV_ns(iCha);
                sumltemp += calData->GetTimeL_ns(iCha);
                sumttemp += calData->GetTimeT_ns(iCha);
            }
            if (!(IS_NAN(sumvtemp)) && !(IS_NAN(sumltemp)) && !(IS_NAN(sumltemp)))
            {
                nPartc[iDet - 1]++;
            }
            else
            {
                continue;
            }

            for (Int_t iCha = 0; iCha < 8; iCha++)
            {
                if (iDet == 1)
                {
                    time_V_LOS1[nPartc[iDet - 1] - 1][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeV_ns(iCha))))
                    { // VFTX
                        time_V_LOS1[nPartc[iDet - 1] - 1][iCha] = calData->GetTimeV_ns(iCha);
                    }
                    if (!(IS_NAN(calData->GetTimeL_ns(iCha))))
                    { // TAMEX leading
                        time_L[iDet - 1][nPartc[iDet - 1] - 1][iCha] = calData->GetTimeL_ns(iCha);
                    }

                    if (!(IS_NAN(calData->GetTimeT_ns(iCha))))
                    { // TAMEX trailing
                        time_T[iDet - 1][nPartc[iDet - 1] - 1][iCha] = calData->GetTimeT_ns(iCha);
                    }
                    time_MTDC[nPartc[iDet - 1]][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeM_ns(iCha))))
                    { // MTDC
                        time_MTDC[nPartc[iDet - 1] - 1][iCha] = calData->GetTimeM_ns(iCha);
                    }
                }
                if (iDet == 2)
                {
                    time_V_LOS2[nPartc[iDet - 1] - 1][iCha] = 0. / 0.;
                    if (!(IS_NAN(calData->GetTimeV_ns(iCha))))
                    { // VFTX
                        time_V_LOS2[nPartc[iDet - 1] - 1][iCha] = calData->GetTimeV_ns(iCha);
                    }
                    if (!(IS_NAN(calData->GetTimeL_ns(iCha))))
                    { // TAMEX leading
                        time_L[iDet - 1][nPartc[iDet - 1] - 1][iCha] = calData->GetTimeL_ns(iCha);
                    }
                    if (!(IS_NAN(calData->GetTimeT_ns(iCha))))
                    { // TAMEX trailing
                        time_T[iDet - 1][nPartc[iDet - 1] - 1][iCha] = calData->GetTimeT_ns(iCha);
                    }
                }
            }

            if (!calData)
            {
                cout << "LOS !calData" << endl;
                continue; // can this happen?
            }
        }

        // Sorting VFTX data:

        // detector 1
        if (nPartc[0] > 0)
        {
            std::qsort(time_V_LOS1,
                       nPartc[0],
                       sizeof(*time_V_LOS1),
                       [](const void* arg1, const void* arg2) -> int
                       {
                           double const* lhs = static_cast<double const*>(arg1);
                           double const* rhs = static_cast<double const*>(arg2);

                           return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
                       });
            for (Int_t iPart = 0; iPart < nPartc[0]; iPart++)
            {
                for (int ipm = 0; ipm < 8; ipm++)
                {
                    time_V[0][iPart][ipm] = time_V_LOS1[iPart][ipm];
                }
            }
        }

        // detector 2
        if (fNofLosDetectors > 1 && nPartc[1] > 0)
        {
            std::qsort(time_V_LOS2,
                       nPartc[1],
                       sizeof(*time_V_LOS2),
                       [](const void* arg1, const void* arg2) -> int
                       {
                           double const* lhs = static_cast<double const*>(arg1);
                           double const* rhs = static_cast<double const*>(arg2);

                           return (lhs[0] < rhs[0]) ? -1 : ((rhs[0] < lhs[0]) ? 1 : 0);
                       });
            for (Int_t iPart = 0; iPart < nPartc[1]; iPart++)
            {
                for (int ipm = 0; ipm < 8; ipm++)
                {
                    time_V[1][iPart][ipm] = time_V_LOS2[iPart][ipm];
                }
            }
        }

        // End sorting

        std::vector<double_t> time_first(fNofLosDetectors, -1.);
        std::vector<double_t> time0(fNofLosDetectors, -1.);
        std::vector<double_t> time1(fNofLosDetectors, -1.);
        std::vector<double_t> time_abs(fNofLosDetectors, -1.);

        for (iDet = 1; iDet <= fNofLosDetectors; iDet++)
        {
            cpt = 0;
            for (Int_t iPart = 0; iPart < nPartc[iDet - 1]; iPart++)
            {
                Bool_t iLOSTypeMCFD = false;
                Bool_t iLOSTypeTAMEX = false;

                if (time_V[iDet - 1][iPart][0] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][0])) &&
                    time_V[iDet - 1][iPart][1] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][1])) &&
                    time_V[iDet - 1][iPart][2] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][2])) &&
                    time_V[iDet - 1][iPart][3] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][3])) &&
                    time_V[iDet - 1][iPart][4] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][4])) &&
                    time_V[iDet - 1][iPart][5] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][5])) &&
                    time_V[iDet - 1][iPart][6] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][6])) &&
                    time_V[iDet - 1][iPart][7] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][7])))
                {
                    iLOSTypeMCFD = true; // all 8 MCFD times
                }

                if (time_L[iDet - 1][iPart][0] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][0])) &&
                    time_L[iDet - 1][iPart][1] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][1])) &&
                    time_L[iDet - 1][iPart][2] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][2])) &&
                    time_L[iDet - 1][iPart][3] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][3])) &&
                    time_L[iDet - 1][iPart][4] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][4])) &&
                    time_L[iDet - 1][iPart][5] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][5])) &&
                    time_L[iDet - 1][iPart][6] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][6])) &&
                    time_L[iDet - 1][iPart][7] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][7])) &&

                    time_T[iDet - 1][iPart][0] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][0])) &&
                    time_T[iDet - 1][iPart][1] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][1])) &&
                    time_T[iDet - 1][iPart][2] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][2])) &&
                    time_T[iDet - 1][iPart][3] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][3])) &&
                    time_T[iDet - 1][iPart][4] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][4])) &&
                    time_T[iDet - 1][iPart][5] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][5])) &&
                    time_T[iDet - 1][iPart][6] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][6])) &&
                    time_T[iDet - 1][iPart][7] > 0. && !(IS_NAN(time_T[iDet - 1][iPart][7])))
                {
                    iLOSTypeTAMEX = true; // all 8 leading and trailing times
                }

                // We will consider only events in which booth MCFD and TAMEX see same number of channels:
                if (iLOSTypeTAMEX && iLOSTypeMCFD)
                    iLOSType[iDet - 1][iPart] = true;

                // if(iDet < fNofLosDetectors+1)
                // {

                if (iLOSType[iDet - 1][iPart])
                {
                    int nPMT = 0;
                    int nPMV = 0;

                    for (int ipm = 0; ipm < 8; ipm++)
                    {

                        if (time_T[iDet - 1][iPart][ipm] > 0. && time_L[iDet - 1][iPart][ipm] > 0. &&
                            !(IS_NAN(time_T[iDet - 1][iPart][ipm])) && !(IS_NAN(time_L[iDet - 1][iPart][ipm])))
                        {
                            while (time_T[iDet - 1][iPart][ipm] - time_L[iDet - 1][iPart][ipm] <= 0.)
                            {
                                time_T[iDet - 1][iPart][ipm] = time_T[iDet - 1][iPart][ipm] + 2048. * fClockFreq;
                            }

                            nPMT = nPMT + 1;
                            tot[iDet - 1][iPart][ipm] = time_T[iDet - 1][iPart][ipm] - time_L[iDet - 1][iPart][ipm];

                            // pileup rejection
                            if (tot[iDet - 1][iPart][ipm] > fEpileup)
                                iLOSPileUp[iDet - 1][iPart] = true;
                        }

                        if (tot[iDet - 1][iPart][ipm] != 0. && !(IS_NAN(tot[iDet - 1][iPart][ipm])))
                            totsum[iDet - 1][iPart] += tot[iDet - 1][iPart][ipm];

                        if (time_L[iDet - 1][iPart][ipm] > 0. && !(IS_NAN(time_L[iDet - 1][iPart][ipm])))
                            timeLosT[iDet - 1][iPart] += time_L[iDet - 1][iPart][ipm];

                        // Calculate detector time
                        if (time_V[iDet - 1][iPart][ipm] > 0. && !(IS_NAN(time_V[iDet - 1][iPart][ipm])))
                        {
                            timeLosV[iDet - 1][iPart] += time_V[iDet - 1][iPart][ipm];
                            nPMV = nPMV + 1;
                        }
                    }

                    totsum[iDet - 1][iPart] = totsum[iDet - 1][iPart] / nPMT;

                    timeLosV[iDet - 1][iPart] = timeLosV[iDet - 1][iPart] / nPMV;

                    timeLosT[iDet - 1][iPart] = timeLosT[iDet - 1][iPart] / nPMT;

                    timeLos[iDet - 1][iPart] = timeLosV[iDet - 1][iPart];

                    // cout<<"LOS: "<<iDet<<", "<<iPart<<"; "<<timeLos[iDet-1][iPart]<<endl;

                    LosTresV[iDet - 1][iPart] = ((time_V[iDet - 1][iPart][0] + time_V[iDet - 1][iPart][2] +
                                                  time_V[iDet - 1][iPart][4] + time_V[iDet - 1][iPart][6]) -
                                                 (time_V[iDet - 1][iPart][1] + time_V[iDet - 1][iPart][3] +
                                                  time_V[iDet - 1][iPart][5] + time_V[iDet - 1][iPart][7])) /
                                                4.;

                    LosTresT[iDet - 1][iPart] = ((time_L[iDet - 1][iPart][0] + time_L[iDet - 1][iPart][2] +
                                                  time_L[iDet - 1][iPart][4] + time_L[iDet - 1][iPart][6]) -
                                                 (time_L[iDet - 1][iPart][1] + time_L[iDet - 1][iPart][3] +
                                                  time_L[iDet - 1][iPart][5] + time_L[iDet - 1][iPart][7])) /
                                                4.;

                    LosTresT[iDet - 1][iPart] = ((time_V[iDet - 1][iPart][0] * 0. + time_V[iDet - 1][iPart][2] +
                                                  time_V[iDet - 1][iPart][4] + time_V[iDet - 1][iPart][6]) /
                                                     3. -
                                                 (time_V[iDet - 1][iPart][1] + time_V[iDet - 1][iPart][3] +
                                                  time_V[iDet - 1][iPart][5] + time_V[iDet - 1][iPart][7]) /
                                                     4.);

                    // right koord.-system, Z-axis beam direction:
                    // Position from tamex:
                    xT_cm[iDet - 1][iPart] = (time_L[iDet - 1][iPart][1] + time_L[iDet - 1][iPart][2]) / 2. -
                                             (time_L[iDet - 1][iPart][5] + time_L[iDet - 1][iPart][6]) / 2.;
                    yT_cm[iDet - 1][iPart] = (time_L[iDet - 1][iPart][3] + time_L[iDet - 1][iPart][4]) / 2. -
                                             (time_L[iDet - 1][iPart][7] + time_L[iDet - 1][iPart][0]) / 2.;
                    xT_cm[iDet - 1][iPart] = (xT_cm[iDet - 1][iPart] - flosOffsetXT[iDet - 1]) * flosVeffXT[iDet - 1];
                    yT_cm[iDet - 1][iPart] = (yT_cm[iDet - 1][iPart] - flosOffsetYT[iDet - 1]) * flosVeffYT[iDet - 1];

                    // Position from VFTX:
                    xV_cm[iDet - 1][iPart] = (time_V[iDet - 1][iPart][1] + time_V[iDet - 1][iPart][2]) / 2. -
                                             (time_V[iDet - 1][iPart][5] + time_V[iDet - 1][iPart][6]) / 2.;
                    yV_cm[iDet - 1][iPart] = (time_V[iDet - 1][iPart][3] + time_V[iDet - 1][iPart][4]) / 2. -
                                             (time_V[iDet - 1][iPart][7] + time_V[iDet - 1][iPart][0]) / 2.;
                    xV_cm[iDet - 1][iPart] = (xV_cm[iDet - 1][iPart] - flosOffsetXV[iDet - 1]) * flosVeffXV[iDet - 1];
                    yV_cm[iDet - 1][iPart] = (yV_cm[iDet - 1][iPart] - flosOffsetYV[iDet - 1]) * flosVeffYV[iDet - 1];

                    // Position from ToT:
                    if (tot[iDet - 1][iPart][1] > 0. && tot[iDet - 1][iPart][2] > 0. && tot[iDet - 1][iPart][5] > 0. &&
                        tot[iDet - 1][iPart][6] > 0. && tot[iDet - 1][iPart][0] > 0. && tot[iDet - 1][iPart][3] > 0. &&
                        tot[iDet - 1][iPart][4] > 0. && tot[iDet - 1][iPart][7] > 0.)
                    {
                        xToT_cm[iDet - 1][iPart] = (((tot[iDet - 1][iPart][5] + tot[iDet - 1][iPart][6]) / 2. -
                                                     (tot[iDet - 1][iPart][1] + tot[iDet - 1][iPart][2]) / 2.) /
                                                    ((tot[iDet - 1][iPart][1] + tot[iDet - 1][iPart][2] +
                                                      tot[iDet - 1][iPart][5] + tot[iDet - 1][iPart][6]) /
                                                     4.));

                        yToT_cm[iDet - 1][iPart] = (((tot[iDet - 1][iPart][0] + tot[iDet - 1][iPart][7]) / 2. -
                                                     (tot[iDet - 1][iPart][3] + tot[iDet - 1][iPart][4]) / 2.) /
                                                    ((tot[iDet - 1][iPart][7] + tot[iDet - 1][iPart][0] +
                                                      tot[iDet - 1][iPart][3] + tot[iDet - 1][iPart][4]) /
                                                     4.));
                    }

                    xToT_cm[iDet - 1][iPart] =
                        (xToT_cm[iDet - 1][iPart] - flosOffsetXQ[iDet - 1]) * flosVeffXQ[iDet - 1];
                    yToT_cm[iDet - 1][iPart] =
                        (yToT_cm[iDet - 1][iPart] - flosOffsetYQ[iDet - 1]) * flosVeffYQ[iDet - 1];

                    if (timeLosV[iDet - 1][iPart] > 0. && timeLosV[iDet - 1][iPart] < 8192. * 5. &&
                        !(IS_NAN(timeLosV[iDet - 1][iPart])))
                    {
                        while (timeLosT[iDet - 1][iPart] - timeLosV[iDet - 1][iPart] < 2048. * 5. / 2.)
                        {
                            timeLosT[iDet - 1][iPart] = timeLosT[iDet - 1][iPart] + 2048. * 5.;
                        }
                        while (timeLosT[iDet - 1][iPart] - timeLosV[iDet - 1][iPart] > 2048. * 5. / 2.)
                        {
                            timeLosT[iDet - 1][iPart] = timeLosT[iDet - 1][iPart] - 2048. * 5.;
                        }
                    }

                } // if iLosType

            } // for iPart

            // select multiplicity == 1 at Cave C and S2
            if (iDet == 1 && nPartc[0] == 1 && multTcal[0] == 1 && multTcal[1] == 1)
            {

                // --- -----------------------------
                // --- secondary beam identification
                // --- -----------------------------
                // if X is increasing from left to right:
                //    Brho = fBhro0 * (1 - xMwpc0/fDCC + xS2/fDS2)
                // in R3BRoot, X is increasing from right to left
                //    Bro = fBrho0 * (1 + xMwpc0/fDCC - xS2/fDS2)

                ToFraw_m1 = timeLosV[0][0] - 0.5 * (iRawTimeNs[0][0] + iRawTimeNs[1][0]);
                Velo_m1 = 1. / (fTof2InvV_p0 + fTof2InvV_p1 * (fToFoffset + ToFraw_m1)); // [m/ns]
                Beta_m1 = Velo_m1 / 0.299792458;
                Gamma_m1 = 1. / (TMath::Sqrt(1. - TMath::Power(Beta_m1, 2)));
                PosRaw_m1 = iRawTimeNs[0][0] - iRawTimeNs[1][0]; // [ns]
                PosCal_m1 = fPos_p0 + fPos_p1 * PosRaw_m1;       // [mm] at S2
                Brho_m1 = fBrho0_S2toCC * (1. - PosCal_m1 / fDispersionS2);
                AoQ_m1 = Brho_m1 / (3.10716 * Beta_m1 * Gamma_m1);
                if (Zmusic > 0.)
                {
                    double Emus = ((Zmusic + 4.7) / 0.28) * ((Zmusic + 4.7) / 0.28);
                    double zcor = sqrt(Emus * Beta_m1) * 0.277;
                    double x0_point = 0.;
                    double y0_point = 50.39;
                    double rot_ang = 0.0375;
                    double zcorang =
                        y0_point + (Music_ang - x0_point) * sin(rot_ang) + (zcor - y0_point) * cos(rot_ang) + 0.2;
                    AddData(1, 2, zcorang, AoQ_m1, Beta_m1, Brho_m1, PosCal_m1, 0.);
                }
            }
        } // for iDet

    } // if fCallItems

    fNEvents += 1;
}

void R3BAnalysisIncomingFrs::FinishEvent()
{

    for (Int_t det = 0; det < DET_MAX; det++)
    {
        if (fMappedItems.at(det))
        {
            fMappedItems.at(det)->Clear();
        }
        if (fCalItems.at(det))
        {
            fCalItems.at(det)->Clear();
        }
    }
    if (fTcalSci2)
        fTcalSci2->Clear();
    if (fHitItemsMus)
        fHitItemsMus->Clear();
}

void R3BAnalysisIncomingFrs::FinishTask() {}

// -----   Private method AddData  --------------------------------------------
R3BFrsData* R3BAnalysisIncomingFrs::AddData(Int_t StaId,
                                            Int_t StoId,
                                            Double_t z,
                                            Double_t aq,
                                            Double_t beta,
                                            Double_t brho,
                                            Double_t xs2,
                                            Double_t xc)
{
    // It fills the R3BSofFrsData
    TClonesArray& clref = *fFrsDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BFrsData(StaId, StoId, z, aq, beta, brho, xs2, xc);
}

ClassImp(R3BAnalysisIncomingFrs)
