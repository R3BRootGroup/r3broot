/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
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
// -----               R3BIncomingBeta                    -----
// -----    Created 10/05/22 by J.L. Rodriguez-Sanchez    -----
// ------------------------------------------------------------
#define MAXMULT 64

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "R3BCoarseTimeStitch.h"
#include "R3BEventHeader.h"
#include "R3BIncomingBeta.h"
#include "R3BIncomingIDPar.h"
#include "R3BLogger.h"
#include "R3BLosCalData.h"
#include "R3BLosHitData.h"
#include "R3BSci2HitData.h"
#include "R3BSci2TcalData.h"

#include "TClonesArray.h"
#include "TMath.h"

R3BIncomingBeta::R3BIncomingBeta()
    : R3BIncomingBeta("IncomingBeta", 1)
{
}

R3BIncomingBeta::R3BIncomingBeta(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTimeStitch(nullptr)
    , fIncomingID_Par(NULL)
    , fFrsDataCA(NULL)
    , fHitSci2(NULL)
    , fPosCalFrsSci(NULL)
    , fTcalSci2(NULL)
    , fOnline(kFALSE)
    , fStaId(1)
    , fStoId(2)
    , fLosRefCh(0)
    , fLosCalTrig_Low(3160.)
    , fLosCalTrig_High(3220.)
    , fNumDet(1)
    , fUseTref(kFALSE)
    , fUseMultHit(kFALSE)
    , fUseFrsSci(kFALSE)
{
    fToFoffset = new TArrayF(fNumDet);
    fPosS2Left = new TArrayF(fNumDet);
    fPosS2Right = new TArrayF(fNumDet);
    fTof2InvV_p0 = new TArrayF(fNumDet);
    fTof2InvV_p1 = new TArrayF(fNumDet);
}

R3BIncomingBeta::~R3BIncomingBeta()
{
    R3BLOG(debug1, "");
    if (fFrsDataCA)
        delete fFrsDataCA;
}

void R3BIncomingBeta::SetParContainers()
{
    R3BLOG(info, "");
    // Reading IncomingIDPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");

    fIncomingID_Par = dynamic_cast<R3BIncomingIDPar*>(rtdb->getContainer("IncomingIDPar"));
    R3BLOG_IF(fatal, !fIncomingID_Par, "Couldn't get handle on IncomingIDPar container");
    R3BLOG_IF(info, fIncomingID_Par, "IncomingIDPar container was found");

    return;
}

void R3BIncomingBeta::SetParameter()
{
    //--- Parameter Container ---
    fBeta_min = fIncomingID_Par->GetBeta_min();
    fBeta_max = fIncomingID_Par->GetBeta_max();

    for (Int_t i = 1; i < fNumDet + 1; i++)
    {
        fToFoffset->AddAt(fIncomingID_Par->GetToFoffset(i), i - 1);
        fPosS2Left->AddAt(fIncomingID_Par->GetPosS2Left(i), i - 1);
        fPosS2Right->AddAt(fIncomingID_Par->GetPosS2Right(i), i - 1);
        fTof2InvV_p0->AddAt(fIncomingID_Par->GetTof2InvV_p0(i), i - 1);
        fTof2InvV_p1->AddAt(fIncomingID_Par->GetTof2InvV_p1(i), i - 1);
    }
    return;
}

InitStatus R3BIncomingBeta::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();
    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    fHeader = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    if (!fUseFrsSci)
    {
        // Get access to Sci2 data at hit level
        fHitSci2 = dynamic_cast<TClonesArray*>(mgr->GetObject("Sci2Hit"));
        R3BLOG_IF(fatal, !fHitSci2, "Could not find Sci2Hit");
    }
    else
    {
        // Get access to Sci2 data at hit level
        fPosCalFrsSci = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsSciPosCalData"));
        R3BLOG_IF(fatal, !fPosCalFrsSci, "Could not find FrsSciPosCalData");
    }

    // Get access to hit data of the LOS
    fHitLos = dynamic_cast<TClonesArray*>(mgr->GetObject("LosHit"));
    R3BLOG_IF(fatal, !fHitLos, "LosHit not found");

    if (fHeader->GetExpId() == 91 || fHeader->GetExpId() == 118)
    {
        fCalLos = dynamic_cast<TClonesArray*>(mgr->GetObject("LosCal"));
        R3BLOG_IF(fatal, !fCalLos, "LosCal not found");

        fCalLosTrig = dynamic_cast<TClonesArray*>(mgr->GetObject("LosTriggerCal"));
        R3BLOG_IF(fatal, !fCalLosTrig, "LosTriggerCal not found");
    }

    // Get access to Sci2 data at Tcal level
    if (fHeader->GetExpId() == 509)
    {
        fTcalSci2 = dynamic_cast<TClonesArray*>(mgr->GetObject("Sci2Tcal"));
        R3BLOG_IF(fatal, !fTcalSci2, "Could not find Sci2Tcal");
    }

    // Output data
    fFrsDataCA = dynamic_cast<TClonesArray*>(mgr->GetObject("FrsData"));
    if (fFrsDataCA == NULL)
    {
        fFrsDataCA = new TClonesArray("R3BFrsData");
        mgr->Register("FrsData", "Analysis FRS", fFrsDataCA, !fOnline);
        R3BLOG(info, "FrsData register done, it is OK");
    }

    SetParameter();

    // Definition of a time stich object to correlate times coming from different systems
    fTimeStitch = new R3BCoarseTimeStitch();

    return kSUCCESS;
}

InitStatus R3BIncomingBeta::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BIncomingBeta::Exec(Option_t* option)
{
    // Reset entries in output arrays
    Reset();

    // --- local variables --- //
    Double_t timeLosV[fNumDet][MAXMULT];
    Double_t posLosX_cm[fNumDet][MAXMULT];
    Double_t TimeSci2_m1[fNumDet][MAXMULT];
    Double_t TimeSci2wTref_m1[fNumDet][MAXMULT];
    Double_t TimeSci2_tcal[fNumDet * 3][MAXMULT];
    Double_t PosSci2_m1[fNumDet][MAXMULT];
    UInt_t nHits = 0;
    Double_t ToFraw_m1 = 0.;
    Double_t Velo_m1 = 0., Beta_m1 = 0., Gamma_m1 = 0.;

    Int_t multSci2[fNumDet], multSci2Tcal[fNumDet * 3];
    Int_t multLos[fNumDet];

    for (Int_t i = 0; i < fNumDet; i++)
    {
        multSci2[i] = 0;
        multLos[i] = 0;
        for (Int_t m = 0; m < MAXMULT; m++)
        {
            PosSci2_m1[i][m] = 0.;
            TimeSci2_m1[i][m] = 0.;
            TimeSci2wTref_m1[i][m] = 0.;
            posLosX_cm[i][m] = 0.;
            timeLosV[i][m] = 0.;
            for (Int_t k = 0; k < 3; k++)
            {
                multSci2Tcal[i * 3 + k] = 0;
                TimeSci2_tcal[i * 3 + k][m] = 0;
            }
        }
    }

    // --- read Tcal from Sci2 data --- //
    if (fTcalSci2 && fTcalSci2->GetEntriesFast() > 0)
    {
        nHits = fTcalSci2->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSci2TcalData* hittcal = dynamic_cast<R3BSci2TcalData*>(fTcalSci2->At(ihit));
            UInt_t numDet = hittcal->GetDetector();
            UInt_t ch = hittcal->GetChannel() - 1;
            if (numDet > fNumDet)
            {
                R3BLOG(warn, "Sci2 detector id:" << numDet << " is out of range!");
                continue;
            }
            if (multSci2Tcal[(numDet - 1) * 3 + ch] >= MAXMULT)
                continue;
            TimeSci2_tcal[(numDet - 1) * 3 + ch][multSci2Tcal[(numDet - 1) * 3 + ch]] = hittcal->GetRawTimeNs();
            multSci2Tcal[(numDet - 1) * 3 + ch]++;
        } // --- end of loop over Sci2 Tcal data --- //
    }
    // --- read hit from Sci2 data --- //
    if (fHitSci2 && fHitSci2->GetEntriesFast() > 0)
    {
        Int_t numDet = 1;
        nHits = fHitSci2->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BSci2HitData* hittcal = dynamic_cast<R3BSci2HitData*>(fHitSci2->At(ihit));
            numDet = hittcal->GetSciId();
            if (numDet > fNumDet)
            {
                R3BLOG(warn, "Sci2 detector id:" << numDet << " is out of range!");
                continue;
            }
            if (multSci2[numDet - 1] >= MAXMULT)
                continue;
            PosSci2_m1[numDet - 1][multSci2[numDet - 1]] = hittcal->GetX();
            TimeSci2_m1[numDet - 1][multSci2[numDet - 1]] = hittcal->GetTime();
            TimeSci2wTref_m1[numDet - 1][multSci2[numDet - 1]] = hittcal->GetTimeWithTref();
            multSci2[numDet - 1]++;
        } // --- end of loop over hit data --- //
    }

    // --- read hit from Sci2 data --- //
    if (fPosCalFrsSci && fPosCalFrsSci->GetEntriesFast() > 0)
    {
        Int_t numDet = -1;
        nHits = fPosCalFrsSci->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BFrsSciPosCalData* hittcal = dynamic_cast<R3BFrsSciPosCalData*>(fPosCalFrsSci->At(ihit));
            numDet = hittcal->GetDetector();
            if (numDet > fNumDet)
            {
                R3BLOG(warn, "FrsSci detector id:" << numDet << " is out of range!");
                continue;
            }
            if (multSci2[numDet - 1] >= MAXMULT)
                continue;
            PosSci2_m1[numDet - 1][multSci2[numDet - 1]] = hittcal->GetCalPosMm();
            TimeSci2_m1[numDet - 1][multSci2[numDet - 1]] = hittcal->GetRawTimeNs();
            TimeSci2wTref_m1[numDet - 1][multSci2[numDet - 1]] = hittcal->GetRawTimeNsWithTref();
            multSci2[numDet - 1]++;
        } // --- end of loop over hit data --- //
    }

    if (fHeader->GetExpId() == 91 || fHeader->GetExpId() == 118)
    {
        // --- read hit from LOSCal data --- //
        if (fCalLos && fCalLos->GetEntriesFast() > 0)
        {
            nHits = fCalLos->GetEntriesFast();
            if (nHits == 1 || (fCalLosTrig && fCalLosTrig->GetEntriesFast() == 1))
            {
                for (Int_t ihit = 0; ihit < nHits; ihit++)
                {
                    auto* hittcal = dynamic_cast<R3BLosCalData*>(fCalLos->At(ihit));
                    auto* hittcaltrig = dynamic_cast<R3BLosCalData*>(fCalLosTrig->At(0));
                    auto numDet = hittcal->GetDetector();
                    if (multLos[numDet - 1] >= MAXMULT)
                        break;
                    auto los_tdiff = fTimeStitch->GetTime(
                        hittcaltrig->GetTimeV_ns(0) - hittcal->GetTimeV_ns(fLosRefCh), "vftx", "vftx");
                    if (los_tdiff < fLosCalTrig_Low || los_tdiff > fLosCalTrig_High)
                    {
                        continue;
                    }
                    timeLosV[numDet - 1][multLos[numDet - 1]] =
                        hittcal->GetMeanTimeVFTX() - hittcal->GetTimeV_ns(fLosRefCh);
                    // posLosX_cm[numDet - 1][multLos[numDet - 1]] = hittcal->GetX_cm();
                    multLos[numDet - 1]++;
                } // --- end of loop over hit data --- //
            }
        }
    }
    else
    {
        // --- read hit from LOS data --- //
        if (fHitLos && fHitLos->GetEntriesFast() > 0)
        {
            Int_t numDet = -1;

            nHits = fHitLos->GetEntriesFast();

            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BLosHitData* hittcal = dynamic_cast<R3BLosHitData*>(fHitLos->At(ihit));
                numDet = hittcal->GetDetector();

                if (multLos[numDet - 1] >= MAXMULT)
                    break;
                timeLosV[numDet - 1][multLos[numDet - 1]] = hittcal->GetTime();
                posLosX_cm[numDet - 1][multLos[numDet - 1]] = hittcal->GetX_cm();
                multLos[numDet - 1]++;
            } // --- end of loop over hit data --- //
        }
    }

    // Note: If the objective is to use only Los Z, FRSdata can
    // now easily be made multihit capable from this code itself,
    // just save the Los Z and calculate the Brho and A/Q value here
    // no need to run the separate R3BAnalysisIncomingID task which
    // can currently take only single hits because the Z information
    // is generally taken from other detectors
    // -Nikhil
    Double_t good_beta = NAN, good_pos_s2 = NAN, good_pos_los = NAN, good_tof = NAN;

    for (int i = 0; i < fNumDet; i++)
    {
        if (fUseFrsSci && (fStoId != 0 || fStaId != i + 1))
        {
            continue;
        }
        for (Int_t i_L = 0; i_L < multLos[0]; i_L++) // Here we assume there is only one LOS detector
        {
            Int_t num_tof_candidates = 0;
            for (Int_t i_2 = 0; i_2 < multSci2[i]; i_2++)
            {
                if (fUseTref)
                {
                    ToFraw_m1 = fTimeStitch->GetTime(fHeader->GetTStart() - TimeSci2wTref_m1[i][i_2], "vftx", "vftx");
                }
                else if (fHeader->GetExpId() == 91 || fHeader->GetExpId() == 118)
                {
                    ToFraw_m1 =
                        fTimeStitch->GetTime(timeLosV[0][i_L] - TimeSci2wTref_m1[i][i_2],
                                             "vftx",
                                             "vftx"); // Ref time is subtracted already on timeLosV for this case
                }
                else
                {
                    ToFraw_m1 = fTimeStitch->GetTime(timeLosV[0][i_L] - TimeSci2_m1[i][i_2], "vftx", "vftx");
                    if (ToFraw_m1 > 0. && fHeader->GetExpId() == 515)
                        ToFraw_m1 = ToFraw_m1 - 40960.;
                }
                Velo_m1 = 1. / (fTof2InvV_p0->GetAt(i) +
                                fTof2InvV_p1->GetAt(i) * (fToFoffset->GetAt(i) + ToFraw_m1)); // [m/ns]
                Beta_m1 = Velo_m1 / (TMath::C() / pow(10, 9));
                // Select good ToF hit with gating beta
                if (Beta_m1 < fBeta_max && Beta_m1 > fBeta_min)
                {
                    good_beta = Beta_m1;
                    good_tof = ToFraw_m1;
                    good_pos_s2 = PosSci2_m1[i][i_2];
                    good_pos_los = posLosX_cm[0][i_L];
                    num_tof_candidates++;
                }
            }
            if (num_tof_candidates == 1)
            {
                AddData(fStaId,
                        fStoId,
                        0.,
                        0.,
                        good_beta,
                        0.,
                        good_pos_s2,
                        good_pos_los,
                        good_tof); // NaN indicator of only one
                if (!fUseMultHit)
                    break;
            }
            if (num_tof_candidates == 0 && fHeader->GetExpId() == 509)
            {
                Int_t num_tof_ch = 0;
                for (Int_t nCh = 0; nCh < 2; nCh++)
                {
                    if (!(multSci2Tcal[i * 3 + 2] == 1))
                        break;
                    for (Int_t i0 = 0; i0 < multSci2Tcal[i * 3 + nCh]; i0++)
                    {
                        if (fUseTref)
                        {
                            ToFraw_m1 = fTimeStitch->GetTime(
                                fHeader->GetTStart() - (TimeSci2_tcal[i * 3 + nCh][i0] - TimeSci2_tcal[i * 3 + 2][0]),
                                "vftx",
                                "vftx");
                        }
                        else
                        {
                            ToFraw_m1 = fTimeStitch->GetTime(
                                timeLosV[i][i_L] - (TimeSci2_tcal[i * 3 + nCh][i0]), "vftx", "vftx");
                        }
                        Velo_m1 = 1. / (fTof2InvV_p0->GetAt(i) +
                                        fTof2InvV_p1->GetAt(i) * (fToFoffset->GetAt(i) + ToFraw_m1)); // [m/ns]
                        Beta_m1 = Velo_m1 / (TMath::C() / pow(10, 9));
                        // Select good ToF hit with gating beta
                        if (Beta_m1 < fBeta_max && Beta_m1 > fBeta_min)
                        {
                            good_beta = Beta_m1;
                            good_tof = ToFraw_m1;
                            good_pos_los = posLosX_cm[i][i_L];
                            num_tof_ch++;
                        }
                    }
                    if (num_tof_ch > 0)
                        break;
                }
                if (num_tof_ch == 1)
                {
                    AddData(1, 2, 0., 0., good_beta, 0., NAN, good_pos_los,
                            good_tof); // NaN indicator of only one
                    // S2 pmt present
                    if (!fUseMultHit)
                        break;
                }
            }             // end: if (num_tof_candidates == 0 && fHeader->GetExpId() == 509)
            if (fUseTref) // Tref defines when the trigger time is. And only one hit is recorded.
                break;
        } // End of loop with LOS mult
    }

    // This part can be done by FrsSciTCal2Cal and resulting as FrsSciTofCalData.
    /*
    if (fUseFrsSci && fStoId > 0 && fStaId > 0)
    { // Use FrsSci for both Start and Stop
        for (Int_t i_stop = 0; i_stop < multSci2[fStoId - 1]; i_stop++)
        {
            for (Int_t i_2 = 0; i_2 < multSci2[fStaId - 1]; i_2++)
            {
                if (fUseTref)
                {
                    ToFraw_m1 =
                        fTimeStitch->GetTime(fHeader->GetTStart() - TimeSci2wTref_m1[fStaId - 1][i_2], "vftx", "vftx");
                } // Enable only when the FrsSciProvideTStart is ready
                else
                {
                    ToFraw_m1 = fTimeStitch->GetTime(
                        TimeSci2_m1[fStoId - 1][i_stop] - TimeSci2_m1[fStaId - 1][i_2], "vftx", "vftx");
                }
                Velo_m1 =
                    1. / (fTof2InvV_p0->GetAt(fStaId - 1) +
                          fTof2InvV_p1->GetAt(fStaId - 1) * (fToFoffset->GetAt(fStaId - 1) + ToFraw_m1)); // [m/ns]
                Beta_m1 = Velo_m1 / (TMath::C() / pow(10, 9));

                // Select good ToF hit with gating beta
                if (Beta_m1 < fBeta_max && Beta_m1 > fBeta_min)
                {
                    good_beta = Beta_m1;
                    good_tof = ToFraw_m1;
                    good_pos_s2 = PosSci2_m1[fStaId - 1][i_2];
                    good_pos_los = PosSci2_m1[fStoId - 1][i_2]; // Let's keep calling stop det as los
                    num_tof_candidates++;
                    AddData(fStaId, fStoId, 0., 0., good_beta, 0., good_pos_s2, good_pos_los, good_tof);
                    if (!fUseMultHit)
                        break;
                }
            }
        }
    }
    */
}

void R3BIncomingBeta::FinishEvent()
{
    R3BLOG(debug1, "Clearing Loaded Structures");
    if (fHitSci2)
    {
        fHitSci2->Clear();
    }
    if (fHeader->GetExpId() == 91 || fHeader->GetExpId() == 118)
    {
        if (fCalLos)
        {
            fCalLos->Clear();
        }
        if (fCalLosTrig)
        {
            fCalLosTrig->Clear();
        }
    }
    if (fHitLos)
    {
        fHitLos->Clear();
    }
    if (fFrsDataCA)
    {
        fFrsDataCA->Clear();
    }
}

void R3BIncomingBeta::Reset() {}

// -----   Private method AddData  --------------------------------------------
R3BFrsData* R3BIncomingBeta::AddData(Int_t StaId,
                                     Int_t StoId,
                                     Double_t z,
                                     Double_t aq,
                                     Double_t beta,
                                     Double_t brho,
                                     Double_t xs2,
                                     Double_t xc,
                                     Double_t tof)
{
    // It fills the R3BFrsData
    TClonesArray& clref = *fFrsDataCA;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BFrsData(StaId, StoId, z, aq, beta, brho, xs2, xc, tof);
}

ClassImp(R3BIncomingBeta);
