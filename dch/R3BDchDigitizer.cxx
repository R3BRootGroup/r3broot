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

#include "R3BDchDigitizer.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"

// includes for modeling
#include "TGeoManager.h"
#include "TParticle.h"
#include "TVirtualMC.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <iostream>
#include <string>

#include "R3BDchPoint.h"
#include "R3BMCTrack.h"

// hardcoded detector centres wrt target, from Ralf's Tracker: world coords (wrt target) of detector centres, in
// centimeters. PDC1:    x0=-123.219446 y0=3.597104 z0=444.126271 PDC2:    x0=-167.015888 y0=1.016917 z0=535.093884 PDC1
// angle:    x0=0.000000 y0=31.000000 z0=8.880000 PDC2 angle:    x0=0.000000 y0=31.000000 z0=-9.350000
#define PDC1_X0 -123.219446
#define PDC1_Y0 3.597104
#define PDC1_Z0 444.126271
#define PDC2_X0 -167.015888
#define PDC2_Y0 1.016917
#define PDC2_Z0 535.093884
#define PDC1_Aparm 31.000000
#define PDC1_Atilt 8.880000
// #define PDC1_Atilt	-8.880000
#define PDC2_Aparm 31.000000
#define PDC2_Atilt -9.350000
// #define PDC2_Atilt	9.350000

// Dead wire areas, determined from s318 DHIT data:
#define PDC1_X_DEAD_1_LO 22.16 // cm, DHIT level, >
#define PDC1_X_DEAD_1_HI 23.55 // cm, DHIT level, <=
// #define PDC1_Y_DEAD_1_LO	13.85	//cm, DHIT level, >
// #define PDC1_Y_DEAD_1_HI	15.24	//cm, DHIT level, <=
// #define PDC1_Y_DEAD_2_LO	54.73	//cm, DHIT level, >
// #define PDC1_Y_DEAD_2_HI	55.42	//cm, DHIT level, <=
// #define PDC1_Y_DEAD_3_LO	58.88	//cm, DHIT level, >
// #define PDC1_Y_DEAD_3_HI	60.97	//cm, DHIT level, <=
// #define PDC1_Y_DEAD_4_LO	64.41	//cm, DHIT level, >
// #define PDC1_Y_DEAD_4_HI	66.50	//cm, DHIT level, <=
//  update: HIT level values, translated to dhit
#define PDC1_Y_DEAD_1_LO 14.3 // cm, DHIT level, >
#define PDC1_Y_DEAD_1_HI 15.5 // cm, DHIT level, <=
#define PDC1_Y_DEAD_2_LO 54.8 // cm, DHIT level, >
#define PDC1_Y_DEAD_2_HI 56.2 // cm, DHIT level, <=
#define PDC1_Y_DEAD_3_LO 56.5 // cm, DHIT level, >
#define PDC1_Y_DEAD_3_HI 57.1 // cm, DHIT level, <=
#define PDC1_Y_DEAD_4_LO 60.7 // cm, DHIT level, >
#define PDC1_Y_DEAD_4_HI 62.7 // cm, DHIT level, <=
#define PDC1_Y_DEAD_5_LO 66.2 // cm, DHIT level, >
#define PDC1_Y_DEAD_5_HI 66.8 // cm, DHIT level, <=

#define PDC2_X_DEAD_1_LO 13.84 // cm, DHIT level, >
#define PDC2_X_DEAD_1_HI 15.26 // cm, DHIT level, <=

// Size of active area
// #define PDC_L_X		100.0	//length of active area
// #define PDC_L_Y		80.0	//length of active area
#define PDC_L_X 99.07 // x length of active area
#define PDC_L_Y 77.60 // y length of active area

R3BDchDigitizer::R3BDchDigitizer()
    : FairTask("R3B Dch Digitization scheme ")
{
    detectorResolution = 0;
}
R3BDchDigitizer::R3BDchDigitizer(int detres)
    : FairTask("R3B Dch Digitization scheme ")
{
    detectorResolution = detres; // 0=all off, 1= all on, ...
}

R3BDchDigitizer::~R3BDchDigitizer() {}

void R3BDchDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fDchDigiPar = dynamic_cast<R3BDchDigiPar*>(rtdb->getContainer("R3BDchDigiPar"));

    if (fDchDigiPar)
    {
        LOG(info) << "R3BDchDigitizer::SetParContainers() ";
        LOG(info) << "Container R3BDchDigiPar  loaded ";
    }
}

InitStatus R3BDchDigitizer::Init()
{

    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fDchPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("DCHPoint"));
    fDchMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array DchDigi
    fDchDigi = new TClonesArray("R3BDchDigi", 1000);
    ioman->Register("DchDigi", "Digital response in Dch", fDchDigi, kTRUE);

    eventNoDch = 0;

    // Initialise control histograms
    DCH1Px = new TH1F("DCH1Px", "DCH1Px", 500, -1., 0.2);
    DCH1Px->GetXaxis()->SetTitle("Momentum");
    DCH1Px->GetYaxis()->SetTitle("Counts");

    DCH2Px = new TH1F("DCH2Px", "DCH2Px", 500, -1., 0.2);
    DCH2Px->GetXaxis()->SetTitle("Momentum");
    DCH2Px->GetYaxis()->SetTitle("Counts");

    DCH1Py = new TH1F("DCH1Py", "DCH1Py", 500, -0.15, 0.15);
    DCH1Py->GetXaxis()->SetTitle("Momentum");
    DCH1Py->GetYaxis()->SetTitle("Counts");

    DCH2Py = new TH1F("DCH2Py", "DCH2Py", 500, -0.15, 0.15);
    DCH2Py->GetXaxis()->SetTitle("Momentum");
    DCH2Py->GetYaxis()->SetTitle("Counts");

    DCH1Pz = new TH1F("DCH1Pz", "DCH1Pz", 500, -0.2, 1.2);
    DCH1Pz->GetXaxis()->SetTitle("Momentum");
    DCH1Pz->GetYaxis()->SetTitle("Counts");

    DCH2Pz = new TH1F("DCH2Pz", "DCH2Pz", 500, -0.2, 1.2);
    DCH2Pz->GetXaxis()->SetTitle("Momentum");
    DCH2Pz->GetYaxis()->SetTitle("Counts");

    DCH1X = new TH1F("DCH1X", "DCH1X", 2400, -120, 120);
    DCH1X->GetXaxis()->SetTitle("Position");
    DCH1X->GetYaxis()->SetTitle("Counts");

    DCH2X = new TH1F("DCH2X", "DCH2X", 2400, -120, 120);
    DCH2X->GetXaxis()->SetTitle("Position");
    DCH2X->GetYaxis()->SetTitle("Counts");

    DCH1Y = new TH1F("DCH1Y", "DCH1Y", 2000, -100, 100);
    DCH1Y->GetXaxis()->SetTitle("Position");
    DCH1Y->GetYaxis()->SetTitle("Counts");

    DCH2Y = new TH1F("DCH2Y", "DCH2Y", 2000, -100, 100);
    DCH2Y->GetXaxis()->SetTitle("Position");
    DCH2Y->GetYaxis()->SetTitle("Counts");

    DCH1elosshis = new TH1F("DCH1elosshis", "DCH1elosshis", 500, 0., 0.1);

    DCH2elosshis = new TH1F("DCH2elosshis", "DCH2elosshis", 500, 0., 0.1);

    TrackPx = new TH1F("TrackPx", "TrackPx", 500, -0.15, 0.15);
    TrackPx->GetXaxis()->SetTitle("Momentum");
    TrackPx->GetYaxis()->SetTitle("Counts");

    TrackPy = new TH1F("TrackPy", "TrackPy", 500, -0.15, 0.15);
    TrackPy->GetXaxis()->SetTitle("Momentum");
    TrackPy->GetYaxis()->SetTitle("Counts");

    TrackPz = new TH1F("TrackPz", "TrackPz", 500, 0.5, 1.5);
    TrackPz->GetXaxis()->SetTitle("Momentum");
    TrackPz->GetYaxis()->SetTitle("Counts");

    TrackPxVSDCH1Px = new TH2F("TrackPxVSDCH1Px", "TrackPxVSDCH1Px", 500, -0.8, -0.2, 500, -0.15, 0.15);
    TrackPxVSDCH1Px->GetXaxis()->SetTitle("DCH1Px");
    TrackPxVSDCH1Px->GetYaxis()->SetTitle("TrackPx");

    return kSUCCESS;
}

void R3BDchDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoDch += 1;

    Int_t nentriesDch = fDchPoints->GetEntries();
    LOG(info) << nentriesDch;

    Int_t TrackId = 0;

    Double_t total_energy_dch1 = 0.;
    Double_t total_energy_dch2 = 0.;

    Double_t DCH1eloss;
    Double_t DCH2eloss;

    Int_t pd1mul = 0;
    Int_t pd2mul = 0;
    Double_t Pdx1_p1 = -1000.; // so they stay at that when not set otherwise. -> Easy to disentangle
    Double_t Pdx2_p1 = -1000.;
    Double_t Pdy1_p1 = -1000.;
    Double_t Pdy2_p1 = -1000.;

    //********************** Energy loss ***************************//

    for (Int_t l = 0; l < nentriesDch; l++)
    {

        R3BDchPoint* dch_obj = dynamic_cast<R3BDchPoint*>(fDchPoints->At(l));

        Int_t DetID = dch_obj->GetDetectorID();

        if (DetID == 0)
        {
            DCH1eloss = dch_obj->GetEnergyLoss() * 1000.;
            total_energy_dch1 += DCH1eloss;
            DCH1elosshis->Fill(total_energy_dch1);
        }

        if (DetID == 1)
        {
            DCH2eloss = dch_obj->GetEnergyLoss() * 1000.;
            total_energy_dch2 += DCH2eloss;
            DCH2elosshis->Fill(total_energy_dch2);
        }
    }

    //******************** DCH **************************//

    for (Int_t l = 0; l < nentriesDch; l++)
    {

        R3BDchPoint* dch_obj = dynamic_cast<R3BDchPoint*>(fDchPoints->At(l));

        Int_t DetID = dch_obj->GetDetectorID();
        Double_t fPx = dch_obj->GetPx();
        Double_t fPy = dch_obj->GetPy();
        Double_t fPz = dch_obj->GetPz();

        Double_t fX_Local_In = dch_obj->GetXLocalIn();
        Double_t fY_Local_In = dch_obj->GetYLocalIn();
        Double_t fX_Local_Out = dch_obj->GetXLocalOut();
        Double_t fY_Local_Out = dch_obj->GetYLocalOut();

        //     Double_t fX_Global_In = dch_obj->GetXIn();
        //     Double_t fY_Global_In = dch_obj->GetYIn();
        //     Double_t fZ_Global_In = dch_obj->GetZIn();
        //     Double_t fX_Global_Out = dch_obj->GetXOut();
        //     Double_t fY_Global_Out = dch_obj->GetYOut();
        //     Double_t fZ_Global_Out = dch_obj->GetZOut();
        TrackId = dch_obj->GetTrackID();
        R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fDchMCTrack->At(TrackId));
        Int_t PID = aTrack->GetPdgCode();
        Int_t mother = aTrack->GetMotherId();
        Double_t fPx_track = aTrack->GetPx();
        Double_t fPy_track = aTrack->GetPy();
        Double_t fPz_track = aTrack->GetPz();

        Double_t fX_Local_sim = ((fX_Local_In + fX_Local_Out) / 2);
        Double_t fY_Local_sim = ((fY_Local_In + fY_Local_Out) / 2);

        // Need modifications
        //     Double_t fX_Global = ((fX_Global_In + fX_Global_Out)/2);
        //     Double_t fY_Global = ((fY_Global_In + fY_Global_Out)/2);
        //     Double_t fZ_Global = ((fZ_Global_In + fZ_Global_Out)/2);

        // hardcoded detector centres wrt target, from Ralf's Tracker:
        // PDC1:    x0=-123.219446 y0=3.597104 z0=444.126271
        // PDC2:    x0=-167.015888 y0=1.016917 z0=535.093884

        Double_t fX_Local = 0.; // initialisation, values will be PDC-specific
        Double_t fY_Local = 0.;
        //     Double_t fZ_Local = 0.;

        // if(PID==2212 && mother<0) {
        if (PID == 1000010010 && mother < 0)
        {
            if (DetID == 0)
            {
                // Check tilt direction, positive or negative angle?! Make it consistent with rotations in R3BDch.cxx
                // !!! using this manual calculation (based on the manual offsets above), a consistency with the s318
                // tracker is achieved. F. Wamers. fX_Local	=
                // (fX_Global-PDC1_X0)*cos(PDC1_Atilt*TMath::Pi()/180.)/cos(PDC1_Aparm*TMath::Pi()/180.) -
                // (fY_Global-PDC1_Y0)*sin(PDC1_Atilt*TMath::Pi()/180.); fY_Local	=
                // (fX_Global-PDC1_X0)*sin(PDC1_Atilt*TMath::Pi()/180.)/cos(PDC1_Aparm*TMath::Pi()/180.) +
                // (fY_Global-PDC1_Y0)*cos(PDC1_Atilt*TMath::Pi()/180.);
                fX_Local = fX_Local_sim; // values identical to manual ones. both correct.
                fY_Local = fY_Local_sim;

                // Here, discard events where areas of dead wires were hit:
                // if ( (fX_Local+PDC_L_X/2>PDC1_X_DEAD_1_LO && fX_Local+PDC_L_X/2<=PDC1_X_DEAD_1_HI)
                if ((-fX_Local + PDC_L_X / 2 > PDC1_X_DEAD_1_LO &&
                     -fX_Local + PDC_L_X / 2 <= PDC1_X_DEAD_1_HI) //- for dhit-hit x-swap
                    || (fY_Local + PDC_L_Y / 2 > PDC1_Y_DEAD_1_LO && fY_Local + PDC_L_Y / 2 <= PDC1_Y_DEAD_1_HI) ||
                    (fY_Local + PDC_L_Y / 2 > PDC1_Y_DEAD_2_LO && fY_Local + PDC_L_Y / 2 <= PDC1_Y_DEAD_2_HI) ||
                    (fY_Local + PDC_L_Y / 2 > PDC1_Y_DEAD_3_LO && fY_Local + PDC_L_Y / 2 <= PDC1_Y_DEAD_3_HI) ||
                    (fY_Local + PDC_L_Y / 2 > PDC1_Y_DEAD_4_LO && fY_Local + PDC_L_Y / 2 <= PDC1_Y_DEAD_4_HI) ||
                    (fY_Local + PDC_L_Y / 2 > PDC1_Y_DEAD_5_LO && fY_Local + PDC_L_Y / 2 <= PDC1_Y_DEAD_5_HI))
                {
                    continue;
                }

                Pdx1_p1 = fX_Local; // not += !!!
                Pdy1_p1 = fY_Local; // not += !!!
                DCH1Px->Fill(fPx);
                DCH1Py->Fill(fPy);
                DCH1Pz->Fill(fPz);
                DCH1X->Fill(fX_Local + PDC_L_X / 2); // modified by additional addon
                DCH1Y->Fill(fY_Local + PDC_L_Y / 2);
                TrackPxVSDCH1Px->Fill(fPx, fPx_track);

                // resolutions:
                if (detectorResolution == 1)
                {
                    Pdx1_p1 = gRandom->Gaus(Pdx1_p1, 0.0278); // cm;
                    Pdy1_p1 = gRandom->Gaus(Pdy1_p1, 0.0215); // cm;
                    // Pdx1_p1 = gRandom->Gaus(Pdx1_p1,0.02);	//cm, Justynas default value;
                    // Pdy1_p1 = gRandom->Gaus(Pdy1_p1,0.02); 	//cm, Justynas default value;
                }

                pd1mul++;
            }

            if (DetID == 1)
            {
                // Check tilt direction, positive or negative angle?! Make it consistent with rotations in R3BDch.cxx
                // !!!
                // using this manual calculation (based on the manual offsets above), a consistency with the s318
                // tracker is achieved. F. Wamers.
                // fX_Local	= (fX_Global-PDC2_X0)*cos(PDC2_Atilt*TMath::Pi()/180.)/cos(PDC2_Aparm*TMath::Pi()/180.) -
                // (fY_Global-PDC2_Y0)*sin(PDC2_Atilt*TMath::Pi()/180.); fY_Local	=
                // (fX_Global-PDC2_X0)*sin(PDC2_Atilt*TMath::Pi()/180.)/cos(PDC2_Aparm*TMath::Pi()/180.) +
                // (fY_Global-PDC2_Y0)*cos(PDC2_Atilt*TMath::Pi()/180.);
                fX_Local = fX_Local_sim;
                fY_Local = fY_Local_sim;

                // Here, discard events where areas of dead wires were hit:
                // if ( (fX_Local+PDC_L_X/2>PDC2_X_DEAD_1_LO && fX_Local+PDC_L_X/2<=PDC2_X_DEAD_1_HI)
                if ((-fX_Local + PDC_L_X / 2 > PDC2_X_DEAD_1_LO &&
                     -fX_Local + PDC_L_X / 2 <= PDC2_X_DEAD_1_HI) //- for dhit-hit x-swap
                )
                    continue;

                Pdx2_p1 = fX_Local; // not += !!!
                Pdy2_p1 = fY_Local; // not += !!!
                DCH2Px->Fill(fPx);
                DCH2Py->Fill(fPy);
                DCH2Pz->Fill(fPz);
                DCH2X->Fill(fX_Local + PDC_L_X / 2);
                DCH2Y->Fill(fY_Local + PDC_L_Y / 2);

                // resolutions:
                if (detectorResolution == 1)
                {
                    Pdx2_p1 = gRandom->Gaus(Pdx2_p1, 0.0268); // cm;
                    Pdy2_p1 = gRandom->Gaus(Pdy2_p1, 0.0213); // cm;
                    // Pdx2_p1 = gRandom->Gaus(Pdx2_p1,0.02);	//cm, Justynas default value;
                    // Pdy2_p1 = gRandom->Gaus(Pdy2_p1,0.02); 	//cm, Justynas default value;
                }

                pd2mul++;
            }

            if (mother < 0)
            {
                TrackPx->Fill(fPx_track);
                TrackPy->Fill(fPy_track);
                TrackPz->Fill(fPz_track);
            }

        } // if primary proton

    } // loop over DCH points

    AddHit(pd1mul, Pdx1_p1, Pdy1_p1, pd2mul, Pdx2_p1, Pdy2_p1);
}
// -------------------------------------------------------------------------

void R3BDchDigitizer::Reset()
{
    if (fDchDigi)
        fDchDigi->Clear();
}

void R3BDchDigitizer::Finish()
{
    // Write control histograms
    DCH1Px->Write();
    DCH2Px->Write();

    DCH1Py->Write();
    DCH2Py->Write();

    DCH1Pz->Write();
    DCH2Pz->Write();

    DCH1X->Write();
    DCH2X->Write();

    DCH1Y->Write();
    DCH2Y->Write();

    DCH1elosshis->Write();
    DCH2elosshis->Write();

    TrackPx->Write();
    TrackPy->Write();
    TrackPz->Write();

    TrackPxVSDCH1Px->Write();
}

R3BDchDigi* R3BDchDigitizer::AddHit(Int_t pd1mul,
                                    Double_t Pdx1_p1,
                                    Double_t Pdy1_p1,
                                    Int_t pd2mul,
                                    Double_t Pdx2_p1,
                                    Double_t Pdy2_p1)
{
    TClonesArray& clref = *fDchDigi;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BDchDigi(pd1mul, Pdx1_p1, Pdy1_p1, pd2mul, Pdx2_p1, Pdy2_p1);
}

// R3BDchDigi* R3BDchDigitizer::AddHit(
// return new(clref[size]) R3BDchDigi();
//}

ClassImp(R3BDchDigitizer)
