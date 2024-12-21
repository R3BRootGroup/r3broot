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

#include "R3BTargetDigitizer.h"
#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRuntimeDb.h"
#include "TClonesArray.h"

// includes for modeling
#include "TGeoBBox.h"
#include "TGeoCompositeShape.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMatrix.h"
#include "TGeoMedium.h"
#include "TGeoShapeAssembly.h"
#include "TParticle.h"
#include "TVirtualMC.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TRandom.h"
#include "TVector3.h"
#include <string>

#include "R3BMCTrack.h"
#include "R3BTraPoint.h"

#define SQR(x) ((x) * (x))
// #define U_MEV_C2 931.494
// #define MASS_PROTON_MEV_C2 938.272046 //938.279
// #define MASS_17NE_MEV_C2   29796.65565 //Ar //15846.9
// #define MASS_15O_MEV_C2    28860.78 //Cl //13971.3

// Updated values by Felix: wikipedia and http://ie.lbl.gov/toi2003/MassSearch.asp
#define U_MEV_C2 931.494061
#define MASS_PROTON_MEV_C2 938.272032973
#define MASS_17NE_MEV_C2 15846.7485478
#define MASS_16F_MEV_C2 14909.9864966 // GS resonance 535keV above 15O+p
#define MASS_15O_MEV_C2 13971.1785118

R3BTargetDigitizer::R3BTargetDigitizer()
    : FairTask("R3B Target Digitization scheme ")
{
}

R3BTargetDigitizer::~R3BTargetDigitizer() {}

void R3BTargetDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fTargetDigiPar = dynamic_cast<R3BTargetDigiPar*>(rtdb->getContainer("R3BTargetDigiPar"));

    if (fTargetDigiPar)
    {
        LOG(info) << "R3BTargetDigitizer::SetParContainers() ";
        LOG(info) << "Container R3BTargetDigiPar  loaded ";
    }
}

InitStatus R3BTargetDigitizer::Init()
{

    //  LOG(info)<<"Init ";
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fTargetPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("TraPoint"));
    fTargetMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    fMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array TraDigi
    fTargetDigi = new TClonesArray("R3BTargetDigi", 1000);
    ioman->Register("TargetDigi", "Digital response in Target", fTargetDigi, kTRUE);

    eventNoTra = 0;

    // Initialise control histograms
    ExEnIn_his = new TH1F("ExEnIn_his", "ExEnIn_his", 110, 0., 22.);

    return kSUCCESS;
}

void R3BTargetDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoTra += 1;
    if (eventNoTra / 1000. == (int)eventNoTra / 1000.)
        LOG(info) << "Event #: " << eventNoTra - 1;

    Int_t nentriesTra = fTargetPoints->GetEntries();
    Int_t nentries = fMCTrack->GetEntries();

    Int_t TrackIdTra = 0;
    //     Int_t TrackId=0;

    Double_t x0 = nan("");
    Double_t y0 = nan("");
    Double_t t0 = nan("");

    Int_t ss03_smul;
    Int_t ss03_kmul;
    Int_t ss06_smul;
    Int_t ss06_kmul;

    Double_t Pxf = 0;
    Double_t Pyf = 0;
    Double_t Pzf = 0;

    Double_t Pxp1 = 0;
    Double_t Pyp1 = 0;
    Double_t Pzp1 = 0;

    Double_t Pf_tot = nan("");
    Double_t Pp1_tot = nan("");

    Double_t f_beta = nan("");
    Double_t p1_beta = nan("");

    Double_t p2 = nan("");
    Double_t E2 = nan("");
    Double_t estar = nan("");

    //******************** Target **************************//

    for (Int_t l = 0; l < nentries; l++)
    {
        //   LOG(info)<<"entries "<<l;

        R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fMCTrack->At(l));

        Int_t PID = aTrack->GetPdgCode();
        Int_t mother = aTrack->GetMotherId();

        Double_t fX = aTrack->GetStartX();
        Double_t fY = aTrack->GetStartY();
        Double_t fZ = aTrack->GetStartZ();
        Double_t fT = aTrack->GetStartT();

        Double_t Px = aTrack->GetPx();
        Double_t Py = aTrack->GetPy();
        Double_t Pz = aTrack->GetPz();

        if (mother < 0 && PID == 1000080150)
        {
            Pxf = (Px * 1000);
            Pyf = (Py * 1000);
            Pzf = (Pz * 1000);
            //   LOG(info)<<"In "<<"Pxf "<<Pxf<<" Pyf "<<Pyf<<" Pzf "<<Pzf;
            Pf_tot = sqrt(SQR(Pxf) + SQR(Pyf) + SQR(Pzf));
            //   f_beta=0.7579865;
            f_beta = sqrt((SQR(Pf_tot)) / ((SQR(MASS_15O_MEV_C2)) + (SQR(Pf_tot))));
            //   LOG(info)<<"In "<<"Pf_tot "<<Pf_tot<<" f_beta "<<f_beta;
        }

        if (mother < 0 && PID == 2212 /*&& l==1*/)
        {
            Pxp1 = (Px * 1000);
            Pyp1 = (Py * 1000);
            Pzp1 = (Pz * 1000);
            //   LOG(info)<<"In "<<"Pxp1 "<<Pxp1<<" Pyp1 "<<Pyp1<<" Pzp1 "<<Pzp1;
            Pp1_tot = sqrt(SQR(Pxp1) + SQR(Pyp1) + SQR(Pzp1));
            //   p1_beta=0.7579865;
            p1_beta = sqrt((SQR(Pp1_tot)) / ((SQR(MASS_PROTON_MEV_C2)) + (SQR(Pp1_tot))));
            //   LOG(info)<<"In "<<"Pp1_tot "<<Pp1_tot<<" p1_beta "<<p1_beta;
        }

        if (mother < 0)
        {
            // x0=fX;
            // y0=fY;

            // x0=(((fX - 0.0)*1) - ((fZ - 0.0)*0.0));
            // y0=(fY - 0.0);

            // Change to have lab coordinates instead (using the tracker offsets of the target):
            x0 = (((fX + 0.202437) * 1) - ((fZ - 0.0) * 0.0));
            y0 = (fY + 0.077698);

            t0 = fT;
        }
    }

    Double_t px = SQR(Pxp1 + Pxf);
    Double_t py = SQR(Pyp1 + Pyf);
    Double_t pz = SQR(Pzp1 + Pzf);
    p2 = px + py + pz; // MeV^2/c^2
                       //   LOG(info)<<"In "<<"px "<<px<<" py "<<py<<" pz "<<pz<<" sqrt_p2 "<<sqrt(p2);

    Double_t E_f = 1.0 / sqrt(1 - f_beta * f_beta) * MASS_15O_MEV_C2;       // *c
    Double_t E_p1 = 1.0 / sqrt(1 - p1_beta * p1_beta) * MASS_PROTON_MEV_C2; // *c
    E2 = (E_f + E_p1) * (E_f + E_p1);                                       // MeV^2/c^2
    //   LOG(info)<<"In "<<"E_f "<<E_f<<" E_p1 "<<E_p1<<" sqrt_E2 "<<sqrt(E2);

    //  sqrt(MeV^2/c^2 - MeV^2/c^2)=MeV/c
    // estar=sqrt(E2-p2)-MASS_17NE_MEV_C2; // *c2
    estar = sqrt(E2 - p2) - (MASS_15O_MEV_C2 + MASS_PROTON_MEV_C2); // *c2

    //   LOG(info)<<"Estar In "<<Estar;

    ExEnIn_his->Fill(estar);

    //******************** SSTs **************************//

    ss03_smul = 0;
    ss03_kmul = 0;
    ss06_smul = 0;
    ss06_kmul = 0;

    for (Int_t l = 0; l < nentriesTra; l++)
    {
        //   LOG(info)<<"entries "<<l;

        R3BTraPoint* Tra_obj = dynamic_cast<R3BTraPoint*>(fTargetPoints->At(l));

        //     Int_t DetID = Tra_obj->GetDetectorID();

        TrackIdTra = Tra_obj->GetTrackID();
        R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fTargetMCTrack->At(TrackIdTra));
        Int_t PID = aTrack->GetPdgCode();
        Int_t mother = aTrack->GetMotherId();

        Double_t fZ_In = Tra_obj->GetZIn();
        Double_t fZ_Out = Tra_obj->GetZOut();
        Double_t fZ = ((fZ_In + fZ_Out) / 2);

        if (PID == 1000080150 && mother < 0)
        {

            // if (DetID==15)
            // if (DetID==21)
            // if (DetID==23)
            if (fZ < 12)
            {
                ss03_smul++;
                ss03_kmul++;
            }
            // if (DetID==16)
            // if (DetID==22)
            // if (DetID==24)
            if (fZ > 12)
            {
                ss06_smul++;
                ss06_kmul++;
            }
        }

        if (PID == 2212 && mother < 0)
        {

            // if (DetID==15)
            // if (DetID==21)
            // if (DetID==23)
            if (fZ < 12)
            {
                ss03_smul++;
                ss03_kmul++;
            }
            // if (DetID==16)
            // if (DetID==22)
            // if (DetID==24)
            if (fZ > 12)
            {
                ss06_smul++;
                ss06_kmul++;
            }
        }
    }

    AddHit(ss03_smul,
           ss03_kmul,
           ss06_smul,
           ss06_kmul,
           x0,
           y0,
           t0,
           estar,
           Pxf,
           Pyf,
           Pzf,
           Pf_tot,
           f_beta,
           Pxp1,
           Pyp1,
           Pzp1,
           Pp1_tot,
           p1_beta);
}
// -------------------------------------------------------------------------

void R3BTargetDigitizer::Reset()
{
    // Clear the structure
    //   LOG(info) << " -I- Digit Reset() called ";

    if (fTargetDigi)
        fTargetDigi->Clear();
}

void R3BTargetDigitizer::Finish()
{
    // Write control histograms

    ExEnIn_his->Write();
}

R3BTargetDigi* R3BTargetDigitizer::AddHit(Double_t ss03_smul,
                                          Double_t ss03_kmul,
                                          Double_t ss06_smul,
                                          Double_t ss06_kmul,
                                          Double_t x0,
                                          Double_t y0,
                                          Double_t t0,
                                          Double_t estar,
                                          Double_t Pxf,
                                          Double_t Pyf,
                                          Double_t Pzf,
                                          Double_t Pf_tot,
                                          Double_t f_beta,
                                          Double_t Pxp1,
                                          Double_t Pyp1,
                                          Double_t Pzp1,
                                          Double_t Pp1_tot,
                                          Double_t p1_beta)
{
    TClonesArray& clref = *fTargetDigi;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTargetDigi(ss03_smul,
                                           ss03_kmul,
                                           ss06_smul,
                                           ss06_kmul,
                                           x0,
                                           y0,
                                           t0,
                                           estar,
                                           Pxf,
                                           Pyf,
                                           Pzf,
                                           Pf_tot,
                                           f_beta,
                                           Pxp1,
                                           Pyp1,
                                           Pzp1,
                                           Pp1_tot,
                                           p1_beta);
}

// R3BTraDigi* R3BTargetDigitizer::AddHit(
// return new(clref[size]) R3BTargetDigi();
//}

ClassImp(R3BTargetDigitizer)
