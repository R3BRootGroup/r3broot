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

#include "R3BTarget2pDigitizer.h"
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
// #define MASS_PROTON_MEV_C2 938.279
// #define MASS_17NE_MEV_C2   15846.9
// #define MASS_15O_MEV_C2    13971.3
////#define MASS_PROTON_MEV_C2 938.272046
////#define MASS_17NE_MEV_C2   29796.65565  //Ar
////#define MASS_15O_MEV_C2    27922.79434  //S

// Updated values by Felix: wikipedia and http://ie.lbl.gov/toi2003/MassSearch.asp
#define U_MEV_C2 931.494061
#define MASS_PROTON_MEV_C2 938.272032973
#define MASS_17NE_MEV_C2 15846.7485478
#define MASS_16F_MEV_C2 14909.9864966 // GS resonance 535keV above 15O+p
#define MASS_15O_MEV_C2 13971.1785118

R3BTarget2pDigitizer::R3BTarget2pDigitizer()
    : FairTask("R3B Target2p Digitization scheme ")
{
}

R3BTarget2pDigitizer::~R3BTarget2pDigitizer() {}

void R3BTarget2pDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fTarget2pDigiPar = dynamic_cast<R3BTarget2pDigiPar*>(rtdb->getContainer("R3BTarget2pDigiPar"));

    if (fTarget2pDigiPar)
    {
        LOG(info) << "-I- R3BTarget2pDigitizer::SetParContainers() ";
        LOG(info) << "-I- Container R3BTarget2pDigiPar  loaded ";
    }
}

InitStatus R3BTarget2pDigitizer::Init()
{

    //  LOG(info)<<"Init ";
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fTarget2pPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("TraPoint"));
    fTarget2pMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    fMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array TraDigi
    fTarget2pDigi = new TClonesArray("R3BTarget2pDigi", 1000);
    ioman->Register("Target2pDigi", "Digital response in Target2p", fTarget2pDigi, kTRUE);

    eventNoTra = 0;

    // Initialise control histograms
    ExEnIn_his = new TH1F("ExEnIn_his", "ExEnIn_his", 300, 0., 30.);

    return kSUCCESS;
}

void R3BTarget2pDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoTra += 1;
    if (eventNoTra / 1000. == (int)eventNoTra / 1000.)
        LOG(info) << "Event #: " << eventNoTra - 1;

    Int_t nentriesTra = fTarget2pPoints->GetEntries();
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

    Double_t Pxp2 = 0;
    Double_t Pyp2 = 0;
    Double_t Pzp2 = 0;

    Double_t Pf_tot = nan("");
    Double_t Pp1_tot = nan("");
    Double_t Pp2_tot = nan("");

    Double_t f_beta = nan("");
    Double_t p1_beta = nan("");
    Double_t p2_beta = nan("");

    Double_t p2 = nan("");
    Double_t E2 = nan("");
    Double_t estar = nan("");

    // Todo: do those properly, maybe 15O+2p mass
    Double_t in_beta = 0.7579865;
    Double_t Ne17mass = 15851.85984; // MeV/c2

    Int_t ppmul;

    //******************** Target **************************//

    ppmul = 0;

    for (Int_t l = 0; l < nentries; l++)
    {

        // if (l<4){
        //   LOG(info)<<"entries "<<l;
        //}
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
            //   if (mother<0 && PID==1000160300){   //Christoph 2p
            Pxf = (Px * 1000);
            Pyf = (Py * 1000);
            Pzf = (Pz * 1000);
            //   LOG(info)<<"In "<<"Pxf "<<Pxf<<" Pyf "<<Pyf<<" Pzf "<<Pzf;
            Pf_tot = sqrt(SQR(Pxf) + SQR(Pyf) + SQR(Pzf));
            //   f_beta=0.7579865;
            f_beta = sqrt((SQR(Pf_tot)) / ((SQR(MASS_15O_MEV_C2)) + (SQR(Pf_tot))));
            //   LOG(info)<<"In "<<"Pf_tot "<<Pf_tot<<" f_beta "<<f_beta;
        }

        if (mother < 0 && PID == 2212)
        {

            if (l == 1)
            {
                Pxp1 = (Px * 1000);
                Pyp1 = (Py * 1000);
                Pzp1 = (Pz * 1000);

                Pp1_tot = sqrt(SQR(Pxp1) + SQR(Pyp1) + SQR(Pzp1));
                p1_beta = sqrt((SQR(Pp1_tot)) / ((SQR(MASS_PROTON_MEV_C2)) + (SQR(Pp1_tot))));
            }

            if (l == 2)
            {
                Pxp2 = (Px * 1000);
                Pyp2 = (Py * 1000);
                Pzp2 = (Pz * 1000);

                Pp2_tot = sqrt(SQR(Pxp2) + SQR(Pyp2) + SQR(Pzp2));
                p2_beta = sqrt((SQR(Pp2_tot)) / ((SQR(MASS_PROTON_MEV_C2)) + (SQR(Pp2_tot))));
            }
            ppmul++;
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

    Double_t px = SQR(Pxp1 + Pxp2 + Pxf);
    Double_t py = SQR(Pyp1 + Pyp2 + Pyf);
    Double_t pz = SQR(Pzp1 + Pzp2 + Pzf);
    p2 = px + py + pz; // MeV^2/c^2
                       //   LOG(info)<<"In "<<"px "<<px<<" py "<<py<<" pz "<<pz<<" sqrt_p2 "<<sqrt(p2);

    Double_t E_f = 1.0 / sqrt(1 - f_beta * f_beta) * MASS_15O_MEV_C2;       // *c
    Double_t E_p1 = 1.0 / sqrt(1 - p1_beta * p1_beta) * MASS_PROTON_MEV_C2; // *c
    Double_t E_p2 = 1.0 / sqrt(1 - p2_beta * p2_beta) * MASS_PROTON_MEV_C2; // *c
    E2 = (E_f + E_p1 + E_p2) * (E_f + E_p1 + E_p2);                         // MeV^2/c^2
    //   LOG(info)<<"In "<<"E_f "<<E_f<<" E_p1 "<<E_p1<<" E_p2 "<<E_p2<<" sqrt_E2 "<<sqrt(E2);

    //  sqrt(MeV^2/c^2 - MeV^2/c^2)=MeV/c
    // estar=sqrt(E2-p2)-MASS_17NE_MEV_C2; // *c2
    estar = sqrt(E2 - p2) - (MASS_15O_MEV_C2 + 2.0 * MASS_PROTON_MEV_C2); // *c2

    //   LOG(info)<<"Estar In "<<Estar;

    ExEnIn_his->Fill(estar);

    Double_t Ne17_px = (Pxp1 + Pxp2 + Pxf);
    Double_t Ne17_py = (Pyp1 + Pyp2 + Pyf);
    Double_t Ne17_pz = (Pzp1 + Pzp2 + Pzf);
    Double_t Ne17_p = sqrt(p2);

    Double_t Indx = 0;
    Double_t Indy = 0;
    Double_t Indz = 1;

    Double_t In_vector = (Ne17mass * in_beta) / sqrt(1 - in_beta * in_beta);
    Double_t In_px = Indx * In_vector;
    Double_t In_py = Indy * In_vector;
    Double_t In_pz = Indz * In_vector;

    Double_t cosalpha = ((Ne17_px * In_px) + (Ne17_py * In_py) + (Ne17_pz * In_pz)) / (Ne17_p * In_vector);
    Double_t alpha = acos(cosalpha);

    //******************** SSTs **************************//

    ss03_smul = 0;
    ss03_kmul = 0;
    ss06_smul = 0;
    ss06_kmul = 0;

    for (Int_t l = 0; l < nentriesTra; l++)
    {
        //   LOG(info)<<"entries "<<l;

        R3BTraPoint* Tra_obj = dynamic_cast<R3BTraPoint*>(fTarget2pPoints->At(l));

        //     Int_t DetID = Tra_obj->GetDetectorID();

        TrackIdTra = Tra_obj->GetTrackID();
        R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fTarget2pMCTrack->At(TrackIdTra));
        Int_t PID = aTrack->GetPdgCode();
        Int_t mother = aTrack->GetMotherId();

        Double_t fZ_In = Tra_obj->GetZIn();
        Double_t fZ_Out = Tra_obj->GetZOut();
        Double_t fZ = ((fZ_In + fZ_Out) / 2);

        if (PID == 1000080150 && mother < 0)
        {
            //    if(PID==1000160300 && mother<0){  //Christoph 2p

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

            // if (DetID==15 && l==0)
            // if (DetID==21 && l==0)
            // if (DetID==23 && l==0)
            if (fZ < 12 && l == 0)
            {
                ss03_smul++;
                ss03_kmul++;
            }
            // if (DetID==15 && l>1)
            // if (DetID==21 && l>1)
            // if (DetID==23 && l>1)
            if (fZ < 12 && l > 0)
            {
                ss03_smul++;
                ss03_kmul++;
            }

            // if (DetID==16 && l==1)
            // if (DetID==22 && l==1)
            // if (DetID==24 && l==1)
            if (fZ > 12 && l == 1)
            {
                ss06_smul++;
                ss06_kmul++;
            }
            // if (DetID==16 && l>1)
            // if (DetID==22 && l>1)
            // if (DetID==24 && l>1)
            if (fZ > 12 && l > 1)
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
           p1_beta,
           // Pxp2,Pyp2,Pzp2,Pp2_tot,p2_beta);
           Pxp2,
           Pyp2,
           Pzp2,
           Pp2_tot,
           p2_beta,
           alpha,
           ppmul);
}
// -------------------------------------------------------------------------

void R3BTarget2pDigitizer::Reset()
{
    // Clear the structure
    //   LOG(info) << " -I- Digit Reset() called ";

    if (fTarget2pDigi)
        fTarget2pDigi->Clear();
}

void R3BTarget2pDigitizer::Finish()
{
    // Write control histograms

    ExEnIn_his->Write();
}

R3BTarget2pDigi* R3BTarget2pDigitizer::AddHit(Double_t ss03_smul,
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
                                              Double_t p1_beta,
                                              Double_t Pxp2,
                                              Double_t Pyp2,
                                              Double_t Pzp2,
                                              // Double_t Pp2_tot,Double_t p2_beta){
                                              Double_t Pp2_tot,
                                              Double_t p2_beta,
                                              Double_t alpha,
                                              Int_t ppmul)
{
    TClonesArray& clref = *fTarget2pDigi;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTarget2pDigi(ss03_smul,
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
                                             // Pxp1,Pyp1,Pzp1,Pp1_tot,p1_beta,Pxp2,Pyp2,Pzp2,Pp2_tot,p2_beta);
                                             Pxp1,
                                             Pyp1,
                                             Pzp1,
                                             Pp1_tot,
                                             p1_beta,
                                             Pxp2,
                                             Pyp2,
                                             Pzp2,
                                             Pp2_tot,
                                             p2_beta,
                                             alpha,
                                             ppmul);
}

// R3BTarget2pDigi* R3BTarget2pDigitizer::AddHit(
// return new(clref[size]) R3BTarget2pDigi();
//}

ClassImp(R3BTarget2pDigitizer)
