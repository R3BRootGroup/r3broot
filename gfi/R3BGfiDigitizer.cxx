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

#include "R3BGfiDigitizer.h"
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

#include "R3BGfiPoint.h"
#include "R3BMCTrack.h"

using std::cout;
using std::endl;

R3BGfiDigitizer::R3BGfiDigitizer()
    : FairTask("R3B Gfi Digitization scheme ")
{
}

R3BGfiDigitizer::~R3BGfiDigitizer() {}

void R3BGfiDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fGfiDigiPar = dynamic_cast<R3BGfiDigiPar*>(rtdb->getContainer("R3BGfiDigiPar"));

    if (fGfiDigiPar)
    {
        cout << "-I- R3BGfiDigitizer::SetParContainers() " << endl;
        cout << "-I- Container R3BGfiDigiPar  loaded " << endl;
    }
}

InitStatus R3BGfiDigitizer::Init()
{

    //  cout<<"Init "<<endl;
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fGfiPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("GFIPoint"));
    fGfiMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array GfiDigi
    fGfiDigi = new TClonesArray("R3BGfiDigi", 1000);
    ioman->Register("GfiDigi", "Digital response in Gfi", fGfiDigi, kTRUE);

    eventNoGfi = 0;

    // Initialise control histograms

    return kSUCCESS;
}

void R3BGfiDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoGfi += 1;
    //     if(eventNoGfi/1000. == (int)eventNoGfi/1000.) cout<<"Event #: "<<eventNoGfi-1<<endl;

    Int_t nentriesGfi = fGfiPoints->GetEntries();

    Int_t TrackIdGfi = 0;

    Int_t gf1mul;
    Int_t gf2mul;
    Double_t gf1_01x = nan("");
    Double_t gf2_01x = nan("");

    //******************** GFI **************************//
    gf1mul = 0;
    gf2mul = 0;

    for (Int_t l = 0; l < nentriesGfi; l++)
    {
        //   cout<<"entries "<<l<<endl;

        R3BGfiPoint* Gfi_obj = dynamic_cast<R3BGfiPoint*>(fGfiPoints->At(l));

        //     Int_t DetID = Gfi_obj->GetDetectorID();
        Double_t fX_In = Gfi_obj->GetXIn();
        Double_t fX_Out = Gfi_obj->GetXOut();
        Double_t fZ_In = Gfi_obj->GetZIn();
        Double_t fZ_Out = Gfi_obj->GetZOut();
        TrackIdGfi = Gfi_obj->GetTrackID();
        R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fGfiMCTrack->At(TrackIdGfi));
        Int_t PID = aTrack->GetPdgCode();
        Int_t mother = aTrack->GetMotherId();

        Double_t fX = ((fX_In + fX_Out) / 2);
        Double_t fZ = ((fZ_In + fZ_Out) / 2);

        // offset - position of GFI1 detector x: -71.97331; z: +513.967775
        // offset - position of GFI2 detector x: -145.270039; z: +730.318956
        // angle = -16.7 degrees
        // cos(angle) = 0.957822495
        // sin(angle) = -0.28736052

        if (PID == 1000080150 && mother < 0)
        {
            //     if(PID==1000170310 && mother<0){  //Christoph 1p
            //     if(PID==1000160300 && mother<0){  //Christoph 2p

            if (fZ < 600) // no DetectorID for GFIs  this is GFI1
            {
                gf1_01x = (-(((fX + 72.164874) * 0.957822495) -
                             ((fZ - 513.910302) *
                              (-0.28736052)))); // x position swop (with -) - requirement for tracker //justyna new
                // gf1_01x=(-(((fX+71.97331) * 0.957822495) - ((fZ- 513.967775) * (-0.28736052))));//x position swop
                // (with -) - requirement for tracker
                //        gf1_01x=(-(((fX+88.55368) * 0.957822495) - ((fZ- 515.015418) * (-0.28736052)))); //Christoph
                //	cout<<"GFI1 - fragment detector "<<PID<<endl;
                gf1mul++;
            }
            if (fZ > 600) // GFI2
            {
                gf2_01x = (-(((fX + 145.212570) * 0.957822495) -
                             ((fZ - 730.336197) *
                              (-0.28736052)))); // x position swop (with -) - requirement for tracker //justyna new
                // gf2_01x=(-(((fX+145.270039) * 0.957822495) - ((fZ - 730.318956) * (-0.28736052))));//x position swop
                // (with -) - requirement for tracker
                //        gf2_01x=(-(((fX+149.791168) * 0.957822495) - ((fZ - 744.623074) * (-0.28736052))));
                //        //Christoph
                //	cout<<"GFI2 - fragment detector "<<PID<<endl;
                gf2mul++;
            }
        }
    }

    AddHit(gf1mul, gf1_01x, gf2mul, gf2_01x);
}
// -------------------------------------------------------------------------

void R3BGfiDigitizer::Reset()
{
    // Clear the structure
    //   cout << " -I- Digit Reset() called " << endl;

    if (fGfiDigi)
        fGfiDigi->Clear();
}

void R3BGfiDigitizer::Finish()
{
    // Write control histograms
}

R3BGfiDigi* R3BGfiDigitizer::AddHit(Int_t gf1mul, Double_t gf1_01x, Int_t gf2mul, Double_t gf2_01x)
{
    TClonesArray& clref = *fGfiDigi;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BGfiDigi(gf1mul, gf1_01x, gf2mul, gf2_01x);
}

// R3BGfiDigi* R3BGfiDigitizer::AddHit(
// return new(clref[size]) R3BGfiDigi();
//}

ClassImp(R3BGfiDigitizer)
