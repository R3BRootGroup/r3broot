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

#include "R3BTof2pDigitizer.h"
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
#include <iostream>
#include <string>

#include "R3BMCTrack.h"
#include "R3BTofPoint.h"

using std::cout;
using std::endl;

R3BTof2pDigitizer::R3BTof2pDigitizer()
    : FairTask("R3B Tof2p Digitization scheme ")
{
}

R3BTof2pDigitizer::~R3BTof2pDigitizer() {}

void R3BTof2pDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fTof2pDigiPar = dynamic_cast<R3BTof2pDigiPar*>(rtdb->getContainer("R3BTof2pDigiPar"));

    if (fTof2pDigiPar)
    {
        cout << "-I- R3BTof2pDigitizer::SetParContainers() " << endl;
        cout << "-I- Container R3BTofDigiPar  loaded " << endl;
    }
}

InitStatus R3BTof2pDigitizer::Init()
{

    //  cout<<"Init "<<endl;
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fTofPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("TOFPoint"));
    fTofMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array TofDigi
    fTof2pDigi = new TClonesArray("R3BTof2pDigi", 1000);
    ioman->Register("Tof2pDigi", "Digital response in Tof", fTof2pDigi, kTRUE);

    eventNoTof = 0;

    // Initialise control histograms

    return kSUCCESS;
}

void R3BTof2pDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoTof += 1;
    //     if(eventNoTof/1000. == (int)eventNoTof/1000.) cout<<"Event #: "<<eventNoTof-1<<endl;

    Int_t nentriesTof = fTofPoints->GetEntries();
    Int_t TrackIdTof = 0;

    // Double_t total_energy_tof=0.;
    // Double_t TOFeloss;

    Int_t tfmul;

    Double_t tfwx_p1;
    Double_t tfwy_p1;
    Double_t tfwt_p1;
    Double_t tfwpath_p1;
    Double_t tfwpx_p1;
    Double_t tfwpy_p1;
    Double_t tfwpz_p1;
    Double_t tfwe_p1;

    Double_t tfwx_p2;
    Double_t tfwy_p2;
    Double_t tfwt_p2;
    Double_t tfwpath_p2;
    Double_t tfwpx_p2;
    Double_t tfwpy_p2;
    Double_t tfwpz_p2;
    Double_t tfwe_p2;

    //********************* energy looping *********************//
    /*
       for (Int_t l=0;l<nentriesTof;l++){

         R3BTofPoint *tof_obj = (R3BTofPoint*) fTofPoints->At(l);

         TOFeloss = tof_obj->GetEnergyLoss()*1000;

         TrackIdTof = tof_obj->GetTrackID();
         R3BMCTrack *aTrack = (R3BMCTrack*) fTofMCTrack->At(TrackIdTof);
         Int_t PID = aTrack->GetPdgCode();


         if (PID==2212){
         total_energy_tof += TOFeloss;
         }


       }
    */

    //******************** TFW ********************//
    tfmul = 0;
    tfwx_p1 = 0;
    tfwy_p1 = 0;
    tfwt_p1 = 0;
    tfwpath_p1 = 0;
    tfwpx_p1 = 0;
    tfwpy_p1 = 0;
    tfwpz_p1 = 0;
    tfwe_p1 = 0;

    tfwx_p2 = 0;
    tfwy_p2 = 0;
    tfwt_p2 = 0;
    tfwpath_p2 = 0;
    tfwpx_p2 = 0;
    tfwpy_p2 = 0;
    tfwpz_p2 = 0;
    tfwe_p2 = 0;

    for (Int_t l = 0; l < nentriesTof; l++)
    {

        R3BTofPoint* tof_obj = dynamic_cast<R3BTofPoint*>(fTofPoints->At(l));

        TrackIdTof = tof_obj->GetTrackID();
        R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fTofMCTrack->At(TrackIdTof));
        Int_t PID = aTrack->GetPdgCode();
        Int_t mother = aTrack->GetMotherId();

        Double_t fX_in = tof_obj->GetXIn();
        Double_t fY_in = tof_obj->GetYIn();
        Double_t fZ_in = tof_obj->GetZIn();
        Double_t fX_out = tof_obj->GetXOut();
        Double_t fY_out = tof_obj->GetYOut();
        Double_t fZ_out = tof_obj->GetZOut();
        Double_t PX_in = tof_obj->GetPx();
        Double_t PY_in = tof_obj->GetPy();
        Double_t PZ_in = tof_obj->GetPz();
        Double_t PX_out = tof_obj->GetPxOut();
        Double_t PY_out = tof_obj->GetPyOut();
        Double_t PZ_out = tof_obj->GetPzOut();
        Double_t ftime = tof_obj->GetTime();
        Double_t flength = tof_obj->GetLength();
        Double_t TOFeloss = tof_obj->GetEnergyLoss() * 1000;

        Double_t fX = ((fX_in + fX_out) / 2);
        Double_t fY = ((fY_in + fY_out) / 2);
        Double_t fZ = ((fZ_in + fZ_out) / 2);

        Double_t PX = ((PX_in + PX_out) / 2);
        Double_t PY = ((PY_in + PY_out) / 2);
        Double_t PZ = ((PZ_in + PZ_out) / 2);

        // offset - position of detector x: -417.359574; y: +2.4; z: +960.777114
        // angle = -31 degrees
        // cos(angle) = 0.857167301
        // sin(angle) = -0.515038075

        if (PID == 2212 && mother < 0)
        {
            if (l == 0)
            {
                tfwx_p1 = (((fX + 417.359574) * 0.857167301) - ((fZ - 960.777114) * (-0.515038075)));
                //       tfwx_p1=(((fX+421.33683) * 0.857167301) - ((fZ - 958.387337) * (-0.515038075)));  //Christoph
                tfwy_p1 = (fY - 2.4);
                //       tfwy_p1=(fY-2.12);  //Christoph
                tfwt_p1 = ftime;
                tfwpath_p1 = flength;
                tfwpx_p1 = PX;
                tfwpy_p1 = PY;
                tfwpz_p1 = PZ;
                tfwe_p1 += TOFeloss;

                //       cout<<"TFW 1p X "<<fX<<" Y "<<fY<<" Z "<<fZ<<endl;
                //       cout<<"TFW - first proton"<<endl;
            }
            if (l > 0)
            {
                tfwx_p2 = (((fX + 417.359574) * 0.857167301) - ((fZ - 960.777114) * (-0.515038075)));
                //       tfwx_p2=(((fX+421.33683) * 0.857167301) - ((fZ - 958.387337) * (-0.515038075)));  //Christoph
                tfwy_p2 = (fY - 2.4);
                //       tfwy_p2=(fY-2.12);  //Christoph
                tfwt_p2 = ftime;
                tfwpath_p2 = flength;
                tfwpx_p2 = PX;
                tfwpy_p2 = PY;
                tfwpz_p2 = PZ;
                tfwe_p2 += TOFeloss;

                //       cout<<"TFW 2p X "<<fX<<" Y "<<fY<<" Z "<<fZ<<endl;
                //       cout<<"TFW - second proton "<<endl;
            }

            tfmul++;
        }
    }

    AddHit(tfmul,
           tfwx_p1,
           tfwy_p1,
           tfwt_p1,
           tfwx_p2,
           tfwy_p2,
           tfwt_p2,
           tfwpath_p1,
           tfwpath_p2,
           tfwpx_p1,
           tfwpy_p1,
           tfwpz_p1,
           tfwpx_p2,
           // tfwpy_p2,tfwpz_p2);
           tfwpy_p2,
           tfwpz_p2,
           tfwe_p1,
           tfwe_p2);

    //     cout<<"addhit"<<nentriesTof<<" tfmul "<<tfmul<<" tfwx_p1 "<<tfwx_p1<<" tfwy_p1 "<<tfwy_p1<<endl;
}
// -------------------------------------------------------------------------

void R3BTof2pDigitizer::Reset()
{
    // Clear the structure
    //   cout << " -I- Digit Reset() called " << endl;

    if (fTof2pDigi)
        fTof2pDigi->Clear();
}

void R3BTof2pDigitizer::Finish()
{
    // Write control histograms
}

R3BTof2pDigi* R3BTof2pDigitizer::AddHit(
    Int_t tfmul,
    Double_t tfwx_p1,
    Double_t tfwy_p1,
    Double_t tfwt_p1,
    Double_t tfwx_p2,
    Double_t tfwy_p2,
    Double_t tfwt_p2,
    Double_t tfwpath_p1,
    Double_t tfwpath_p2,
    Double_t tfwpx_p1,
    // Double_t tfwpy_p1,Double_t tfwpz_p1,Double_t tfwpx_p2,Double_t tfwpy_p2,Double_t tfwpz_p2){
    Double_t tfwpy_p1,
    Double_t tfwpz_p1,
    Double_t tfwpx_p2,
    Double_t tfwpy_p2,
    Double_t tfwpz_p2,
    Double_t tfwe_p1,
    Double_t tfwe_p2)
{
    TClonesArray& clref = *fTof2pDigi;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTof2pDigi(tfmul,
                                          tfwx_p1,
                                          tfwy_p1,
                                          tfwt_p1,
                                          tfwx_p2,
                                          tfwy_p2,
                                          tfwt_p2,
                                          tfwpath_p1,
                                          tfwpath_p2,
                                          // tfwpx_p1,tfwpy_p1,tfwpz_p1,tfwpx_p2,tfwpy_p2,tfwpz_p2);
                                          tfwpx_p1,
                                          tfwpy_p1,
                                          tfwpz_p1,
                                          tfwpx_p2,
                                          tfwpy_p2,
                                          tfwpz_p2,
                                          tfwe_p1,
                                          tfwe_p2);
}

// R3BTofDigi* R3BTofDigitizer::AddHit(
// return new(clref[size]) R3BTofDigi();
//}

ClassImp(R3BTof2pDigitizer)
