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

#include "R3BTofDigitizer.h"
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

R3BTofDigitizer::R3BTofDigitizer()
    : FairTask("R3B Tof Digitization scheme ")
{
}

R3BTofDigitizer::~R3BTofDigitizer() {}

void R3BTofDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fTofDigiPar = dynamic_cast<R3BTofDigiPar*>(rtdb->getContainer("R3BTofDigiPar"));

    if (fTofDigiPar)
    {
        cout << "-I- R3BTofDigitizer::SetParContainers() " << endl;
        cout << "-I- Container R3BTofDigiPar  loaded " << endl;
    }
}

InitStatus R3BTofDigitizer::Init()
{

    //  cout<<"Init "<<endl;
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fTofPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("TOFPoint"));
    fTofMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array TofDigi
    fTofDigi = new TClonesArray("R3BTofDigi", 1000);
    ioman->Register("TofDigi", "Digital response in Tof", fTofDigi, kTRUE);

    eventNoTof = 0;

    // Initialise control histograms

    TOFelosshis = new TH1F("TOFelosshis", "TOFelosshis", 500, 0., 80.);

    TfwXhis = new TH1F("TfwXhis", "TfwXhis", 1000, 0., 100.);
    TfwXhis->GetXaxis()->SetTitle("Position");
    TfwXhis->GetYaxis()->SetTitle("Counts");

    TfwYhis = new TH1F("TfwYhis", "TfwYhis", 500, -80., 80.);
    TfwYhis->GetXaxis()->SetTitle("Position");
    TfwYhis->GetYaxis()->SetTitle("Counts");

    TfwThis = new TH1F("TfwThis", "TfwThis", 500, -10., 4000.);
    TfwThis->GetXaxis()->SetTitle("Time");
    TfwThis->GetYaxis()->SetTitle("Counts");

    return kSUCCESS;
}

void R3BTofDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoTof += 1;
    if (eventNoTof / 1000. == (int)eventNoTof / 1000)
        cout << "Event #: " << eventNoTof << endl;

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

    tfmul = 0;

    Double_t x4He = 0.;
    Double_t x12C = 0.;

    for (Int_t l = 0; l < nentriesTof; l++)
    {
        R3BTofPoint* tof_obj = dynamic_cast<R3BTofPoint*>(fTofPoints->At(l));
        //        if (tof_obj==NULL) continue;
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

        if (PID == 1000020040 && mother < 0)
        { // alpha
            cout << "x4He: " << fX << endl;
            x4He = fX;
        }

        if (PID == 1000060120 && mother < 0)
        { // 12C
            cout << "x12C: " << fX << endl;
            x12C = fX;
        }

        // offset - position of detector x: -417.359574; y: +2.4; z: +960.777114
        // angle = -31 degrees
        // cos(angle) = 0.857167301
        // sin(angle) = -0.515038075

        tfmul++;
    }

    if (x4He != 0. && x12C != 0)
    {
        cout << "dX: " << fabs(x4He - x12C) << endl;
        TfwXhis->Fill(fabs(x4He - x12C));
    }
    //    TfwYhis->Fill(tfwy_p1);
    //    TfwThis->Fill(tfwt_p1);

    // AddHit(tfmul,tfwx_p1,tfwy_p1,tfwt_p1,tfwpath_p1,tfwpx_p1,tfwpy_p1,tfwpz_p1);
    // AddHit(tfmul,tfwx_p1,tfwy_p1,tfwt_p1,tfwpath_p1,tfwpx_p1,tfwpy_p1,tfwpz_p1,tfwe_p1);

    //     cout<<"addhit"<<nentriesTof<<" tfmul "<<tfmul<<" tfwx_p1 "<<tfwx_p1<<" tfwy_p1 "<<tfwy_p1<<endl;
}
// -------------------------------------------------------------------------

void R3BTofDigitizer::Reset()
{
    // Clear the structure
    //   cout << " -I- Digit Reset() called " << endl;

    if (fTofDigi)
        fTofDigi->Clear();
}

void R3BTofDigitizer::Finish()
{
    // Write control histograms

    TOFelosshis->Write();

    TfwXhis->Write();
    TfwYhis->Write();
    TfwThis->Write();
}

R3BTofDigi* R3BTofDigitizer::AddHit(Int_t tfmul,
                                    Double_t tfwx_p1,
                                    Double_t tfwy_p1,
                                    Double_t tfwt_p1,
                                    Double_t tfwpath_p1,
                                    // Double_t tfwpx_p1,Double_t tfwpy_p1,Double_t tfwpz_p1){
                                    Double_t tfwpx_p1,
                                    Double_t tfwpy_p1,
                                    Double_t tfwpz_p1,
                                    Double_t tfwe_p1)
{
    TClonesArray& clref = *fTofDigi;
    Int_t size = clref.GetEntriesFast();
    // return new(clref[size]) R3BTofDigi(tfmul,tfwx_p1,tfwy_p1,tfwt_p1,tfwpath_p1,tfwpx_p1,tfwpy_p1,tfwpz_p1);
    return new (clref[size])
        R3BTofDigi(tfmul, tfwx_p1, tfwy_p1, tfwt_p1, tfwpath_p1, tfwpx_p1, tfwpy_p1, tfwpz_p1, tfwe_p1);
}

// R3BTofDigi* R3BTofDigitizer::AddHit(
// return new(clref[size]) R3BTofDigi();
//}

ClassImp(R3BTofDigitizer)
