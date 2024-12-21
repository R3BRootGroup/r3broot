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

#include "R3BTraFraDigitizer.h"
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

#define SST_HALF_WIDTH_CM 3.5200 // from Ralf's tracker
#define SST_HALF_HEIGHT_CM 1.9968

R3BTraFraDigitizer::R3BTraFraDigitizer()
    : FairTask("R3B TraFra Digitization scheme ")
{
}

R3BTraFraDigitizer::~R3BTraFraDigitizer() {}

void R3BTraFraDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fTraFraDigiPar = dynamic_cast<R3BTraFraDigiPar*>(rtdb->getContainer("R3BTraFraDigiPar"));

    if (fTraFraDigiPar)
    {
        LOG(info) << "-I- R3BTraFraDigitizer::SetParContainers() ";
        LOG(info) << "-I- Container R3BTraFraDigiPar  loaded ";
    }
}

InitStatus R3BTraFraDigitizer::Init()
{

    //  LOG(info)<<"Init ";
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fTraFraPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("TraPoint"));
    fTraFraMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    fMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array TraFraDigi
    fTraFraDigi = new TClonesArray("R3BTraFraDigi", 1000);
    ioman->Register("TraFraDigi", "Digital response in TraFra", fTraFraDigi, kTRUE);

    eventNoTra = 0;

    // Initialise control histograms

    return kSUCCESS;
}

void R3BTraFraDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoTra += 1;
    //     if(eventNoTra/1000. == (int)eventNoTra/1000.) LOG(info)<<"Event #: "<<eventNoTra-1;

    Int_t nentriesTra = fTraFraPoints->GetEntries();
    // Int_t nentries = fMCTrack->GetEntries();

    Int_t TrackIdTra = 0;
    // Int_t TrackId=0;

    Double_t ss03_se_f;
    Double_t ss03_spos_f = nan("");
    Double_t ss03_sbw_f;
    Double_t ss03_sarea_f;
    Double_t ss03_seta_f;

    Double_t ss03_ke_f;
    Double_t ss03_kpos_f = nan("");
    Double_t ss03_kbw_f;
    Double_t ss03_karea_f;
    Double_t ss03_keta_f;

    Double_t ss06_se_f;
    Double_t ss06_spos_f = nan("");
    Double_t ss06_sbw_f;
    Double_t ss06_sarea_f;
    Double_t ss06_seta_f;

    Double_t ss06_ke_f;
    Double_t ss06_kpos_f = nan("");
    Double_t ss06_kbw_f;
    Double_t ss06_karea_f;
    Double_t ss06_keta_f;

    //******************** SSTs **************************//

    for (Int_t l = 0; l < nentriesTra; l++)
    {
        //   LOG(info)<<"entries "<<l;

        R3BTraPoint* TraFra_obj = dynamic_cast<R3BTraPoint*>(fTraFraPoints->At(l));

        // Int_t DetID = TraFra_obj->GetDetectorID();
        Double_t fX_In = TraFra_obj->GetXIn();
        Double_t fY_In = TraFra_obj->GetYIn();
        Double_t fZ_In = TraFra_obj->GetZIn();
        Double_t fX_Out = TraFra_obj->GetXOut();
        Double_t fY_Out = TraFra_obj->GetYOut();
        Double_t fZ_Out = TraFra_obj->GetZOut();
        TrackIdTra = TraFra_obj->GetTrackID();
        R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fTraFraMCTrack->At(TrackIdTra));
        Int_t PID = aTrack->GetPdgCode();
        Int_t mother = aTrack->GetMotherId();

        Double_t fX = ((fX_In + fX_Out) / 2);
        Double_t fY = ((fY_In + fY_Out) / 2);
        Double_t fZ = ((fZ_In + fZ_Out) / 2);

        if (PID == 1000080150 && mother < 0)
        {

            // if (DetID==15)
            // if (DetID==21)
            // if (DetID==23)
            if (fZ < 12)
            {
                ss03_se_f = 2000;
                // ss03_spos_f = (-0.038230 + (3.5 + fX)); //detectors couldn't be shifted in simulation, they are
                // shifted here 0.038230
                ss03_spos_f =
                    (-0.038230 + (SST_HALF_WIDTH_CM +
                                  fX)); // detectors couldn't be shifted in simulation, they are shifted here 0.038230
                // ss03_spos_f = (3.5 + fX);
                ss03_sbw_f = 2;
                ss03_sarea_f = 0;
                ss03_seta_f = 0;

                ss03_ke_f = 2000;
                // ss03_kpos_f = (3.9936-((2 + fY)+0.006402 )); //detectors couldn't be shifted in simulation, they are
                // shifted here -0.006402
                ss03_kpos_f =
                    (2. * SST_HALF_HEIGHT_CM -
                     ((SST_HALF_HEIGHT_CM + fY) +
                      0.006402)); // detectors couldn't be shifted in simulation, they are shifted here -0.006402
                // ss03_kpos_f = (2 + fY);
                ss03_kbw_f = 2; // kpos for 1 SST swop (with 3.9936-) - requirement for tracker
                ss03_karea_f = 0;
                ss03_keta_f = 0;
                //	LOG(info)<<"ss03_se_f "<<ss03_se_f;
            }
            // if (DetID==16)
            // if (DetID==22)
            // if (DetID==24)
            if (fZ > 12)
            {
                ss06_se_f = 2000;
                // ss06_spos_f = (0.038495 + (3.5 + fX));//detectors couldn't be shifted in simulation, they are shifted
                // here -0.038495
                ss06_spos_f =
                    (0.038495 + (SST_HALF_WIDTH_CM +
                                 fX)); // detectors couldn't be shifted in simulation, they are shifted here -0.038495
                // ss06_spos_f = (3.5 + fX);
                ss06_sbw_f = 2;
                ss06_sarea_f = 0;
                ss06_seta_f = 0;

                ss06_ke_f = 2000;
                // ss06_kpos_f = (-0.00798 + (2 + fY));//detectors couldn't be shifted in simulation, they are shifted
                // here 0.00798
                ss06_kpos_f =
                    (-0.00798 + (SST_HALF_HEIGHT_CM +
                                 fY)); // detectors couldn't be shifted in simulation, they are shifted here 0.00798
                // ss06_kpos_f = (2 + fY);
                ss06_kbw_f = 2;
                ss06_karea_f = 0;
                ss06_keta_f = 0;
            }
        }
    }

    AddHit(ss03_se_f,
           ss03_spos_f,
           ss03_sbw_f,
           ss03_sarea_f,
           ss03_seta_f,
           ss03_ke_f,
           ss03_kpos_f,
           ss03_kbw_f,
           ss03_karea_f,
           ss03_keta_f,
           ss06_se_f,
           ss06_spos_f,
           ss06_sbw_f,
           ss06_sarea_f,
           ss06_seta_f,
           ss06_ke_f,
           ss06_kpos_f,
           ss06_kbw_f,
           ss06_karea_f,
           ss06_keta_f);
}
// -------------------------------------------------------------------------

void R3BTraFraDigitizer::Reset()
{
    // Clear the structure
    //   LOG(info) << " -I- Digit Reset() called ";

    if (fTraFraDigi)
        fTraFraDigi->Clear();
}

void R3BTraFraDigitizer::Finish()
{
    // Write control histograms
}

R3BTraFraDigi* R3BTraFraDigitizer::AddHit(Double_t ss03_se_f,
                                          Double_t ss03_spos_f,
                                          Double_t ss03_sbw_f,
                                          Double_t ss03_sarea_f,
                                          Double_t ss03_seta_f,
                                          Double_t ss03_ke_f,
                                          Double_t ss03_kpos_f,
                                          Double_t ss03_kbw_f,
                                          Double_t ss03_karea_f,
                                          Double_t ss03_keta_f,
                                          Double_t ss06_se_f,
                                          Double_t ss06_spos_f,
                                          Double_t ss06_sbw_f,
                                          Double_t ss06_sarea_f,
                                          Double_t ss06_seta_f,
                                          Double_t ss06_ke_f,
                                          Double_t ss06_kpos_f,
                                          Double_t ss06_kbw_f,
                                          Double_t ss06_karea_f,
                                          Double_t ss06_keta_f)
{
    TClonesArray& clref = *fTraFraDigi;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTraFraDigi(ss03_se_f,
                                           ss03_spos_f,
                                           ss03_sbw_f,
                                           ss03_sarea_f,
                                           ss03_seta_f,
                                           ss03_ke_f,
                                           ss03_kpos_f,
                                           ss03_kbw_f,
                                           ss03_karea_f,
                                           ss03_keta_f,
                                           ss06_se_f,
                                           ss06_spos_f,
                                           ss06_sbw_f,
                                           ss06_sarea_f,
                                           ss06_seta_f,
                                           ss06_ke_f,
                                           ss06_kpos_f,
                                           ss06_kbw_f,
                                           ss06_karea_f,
                                           ss06_keta_f);
}

// R3BTraFraDigi* R3BTraFraDigitizer::AddHit(
// return new(clref[size]) R3BTraFraDigi();
//}

ClassImp(R3BTraFraDigitizer)
