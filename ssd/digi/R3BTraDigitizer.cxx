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

#include "R3BTraDigitizer.h"
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

#define SST_HALF_WIDTH_CM 3.5200
#define SST_HALF_HEIGHT_CM 1.9968

R3BTraDigitizer::R3BTraDigitizer()
    : FairTask("R3B Tra Digitization scheme ")
{
}

R3BTraDigitizer::~R3BTraDigitizer() {}

void R3BTraDigitizer::SetParContainers()
{

    // Get run and runtime database
    FairRunAna* run = FairRunAna::Instance();
    if (!run)
        LOG(fatal) << "SetParContainers: No analysis run";

    FairRuntimeDb* rtdb = run->GetRuntimeDb();
    if (!rtdb)
        LOG(fatal) << "SetParContainers: No runtime database";

    fTraDigiPar = dynamic_cast<R3BTraDigiPar*>(rtdb->getContainer("R3BTraDigiPar"));

    if (fTraDigiPar)
    {
        LOG(info) << "-I- R3BTraDigitizer::SetParContainers() ";
        LOG(info) << "-I- Container R3BTraDigiPar  loaded ";
    }
}

InitStatus R3BTraDigitizer::Init()
{
    //  LOG(info)<<"Init ";
    // Get input array
    FairRootManager* ioman = FairRootManager::Instance();
    if (!ioman)
        LOG(fatal) << "Init: No FairRootManager";
    fTraPoints = dynamic_cast<TClonesArray*>(ioman->GetObject("TraPoint"));
    fTraMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));
    fMCTrack = dynamic_cast<TClonesArray*>(ioman->GetObject("MCTrack"));

    // Register output array TraDigi
    fTraDigi = new TClonesArray("R3BTraDigi", 1000);
    ioman->Register("TraDigi", "Digital response in Tra", fTraDigi, kTRUE);

    eventNoTra = 0;

    // Initialise control histograms

    return kSUCCESS;
}

void R3BTraDigitizer::Exec(Option_t* opt)
{

    Reset();
    eventNoTra += 1;
    //     if(eventNoTra/1000. == (int)eventNoTra/1000.) LOG(info)<<"Event #: "<<eventNoTra-1;

    Int_t nentriesTra = fTraPoints->GetEntries();
    //     Int_t nentries = fMCTrack->GetEntries();

    Int_t TrackIdTra = 0;
    //     Int_t TrackId=0;

    Double_t ss03_se_f = nan("");
    Double_t ss03_spos_f = nan("");
    Double_t ss03_sbw_f = nan("");
    Double_t ss03_sarea_f = nan("");
    Double_t ss03_seta_f = nan("");

    Double_t ss03_ke_f = nan("");
    Double_t ss03_kpos_f = nan("");
    Double_t ss03_kbw_f = nan("");
    Double_t ss03_karea_f = nan("");
    Double_t ss03_keta_f = nan("");

    Double_t ss06_se_f = nan("");
    Double_t ss06_spos_f = nan("");
    Double_t ss06_sbw_f = nan("");
    Double_t ss06_sarea_f = nan("");
    Double_t ss06_seta_f = nan("");

    Double_t ss06_ke_f = nan("");
    Double_t ss06_kpos_f = nan("");
    Double_t ss06_kbw_f = nan("");
    Double_t ss06_karea_f = nan("");
    Double_t ss06_keta_f = nan("");

    Double_t ss03_se_p1 = nan("");
    Double_t ss03_spos_p1 = nan("");
    Double_t ss03_sbw_p1 = nan("");
    Double_t ss03_sarea_p1 = nan("");
    Double_t ss03_seta_p1 = nan("");

    Double_t ss03_ke_p1 = nan("");
    Double_t ss03_kpos_p1 = nan("");
    Double_t ss03_kbw_p1 = nan("");
    Double_t ss03_karea_p1 = nan("");
    Double_t ss03_keta_p1 = nan("");

    Double_t ss06_se_p1 = nan("");
    Double_t ss06_spos_p1 = nan("");
    Double_t ss06_sbw_p1 = nan("");
    Double_t ss06_sarea_p1 = nan("");
    Double_t ss06_seta_p1 = nan("");

    Double_t ss06_ke_p1 = nan("");
    Double_t ss06_kpos_p1 = nan("");
    Double_t ss06_kbw_p1 = nan("");
    Double_t ss06_karea_p1 = nan("");
    Double_t ss06_keta_p1 = nan("");

    //******************** SSTs **************************//

    for (Int_t l = 0; l < nentriesTra; l++)
    {
        //   LOG(info)<<"entries "<<l;

        R3BTraPoint* Tra_obj = dynamic_cast<R3BTraPoint*>(fTraPoints->At(l));

        //     Int_t DetID = Tra_obj->GetDetectorID();
        Double_t fX_In = Tra_obj->GetXIn();
        Double_t fY_In = Tra_obj->GetYIn();
        Double_t fZ_In = Tra_obj->GetZIn();
        Double_t fX_Out = Tra_obj->GetXOut();
        Double_t fY_Out = Tra_obj->GetYOut();
        Double_t fZ_Out = Tra_obj->GetZOut();
        TrackIdTra = Tra_obj->GetTrackID();
        R3BMCTrack* aTrack = dynamic_cast<R3BMCTrack*>(fTraMCTrack->At(TrackIdTra));
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
                // ss03_kpos_f = (3.9936-(0.006402 + (2 + fY)));//detectors couldn't be shifted in simulation, they are
                // shifted here -0.006402
                ss03_kpos_f =
                    (2. * SST_HALF_HEIGHT_CM -
                     (0.006402 + (SST_HALF_HEIGHT_CM +
                                  fY))); // detectors couldn't be shifted in simulation, they are shifted here -0.006402
                // ss03_kpos_f = (2 + fY);
                ss03_kbw_f = 2; // kpos for 1 SST swop (with 3.9936-) - requirement for tracker
                ss03_karea_f = 0;
                ss03_keta_f = 0;
                //	LOG(info)<<"SST1 - fragment "<<PID;
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
                //	LOG(info)<<"SST2 - fragment "<<PID;
            }
        }

        if (PID == 2212 && mother < 0)
        {

            // if (DetID==15)
            // if (DetID==21)
            // if (DetID==23)
            if (fZ < 12)
            {
                ss03_se_p1 = 50;
                // ss03_spos_p1 = (-0.038230 + (3.5 + fX)); //detectors couldn't be shifted in simulation, they are
                // shifted here 0.038230
                ss03_spos_p1 =
                    (-0.038230 + (SST_HALF_WIDTH_CM +
                                  fX)); // detectors couldn't be shifted in simulation, they are shifted here 0.038230
                // ss03_spos_p1 = (3.5 + fX);
                ss03_sbw_p1 = 2;
                ss03_sarea_p1 = 0;
                ss03_seta_p1 = 0;

                ss03_ke_p1 = 50;
                // ss03_kpos_p1 = (3.9936-(0.006402 + (2 + fY))); //detectors couldn't be shifted in simulation, they
                // are shifted here -0.006402
                ss03_kpos_p1 =
                    (2. * SST_HALF_HEIGHT_CM -
                     (0.006402 + (SST_HALF_HEIGHT_CM +
                                  fY))); // detectors couldn't be shifted in simulation, they are shifted here -0.006402
                // ss03_kpos_p1 = (2 + fY);
                ss03_kbw_p1 = 2; // kpos for 1 SST swop (with 3.9936-) - requirement for tracker
                ss03_karea_p1 = 0;
                ss03_keta_p1 = 0;
                //	LOG(info)<<"SST1 - proton "<<PID;
                //	LOG(info)<<"l "<<l<<" DetID "<<DetID<<" PID "<<PID;
            }
            // if (DetID==16)
            // if (DetID==22)
            // if (DetID==24)
            if (fZ > 12)
            {
                ss06_se_p1 = 50;
                // ss06_spos_p1 = (0.038495 + (3.5 + fX));//detectors couldn't be shifted in simulation, they are
                // shifted here -0.038495
                ss06_spos_p1 =
                    (0.038495 + (SST_HALF_WIDTH_CM +
                                 fX)); // detectors couldn't be shifted in simulation, they are shifted here -0.038495
                // ss06_spos_p1 = (3.5 + fX);
                ss06_sbw_p1 = 2;
                ss06_sarea_p1 = 0;
                ss06_seta_p1 = 0;

                ss06_ke_p1 = 50;
                // ss06_kpos_p1 = (-0.00798 + (2 + fY));//detectors couldn't be shifted in simulation, they are shifted
                // here 0.00798
                ss06_kpos_p1 =
                    (-0.00798 + (SST_HALF_HEIGHT_CM +
                                 fY)); // detectors couldn't be shifted in simulation, they are shifted here 0.00798
                // ss06_kpos_p1 = (2 + fY);
                ss06_kbw_p1 = 2;
                ss06_karea_p1 = 0;
                ss06_keta_p1 = 0;
                //	LOG(info)<<"SST2 - proton "<<PID;
                //	LOG(info)<<"l "<<l<<" DetID "<<DetID<<" PID "<<PID;
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
           ss06_keta_f,
           ss03_se_p1,
           ss03_spos_p1,
           ss03_sbw_p1,
           ss03_sarea_p1,
           ss03_seta_p1,
           ss03_ke_p1,
           ss03_kpos_p1,
           ss03_kbw_p1,
           ss03_karea_p1,
           ss03_keta_p1,
           ss06_se_p1,
           ss06_spos_p1,
           ss06_sbw_p1,
           ss06_sarea_p1,
           ss06_seta_p1,
           ss06_ke_p1,
           ss06_kpos_p1,
           ss06_kbw_p1,
           ss06_karea_p1,
           ss06_keta_p1);
}
// -------------------------------------------------------------------------

void R3BTraDigitizer::Reset()
{
    // Clear the structure
    //   LOG(info) << " -I- Digit Reset() called ";

    if (fTraDigi)
        fTraDigi->Clear();
}

void R3BTraDigitizer::Finish()
{
    // Write control histograms
}

R3BTraDigi* R3BTraDigitizer::AddHit(Double_t ss03_se_f,
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
                                    Double_t ss06_keta_f,
                                    Double_t ss03_se_p1,
                                    Double_t ss03_spos_p1,
                                    Double_t ss03_sbw_p1,
                                    Double_t ss03_sarea_p1,
                                    Double_t ss03_seta_p1,
                                    Double_t ss03_ke_p1,
                                    Double_t ss03_kpos_p1,
                                    Double_t ss03_kbw_p1,
                                    Double_t ss03_karea_p1,
                                    Double_t ss03_keta_p1,
                                    Double_t ss06_se_p1,
                                    Double_t ss06_spos_p1,
                                    Double_t ss06_sbw_p1,
                                    Double_t ss06_sarea_p1,
                                    Double_t ss06_seta_p1,
                                    Double_t ss06_ke_p1,
                                    Double_t ss06_kpos_p1,
                                    Double_t ss06_kbw_p1,
                                    Double_t ss06_karea_p1,
                                    Double_t ss06_keta_p1)
{
    TClonesArray& clref = *fTraDigi;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BTraDigi(ss03_se_f,
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
                                        ss06_keta_f,
                                        ss03_se_p1,
                                        ss03_spos_p1,
                                        ss03_sbw_p1,
                                        ss03_sarea_p1,
                                        ss03_seta_p1,
                                        ss03_ke_p1,
                                        ss03_kpos_p1,
                                        ss03_kbw_p1,
                                        ss03_karea_p1,
                                        ss03_keta_p1,
                                        ss06_se_p1,
                                        ss06_spos_p1,
                                        ss06_sbw_p1,
                                        ss06_sarea_p1,
                                        ss06_seta_p1,
                                        ss06_ke_p1,
                                        ss06_kpos_p1,
                                        ss06_kbw_p1,
                                        ss06_karea_p1,
                                        ss06_keta_p1);
    LOG(info) << "Point 19";
}

// R3BTraDigi* R3BTraDigitizer::AddHit(
// return new(clref[size]) R3BTraDigi();
//}

ClassImp(R3BTraDigitizer)
