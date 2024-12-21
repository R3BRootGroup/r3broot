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

#include "R3BCalifavsFootOnlineSpectra.h"
#include "R3BCalifaClusterData.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BEventHeader.h"
#include "R3BFootCalData.h"
#include "R3BFootHitData.h"
#include "R3BFootMappedData.h"
#include "R3BLogger.h"
#include "R3BWRData.h"

#include "FairLogger.h"
#include "FairRootManager.h"
#include "FairRunAna.h"
#include "FairRunOnline.h"
#include "FairRuntimeDb.h"

#include "TCanvas.h"
#include "TFolder.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2F.h"
#include "THttpServer.h"
#include "TVector3.h"

#include "TClonesArray.h"
#include "TMath.h"
#include "TRandom.h"
#include <array>

using namespace std;

R3BCalifavsFootOnlineSpectra::R3BCalifavsFootOnlineSpectra()
    : R3BCalifavsFootOnlineSpectra("CALIFAvsFootOnlineSpectra", 1)
{
}

R3BCalifavsFootOnlineSpectra::R3BCalifavsFootOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fCalItemsCalifa(NULL)
    , fHitItemsCalifa(NULL)
    , fMappedItemsFoot(NULL)
    , fCalItemsFoot(NULL)
    , fHitItemsFoot(NULL)
    , fNEvents(0)
    , fTpat(-1)
    , fNbCalifaCrystals(4864)
    , fNumSides(Nb_Sides)
    , fMapHistos_bins(500)
    , fMapHistos_max(4000)
    , fBinsChannelFebex(5000)
    , fMaxBinChannelFebex(65535)
    , fMaxEnergyBarrel(10)
    , fMaxEnergyIphos(30)
    , fMinProtonE(50000.)
    , fRaw2Cal(kFALSE)
{
}

R3BCalifavsFootOnlineSpectra::~R3BCalifavsFootOnlineSpectra()
{
    R3BLOG(debug1, "");
    if (fHitItemsCalifa)
        delete fHitItemsCalifa;
    if (fMappedItemsFoot)
        delete fMappedItemsFoot;
}

void R3BCalifavsFootOnlineSpectra::SetParContainers()
{
    // Parameter Container
    // Reading amsStripCalPar from FairRuntimeDb
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");
}

InitStatus R3BCalifavsFootOnlineSpectra::Init()
{
    R3BLOG(info, "");
    FairRootManager* mgr = FairRootManager::Instance();

    R3BLOG_IF(fatal, NULL == mgr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // get access to Cal data
    fCalItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaCrystalCalData"));
    R3BLOG_IF(fatal, !fCalItemsCalifa, "CalifaCrystalCalData not found");

    // get access to Hit data
    fHitItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaClusterData"));
    R3BLOG_IF(fatal, !fHitItemsCalifa, "CalifaClusterData not found");

    fMappedItemsFoot = dynamic_cast<TClonesArray*>(mgr->GetObject("FootMappedData"));
    R3BLOG_IF(fatal, !fMappedItemsFoot, "FootMappedData not found");

    fCalItemsFoot = dynamic_cast<TClonesArray*>(mgr->GetObject("FootCalData"));
    R3BLOG_IF(warn, !fCalItemsFoot, "FootCalData not found");

    fHitItemsFoot = dynamic_cast<TClonesArray*>(mgr->GetObject("FootHitData"));
    R3BLOG_IF(warn, !fHitItemsFoot, "FootHitData not found");

    // Create histograms for detectors
    char Name1[255];
    char Name2[255];
    char Name3[255];

    // CANVAS Theta vs Phi
    cCalifa_angles = new TCanvas("Califa_Theta_vs_Phi_for_protons", "Theta vs Phi", 10, 10, 500, 500);
    cCalifa_angles->Divide(2, 1);
    for (int i = 0; i < 2; i++)
    {
        cCalifa_angles->cd(i + 1);
        char buf[512];
        snprintf(buf, 512, "%s%s", "fh2_Califa_theta_vs_phi", (i == 0) ? "" : "_withFoot");
        fh2_Califa_theta_phi[i] = new TH2F(buf, buf, 50, 0, 90, 180, -180, 180);
        fh2_Califa_theta_phi[i]->GetXaxis()->SetTitle("Theta [degrees]");
        fh2_Califa_theta_phi[i]->GetYaxis()->SetTitle("Phi [degrees]");
        fh2_Califa_theta_phi[i]->GetYaxis()->SetTitleOffset(1.2);
        fh2_Califa_theta_phi[i]->GetXaxis()->CenterTitle(true);
        fh2_Califa_theta_phi[i]->GetYaxis()->CenterTitle(true);
        fh2_Califa_theta_phi[i]->Draw("COLZ");
    }

    // MAIN FOLDER-Califa
    TFolder* mainfolCalifa = new TFolder("CALIFAvsFoot", "CALIFA vs Foot info");

    if (fHitItemsCalifa && fMappedItemsFoot)
    {
        mainfolCalifa->Add(cCalifa_angles);
    }
    run->AddObject(mainfolCalifa);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_CalifavsFoot",
                                          Form("/Objects/%s/->Reset_CALIFAFOOT_Histo()", GetName()));

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BCalifavsFootOnlineSpectra::ReInit()
{
    SetParContainers();
    // SetParameter();
    return kSUCCESS;
}

void R3BCalifavsFootOnlineSpectra::Reset_CALIFAFOOT_Histo()
{

    LOG(info) << "R3BCalifavsFootOnlineSpectra::Reset_CALIFAFOOT_Histo";

    for (int i = 0; i < 2; i++)
    {
        fh2_Califa_theta_phi[i]->Reset();
    }
}

void R3BCalifavsFootOnlineSpectra::Exec(Option_t* option)
{
    if ((fTpat >= 0) && (header) && ((header->GetTpat() & fTpat) != fTpat))
        return;

    bool ffoot = false;
    for (Int_t ihit = 0; ihit < fMappedItemsFoot->GetEntriesFast(); ihit++)
    {
        auto hit = dynamic_cast<R3BFootMappedData*>(fMappedItemsFoot->At(ihit));
        if (!hit)
            continue;
        if (hit->GetDetId() == 11 && hit->GetEnergy() > 500)
            ffoot = true;
    }

    Int_t nHits = fHitItemsCalifa->GetEntriesFast();

    for (Int_t ihit = 0; ihit < nHits; ihit++)
    {
        auto hit = dynamic_cast<R3BCalifaClusterData*>(fHitItemsCalifa->At(ihit));
        if (hit->GetEnergy() < 50e3) // 50MeV
            continue;

        double theta = hit->GetTheta() * TMath::RadToDeg();
        double phi = hit->GetPhi() * TMath::RadToDeg();
        fh2_Califa_theta_phi[0]->Fill(theta, phi); // always
        if (ffoot)
            fh2_Califa_theta_phi[1]->Fill(theta, phi); // only with foot
    }

    fNEvents += 1;
}

void R3BCalifavsFootOnlineSpectra::FinishEvent()
{

    if (fCalItemsCalifa)
    {
        fCalItemsCalifa->Clear();
    }
    if (fHitItemsCalifa)
    {
        fHitItemsCalifa->Clear();
    }
    if (fMappedItemsFoot)
    {
        fMappedItemsFoot->Clear();
    }
    if (fCalItemsFoot)
    {
        fCalItemsFoot->Clear();
    }
    if (fHitItemsFoot)
    {
        fHitItemsFoot->Clear();
    }
}

void R3BCalifavsFootOnlineSpectra::FinishTask()
{
    // Write canvas for Hit data
    if (fHitItemsCalifa)
    {
        fh2_Califa_theta_phi[0]->Write();

        if (fMappedItemsFoot)
            fh2_Califa_theta_phi[1]->Write();
    }
}

ClassImp(R3BCalifavsFootOnlineSpectra);
