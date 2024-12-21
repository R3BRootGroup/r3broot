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

#include "R3BCalifaOnlineSpectra.h"
#include "R3BCalifaClusterData.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaMappedData.h"
#include "R3BCalifaMappingPar.h"
#include "R3BEventHeader.h"
#include "R3BLogger.h"
#include "R3BWRData.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRunOnline.h>
#include <FairRuntimeDb.h>

#include <TCanvas.h>
#include <TClonesArray.h>
#include <TFolder.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <THttpServer.h>
#include <TMath.h>
#include <TRandom.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

R3BCalifaOnlineSpectra::R3BCalifaOnlineSpectra()
    : R3BCalifaOnlineSpectra("CALIFAOnlineSpectra", 1)
{
}

R3BCalifaOnlineSpectra::R3BCalifaOnlineSpectra(const TString& name, Int_t iVerbose)
    : FairTask(name, iVerbose)
{
}

void R3BCalifaOnlineSpectra::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");

    fMap_Par = dynamic_cast<R3BCalifaMappingPar*>(rtdb->getContainer("califaMappingPar"));
    if (!fMap_Par)
    {
        R3BLOG(error, "Couldn't get handle on califaMappingPar container");
    }
    else
    {
        R3BLOG(info, "CalifaMappingPar container open");
    }
}

void R3BCalifaOnlineSpectra::SetParameter()
{
    R3BLOG_IF(error, !fMap_Par, "CalifaMappingPar container not found");

    //--- Parameter Container ---
    fNbCalifaCrystals = fMap_Par->GetNumCrystals(); // Number of crystals
    assert(std::clamp(fNbCalifaCrystals, 0, MaxNbCrystals) == fNbCalifaCrystals && "Number of crystals out of range");

    R3BLOG(info, "Number of crystals (gamma+proton): " << fNbCalifaCrystals);
    // fMap_Par->printParams();

    fFebexInfo.resize(extents[Nb_Sides][Nb_Rings][Nb_Preamps][Nb_SlotandModule]);
    std::fill(fFebexInfo.data(), fFebexInfo.data() + fFebexInfo.num_elements(), -1);

    for (Int_t c = 1; c <= fNbCalifaCrystals; c++)
    {
        if (c <= fNbCalifaCrystals / 2)
        {
            if (fMap_Par->GetInUse(c) == 1 && fMap_Par->GetHalf(c) > 0)
            { // only for installed crystals (see issue 681)
                fFebexInfo[fMap_Par->GetHalf(c) - 1][fMap_Par->GetRing(c) - 1][fMap_Par->GetPreamp(c) - 1][0] =
                    fMap_Par->GetFebexSlot(c);
                fFebexInfo[fMap_Par->GetHalf(c) - 1][fMap_Par->GetRing(c) - 1][fMap_Par->GetPreamp(c) - 1][1] =
                    fMap_Par->GetFebexMod(c);
            }
        }
        else
        {
            if (fMap_Par->GetInUse(c) == 1 && fMap_Par->GetHalf(c) > 0)
            {
                fFebexInfo[fMap_Par->GetHalf(c) - 1][fMap_Par->GetRing(c) - 1][fMap_Par->GetPreamp(c) - 1][2] =
                    fMap_Par->GetFebexSlot(c);
                fFebexInfo[fMap_Par->GetHalf(c) - 1][fMap_Par->GetRing(c) - 1][fMap_Par->GetPreamp(c) - 1][3] =
                    fMap_Par->GetFebexMod(c);
            }
        }
    }
}

InitStatus R3BCalifaOnlineSpectra::Init()
{
    R3BLOG(info, "");

    FairRootManager* mgr = FairRootManager::Instance();

    R3BLOG_IF(fatal, mgr == nullptr, "FairRootManager not found");

    header = dynamic_cast<R3BEventHeader*>(mgr->GetObject("EventHeader."));

    FairRunOnline* run = FairRunOnline::Instance();
    run->GetHttpServer()->Register("", this);

    // get access to Mapped data
    fMappedItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaMappedData"));
    R3BLOG_IF(fatal, fMappedItemsCalifa == nullptr, "CalifaMappedData not found");

    // get access to trigger Mapped data
    fTrigMappedItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaMappedtrigData"));
    R3BLOG_IF(warn, !fTrigMappedItemsCalifa, "CalifaMappedtrigData not found");

    // get access to Cal data
    fCalItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaCrystalCalData"));
    R3BLOG_IF(warn, !fCalItemsCalifa, "CalifaCrystalCalData not found");

    // get access to Hit data
    fHitItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("CalifaClusterData"));
    R3BLOG_IF(warn, !fHitItemsCalifa, "CalifaClusterData not found");

    // get access to WR-Califa data
    fWRItemsCalifa = dynamic_cast<TClonesArray*>(mgr->GetObject("WRCalifaData"));
    R3BLOG_IF(warn, !fWRItemsCalifa, "WRCalifaData not found");

    // get access to WR-Master data
    fWRItemsMaster = dynamic_cast<TClonesArray*>(mgr->GetObject("WRMasterData"));
    R3BLOG_IF(warn, !fWRItemsMaster, "WRMasterData not found");

    SetParameter();

    // reading the file
    ifstream* FileHistos = new ifstream(fCalifaFile);
    if (FileHistos->is_open() == false)
    {
        R3BLOG(warn, "No Histogram definition file");
    }

    Double_t arry_bins[fNumSides][fNumRings][fNumPreamps][fNumCrystalPreamp];
    Double_t arry_maxE[fNumSides][fNumRings][fNumPreamps][fNumCrystalPreamp];
    Double_t arry_minE[fNumSides][fNumRings][fNumPreamps][fNumCrystalPreamp];

    if (FileHistos->is_open())
    {
        for (int s = 0; s < fNumSides; s++)
            for (int r = 0; r < fNumRings; r++)
                for (int p = 0; p < fNumPreamps; p++)
                    for (int ch = 0; ch < fNumCrystalPreamp; ch++)
                        *FileHistos >> arry_bins[s][r][p][ch] >> arry_minE[s][r][p][ch] >> arry_maxE[s][r][p][ch];
    }
    else
    {
        for (int s = 0; s < fNumSides; s++)
            for (int r = 0; r < fNumRings; r++)
                for (int p = 0; p < fNumPreamps; p++)
                    for (int ch = 0; ch < fNumCrystalPreamp; ch++)
                    {
                        arry_bins[s][r][p][ch] = fMapHistos_bins;
                        arry_minE[s][r][p][ch] = 0;
                        arry_maxE[s][r][p][ch] = fMapHistos_max;
                    }
    }

    // Create histograms for detectors
    std::string Name1;
    std::string Name2;
    std::string Name3;
    std::string Xaxis;
    double bins = fMapHistos_bins;
    double maxE = fMapHistos_max;
    double minE = 0.;

    // CANVAS Crystal_ID vs energy
    cCalifa_cry_energy = new TCanvas("Califa_Map_cryID_energy", "Califa_Map energy vs cryID", 10, 10, 500, 500);
    fh2_Califa_cryId_energy = R3B::root_owned<TH2F>("fh2_Califa_Map_cryID_energy",
                                                    "Califa: CryID vs energy",
                                                    fMap_Par->GetNumCrystals(),
                                                    0.5,
                                                    fMap_Par->GetNumCrystals() + 0.5,
                                                    fBinsChannelFebex,
                                                    0.0,
                                                    fMaxBinChannelFebex);
    fh2_Califa_cryId_energy->GetXaxis()->SetTitle("Crystal ID");
    fh2_Califa_cryId_energy->GetYaxis()->SetTitle("Energy [Channels]");
    fh2_Califa_cryId_energy->GetYaxis()->SetTitleOffset(1.4);
    fh2_Califa_cryId_energy->GetXaxis()->CenterTitle(true);
    fh2_Califa_cryId_energy->GetYaxis()->CenterTitle(true);
    gPad->SetLogz();
    fh2_Califa_cryId_energy->Draw("COLZ");

    // CANVAS Crystal_ID vs energy
    cCalifa_cry_energy_cal = new TCanvas("Califa_Cal_cryID_energy", "Califa_Cal energy vs cryID", 10, 10, 500, 500);
    fh2_Califa_cryId_energy_cal = R3B::root_owned<TH2F>("fh2_Califa_Cal_cryID_energy",
                                                        "Califa: CryID vs calibrated energy",
                                                        fMap_Par->GetNumCrystals(),
                                                        0.5,
                                                        fMap_Par->GetNumCrystals() + 0.5,
                                                        bins,
                                                        minE,
                                                        maxE);
    fh2_Califa_cryId_energy_cal->GetXaxis()->SetTitle("Crystal ID");
    fh2_Califa_cryId_energy_cal->GetYaxis()->SetTitle("Energy [keV]");
    fh2_Califa_cryId_energy_cal->GetYaxis()->SetTitleOffset(1.4);
    fh2_Califa_cryId_energy_cal->GetXaxis()->CenterTitle(true);
    fh2_Califa_cryId_energy_cal->GetYaxis()->CenterTitle(true);
    gPad->SetLogz();
    fh2_Califa_cryId_energy_cal->Draw("COLZ");

    cCalifa_NsNf = new TCanvas("Califa_Cal_NsNf", "Califa_Cal Ns vs Nf", 10, 10, 500, 500);
    cCalifa_NsNf->Divide(1, 3);
    fh2_Califa_NsNf.resize(3);

    std::vector<std::string> region = { "Barrel", "iPhos", "CEPA" };
    for (Int_t i = 0; i < fh2_Califa_NsNf.size(); i++)
    {
        std::stringstream ss1;
        ss1 << "PID_" << region[i];
        fh2_Califa_NsNf[i] =
            R3B::root_owned<TH2F>(ss1.str().c_str(), "Califa PID: Ns and Nf energies", 500, 0., 500., 600, 0., 1.);
        fh2_Califa_NsNf[i]->GetXaxis()->SetTitle("Ns+Nf Energies [MeV]");
        fh2_Califa_NsNf[i]->GetYaxis()->SetTitle("Nf/(Nf+Ns)");
        fh2_Califa_NsNf[i]->GetYaxis()->SetTitleOffset(1.4);
        fh2_Califa_NsNf[i]->GetXaxis()->CenterTitle(true);
        fh2_Califa_NsNf[i]->GetYaxis()->CenterTitle(true);
        cCalifa_NsNf->cd(i + 1);
        gPad->SetLogz();
        fh2_Califa_NsNf[i]->Draw("COLZ");
    }

    // Map data
    for (Int_t i = 0; i < fNumRings; i++)
    {
        // Left side
        std::stringstream ss1;
        ss1 << "Ring_" << i + 1 << "_Left";
        cMap_RingL[i] = new TCanvas(ss1.str().c_str(), ss1.str().c_str(), 10, 10, 800, 700);
        std::stringstream ss2;
        ss2 << "fh2_Ring_" << i + 1 << "_Left";
        std::stringstream ss3;
        ss3 << "Ring " << i + 1 << " left: Preamp number vs channel number";

        fh2_Preamp_vs_ch_L[i] = R3B::root_owned<TH2F>(ss2.str().c_str(),
                                                      ss3.str().c_str(),
                                                      fNumPreamps,
                                                      0.5,
                                                      fNumPreamps + 0.5,
                                                      fNumCrystalPreamp,
                                                      0.5,
                                                      fNumCrystalPreamp + 0.5);
        fh2_Preamp_vs_ch_L[i]->GetXaxis()->SetTitle("Preamp number [1-16]");
        std::stringstream ss4;
        ss4 << "Channel number [1-" << fNumCrystalPreamp << "]";
        fh2_Preamp_vs_ch_L[i]->GetYaxis()->SetTitle(ss4.str().c_str());
        fh2_Preamp_vs_ch_L[i]->GetYaxis()->SetTitleOffset(1.2);
        fh2_Preamp_vs_ch_L[i]->GetXaxis()->CenterTitle(true);
        fh2_Preamp_vs_ch_L[i]->GetYaxis()->CenterTitle(true);
        fh2_Preamp_vs_ch_L[i]->Draw("colz");

        // Right side
        std::stringstream ss5;
        ss5 << "Ring_" << i + 1 << "_Right";
        cMap_RingR[i] = new TCanvas(ss5.str().c_str(), ss5.str().c_str(), 10, 10, 800, 700);
        std::stringstream ss6;
        ss6 << "fh2_Ring_" << i + 1 << "_Right";
        std::stringstream ss7;
        ss7 << "Ring " << i + 1 << " right: Preamp number vs channel number";

        fh2_Preamp_vs_ch_R[i] = R3B::root_owned<TH2F>(ss6.str().c_str(),
                                                      ss7.str().c_str(),
                                                      fNumPreamps,
                                                      0.5,
                                                      fNumPreamps + 0.5,
                                                      fNumCrystalPreamp,
                                                      0.5,
                                                      fNumCrystalPreamp + 0.5);
        fh2_Preamp_vs_ch_R[i]->GetXaxis()->SetTitle("Preamp number [1-16]");
        fh2_Preamp_vs_ch_R[i]->GetYaxis()->SetTitle(ss4.str().c_str());
        fh2_Preamp_vs_ch_R[i]->GetYaxis()->SetTitleOffset(1.2);
        fh2_Preamp_vs_ch_R[i]->GetXaxis()->CenterTitle(true);
        fh2_Preamp_vs_ch_R[i]->GetYaxis()->CenterTitle(true);
        fh2_Preamp_vs_ch_R[i]->Draw("colz");
    }

    auto cMap_ECor = new TCanvas("Trigger_ECor", "Trigger_ECor", 10, 10, 800, 700);
    cMap_ECor->Divide(2, 2);

    fh2_Califa_EtrigCor[0] = R3B::root_owned<TH2F>(
        "fh2Trigger_ECor", "Correlation of trigger energies (all triggers)", 2000, 0., 4000., 2000, 0., 4000.);
    fh2_Califa_EtrigCor[0]->GetXaxis()->SetTitle("Energy messel side");
    fh2_Califa_EtrigCor[0]->GetYaxis()->SetTitle("Energy wixhausen side");
    fh2_Califa_EtrigCor[0]->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_EtrigCor[0]->GetXaxis()->CenterTitle(true);
    fh2_Califa_EtrigCor[0]->GetYaxis()->CenterTitle(true);
    cMap_ECor->cd(1);
    fh2_Califa_EtrigCor[0]->Draw("colz");

    fh1_Califa_Etrig[0] =
        R3B::root_owned<TH1F>("fh1Trigger_Emessel", "Messel trigger energies (all triggers)", 2000, 0., 4000.);
    fh1_Califa_Etrig[0]->GetXaxis()->SetTitle("Energy messel side");
    fh1_Califa_Etrig[0]->GetYaxis()->SetTitle("Counts");
    fh1_Califa_Etrig[0]->GetYaxis()->SetTitleOffset(1.2);
    fh1_Califa_Etrig[0]->GetXaxis()->CenterTitle(true);
    fh1_Califa_Etrig[0]->GetYaxis()->CenterTitle(true);
    fh1_Califa_Etrig[0]->SetFillColor(31);
    cMap_ECor->cd(2);
    fh1_Califa_Etrig[0]->Draw();

    fh1_Califa_Etrig[1] =
        R3B::root_owned<TH1F>("fh1Trigger_Ewixhausen", "Wixhausen trigger energies (all triggers)", 2000, 0., 4000.);
    fh1_Califa_Etrig[1]->GetXaxis()->SetTitle("Energy wixhausen side");
    fh1_Califa_Etrig[1]->GetYaxis()->SetTitle("Counts");
    fh1_Califa_Etrig[1]->GetYaxis()->SetTitleOffset(1.2);
    fh1_Califa_Etrig[1]->GetXaxis()->CenterTitle(true);
    fh1_Califa_Etrig[1]->GetYaxis()->CenterTitle(true);
    fh1_Califa_Etrig[1]->SetFillColor(31);
    cMap_ECor->cd(3);
    fh1_Califa_Etrig[1]->Draw();

    fh2_Califa_EtrigCor[1] = R3B::root_owned<TH2F>(
        "fh2Trigger_ECor_trg1", "Correlation of trigger energies (trigger 1)", 2000, 0., 4000., 2000, 0., 4000.);
    fh2_Califa_EtrigCor[1]->GetXaxis()->SetTitle("Energy messel side");
    fh2_Califa_EtrigCor[1]->GetYaxis()->SetTitle("Energy wixhausen side");
    fh2_Califa_EtrigCor[1]->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_EtrigCor[1]->GetXaxis()->CenterTitle(true);
    fh2_Califa_EtrigCor[1]->GetYaxis()->CenterTitle(true);
    cMap_ECor->cd(4);
    fh2_Califa_EtrigCor[1]->Draw("colz");

    std::vector<std::string> side = { "Right", "Left" };
    for (int s = 0; s < fNumSides; s++) // Side
    {

        for (Int_t r = 0; r < fNumRings; r++)       // Ring
            for (Int_t p = 0; p < fNumPreamps; p++) // Preamp
            {
                if (fFebexInfo[s][r][p][0] != -1)
                {
                    std::stringstream ss1;
                    ss1 << "Ring_" << r + 1 << "_" << side[s] << "_Preamp_" << p + 1;

                    std::stringstream ss2;
                    ss2 << "Ring " << r + 1 << ", " << side[s] << " side, Preamp " << p + 1;

                    cMapCry[s][r][p] = new TCanvas(ss1.str().c_str(), ss2.str().c_str(), 10, 10, 500, 500);
                    cMapCry[s][r][p]->Divide(4, 4);
                    // for TOT correlations
                    std::stringstream ss3;
                    ss3 << "Ring_" << r + 1 << "_" << side[s] << "_Preamp_" << p + 1 << "_tot";

                    cMapCryTot[s][r][p] = new TCanvas(ss3.str().c_str(), ss2.str().c_str(), 10, 10, 500, 500);
                    cMapCryTot[s][r][p]->Divide(4, 4);
                }
                if (fFebexInfo[s][r][p][2] != -1)
                {
                    std::stringstream ss1;
                    ss1 << "Ring_" << r + 1 << "_" << side[s] << "_Preamp_" << p + 1 << "_pr";

                    std::stringstream ss2;
                    ss2 << "Ring " << r + 1 << ", " << side[s] << " side, Preamp " << p + 1 << " for PR";

                    cMapCryP[s][r][p] = new TCanvas(ss1.str().c_str(), ss2.str().c_str(), 10, 10, 500, 500);
                    cMapCryP[s][r][p]->Divide(4, 4);
                    // for TOT correlations
                    std::stringstream ss3;
                    ss3 << "Ring_" << r + 1 << "_" << side[s] << "_Preamp_" << p + 1 << "_tot_pr";

                    std::stringstream ss4;
                    ss4 << "Ring " << r + 1 << ", " << side[s] << " side, Preamp " << p + 1 << " for ToT-PR";

                    cMapCryPTot[s][r][p] = new TCanvas(ss3.str().c_str(), ss4.str().c_str(), 10, 10, 500, 500);
                    cMapCryPTot[s][r][p]->Divide(4, 4);
                }
                for (Int_t j = 0; j < fNumCrystalPreamp; j++)
                { // Channel
                    Xaxis = "Energy [channels]";

                    if (fFebexInfo[s][r][p][0] != -1)
                    {
                        std::stringstream ss1;
                        ss1 << "fh1_Map_Side_" << side[s] << "_Ring_" << r + 1 << "_Preamp_" << p + 1 << "_Ch_" << j + 1
                            << "_energy";

                        // sprintf(Name1, "fh1_Map_Side_%s_Ring_%d_Preamp_%d_Ch_%d_energy", Side, r + 1, p + 1, j + 1);
                        // sprintf(Name2, "Map level, Side %s, Ring %d, Preamp %d, ch. %d", Side, r + 1, p + 1, j + 1);
                        std::stringstream ss2;
                        ss2 << "Map level, Side " << side[s] << ", Ring " << r + 1 << ", Preamp " << p + 1 << ", ch. "
                            << j + 1;

                        fh1_crystals[s][r][p][j] = R3B::root_owned<TH1F>(
                            ss1.str().c_str(), ss2.str().c_str(), fBinsChannelFebex, 0, fMaxBinChannelFebex);
                        fh1_crystals[s][r][p][j]->SetTitleSize(1.6, "t");
                        fh1_crystals[s][r][p][j]->GetXaxis()->SetTitle(Xaxis.c_str());
                        fh1_crystals[s][r][p][j]->GetXaxis()->SetLabelSize(0.06);
                        fh1_crystals[s][r][p][j]->GetYaxis()->SetLabelSize(0.07);
                        fh1_crystals[s][r][p][j]->GetXaxis()->SetTitleSize(0.05);
                        fh1_crystals[s][r][p][j]->GetXaxis()->CenterTitle(true);
                        fh1_crystals[s][r][p][j]->GetYaxis()->CenterTitle(true);
                        fh1_crystals[s][r][p][j]->GetXaxis()->SetTitleOffset(1.);
                        fh1_crystals[s][r][p][j]->SetFillColor(45);
                        cMapCry[s][r][p]->cd(j + 1);
                        gPad->SetLogy();
                        fh1_crystals[s][r][p][j]->Draw();

                        // for TOT correlations
                        std::stringstream ss3;
                        ss3 << "fh2_Map_Side_" << side[s] << "_Ring_" << r + 1 << "_Preamp_" << p + 1 << "_Ch_" << j + 1
                            << "_evstot";

                        fh2_crystalsETot[s][r][p][j] = R3B::root_owned<TH2F>(ss3.str().c_str(),
                                                                             ss2.str().c_str(),
                                                                             fBinsChannelFebex,
                                                                             0,
                                                                             fMaxBinChannelFebex,
                                                                             fBinsChannelFebex,
                                                                             0,
                                                                             fMaxBinChannelFebex / 3);
                        fh2_crystalsETot[s][r][p][j]->SetTitleSize(1.6, "t");
                        fh2_crystalsETot[s][r][p][j]->GetXaxis()->SetTitle(Xaxis.c_str());
                        fh2_crystalsETot[s][r][p][j]->GetYaxis()->SetTitle("Tot");
                        fh2_crystalsETot[s][r][p][j]->GetXaxis()->SetLabelSize(0.06);
                        fh2_crystalsETot[s][r][p][j]->GetYaxis()->SetLabelSize(0.06);
                        fh2_crystalsETot[s][r][p][j]->GetXaxis()->SetTitleSize(0.05);
                        fh2_crystalsETot[s][r][p][j]->GetXaxis()->CenterTitle(true);
                        fh2_crystalsETot[s][r][p][j]->GetYaxis()->CenterTitle(true);
                        fh2_crystalsETot[s][r][p][j]->GetXaxis()->SetTitleOffset(1.);
                        cMapCryTot[s][r][p]->cd(j + 1);
                        fh2_crystalsETot[s][r][p][j]->Draw("col");
                    }
                    if (fFebexInfo[s][r][p][2] != -1)
                    {
                        // histograms for proton range
                        std::stringstream ss1;
                        ss1 << "fh1_Map_Side_" << side[s] << "_Ring_" << r + 1 << "_Preamp_" << p + 1 << "_Ch_" << j + 1
                            << "_energy_pr";

                        std::stringstream ss2;
                        ss2 << "Map level (PR), Side " << side[s] << ", Ring " << r + 1 << ", Preamp " << p + 1
                            << ", ch. " << j + 1;

                        fh1_crystals_p[s][r][p][j] = R3B::root_owned<TH1F>(
                            ss1.str().c_str(), ss2.str().c_str(), fBinsChannelFebex, 0, fMaxBinChannelFebex);
                        fh1_crystals_p[s][r][p][j]->SetTitleSize(1.6, "t");
                        fh1_crystals_p[s][r][p][j]->GetXaxis()->SetTitle(Xaxis.c_str());
                        fh1_crystals_p[s][r][p][j]->GetXaxis()->SetLabelSize(0.06);
                        fh1_crystals_p[s][r][p][j]->GetYaxis()->SetLabelSize(0.07);
                        fh1_crystals_p[s][r][p][j]->GetXaxis()->SetTitleSize(0.05);
                        fh1_crystals_p[s][r][p][j]->GetXaxis()->CenterTitle(true);
                        fh1_crystals_p[s][r][p][j]->GetYaxis()->CenterTitle(true);
                        fh1_crystals_p[s][r][p][j]->GetXaxis()->SetTitleOffset(1.);
                        fh1_crystals_p[s][r][p][j]->SetFillColor(45);
                        cMapCryP[s][r][p]->cd(j + 1);
                        gPad->SetLogy();
                        fh1_crystals_p[s][r][p][j]->Draw();

                        // for TOT correlations
                        std::stringstream ss3;
                        ss3 << "fh2_Map_Side_" << side[s] << "_Ring_" << r + 1 << "_Preamp_" << p + 1 << "_Ch_" << j + 1
                            << "_pr_evstot";
                        fh2_crystalsETot_p[s][r][p][j] = R3B::root_owned<TH2F>(ss3.str().c_str(),
                                                                               ss2.str().c_str(),
                                                                               fBinsChannelFebex,
                                                                               0,
                                                                               fMaxBinChannelFebex,
                                                                               fBinsChannelFebex,
                                                                               0,
                                                                               fMaxBinChannelFebex / 3);
                        fh2_crystalsETot_p[s][r][p][j]->SetTitleSize(1.6, "t");
                        fh2_crystalsETot_p[s][r][p][j]->GetXaxis()->SetTitle(Xaxis.c_str());
                        fh2_crystalsETot_p[s][r][p][j]->GetYaxis()->SetTitle("Tot");
                        fh2_crystalsETot_p[s][r][p][j]->GetXaxis()->SetLabelSize(0.06);
                        fh2_crystalsETot_p[s][r][p][j]->GetYaxis()->SetLabelSize(0.06);
                        fh2_crystalsETot_p[s][r][p][j]->GetXaxis()->SetTitleSize(0.05);
                        fh2_crystalsETot_p[s][r][p][j]->GetXaxis()->CenterTitle(true);
                        fh2_crystalsETot_p[s][r][p][j]->GetYaxis()->CenterTitle(true);
                        fh2_crystalsETot_p[s][r][p][j]->GetXaxis()->SetTitleOffset(1.);
                        cMapCryPTot[s][r][p]->cd(j + 1);
                        fh2_crystalsETot_p[s][r][p][j]->Draw("col");
                    }
                }
            }
    }

    // Cal data
    for (Int_t s = 0; s < fNumSides; s++) // Side
    {
        for (Int_t r = 0; r < fNumRings; r++) // Ring

            for (Int_t p = 0; p < fNumPreamps; p++) // Preamp
            {
                if (fFebexInfo[s][r][p][0] != -1)
                {
                    std::stringstream ss1;
                    ss1 << "Ring_" << r + 1 << "_" << side[s] << "_Preamp_" << p + 1 << "_cal";

                    std::stringstream ss2;
                    ss2 << "Ring " << r + 1 << ", " << side[s] << " side, Preamp " << p + 1 << " for Cal";

                    cMapCryCal[s][r][p] = new TCanvas(ss1.str().c_str(), ss2.str().c_str(), 10, 10, 500, 500);
                    cMapCryCal[s][r][p]->Divide(4, 4);
                }
                if (fFebexInfo[s][r][p][2] != -1)
                {
                    std::stringstream ss1;
                    ss1 << "Ring_" << r + 1 << "_" << side[s] << "_Preamp_" << p + 1 << "_calpr";

                    std::stringstream ss2;
                    ss2 << "Ring " << r + 1 << ", " << side[s] << " side, Preamp " << p + 1 << "for PR Cal";

                    cMapCryPCal[s][r][p] = new TCanvas(ss1.str().c_str(), ss2.str().c_str(), 10, 10, 500, 500);
                    cMapCryPCal[s][r][p]->Divide(4, 4);
                }
                for (Int_t j = 0; j < fNumCrystalPreamp; j++)
                { // Channel
                    Xaxis = "Energy [keV]";

                    if (fFebexInfo[s][r][p][0] != -1)
                    {
                        std::stringstream ss1;
                        ss1 << "fh1_Cal_Side_" << side[s] << "_Ring_" << r + 1 << "_Preamp_" << p + 1 << "_Ch_" << j + 1
                            << "_energy";

                        std::stringstream ss2;
                        ss2 << "Cal level, Side " << side[s] << ", Ring " << r + 1 << ", Preamp " << p + 1 << ", ch. "
                            << j + 1;

                        fh1_crystals_cal[s][r][p][j] = R3B::root_owned<TH1F>(ss1.str().c_str(),
                                                                             ss2.str().c_str(),
                                                                             arry_bins[s][r][p][j],
                                                                             arry_minE[s][r][p][j],
                                                                             arry_maxE[s][r][p][j]);
                        fh1_crystals_cal[s][r][p][j]->SetTitleSize(1.6, "t");
                        fh1_crystals_cal[s][r][p][j]->GetXaxis()->SetTitle(Xaxis.c_str());
                        fh1_crystals_cal[s][r][p][j]->GetXaxis()->SetLabelSize(0.06);
                        fh1_crystals_cal[s][r][p][j]->GetYaxis()->SetLabelSize(0.07);
                        fh1_crystals_cal[s][r][p][j]->GetXaxis()->SetTitleSize(0.05);
                        fh1_crystals_cal[s][r][p][j]->GetXaxis()->CenterTitle(true);
                        fh1_crystals_cal[s][r][p][j]->GetYaxis()->CenterTitle(true);
                        fh1_crystals_cal[s][r][p][j]->GetXaxis()->SetTitleOffset(1.);
                        fh1_crystals_cal[s][r][p][j]->SetFillColor(45);
                        cMapCryCal[s][r][p]->cd(j + 1);
                        gPad->SetLogy();
                        fh1_crystals_cal[s][r][p][j]->Draw();
                    }
                    if (fFebexInfo[s][r][p][2] != -1)
                    {
                        // histograms for proton range
                        std::stringstream ss1;
                        ss1 << "fh1_Cal_Side_" << side[s] << "_Ring_" << r + 1 << "_Preamp_" << p + 1 << "_Ch_" << j + 1
                            << "_energy_pr";

                        std::stringstream ss2;
                        ss2 << "Cal level (PR), Side" << side[s] << ", Ring " << r + 1 << ", Preamp " << p + 1
                            << ", ch. " << j + 1;

                        fh1_crystals_p_cal[s][r][p][j] = R3B::root_owned<TH1F>(ss1.str().c_str(),
                                                                               ss2.str().c_str(),
                                                                               arry_bins[s][r][p][j],
                                                                               arry_minE[s][r][p][j],
                                                                               arry_maxE[s][r][p][j]);
                        fh1_crystals_p_cal[s][r][p][j]->SetTitleSize(1.6, "t");
                        fh1_crystals_p_cal[s][r][p][j]->GetXaxis()->SetTitle(Xaxis.c_str());
                        fh1_crystals_p_cal[s][r][p][j]->GetXaxis()->SetLabelSize(0.06);
                        fh1_crystals_p_cal[s][r][p][j]->GetYaxis()->SetLabelSize(0.07);
                        fh1_crystals_p_cal[s][r][p][j]->GetXaxis()->SetTitleSize(0.05);
                        fh1_crystals_p_cal[s][r][p][j]->GetXaxis()->CenterTitle(true);
                        fh1_crystals_p_cal[s][r][p][j]->GetYaxis()->CenterTitle(true);
                        fh1_crystals_p_cal[s][r][p][j]->GetXaxis()->SetTitleOffset(1.);
                        fh1_crystals_p_cal[s][r][p][j]->SetFillColor(45);
                        cMapCryPCal[s][r][p]->cd(j + 1);
                        gPad->SetLogy();
                        fh1_crystals_p_cal[s][r][p][j]->Draw();
                    }
                }
            }
    }

    // CANVAS Multiplicity
    cCalifaMult = new TCanvas("Califa_Multiplicity", "Califa_Multiplicity", 10, 10, 500, 500);
    fh1_Califa_Mult =
        R3B::root_owned<TH1F>("fh1_Califa_Mult", "Califa multiplicity (crystal:blue, cluster:red)", 341, -0.5, 340.5);
    fh1_Califa_MultHit = R3B::root_owned<TH1F>("fh1_Califa_MultHit", "Califa multiplicity", 341, -0.5, 340.5);
    fh1_Califa_Mult->GetXaxis()->SetTitle("Multiplicity");
    fh1_Califa_Mult->GetXaxis()->CenterTitle(true);
    fh1_Califa_Mult->GetYaxis()->CenterTitle(true);
    fh1_Califa_Mult->GetYaxis()->SetTitleOffset(1.2);
    fh1_Califa_Mult->GetXaxis()->SetTitleOffset(1.1);
    gPad->SetLogy();
    fh1_Califa_Mult->Draw();
    fh1_Califa_MultHit->SetLineColor(kRed);
    fh1_Califa_MultHit->Draw("SAME");

    // CANVAS Energy correlations between hits
    cCalifaCoinE = new TCanvas("Califa_energy_correlation_hits", "Energy correlations, hit level", 10, 10, 500, 500);
    cCalifaCoinE->Divide(2, 1);
    cCalifaCoinE->cd(1);

    fh2_Califa_coinE = R3B::root_owned<TH2F>("fh2_Califa_energy_correlations",
                                             "Califa energy correlations",
                                             (maxE - minE) / 1000.,
                                             minE / 1000.,
                                             maxE / 1000.,
                                             (maxE - minE) / 1000.,
                                             minE / 1000.,
                                             maxE / 1000.);
    fh2_Califa_coinE->GetXaxis()->SetTitle("Energy (MeV)");
    fh2_Califa_coinE->GetYaxis()->SetTitle("Energy (MeV)");
    fh2_Califa_coinE->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_coinE->GetXaxis()->CenterTitle(true);
    fh2_Califa_coinE->GetYaxis()->CenterTitle(true);
    fh2_Califa_coinE->Draw("COLZ");

    cCalifaCoinE->cd(2);
    fh2_Califa_coinE_p2p = R3B::root_owned<TH2F>("fh2_Califa_energy_correlations_p2p",
                                                 "Califa energy correlations for p2p",
                                                 (maxE - minE) / 1000.,
                                                 minE / 1000.,
                                                 maxE / 1000.,
                                                 (maxE - minE) / 1000.,
                                                 minE / 1000.,
                                                 maxE / 1000.);
    fh2_Califa_coinE_p2p->GetXaxis()->SetTitle("Energy (MeV)");
    fh2_Califa_coinE_p2p->GetYaxis()->SetTitle("Energy (MeV)");
    fh2_Califa_coinE_p2p->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_coinE_p2p->GetXaxis()->CenterTitle(true);
    fh2_Califa_coinE_p2p->GetYaxis()->CenterTitle(true);
    fh2_Califa_coinE_p2p->Draw("COLZ");

    // CANVAS Theta correlations between hits
    cCalifaCoinTheta = new TCanvas("Califa_theta_correlation_hits", "Theta correlations, hit level", 10, 10, 500, 500);
    cCalifaCoinTheta->Divide(2, 1);
    cCalifaCoinTheta->cd(1);
    fh2_Califa_coinTheta =
        R3B::root_owned<TH2F>("fh2_Califa_theta_correlations", "Califa theta correlations", 500, 0, 100, 500, 0, 100);
    fh2_Califa_coinTheta->GetXaxis()->SetTitle("Theta [deg]");
    fh2_Califa_coinTheta->GetYaxis()->SetTitle("Theta [deg]");
    fh2_Califa_coinTheta->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_coinTheta->GetXaxis()->CenterTitle(true);
    fh2_Califa_coinTheta->GetYaxis()->CenterTitle(true);
    fh2_Califa_coinTheta->Draw("COLZ");

    cCalifaCoinTheta->cd(2);
    fh2_Califa_coinTheta_cutOPA = R3B::root_owned<TH2F>(
        "fh2_Califa_theta_correlations_cutOPA", "Califa theta correlations OPA cut", 500, 0, 100, 500, 0, 100);
    fh2_Califa_coinTheta_cutOPA->GetXaxis()->SetTitle("Theta [deg]");
    fh2_Califa_coinTheta_cutOPA->GetYaxis()->SetTitle("Theta [deg]");
    fh2_Califa_coinTheta_cutOPA->SetTitleOffset(1.2);
    fh2_Califa_coinTheta_cutOPA->GetXaxis()->CenterTitle(true);
    fh2_Califa_coinTheta_cutOPA->GetYaxis()->CenterTitle(true);
    fh2_Califa_coinTheta_cutOPA->Draw("colz");

    // CANVAS Phi correlations between hits
    cCalifaCoinPhi = new TCanvas("Califa_phi_correlation_hits", "Phi correlations, hit level", 10, 10, 500, 500);

    fh2_Califa_coinPhi =
        R3B::root_owned<TH2F>("fh2_Califa_phi_correlations", "Califa phi correlations", 600, -190, 190, 600, -190, 190);
    fh2_Califa_coinPhi->GetXaxis()->SetTitle("Phi [deg]");
    fh2_Califa_coinPhi->GetYaxis()->SetTitle("Phi [deg]");
    fh2_Califa_coinPhi->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_coinPhi->GetXaxis()->CenterTitle(true);
    fh2_Califa_coinPhi->GetYaxis()->CenterTitle(true);
    fh2_Califa_coinPhi->Draw("COLZ");

    // CANVAS Theta vs Phi
    cCalifa_angles = new TCanvas("Califa_Theta_vs_Phi", "Theta vs Phi", 10, 10, 500, 500);
    fh2_Califa_theta_phi =
        R3B::root_owned<TH2F>("fh2_Califa_theta_vs_phi", "Califa theta vs phi", 90, 0, 90, 380, -190, 190);
    fh2_Califa_theta_phi->GetXaxis()->SetTitle("Theta [deg]");
    fh2_Califa_theta_phi->GetYaxis()->SetTitle("Phi [deg]");
    fh2_Califa_theta_phi->GetYaxis()->SetTitleOffset(1.2);
    fh2_Califa_theta_phi->GetXaxis()->CenterTitle(true);
    fh2_Califa_theta_phi->GetYaxis()->CenterTitle(true);
    fh2_Califa_theta_phi->Draw("COLZ");

    // CANVAS Theta vs energy
    Name1 = "Califa_calorimeter_theta_vs_energy";
    Name2 = "fh_Califa_theta_vs_total_energy";
    Name3 = "Califa theta vs energy for full calorimeter";
    cCalifa_theta_energy = new TCanvas(Name1.c_str(), Name1.c_str(), 10, 10, 500, 500);
    fh2_Califa_theta_energy = R3B::root_owned<TH2F>(Name2.c_str(), Name3.c_str(), 500, 0, 90, bins, minE, maxE);
    fh2_Califa_theta_energy->GetXaxis()->SetTitle("Theta [deg]");
    fh2_Califa_theta_energy->GetYaxis()->SetTitle("Energy [keV]");
    fh2_Califa_theta_energy->GetYaxis()->SetTitleOffset(1.4);
    fh2_Califa_theta_energy->GetXaxis()->CenterTitle(true);
    fh2_Califa_theta_energy->GetYaxis()->CenterTitle(true);
    fh2_Califa_theta_energy->Draw("COLZ");

    // CANVAS Total energy
    Name1 = "Califa_calorimeter_total_energy_per_hit";
    Name2 = "fh_Califa_total_energy";
    Name3 = "Califa total energy per hit for the full calorimeter";
    cCalifa_hitenergy = new TCanvas(Name1.c_str(), Name1.c_str(), 10, 10, 500, 500);
    fh1_Califa_total_energy = R3B::root_owned<TH1F>(Name2.c_str(), Name3.c_str(), bins, minE, maxE);
    fh1_Califa_total_energy->GetXaxis()->SetTitle("Energy [keV]");
    fh1_Califa_total_energy->GetYaxis()->SetTitle("Counts");
    fh1_Califa_total_energy->GetYaxis()->SetTitleOffset(1.4);
    fh1_Califa_total_energy->GetXaxis()->CenterTitle(true);
    fh1_Califa_total_energy->GetYaxis()->CenterTitle(true);
    fh1_Califa_total_energy->SetFillColor(29);
    fh1_Califa_total_energy->SetLineColor(1);
    fh1_Califa_total_energy->SetLineWidth(2);
    fh1_Califa_total_energy->Draw("");
    gPad->SetLogy();

    // CANVAS opening angle
    Name1 = "Califa_opening_angle_hit";
    Name2 = "fh1_Califa_opening";
    Name3 = "Califa opening angle";
    cCalifa_opening = new TCanvas(Name1.c_str(), Name1.c_str(), 10, 10, 500, 500);
    fh1_openangle = R3B::root_owned<TH1F>(Name2.c_str(), Name3.c_str(), 160, 10, 170);
    fh1_openangle->GetXaxis()->SetTitle("Opening angle [deg]");
    fh1_openangle->GetYaxis()->SetTitle("Counts");
    fh1_openangle->GetXaxis()->CenterTitle(true);
    fh1_openangle->GetYaxis()->CenterTitle(true);
    fh1_openangle->GetYaxis()->SetTitleOffset(1.2);
    fh1_openangle->GetXaxis()->SetTitleOffset(1.2);
    fh1_openangle->SetFillColor(8);
    fh1_openangle->SetLineColor(1);
    fh1_openangle->SetLineWidth(2);
    fh1_openangle->Draw("");

    // CANVAS opening angle
    Name1 = "Califa_opening_angle_tpats";
    Name2 = "fh1_Califa_opening_tpats";
    Name3 = "Califa opening angle vs tpat";
    cCalifa_opening_tpat = new TCanvas(Name1.c_str(), Name1.c_str(), 10, 10, 500, 500);
    fh2_openangle_tpat = R3B::root_owned<TH2F>(Name2.c_str(), Name3.c_str(), 17, -0.5, 16.5, 160, 10, 170);
    fh2_openangle_tpat->GetXaxis()->SetTitle("Tpat");
    fh2_openangle_tpat->GetYaxis()->SetTitle("Opening angle [deg]");
    fh2_openangle_tpat->GetXaxis()->CenterTitle(true);
    fh2_openangle_tpat->GetYaxis()->CenterTitle(true);
    fh2_openangle_tpat->GetYaxis()->SetTitleOffset(1.2);
    fh2_openangle_tpat->GetXaxis()->SetTitleOffset(1.2);
    fh2_openangle_tpat->Draw("colz");

    // Difference between Califa WRs
    Name1 = "WR_Califa";
    Name2 = "fh1_WR_Califa";
    Name3 = "WR-Wixhausen - WR-Messel";
    cCalifa_wr = new TCanvas(Name1.c_str(), Name1.c_str(), 10, 10, 500, 500);
    fh1_Califa_wr = R3B::root_owned<TH1I>(Name2.c_str(), Name3.c_str(), 4000, -4000, 4000);
    fh1_Califa_wr->GetXaxis()->SetTitle("Difference of Califa WRs");
    fh1_Califa_wr->GetYaxis()->SetTitle("Counts");
    fh1_Califa_wr->GetYaxis()->SetTitleOffset(1.3);
    fh1_Califa_wr->GetXaxis()->CenterTitle(true);
    fh1_Califa_wr->GetYaxis()->CenterTitle(true);
    fh1_Califa_wr->SetFillColor(29);
    fh1_Califa_wr->SetLineColor(1);
    fh1_Califa_wr->SetLineWidth(2);
    fh1_Califa_wr->Draw("");

    // Difference between Califa-Master WRs
    Name1 = "WR_Master_Califa";
    cWrs = new TCanvas(Name1.c_str(), Name1.c_str(), 10, 10, 500, 500);
    Name2 = "fh1_WR_Master_Califa";
    Name3 = "WR-Califa - WR-Master: Messel (blue), Wixhausen (red) - WR-Master";
    fh1_wrs[0] = R3B::root_owned<TH1I>(Name2.c_str(), Name3.c_str(), 4000, -4000, 4000);
    fh1_wrs[0]->SetStats(1);
    fh1_wrs[0]->GetXaxis()->SetTitle("WRTs difference");
    fh1_wrs[0]->GetYaxis()->SetTitle("Counts");
    fh1_wrs[0]->GetYaxis()->SetTitleOffset(1.3);
    fh1_wrs[0]->GetXaxis()->CenterTitle(true);
    fh1_wrs[0]->GetYaxis()->CenterTitle(true);
    fh1_wrs[0]->SetLineColor(2);
    fh1_wrs[0]->SetLineWidth(3);
    fh1_wrs[0]->Draw("");
    fh1_wrs[1] = R3B::root_owned<TH1I>("fh1_WR_Master_Califa_Messel", "", 4000, -4000, 4000);
    fh1_wrs[1]->SetLineColor(4);
    fh1_wrs[1]->SetLineWidth(3);
    fh1_wrs[1]->Draw("same");

    // CANVAS energy vs wrs
    Name1 = "Califa_wr_vs_energy";
    Name2 = "fh2_wr_vs_energy_left";
    Name3 = "Califa WR vs hit-energy left side";
    cCalifa_wr_energy = new TCanvas(Name1.c_str(), Name1.c_str(), 10, 10, 500, 500);
    cCalifa_wr_energy->Divide(1, 2);
    fh2_Cal_wr_energy_l = R3B::root_owned<TH2F>(Name2.c_str(), Name3.c_str(), 700, -4000, 4000, bins, minE, maxE);
    fh2_Cal_wr_energy_l->GetXaxis()->SetTitle("WR difference (Master-Califa)");
    fh2_Cal_wr_energy_l->GetYaxis()->SetTitle("Energy [keV]");
    fh2_Cal_wr_energy_l->GetYaxis()->SetTitleOffset(1.4);
    fh2_Cal_wr_energy_l->GetXaxis()->CenterTitle(true);
    fh2_Cal_wr_energy_l->GetYaxis()->CenterTitle(true);
    cCalifa_wr_energy->cd(1);
    fh2_Cal_wr_energy_l->Draw("COLZ");
    Name2 = "fh2_wr_vs_energy_right";
    Name3 = "Califa WR vs hit-energy right side";
    fh2_Cal_wr_energy_r = R3B::root_owned<TH2F>(Name2.c_str(), Name3.c_str(), 700, -4000, 4000, bins, minE, maxE);
    fh2_Cal_wr_energy_r->GetXaxis()->SetTitle("WR difference (Master-Califa)");
    fh2_Cal_wr_energy_r->GetYaxis()->SetTitle("Energy [keV]");
    fh2_Cal_wr_energy_r->GetYaxis()->SetTitleOffset(1.4);
    fh2_Cal_wr_energy_r->GetXaxis()->CenterTitle(true);
    fh2_Cal_wr_energy_r->GetYaxis()->CenterTitle(true);
    cCalifa_wr_energy->cd(2);
    fh2_Cal_wr_energy_r->Draw("COLZ");

    cCalifaTriggers = new TCanvas("Califa_triggers", "Califa_triggers", 10, 10, 500, 500);
    fh1_CalifaTriggers = R3B::root_owned<TH1F>("califa_triggers", "califa_triggers", 3, -0.5, 2.5);
    fh1_CalifaTriggers->GetXaxis()->SetTitle("Triggers");
    fh1_CalifaTriggers->GetXaxis()->CenterTitle(true);
    fh1_CalifaTriggers->SetFillColor(29);
    fh1_CalifaTriggers->SetLineColor(1);
    fh1_CalifaTriggers->SetLineWidth(2);
    fh1_CalifaTriggers->Draw();

    // FOLDERS for Califa
    TFolder* folder_sta = new TFolder("Statistics_per_ring", "Statistics info");
    for (Int_t i = 1; i < fNumRings; i++)
    { // FIXME in the future
        folder_sta->Add(cMap_RingR[i]);
        folder_sta->Add(cMap_RingL[i]);
    }

    TFolder* folder_el = new TFolder("Energy_Map_per_crystal_Left", "Energy per crystal, left side info");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[1][r][p][0] != -1)
                folder_el->Add(cMapCry[1][r][p]);

    TFolder* folder_etotl = new TFolder("Energy_Tot_per_crystal_Left", "Energy vs Tot per crystal, left side info");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[1][r][p][0] != -1)
                folder_etotl->Add(cMapCryTot[1][r][p]);

    TFolder* folder_eprl = new TFolder("Energy_Map_per_crystal_Left_PR", "Energy per crystal, left side info for PR");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[1][r][p][2] != -1)
                folder_eprl->Add(cMapCryP[1][r][p]);

    TFolder* folder_eprtotl =
        new TFolder("Energy_Tot_per_crystal_Left_PR", "Energy vs Tot per crystal, left side info for PR");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[1][r][p][2] != -1)
                folder_eprtotl->Add(cMapCryPTot[1][r][p]);

    TFolder* folder_er = new TFolder("Energy_Map_per_crystal_Right", "Energy per crystal, right side info");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[0][r][p][0] != -1)
                folder_er->Add(cMapCry[0][r][p]);

    TFolder* folder_etotr = new TFolder("Energy_Tot_per_crystal_Right", "Energy vs Tot per crystal, right side info");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[0][r][p][0] != -1)
                folder_etotr->Add(cMapCryTot[0][r][p]);

    TFolder* folder_eprr = new TFolder("Energy_Map_per_crystal_Right_PR", "Energy per crystal, right side info for PR");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[0][r][p][2] != -1)
                folder_eprr->Add(cMapCryP[0][r][p]);

    TFolder* folder_eprtotr =
        new TFolder("Energy_Tot_per_crystal_Right_PR", "Energy vs Tot per crystal, right side info for PR");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[0][r][p][2] != -1)
                folder_eprtotr->Add(cMapCryPTot[0][r][p]);

    TFolder* folder_ecall = new TFolder("Energy_Cal_per_crystal_Left", "Energy Cal per crystal, left side info");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[1][r][p][0] != -1)
                folder_ecall->Add(cMapCryCal[1][r][p]);

    TFolder* folder_eprcall =
        new TFolder("Energy_Cal_per_crystal_Left_PR", "Energy Cal per crystal, left side info for PR");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[1][r][p][2] != -1)
                folder_eprcall->Add(cMapCryPCal[1][r][p]);

    TFolder* folder_ecalr = new TFolder("Energy_Cal_per_crystal_Right", "Energy Cal per crystal, right side info");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[0][r][p][0] != -1)
                folder_ecalr->Add(cMapCryCal[0][r][p]);

    TFolder* folder_eprcalr =
        new TFolder("Energy_Cal_per_crystal_Right_PR", "Energy Cal per crystal, right side info for PR");
    for (Int_t r = 0; r < fNumRings; r++)
        for (Int_t p = 0; p < fNumPreamps; p++)
            if (fFebexInfo[0][r][p][2] != -1)
                folder_eprcalr->Add(cMapCryPCal[0][r][p]);

    // MAIN FOLDER-Califa
    TFolder* mainfolCalifa = new TFolder("CALIFA", "CALIFA info");
    mainfolCalifa->Add(cCalifaMult);
    mainfolCalifa->Add(cCalifa_cry_energy);

    TFolder* folder_wrs = new TFolder("WRs", "CALIFA white-rabbit info");

    if (fWRItemsCalifa && fWRItemsMaster)
    {
        folder_wrs->Add(cWrs);
    }
    if (fTrigMappedItemsCalifa)
    {
        folder_wrs->Add(cMap_ECor);
    }

    if (fWRItemsCalifa)
    {
        folder_wrs->Add(cCalifa_wr);
    }
    if (fWRItemsCalifa)
        mainfolCalifa->Add(folder_wrs);

    mainfolCalifa->Add(folder_sta);
    mainfolCalifa->Add(folder_el);
    mainfolCalifa->Add(folder_eprl);
    mainfolCalifa->Add(folder_er);
    mainfolCalifa->Add(folder_eprr);

    if (fTotHist)
    {
        mainfolCalifa->Add(folder_etotl);
        mainfolCalifa->Add(folder_eprtotl);
        mainfolCalifa->Add(folder_etotr);
        mainfolCalifa->Add(folder_eprtotr);
    }

    if (fCalItemsCalifa)
    {
        mainfolCalifa->Add(cCalifa_cry_energy_cal);
        mainfolCalifa->Add(cCalifa_NsNf);
        mainfolCalifa->Add(folder_ecall);
        mainfolCalifa->Add(folder_eprcall);
        mainfolCalifa->Add(folder_ecalr);
        mainfolCalifa->Add(folder_eprcalr);
    }
    if (fHitItemsCalifa)
    {
        mainfolCalifa->Add(cCalifaCoinE);
        mainfolCalifa->Add(cCalifaCoinTheta);
        mainfolCalifa->Add(cCalifaCoinPhi);
        mainfolCalifa->Add(cCalifa_angles);
        mainfolCalifa->Add(cCalifa_opening);
        mainfolCalifa->Add(cCalifa_opening_tpat);
        mainfolCalifa->Add(cCalifa_theta_energy);
        mainfolCalifa->Add(cCalifa_hitenergy);
        if (fWRItemsCalifa && fWRItemsMaster)
            mainfolCalifa->Add(cCalifa_wr_energy);
        mainfolCalifa->Add(cCalifaTriggers);
    }
    run->AddObject(mainfolCalifa);

    // Register command to reset histograms
    run->GetHttpServer()->RegisterCommand("Reset_Califa", Form("/Objects/%s/->Reset_CALIFA_Histo()", GetName()));
    // Register command for moving between Febex and Preamp channels
    run->GetHttpServer()->RegisterCommand("Febex2Preamp_Califa",
                                          Form("/Objects/%s/->Febex2Preamp_CALIFA_Histo()", GetName()));
    // Register command to change the histogram scales (Log/Lineal)
    run->GetHttpServer()->RegisterCommand("Log_Califa", Form("/Objects/%s/->Log_CALIFA_Histo()", GetName()));

    return kSUCCESS;
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BCalifaOnlineSpectra::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

void R3BCalifaOnlineSpectra::Reset_CALIFA_Histo()
{
    R3BLOG(info, "");

    if (fWRItemsCalifa)
    {
        fh1_Califa_wr->Reset();
    }

    if (fWRItemsCalifa && fWRItemsMaster)
    {
        fh1_wrs[0]->Reset();
        fh1_wrs[1]->Reset();
        if (fHitItemsCalifa)
        {
            fh2_Cal_wr_energy_r->Reset();
            fh2_Cal_wr_energy_l->Reset();
        }
    }

    if (fTrigMappedItemsCalifa)
    {
        fh2_Califa_EtrigCor[0]->Reset();
        fh2_Califa_EtrigCor[1]->Reset();
        fh1_Califa_Etrig[0]->Reset();
        fh1_Califa_Etrig[1]->Reset();
    }

    if (fMappedItemsCalifa)
    {
        fh1_Califa_Mult->Reset();
        // for (Int_t s = 0; s < 3; s++)
        //   fh1_Califa_sync[s]->Reset();
        fh2_Califa_cryId_energy->Reset();
        for (Int_t i = 0; i < fNumRings; i++)
        {
            fh2_Preamp_vs_ch_R[i]->Reset();
            fh2_Preamp_vs_ch_L[i]->Reset();
        }
        for (Int_t s = 0; s < fNumSides; s++)
            for (Int_t r = 0; r < fNumRings; r++)
                for (Int_t p = 0; p < fNumPreamps; p++)
                    for (Int_t ch = 0; ch < fNumCrystalPreamp; ch++)
                    {
                        if (fFebexInfo[s][r][p][0] != -1)
                        {
                            fh1_crystals[s][r][p][ch]->Reset();
                            fh2_crystalsETot[s][r][p][ch]->Reset();
                        }
                        if (fFebexInfo[s][r][p][2] != -1)
                        {
                            fh1_crystals_p[s][r][p][ch]->Reset();
                            fh2_crystalsETot_p[s][r][p][ch]->Reset();
                        }
                    }
    }

    if (fCalItemsCalifa)
    {
        fh2_Califa_cryId_energy_cal->Reset();
        for (const auto& hist : fh2_Califa_NsNf)
        {
            hist->Reset();
        }
        for (Int_t s = 0; s < fNumSides; s++)
            for (Int_t r = 0; r < fNumRings; r++)
                for (Int_t p = 0; p < fNumPreamps; p++)
                    for (Int_t ch = 0; ch < fNumCrystalPreamp; ch++)
                    {
                        if (fFebexInfo[s][r][p][0] != -1)
                            fh1_crystals_cal[s][r][p][ch]->Reset();
                        if (fFebexInfo[s][r][p][2] != -1)
                            fh1_crystals_p_cal[s][r][p][ch]->Reset();
                    }
    }

    if (fHitItemsCalifa)
    {
        fh1_Califa_MultHit->Reset();
        fh2_Califa_coinE->Reset();
        fh2_Califa_coinE_p2p->Reset();
        fh2_Califa_coinTheta->Reset();
        fh2_Califa_coinTheta_cutOPA->Reset();
        fh2_Califa_coinPhi->Reset();
        fh2_Califa_theta_phi->Reset();
        fh2_Califa_theta_energy->Reset();
        fh1_Califa_total_energy->Reset();
        fh1_openangle->Reset();
        fh2_openangle_tpat->Reset();
        fh1_CalifaTriggers->Reset();
    }
}

void R3BCalifaOnlineSpectra::Log_CALIFA_Histo()
{
    R3BLOG(info, "");
    cCalifa_cry_energy->cd();
    if (fLogScale)
    {
        gPad->SetLogz(0);
    }
    else
    {
        gPad->SetLogz(1);
    }

    cCalifaMult->cd();
    if (fLogScale)
    {
        gPad->SetLogy(0);
    }
    else
    {
        gPad->SetLogy(1);
    }

    for (Int_t s = 0; s < fNumSides; s++)
    {
        for (Int_t r = 0; r < fNumRings; r++)
        {
            for (Int_t p = 0; p < fNumPreamps; p++)
            {

                for (Int_t j = 0; j < fNumCrystalPreamp; j++)
                {
                    if (fFebexInfo[s][r][p][0] != -1)
                    {
                        cMapCry[s][r][p]->cd(j + 1);
                        if (fLogScale)
                        {
                            gPad->SetLogy(0);
                        }
                        else
                        {
                            gPad->SetLogy(1);
                        }
                    }
                    if (fFebexInfo[s][r][p][2] != -1)
                    { // histograms for proton range
                        cMapCryP[s][r][p]->cd(j + 1);
                        if (fLogScale)
                        {
                            gPad->SetLogy(0);
                        }
                        else
                        {
                            gPad->SetLogy(1);
                        }
                    }
                    if (fCalItemsCalifa)
                    {
                        if (fFebexInfo[s][r][p][0] != -1)
                        {
                            cMapCryCal[s][r][p]->cd(p + 1);
                            if (fLogScale)
                            {
                                gPad->SetLogy(0);
                            }
                            else
                            {
                                gPad->SetLogy(1);
                            }
                        }
                        if (fFebexInfo[s][r][p][2] != -1)
                        { // histograms for proton range
                            cMapCryPCal[s][r][p]->cd(j + 1);
                            if (fLogScale)
                            {
                                gPad->SetLogy(0);
                            }
                            else
                            {
                                gPad->SetLogy(1);
                            }
                        }
                    }
                }
            }
        }
    }

    if (fCalItemsCalifa)
    {
        cCalifa_cry_energy_cal->cd();
        if (fLogScale)
        {
            gPad->SetLogz(0);
        }
        else
        {
            gPad->SetLogz(1);
        }
    }

    if (fHitItemsCalifa)
    {
        cCalifa_hitenergy->cd();
        if (fLogScale)
        {
            gPad->SetLogy(0);
        }
        else
        {
            gPad->SetLogy(1);
        }
    }

    if (fLogScale)
        fLogScale = kFALSE;
    else
        fLogScale = kTRUE;
}

void R3BCalifaOnlineSpectra::Febex2Preamp_CALIFA_Histo()
{
    R3BLOG(info, "");

    char Name[255];
    char Side[50];

    if (fFebex2Preamp)
    { // Preamp to Febex sequence
        for (Int_t s = 0; s < fNumSides; s++)
        {
            if (s == 1)
                sprintf(Side, "Left");
            else
                sprintf(Side, "Right");
            for (Int_t r = 0; r < fNumRings; r++)
            {
                for (Int_t p = 0; p < fNumPreamps; p++)
                {
                    for (Int_t j = 0; j < fNumCrystalPreamp; j++)
                    {
                        if (fFebexInfo[s][r][p][0] != -1)
                        {
                            cMapCry[s][r][p]->cd(fOrderFebexPreamp[j] + 1);
                            fh1_crystals[s][r][p][j]->SetFillColor(kGreen);
                            sprintf(Name,
                                    "Map level, Side %s, Ring %d, Slot %d, Febex %d, ch. %d",
                                    Side,
                                    r + 1,
                                    fFebexInfo[s][r][p][0],
                                    fFebexInfo[s][r][p][1],
                                    fOrderFebexPreamp[j]);
                            fh1_crystals[s][r][p][j]->SetTitle(Name);
                            fh1_crystals[s][r][p][j]->Draw();

                            cMapCryTot[s][r][p]->cd(fOrderFebexPreamp[j] + 1);
                            sprintf(Name,
                                    "Map level, Side %s, Ring %d, Slot %d, Febex %d, ch. %d",
                                    Side,
                                    r + 1,
                                    fFebexInfo[s][r][p][0],
                                    fFebexInfo[s][r][p][1],
                                    fOrderFebexPreamp[j]);
                            fh2_crystalsETot[s][r][p][j]->SetTitle(Name);
                            fh2_crystalsETot[s][r][p][j]->Draw();
                        }
                        if (fFebexInfo[s][r][p][2] != -1)
                        {
                            cMapCryP[s][r][p]->cd(fOrderFebexPreamp[j] + 1);
                            fh1_crystals_p[s][r][p][j]->SetFillColor(kGreen);
                            sprintf(Name,
                                    "Map level (PR), Side %s, Ring %d, Slot %d, Febex %d, ch. %d",
                                    Side,
                                    r + 1,
                                    fFebexInfo[s][r][p][2],
                                    fFebexInfo[s][r][p][3],
                                    fOrderFebexPreamp[j]);
                            fh1_crystals_p[s][r][p][j]->SetTitle(Name);
                            fh1_crystals_p[s][r][p][j]->Draw();

                            cMapCryPTot[s][r][p]->cd(fOrderFebexPreamp[j] + 1);
                            sprintf(Name,
                                    "Map level (PR), Side %s, Ring %d, Slot %d, Febex %d, ch. %d",
                                    Side,
                                    r + 1,
                                    fFebexInfo[s][r][p][2],
                                    fFebexInfo[s][r][p][3],
                                    fOrderFebexPreamp[j]);
                            fh2_crystalsETot_p[s][r][p][j]->SetTitle(Name);
                            fh2_crystalsETot_p[s][r][p][j]->Draw();
                        }

                        if (fCalItemsCalifa)
                        {
                            if (fFebexInfo[s][r][p][0] != -1)
                            {
                                cMapCryCal[s][r][p]->cd(fOrderFebexPreamp[j] + 1);
                                fh1_crystals_cal[s][r][p][j]->SetFillColor(kGreen);
                                sprintf(Name,
                                        "Cal level, Side %s, Ring %d, Slot %d, Febex %d, ch. %d",
                                        Side,
                                        r + 1,
                                        fFebexInfo[s][r][p][0],
                                        fFebexInfo[s][r][p][1],
                                        fOrderFebexPreamp[j]);
                                fh1_crystals_cal[s][r][p][j]->SetTitle(Name);
                                fh1_crystals_cal[s][r][p][j]->Draw();
                            }
                            if (fFebexInfo[s][r][p][2] != -1)
                            {
                                cMapCryPCal[s][r][p]->cd(fOrderFebexPreamp[j] + 1);
                                fh1_crystals_p_cal[s][r][p][j]->SetFillColor(kGreen);
                                sprintf(Name,
                                        "Cal level (PR), Side %s, Ring %d, Slot %d, Febex %d, ch. %d",
                                        Side,
                                        r + 1,
                                        fFebexInfo[s][r][p][2],
                                        fFebexInfo[s][r][p][3],
                                        fOrderFebexPreamp[j]);
                                fh1_crystals_p_cal[s][r][p][j]->SetTitle(Name);
                                fh1_crystals_p_cal[s][r][p][j]->Draw();
                            }
                        }
                    }
                }
            }
        }
        fFebex2Preamp = kFALSE;
    }
    else
    { // Febex to Preamp sequence
        for (Int_t s = 0; s < fNumSides; s++)
        {
            if (s == 1)
                sprintf(Side, "Left");
            else
                sprintf(Side, "Right");
            for (Int_t r = 0; r < fNumRings; r++)
            {
                for (Int_t p = 0; p < fNumPreamps; p++)
                {
                    for (Int_t j = 0; j < fNumCrystalPreamp; j++)
                    {
                        if (fFebexInfo[s][r][p][0] != -1)
                        {
                            cMapCry[s][r][p]->cd(j + 1);
                            fh1_crystals[s][r][p][j]->SetFillColor(45);
                            sprintf(Name, "Map level, Side %s, Ring %d, Preamp %d, ch. %d", Side, r + 1, p + 1, j + 1);
                            fh1_crystals[s][r][p][j]->SetTitle(Name);
                            fh1_crystals[s][r][p][j]->Draw();

                            cMapCryTot[s][r][p]->cd(j + 1);
                            sprintf(Name, "Map level, Side %s, Ring %d, Preamp %d, ch. %d", Side, r + 1, p + 1, j + 1);
                            fh2_crystalsETot[s][r][p][j]->SetTitle(Name);
                            fh2_crystalsETot[s][r][p][j]->Draw();
                        }
                        if (fFebexInfo[s][r][p][2] != -1)
                        {
                            cMapCryP[s][r][p]->cd(j + 1);
                            fh1_crystals_p[s][r][p][j]->SetFillColor(45);
                            sprintf(
                                Name, "Map level (PR), Side %s, Ring %d, Preamp %d, ch. %d", Side, r + 1, p + 1, j + 1);
                            fh1_crystals_p[s][r][p][j]->SetTitle(Name);
                            fh1_crystals_p[s][r][p][j]->Draw();

                            cMapCryPTot[s][r][p]->cd(j + 1);
                            sprintf(
                                Name, "Map level (PR), Side %s, Ring %d, Preamp %d, ch. %d", Side, r + 1, p + 1, j + 1);
                            fh2_crystalsETot_p[s][r][p][j]->SetTitle(Name);
                            fh2_crystalsETot_p[s][r][p][j]->Draw();
                        }

                        if (fCalItemsCalifa)
                        {
                            if (fFebexInfo[s][r][p][0] != -1)
                            {
                                cMapCryCal[s][r][p]->cd(j + 1);
                                fh1_crystals_cal[s][r][p][j]->SetFillColor(45);
                                sprintf(
                                    Name, "Cal level, Side %s, Ring %d, Preamp %d, ch. %d", Side, r + 1, p + 1, j + 1);
                                fh1_crystals_cal[s][r][p][j]->SetTitle(Name);
                                fh1_crystals_cal[s][r][p][j]->Draw();
                            }

                            if (fFebexInfo[s][r][p][2] != -1)
                            {
                                cMapCryPCal[s][r][p]->cd(j + 1);
                                fh1_crystals_p_cal[s][r][p][j]->SetFillColor(45);
                                sprintf(Name,
                                        "Cal level (PR), Side %s, Ring %d, Preamp %d, ch. %d",
                                        Side,
                                        r + 1,
                                        p + 1,
                                        j + 1);
                                fh1_crystals_p_cal[s][r][p][j]->SetTitle(Name);
                                fh1_crystals_p_cal[s][r][p][j]->Draw();
                            }
                        }
                    }
                }
            }
        }
        fFebex2Preamp = kTRUE;
    }
}

void R3BCalifaOnlineSpectra::Exec(Option_t* /*option*/)
{
    if ((fTrigger >= 0) && (header) && (header->GetTrigger() != fTrigger))
        return;
    // fTpat = 1-16; fTpat_bit = 0-15
    Int_t fTpat_bit1 = fTpat1 - 1;
    Int_t fTpat_bit2 = fTpat2 - 1;
    Int_t tpatbin = 0;
    std::vector<int> tpatindex;
    if (header && fTpat1 >= 0 && fTpat2 >= 0)
    {
        for (int i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));

            if (tpatbin != 0)
                tpatindex.push_back(i + 1);

            if (tpatbin != 0 && (i < fTpat_bit1 || i > fTpat_bit2))
            {
                return;
            }
        }
    }
    else if (header)
    {
        for (Int_t i = 0; i < 16; i++)
        {
            tpatbin = (header->GetTpat() & (1 << i));
            if (tpatbin != 0)
                tpatindex.push_back(i + 1);
        }
    }

    int64_t wr[2];
    int64_t wrm = 0.0;
    for (int i = 0; i < 2; i++)
        wr[i] = 0;
    // WR data
    if (fWRItemsCalifa && fWRItemsCalifa->GetEntriesFast() > 0)
    {
        // Califa
        Int_t nHits = fWRItemsCalifa->GetEntriesFast();
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsCalifa->At(ihit));
            if (!hit)
                continue;
            wr[ihit] = hit->GetTimeStamp();
        }
        if (nHits == 2)
            fh1_Califa_wr->Fill(wr[1] - wr[0]);

        // Master Ref. (exp. 2020)
        if (fWRItemsMaster && fWRItemsMaster->GetEntriesFast() > 0)
        {
            nHits = fWRItemsMaster->GetEntriesFast();
            for (Int_t ihit = 0; ihit < nHits; ihit++)
            {
                R3BWRData* hit = dynamic_cast<R3BWRData*>(fWRItemsMaster->At(ihit));
                if (!hit)
                    continue;
                wrm = hit->GetTimeStamp();
            }
        }
    }

    // Mapped trigger data
    if (fTrigMappedItemsCalifa && fTrigMappedItemsCalifa->GetEntriesFast() > 0)
    {
        Int_t nHits = fTrigMappedItemsCalifa->GetEntriesFast();
        Double_t e[2];
        e[0] = 0.;
        e[1] = 0.;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hit = dynamic_cast<R3BCalifaMappedData*>(fTrigMappedItemsCalifa->At(ihit));
            if (!hit)
                continue;
            Int_t ch = hit->GetCrystalId() - 1;
            e[ch] = hit->GetEnergy();
        }
        if (e[0] > 0 && e[1] > 0)
            fh2_Califa_EtrigCor[0]->Fill(e[0], e[1]);
        if (e[0] > 0 && e[1] > 0 && header && header->GetTrigger() == 1)
            fh2_Califa_EtrigCor[1]->Fill(e[0], e[1]);

        if (e[0] > 0)
            fh1_Califa_Etrig[0]->Fill(e[0]);

        if (e[1] > 0)
            fh1_Califa_Etrig[1]->Fill(e[1]);
    }

    // Mapped data
    if (fMappedItemsCalifa && fMappedItemsCalifa->GetEntriesFast() > 0)
    {
        Int_t nHits = fMappedItemsCalifa->GetEntriesFast();
        Int_t Crymult = 0;
        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hit = dynamic_cast<R3BCalifaMappedData*>(fMappedItemsCalifa->At(ihit));
            if (!hit)
                continue;

            Int_t cryId = hit->GetCrystalId();

            // compensate slave exploder delays:
            int64_t wrc = hit->GetWrts() + 245 * (fMap_Par->GetPreamp(cryId) > 8);
            if (wrm > 0.)
            {
                bool side = !(fMap_Par->GetHalf(cryId) % 2);
                fh1_wrs[side]->Fill(wrc - wrm);
            }

            if ((fMap_Par->GetInUse(cryId) == 1 && cryId <= fNbCalifaCrystals / 2) ||
                (cryId > fNbCalifaCrystals / 2 &&
                 fMap_Par->GetInUse(cryId) != fMap_Par->GetInUse(cryId - fNbCalifaCrystals / 2)))
                Crymult++;

            fh2_Califa_cryId_energy->Fill(cryId, hit->GetEnergy());

            if (fMap_Par->GetHalf(cryId) == 2)
                fh2_Preamp_vs_ch_L[fMap_Par->GetRing(cryId) - 1]->Fill(fMap_Par->GetPreamp(cryId),
                                                                       fMap_Par->GetChannel(cryId));
            if (fMap_Par->GetHalf(cryId) == 1)
                fh2_Preamp_vs_ch_R[fMap_Par->GetRing(cryId) - 1]->Fill(fMap_Par->GetPreamp(cryId),
                                                                       fMap_Par->GetChannel(cryId));

            if (fMap_Par->GetInUse(cryId) == 1 && cryId <= fNbCalifaCrystals / 2)
            {
                fh1_crystals[fMap_Par->GetHalf(cryId) - 1][fMap_Par->GetRing(cryId) - 1][fMap_Par->GetPreamp(cryId) - 1]
                            [fMap_Par->GetChannel(cryId) - 1]
                                ->Fill(hit->GetEnergy());
                fh2_crystalsETot[fMap_Par->GetHalf(cryId) - 1][fMap_Par->GetRing(cryId) - 1]
                                [fMap_Par->GetPreamp(cryId) - 1][fMap_Par->GetChannel(cryId) - 1]
                                    ->Fill(hit->GetEnergy(), hit->GetTot());
            }

            else if (fMap_Par->GetInUse(cryId) == 1 && cryId > fNbCalifaCrystals / 2)
            {
                fh1_crystals_p[fMap_Par->GetHalf(cryId) - 1][fMap_Par->GetRing(cryId) - 1]
                              [fMap_Par->GetPreamp(cryId) - 1][fMap_Par->GetChannel(cryId) - 1]
                                  ->Fill(hit->GetEnergy());
                fh2_crystalsETot_p[fMap_Par->GetHalf(cryId) - 1][fMap_Par->GetRing(cryId) - 1]
                                  [fMap_Par->GetPreamp(cryId) - 1][fMap_Par->GetChannel(cryId) - 1]
                                      ->Fill(hit->GetEnergy(), hit->GetTot());
            }
        }
        fh1_Califa_Mult->Fill(Crymult);
    }

    // Cal data
    bool cond_califa_veto = false;
    if (fCalItemsCalifa && fCalItemsCalifa->GetEntriesFast() > 0)
    {
        Int_t nHits = fCalItemsCalifa->GetEntriesFast();

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hit = dynamic_cast<R3BCalifaCrystalCalData*>(fCalItemsCalifa->At(ihit));
            if (!hit)
                continue;

            Int_t cryId = hit->GetCrystalId();

            fh2_Califa_cryId_energy_cal->Fill(cryId, hit->GetEnergy());

            if (cryId > 3344 && cryId < 3473 && hit->GetEnergy() > 10000.)
                cond_califa_veto = true;

            auto nf_ns = (hit->GetNs() + hit->GetNf());
            if (nf_ns > 0. && hit->GetNs() > 0 && hit->GetNf() > 0)
            {
                if (cryId > fNbCalifaCrystals / 2)
                {
                    if (cryId <= BarrelCrystals)
                    {
                        fh2_Califa_NsNf[0]->Fill(nf_ns / 1000., hit->GetNf() / nf_ns);
                    }
                    else if (cryId <= iPhosCrystals)
                    {
                        fh2_Califa_NsNf[1]->Fill(nf_ns / 1000., hit->GetNf() / nf_ns);
                    }
                    else
                    {
                        fh2_Califa_NsNf[2]->Fill(nf_ns / 1000., hit->GetNf() / nf_ns);
                    }
                }
            }

            if (fMap_Par->GetInUse(cryId) == 1 && cryId <= fNbCalifaCrystals / 2)
                fh1_crystals_cal[fMap_Par->GetHalf(cryId) - 1][fMap_Par->GetRing(cryId) - 1]
                                [fMap_Par->GetPreamp(cryId) - 1][fMap_Par->GetChannel(cryId) - 1]
                                    ->Fill(hit->GetEnergy());
            else if (fMap_Par->GetInUse(cryId) == 1 && cryId > fNbCalifaCrystals / 2)
                fh1_crystals_p_cal[fMap_Par->GetHalf(cryId) - 1][fMap_Par->GetRing(cryId) - 1]
                                  [fMap_Par->GetPreamp(cryId) - 1][fMap_Par->GetChannel(cryId) - 1]
                                      ->Fill(hit->GetEnergy());
        }
    }

    // Hit data
    if (fHitItemsCalifa && fHitItemsCalifa->GetEntriesFast() > 0)
    {
        Int_t nHits = fHitItemsCalifa->GetEntriesFast();
        fh1_Califa_MultHit->Fill(nHits);

        double theta = 0., phi = 0.;
        std::vector<double> califa_theta;
        std::vector<double> califa_phi;
        std::vector<double> califa_e;

        for (Int_t ihit = 0; ihit < nHits; ihit++)
        {
            auto hit = dynamic_cast<R3BCalifaClusterData*>(fHitItemsCalifa->At(ihit));
            if (!hit)
                continue;
            theta = hit->GetTheta() * TMath::RadToDeg();
            phi = hit->GetPhi() * TMath::RadToDeg();
            if (hit->GetClusterType() == 0) // proton
            {
                califa_theta.push_back(theta);
                califa_phi.push_back(phi);
                califa_e.push_back(hit->GetEnergy());
            }
            fh2_Califa_theta_phi->Fill(theta, phi);
            if (hit->GetClusterType() == 0) // proton
                fh2_Califa_theta_energy->Fill(theta, hit->GetEnergy());
            fh1_Califa_total_energy->Fill(hit->GetEnergy());
        }

        TVector3 master[2];
        Double_t maxEL = 0., maxER = 0.;
        for (Int_t i1 = 0; i1 < califa_theta.size(); i1++)
        {

            if (califa_e[i1] > maxER && TMath::Abs(califa_phi[i1]) > 90.) // wixhausen
            {
                master[0].SetMagThetaPhi(1., califa_theta[i1] * TMath::DegToRad(), califa_phi[i1] * TMath::DegToRad());
                maxER = califa_e[i1];
            }
            if (califa_e[i1] > maxEL && TMath::Abs(califa_phi[i1]) < 90.)
            { // messel
                master[1].SetMagThetaPhi(1., califa_theta[i1] * TMath::DegToRad(), califa_phi[i1] * TMath::DegToRad());
                maxEL = califa_e[i1];
            }
        }
        if (maxEL > fMinProtonE && maxER > fMinProtonE)
        {
            auto opa = master[0].Angle(master[1]) * TMath::RadToDeg();
            fh1_openangle->Fill(opa);
            if (opa < 90. && opa > 68.)
                fh2_Califa_coinE_p2p->Fill(maxEL / 1000., maxER / 1000.);
            for (const auto& itpat : tpatindex)
                fh2_openangle_tpat->Fill(itpat, master[0].Angle(master[1]) * TMath::RadToDeg());
        }

        // 0->veto 1->foot 2->inclusive proton
        // if (master[0].Theta() * TMath::RadToDeg() > 90 || master[1].Theta() * TMath::RadToDeg() > 90)
        if (cond_califa_veto)
        {
            fh1_CalifaTriggers->Fill(0);
        }
        else
        {
            if ((abs(master[0].Phi() * TMath::RadToDeg()) < 45 && master[0].Theta() * TMath::RadToDeg() < 90 &&
                 master[0].Theta() * TMath::RadToDeg() > 45) ||
                (abs(master[1].Phi() * TMath::RadToDeg()) < 45 && master[1].Theta() * TMath::RadToDeg() < 90 &&
                 master[1].Theta() * TMath::RadToDeg() > 45))
                fh1_CalifaTriggers->Fill(1);
            if (master[0].Theta() * TMath::RadToDeg() < 90 && master[1].Theta() * TMath::RadToDeg() < 90)
                fh1_CalifaTriggers->Fill(2);
        }

        // Comparison of hits to get energy, theta and phi correlations between them
        for (Int_t i1 = 0; i1 < califa_theta.size(); i1++)
        {
            for (Int_t i2 = i1 + 1; i2 < califa_theta.size(); i2++)
            {
                if (gRandom->Uniform(0., 1.) < 0.5)
                {
                    fh2_Califa_coinE->Fill(califa_e[i1] / 1000., califa_e[i2] / 1000.);
                    fh2_Califa_coinTheta->Fill(califa_theta[i1], califa_theta[i2]);
                    fh2_Califa_coinPhi->Fill(califa_phi[i1], califa_phi[i2]);

                    if (master[0].Angle(master[1]) * TMath::RadToDeg() > 68 &&
                        master[0].Angle(master[1]) * TMath::RadToDeg() < 90)
                    {
                        fh2_Califa_coinTheta_cutOPA->Fill(califa_theta[i1], califa_theta[i2]);
                    }
                }
                else
                {
                    fh2_Califa_coinE->Fill(califa_e[i2] / 1000., califa_e[i1] / 1000.);
                    fh2_Califa_coinTheta->Fill(califa_theta[i2], califa_theta[i1]);
                    fh2_Califa_coinPhi->Fill(califa_phi[i2], califa_phi[i1]);

                    if (master[0].Angle(master[1]) * TMath::RadToDeg() > 68 &&
                        master[0].Angle(master[1]) * TMath::RadToDeg() < 90)
                    {
                        fh2_Califa_coinTheta_cutOPA->Fill(califa_theta[i2], califa_theta[i1]);
                    }
                }
            }
        }
    }

    fNEvents += 1;
}

void R3BCalifaOnlineSpectra::FinishEvent()
{
    if (fMappedItemsCalifa)
    {
        fMappedItemsCalifa->Clear();
    }
    if (fTrigMappedItemsCalifa)
    {
        fTrigMappedItemsCalifa->Clear();
    }
    if (fCalItemsCalifa)
    {
        fCalItemsCalifa->Clear();
    }
    if (fHitItemsCalifa)
    {
        fHitItemsCalifa->Clear();
    }
    if (fWRItemsCalifa)
    {
        fWRItemsCalifa->Clear();
    }
    if (fWRItemsMaster)
    {
        fWRItemsMaster->Clear();
    }
}

void R3BCalifaOnlineSpectra::FinishTask()
{
    // Write canvas for Califa WR data
    if (fWRItemsCalifa)
    {
        cCalifa_wr->Write();
    }

    // Write canvas for Master-Califa WR data
    if (fWRItemsMaster && fWRItemsCalifa)
    {
        cWrs->Write();
        if (fHitItemsCalifa)
            cCalifa_wr_energy->Write();
    }

    if (fTrigMappedItemsCalifa)
    {
        fh2_Califa_EtrigCor[0]->Write();
        fh2_Califa_EtrigCor[1]->Write();
        fh1_Califa_Etrig[0]->Write();
        fh1_Califa_Etrig[1]->Write();
    }

    // Write canvas for Mapped data
    if (fMappedItemsCalifa)
    {
        cCalifaMult->Write();
        cCalifa_cry_energy->Write();

        for (Int_t i = 0; i < fNumRings; i++)
        {
            cMap_RingR[i]->Write();
            cMap_RingL[i]->Write();
        }
        for (Int_t s = 0; s < fNumSides; s++)
            for (Int_t r = 0; r < fNumRings; r++)
                for (Int_t p = 0; p < fNumPreamps; p++)
                {
                    if (fFebexInfo[s][r][p][0] != -1)
                        cMapCry[s][r][p]->Write();
                    if (fFebexInfo[s][r][p][2] != -1)
                        cMapCryP[s][r][p]->Write();
                }
    }

    // Write canvas for Cal data
    if (fCalItemsCalifa)
    {
        cCalifa_cry_energy_cal->Write();
        for (const auto& hist : fh2_Califa_NsNf)
        {
            hist->Write();
        }
        for (Int_t s = 0; s < fNumSides; s++)
            for (Int_t r = 0; r < fNumRings; r++)
                for (Int_t p = 0; p < fNumPreamps; p++)
                {
                    if (fFebexInfo[s][r][p][0] != -1)
                        cMapCryCal[s][r][p]->Write();
                    if (fFebexInfo[s][r][p][2] != -1)
                        cMapCryPCal[s][r][p]->Write();
                }
    }

    // Write canvas for Hit data
    if (fHitItemsCalifa)
    {
        cCalifaCoinE->Write();
        cCalifaCoinTheta->Write();
        cCalifaCoinPhi->Write();
        cCalifa_angles->Write();
        cCalifa_theta_energy->Write();
        cCalifa_hitenergy->Write();
        cCalifa_opening->Write();
        cCalifa_opening_tpat->Write();
    }
}

ClassImp(R3BCalifaOnlineSpectra)
