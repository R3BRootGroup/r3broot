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

#include "R3BNeulandVisualizer.h"
#include "TF1.h"
#include "TList.h"
#include "TPaletteAxis.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>

/* This function is required to suppress boxes for empty bins - make them transparent.*/
static Double_t gEmptyBinSupressor(const Double_t* x, const Double_t*)
{
    if (x)
    {
        return *x > 0 ? 1. : 0.;
    }
    return 0.;
}

R3BNeulandVisualizer::R3BNeulandVisualizer(const TString& inputFileName, const TString& what)
    : fFile(std::make_shared<TFile>(inputFileName, "read"))
    , fTree(dynamic_cast<TTree*>(fFile->Get("evt")))
    , fh3(nullptr)
    , fCanvas(nullptr)
    , fIndex(0)
{
    fTree->SetBranchAddress(what, &fh3);

    gStyle->SetCanvasPreferGL(kTRUE);
    gStyle->SetPalette(kViridis);
    gStyle->SetOptStat(0);

    // Needs to be after setting gStyle
    fCanvas = std::make_shared<TCanvas>("canvas", "NeuLAND 3D Event Visualization", 0, 0, 800, 800);

    Visualize();
}

void R3BNeulandVisualizer::Visualize()
{
    if (fIndex >= fTree->GetEntries())
    {
        std::cout << "Index larger than number of Events" << std::endl;
        return;
    }

    fTree->GetEntry(fIndex);

    TList* const lof = fh3->GetListOfFunctions();
    lof->Add(new TF1("TransferFunction", gEmptyBinSupressor, 0., 1000., 0));

    fh3->Draw("glcolz");

    fCanvas->Flush();
}

ClassImp(R3BNeulandVisualizer)
