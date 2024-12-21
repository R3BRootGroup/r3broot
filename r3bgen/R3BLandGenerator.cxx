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

// -------------------------------------------------------------------------
// -----                R3BLandGenerator source file                 -------
// -----              Created by M.I.Cherciu -18.06.2010             -------
// -------------------------------------------------------------------------
#include "R3BLandGenerator.h"

#include "FairIon.h"
#include "FairLogger.h"
#include "FairPrimaryGenerator.h"
#include "FairRunSim.h"

// ROOT INCLUDES
#include "TApplication.h"
#include "TBranch.h"
#include "TDatabasePDG.h"
#include "TFile.h"
#include "TLeaf.h"
#include "TParticlePDG.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"

#include <iostream>

using std::cout;
using std::endl;
using std::map;

// -----   Default constructor   ------------------------------------------
R3BLandGenerator::R3BLandGenerator()
    : fInputFile(NULL)
    , tree(NULL)
    , lPDG(NULL)
    , lM(NULL)
    , lZ(NULL)
    , lE(NULL)
    , lp(NULL)
    , ltheta(NULL)
    , lphi(NULL)
    , lpx(NULL)
    , lpy(NULL)
    , lpz(NULL)
    , lx(NULL)
    , ly(NULL)
    , lz(NULL)
    , fFileName("")
    , fPDG(NULL)
    , PDGo(0.)
    , Mo(0.)
    , Zo(0.)
    , Eo(0.)
    , po(0.)
    , thetao(0.)
    , phio(0.)
    , xo(0.)
    , yo(0.)
    , zo(0.)
    , pxo(0.)
    , pyo(0.)
    , pzo(0.)
    , vx(0.)
    , vy(0.)
    , vz(0.)
{
}
// ------------------------------------------------------------------------

// -----   Standard constructor   -----------------------------------------
R3BLandGenerator::R3BLandGenerator(const char* fileName)
    : fInputFile(NULL)
    , tree(NULL)
    , lPDG(NULL)
    , lM(NULL)
    , lZ(NULL)
    , lE(NULL)
    , lp(NULL)
    , ltheta(NULL)
    , lphi(NULL)
    , lpx(NULL)
    , lpy(NULL)
    , lpz(NULL)
    , lx(NULL)
    , ly(NULL)
    , lz(NULL)
    , fFileName(fileName)
    , fPDG(TDatabasePDG::Instance())
    , PDGo(0.)
    , Mo(0.)
    , Zo(0.)
    , Eo(0.)
    , po(0.)
    , thetao(0.)
    , phio(0.)
    , xo(0.)
    , yo(0.)
    , zo(0.)
    , pxo(0.)
    , pyo(0.)
    , pzo(0.)
    , vx(0.)
    , vy(0.)
    , vz(0.)
{
    cout << "-I- R3BLandGenerator: Opening input file " << fileName << endl;
    // Open first the file to read all particles
    fInputFile = new TFile(fFileName, "READ");

    if (!fInputFile->IsOpen())
        LOG(fatal) << "R3BLandGenerator: Cannot open input file.";

    tree = dynamic_cast<TTree*>(fInputFile->Get("Tree"));
    lPDG = (TLeaf*)tree->GetLeaf("PDG");
    lM = (TLeaf*)tree->GetLeaf("M");
    lZ = (TLeaf*)tree->GetLeaf("Z");
    lE = (TLeaf*)tree->GetLeaf("E");
    lp = (TLeaf*)tree->GetLeaf("p");
    ltheta = (TLeaf*)tree->GetLeaf("theta");
    lphi = (TLeaf*)tree->GetLeaf("phi");
    lx = (TLeaf*)tree->GetLeaf("x");
    ly = (TLeaf*)tree->GetLeaf("y");
    lz = (TLeaf*)tree->GetLeaf("z");
    lpx = (TLeaf*)tree->GetLeaf("px");
    lpy = (TLeaf*)tree->GetLeaf("py");
    lpz = (TLeaf*)tree->GetLeaf("pz");
}
// ------------------------------------------------------------------------

R3BLandGenerator::R3BLandGenerator(const R3BLandGenerator& right)
    : fInputFile(right.fInputFile)
    , tree(right.tree)
    , lPDG(right.lPDG)
    , lM(right.lM)
    , lZ(right.lZ)
    , lE(right.lE)
    , lp(right.lp)
    , ltheta(right.ltheta)
    , lphi(right.lphi)
    , lpx(right.lpx)
    , lpy(right.lpy)
    , lpz(right.lpz)
    , lx(right.lx)
    , ly(right.ly)
    , lz(right.lz)
    , fFileName(right.fFileName)
    , fPDG(right.fPDG)
    , PDGo(right.PDGo)
    , Mo(right.Mo)
    , Zo(right.Zo)
    , Eo(right.Eo)
    , po(right.po)
    , thetao(right.thetao)
    , phio(right.phio)
    , xo(right.xo)
    , yo(right.yo)
    , zo(right.zo)
    , pxo(right.pxo)
    , pyo(right.pyo)
    , pzo(right.pzo)
    , vx(right.vx)
    , vy(right.vy)
    , vz(right.vz)
{
}

// -----   Destructor   ---------------------------------------------------
R3BLandGenerator::~R3BLandGenerator() { CloseInput(); }
// ------------------------------------------------------------------------

// -----   Public method ReadEvent   --------------------------------------
Bool_t R3BLandGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{

    // Check for input file
    if (!fInputFile->IsOpen())
    {
        cout << "-E- R3BLandGenerator: Input file not open!" << endl;
        return kFALSE;
    }

    // Define track variables to be read from file
    Int_t nTracks = 0;
    PDGo = 0.;
    Mo = 0.;
    Zo = 0.;
    Eo = 0.;
    po = 0.;
    thetao = 0.;
    phio = 0.;
    xo = 0.;
    yo = 0.;
    zo = 0.;
    pxo = 0.;
    pyo = 0.;
    pzo = 0.;

    // How many tracks?! from the ROOT File
    nTracks = tree->GetEntries();

    // Loop over tracks in the current event
    for (Int_t itrack = 0; itrack < nTracks; itrack++)
    {

        lPDG->GetBranch()->GetEntry(itrack);
        lM->GetBranch()->GetEntry(itrack);
        lZ->GetBranch()->GetEntry(itrack);
        lE->GetBranch()->GetEntry(itrack);
        lp->GetBranch()->GetEntry(itrack);
        ltheta->GetBranch()->GetEntry(itrack);
        lphi->GetBranch()->GetEntry(itrack);
        lx->GetBranch()->GetEntry(itrack);
        ly->GetBranch()->GetEntry(itrack);
        lz->GetBranch()->GetEntry(itrack);
        lpx->GetBranch()->GetEntry(itrack);
        lpy->GetBranch()->GetEntry(itrack);
        lpz->GetBranch()->GetEntry(itrack);

        PDGo = lPDG->GetValue();
        Mo = lM->GetValue();
        Zo = lZ->GetValue();
        Eo = lE->GetValue();
        po = lp->GetValue();
        thetao = ltheta->GetValue();
        phio = lphi->GetValue();
        xo = lx->GetValue();
        yo = ly->GetValue();
        zo = lz->GetValue();
        pxo = lpx->GetValue();
        pyo = lpy->GetValue();
        pzo = lpz->GetValue();
        // Define de vertex points for particles
        vx = 0.;
        vy = 0.;
        vz = 0.;
        cout << "-I- R3BLandGenerator: iPid: " << PDGo << " Z = " << Zo << " px = " << pxo << " py = " << pyo
             << " pz = " << pzo << " vx = " << vx << " vy = " << vy << " vz = " << vz << endl;

        primGen->AddTrack(PDGo, pxo, pyo, pzo, vx, vy, vz);
    }
    return kTRUE;
}
// ------------------------------------------------------------------------
// -----   Private method CloseInput   ------------------------------------
void R3BLandGenerator::CloseInput()
{
    if (fInputFile)
    {
        if (fInputFile->IsOpen())
        {
            cout << "-I- R3BLandGenerator: Closing input file " << fFileName << endl;
            fInputFile->Close();
        }
        delete fInputFile;
        fInputFile = NULL;
    }
}
// ------------------------------------------------------------------------

ClassImp(R3BLandGenerator)
