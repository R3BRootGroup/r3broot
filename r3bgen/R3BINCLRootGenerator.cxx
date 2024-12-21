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

#include "R3BINCLRootGenerator.h"
#include "FairIon.h"
#include "FairLogger.h"
#include "FairPrimaryGenerator.h"
#include "FairRunSim.h"
#include "G4NistManager.hh"
#include "R3BLogger.h"
#include "TMath.h"
#include "TRandom.h"

R3BINCLRootGenerator::R3BINCLRootGenerator()
    : fFileName()
    , fX(0.)
    , fY(0.)
    , fZ(0.)
    , fPointVtxIsSet(kFALSE)
    , fDX(0.)
    , fDY(0.)
    , fDZ(0.)
    , fBoxVtxIsSet(kFALSE)
    , fOnlySpallation(kFALSE)
    , fOnlyP2pFission(kFALSE)
    , fOnlyFission(kFALSE)
    , fEvtRoot(0)
    , fEvt(0)
{
}

R3BINCLRootGenerator::R3BINCLRootGenerator(std::string fileName)
    : fFileName(std::move(fileName))
    , fX(0.)
    , fY(0.)
    , fZ(0.)
    , fPointVtxIsSet(kFALSE)
    , fDX(0.)
    , fDY(0.)
    , fDZ(0.)
    , fBoxVtxIsSet(kFALSE)
    , fOnlySpallation(kFALSE)
    , fOnlyP2pFission(kFALSE)
    , fOnlyFission(kFALSE)
    , fEvtRoot(0)
    , fEvt(0)
{
    RegisterIons();
}

R3BINCLRootGenerator::R3BINCLRootGenerator(const TString& fileName)
    : R3BINCLRootGenerator(std::string(fileName))
{
}

R3BINCLRootGenerator::R3BINCLRootGenerator(const char* fileName)
    : R3BINCLRootGenerator(std::string(fileName))
{
}

R3BINCLRootGenerator::~R3BINCLRootGenerator() { fInput->Close(); }

bool R3BINCLRootGenerator::ReadEvent(FairPrimaryGenerator* primGen)
{
    // Track variables to be read from file
    int iZ = 0;
    int iA = 0;
    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double vx = 0.;
    double vy = 0.;
    double vz = 0.;

    if (fPointVtxIsSet)
    {
        if (fBoxVtxIsSet)
        {
            vx = gRandom->Gaus(fX, fDX);
            vy = gRandom->Gaus(fY, fDY);
            vz = gRandom->Uniform(fZ - fDZ, fZ + fDZ);
        }
        else
        {
            vx = fX;
            vy = fY;
            vz = fZ;
        }
    }

    bool validevent = false;
    while (!validevent)
    {
        Tree->GetEntry(fEvt);

        if (fOnlyFission)
        {
            for (Int_t j = 0; j < fParticles; j++)
                if (fOrigin[j] == 1 || fOrigin[j] == 101)
                    validevent = true;
            if (!validevent)
            {
                fEvt++;
            }
        }
        else if (fOnlySpallation)
        {
            for (Int_t j = 0; j < fParticles; j++)
                if (fOrigin[j] == 0)
                    validevent = true;
            if (!validevent)
            {
                fEvt++;
            }
        }
        else if (fOnlyP2pFission)
        {
            Int_t nbp = 0;
            for (Int_t j = 0; j < fParticles; j++)
                if (fOrigin[j] == 1 || fOrigin[j] == 101)
                    validevent = true;
            if (!validevent)
            {
                fEvt++;
            }

            for (Int_t j = 0; j < fParticles; j++)
                if (fPdgCode[j] == 2212)
                    nbp++;

            if (nbp != 2 && validevent)
            {
                fEvt++;
                validevent = false;
            }
        }
        else
        {
            validevent = true;
        }

        if (fEvt > fEvtRoot)
        {
            LOG(error)
                << "\033[5m\033[31m R3BINCLRootGenerator: Number of simulated events larger than the ones contained "
                   "in the Root file \033[0m ";
            LOG(warn)
                << "\033[5m\033[33m R3BINCLRootGenerator: Please, provide a new Root file with more events \033[0m ";
            return kTRUE;
        }

        R3BLOG(debug, "fParticles: " << fParticles);

        if (validevent)
        {
            for (Int_t j = 0; j < fParticles; j++)
            {
                Int_t pdg = 0;
                if (fMass[j] > 1 && fCharge[j] > 0)
                {
                    iA = fMass[j];
                    iZ = fCharge[j];
                    pdg = GetIonPdgId(iZ, iA);
                    pz = fPzPrime[j] / 1000.;
                    Double_t pt = pz * TMath::Tan(fThetaPrime[j] * TMath::DegToRad());
                    px = pt * TMath::Cos(fPhi[j] * TMath::DegToRad());
                    py = pt * TMath::Sin(fPhi[j] * TMath::DegToRad());
                }
                else
                {
                    pdg = fPdgCode[j];
                    pz = fPzPrime[j] / 1000.;
                    Double_t pt = pz * TMath::Tan(fThetaPrime[j] * TMath::DegToRad());
                    px = pt * TMath::Cos(fPhi[j] * TMath::DegToRad());
                    py = pt * TMath::Sin(fPhi[j] * TMath::DegToRad());
                }
                R3BLOG(debug, "PDG:Px:Py:Pz " << pdg << " " << px << " " << py << " " << pz);
                primGen->AddTrack(pdg, px, py, pz, vx, vy, vz);
            }
        }
    }

    fEvt++;
    return kTRUE;
}

void R3BINCLRootGenerator::RegisterIons()
{
    R3BLOG(info, "Looking for ions ...");

    // Track variables to be read from file
    Int_t iZ = 0;
    Int_t iA = 0;
    // Keep a list of ions to register
    std::map<Int_t, FairIon*> ions;

    fInput = new TFile(fFileName);
    Tree = dynamic_cast<TTree*>(fInput->Get("et"));

    fEvtRoot = Tree->GetEntries();
    R3BLOG(info, "Root file entries " << fEvtRoot);

    // first root file
    Tree->SetBranchAddress("nParticles", &fParticles);
    Tree->SetBranchAddress("origin", fOrigin);
    Tree->SetBranchAddress("A", fMass);
    Tree->SetBranchAddress("Z", fCharge);
    Tree->SetBranchAddress("pzPrime", fPzPrime);
    Tree->SetBranchAddress("thetaPrime", fThetaPrime);
    Tree->SetBranchAddress("phi", fPhi);
    Tree->SetBranchAddress("PDGCode", fPdgCode);

    for (Int_t ie = 0; ie < fEvtRoot; ie++)
    {
        Tree->GetEntry(ie);

        for (Int_t j = 0; j < fParticles; j++)
            if (fMass[j] > 1 && fCharge[j] > 0)
            {
                iA = fMass[j];
                iZ = fCharge[j];
                const Int_t pdg = GetIonPdgId(iZ, iA);
                if (ions.find(pdg) == ions.end())
                {
                    const Double_t mass = G4NistManager::Instance()->GetIsotopeMass(iZ, iA) / CLHEP::GeV;
                    LOG(debug) << "R3BINCLRootGenerator: New ion " << iZ << "\t" << iA << "\t" << mass;
                    ions[pdg] = new FairIon(TString::Format("Ion_%d_%d", iA, iZ), iZ, iA, iZ, 0., mass);
                }
            }
    }

    for (const auto& kv : ions)
    {
        // Note: FairRoot will not register ions known to TDatabasePDG (e.g. alphas)
        FairRunSim::Instance()->AddNewIon(kv.second);
    }
    R3BLOG(info, ions.size() << " ions registered.");
}

void R3BINCLRootGenerator::SetXYZ(Double32_t x, Double32_t y, Double32_t z)
{
    fX = x;
    fY = y;
    fZ = z;
    fPointVtxIsSet = kTRUE;
}

void R3BINCLRootGenerator::SetDxDyDz(Double32_t sx, Double32_t sy, Double32_t sz)
{
    fDX = sx;
    fDY = sy;
    fDZ = sz;
    fBoxVtxIsSet = kTRUE;
}

ClassImp(R3BINCLRootGenerator);
