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

#include "R3BMCTracks.h"
#include "R3BEventManager.h"

#include "FairEventManager.h" // for FairEventManager
#include "FairLogger.h"
#include "FairRootManager.h" // for FairRootManager

#include <TClonesArray.h>        // for TClonesArray
#include <TEveManager.h>         // for TEveManager, gEve
#include <TEvePathMark.h>        // for TEvePathMark
#include <TEveTrack.h>           // for TEveTrackList, TEveTrack
#include <TEveTrackPropagator.h> // for TEveTrackPropagator
#include <TEveVector.h>          // for TEveVector, TEveVectorT
#include <TGeoTrack.h>           // for TGeoTrack
#include <TMathBase.h>           // for Max, Min
#include <TObjArray.h>           // for TObjArray
#include <TParticle.h>           // for TParticle
#include <cstring>               // for strcmp

#include <iostream>
using std::cout;
using std::endl;

// -----   Default constructor   -------------------------------------------
R3BMCTracks::R3BMCTracks()
    : FairTask("R3BMCTracks", 0)
    , fTrackList(nullptr)
    , fTrPr(nullptr)
    , fEventManager(nullptr)
    , fEveTrList(nullptr)
    , fEvent("")
    , fTrList(nullptr)
    , MinEnergyLimit(-1.)
    , MaxEnergyLimit(-1.)
    , PEnergy(-1.)
{
}

R3BMCTracks::R3BMCTracks(const char* name, Int_t iVerbose)
    : FairTask(name, iVerbose)
    , fTrackList(nullptr)
    , fTrPr(nullptr)
    , fEventManager(nullptr)
    , fEveTrList(new TObjArray(16))
    , fEvent("")
    , fTrList(nullptr)
    , MinEnergyLimit(-1.)
    , MaxEnergyLimit(-1.)
    , PEnergy(-1.)
{
}

InitStatus R3BMCTracks::Init()
{
    if (fVerbose > 1)
        cout << "R3BMCTracks::Init()" << endl;

    FairRootManager* fManager = FairRootManager::Instance();
    fTrackList = dynamic_cast<TClonesArray*>(fManager->GetObject("GeoTracks"));
    if (fTrackList == 0)
    {
        cout << "FairMCPointDraw::Init()  branch " << GetName() << " Not found! Task will be deactivated " << endl;
        SetActive(kFALSE);
    }
    if (fVerbose > 2)
        cout << "FairMCTracks::Init() get track list" << fTrackList << endl;
    if (fVerbose > 2)
        cout << "FairMCTracks::Init()  create propagator" << endl;

    fEventManager = R3BEventManager::Instance();
    if (fVerbose > 2)
        cout << "FairMCTracks::Init() get instance of FairEventManager " << endl;
    fEvent = "Current Event";
    MinEnergyLimit = fEventManager->GetEvtMinEnergy();
    MaxEnergyLimit = fEventManager->GetEvtMaxEnergy();
    PEnergy = 0;
    if (IsActive())
        return kSUCCESS;
    else
        return kERROR;
}

void R3BMCTracks::Exec(Option_t* option)
{
    if (IsActive())
    {

        if (fVerbose > 1)
            cout << " FairMCTracks::Exec " << endl;
        TGeoTrack* tr;
        const Double_t* point;

        Reset();

        for (Int_t i = 0; i < fTrackList->GetEntriesFast(); i++)
        {
            if (fVerbose > 2)
                cout << "FairMCTracks::Exec " << i << endl;
            tr = dynamic_cast<TGeoTrack*>(fTrackList->At(i));
            TParticle* P = dynamic_cast<TParticle*>(tr->GetParticle());

            PEnergy = (P->Energy() - P->GetCalcMass()) * 1000; //[MeV]
            MinEnergyLimit = TMath::Min(PEnergy - 10, MinEnergyLimit);
            MinEnergyLimit = TMath::Max(0.0, MinEnergyLimit);
            MaxEnergyLimit = TMath::Max(PEnergy + 10, MaxEnergyLimit);

            fEventManager->SetMaxEnergy(MaxEnergyLimit + 1);

            if (fVerbose > 2)
                cout << "MinEnergyLimit " << MinEnergyLimit << " MaxEnergyLimit " << MaxEnergyLimit << endl;
            if (fEventManager->IsPriOnly() && P->GetMother(0) > -1)
                continue;
            if (fEventManager->GetCurrentPDG() != 0 && fEventManager->GetCurrentPDG() != tr->GetPDG())
                continue;

            if (fVerbose > 2)
                cout << "PEnergy " << PEnergy << " Min " << fEventManager->GetMinEnergy() << " Max "
                     << fEventManager->GetMaxEnergy() << endl;
            if ((PEnergy < fEventManager->GetMinEnergy()) || (PEnergy > fEventManager->GetMaxEnergy()))
                continue;

            (dynamic_cast<R3BEventManager*>(fEventManager))->AddParticlesToPdgDataBase(tr->GetPDG());
            if (fVerbose > 3)
                cout << "Particle with PDG " << tr->GetPDG() << " added to DataBase " << endl;
            if (fVerbose > 3)
                cout << "Particle  " << P << " and propagator " << fTrPr << endl;

            Int_t Np = tr->GetNpoints();
            fTrList = GetTrGroup(P);
            if (fVerbose > 3)
                cout << "Track list: " << fTrList << " - " << fTrList->GetLimP() << " - " << fTrList->GetMaxP() << endl;
            TEveTrack* track = new TEveTrack(P, tr->GetPDG(), fTrPr);
            if (fVerbose > 3)
                cout << "Track: " << track << " - " << track->GetPdg() << " - " << track->GetLabel() << endl;
            if (tr->GetPDG() > 5000000)
            {
                track->SetLineColor(12);
            }
            else
            {
                track->SetLineColor(fEventManager->Color(tr->GetPDG()));
            }
            track->SetLineStyle(9);

            // Set Title / Tooltip
            char title[100];
            sprintf(title,
                    "pdg: %i, name: %s\nTrackID: %i, MotherID: %i\nE: %f MeV\nT: %f ns",
                    tr->GetPDG(),
                    P->GetTitle(),
                    i,
                    P->GetMother(0),
                    PEnergy,
                    P->T());
            track->SetTitle(title);

            // Set the line width depending on energy
            if ((dynamic_cast<R3BEventManager*>(fEventManager))->IsScaleByEnergy())
            {
                Int_t lineWidth =
                    (Int_t)(PEnergy / TMath::Min(fEventManager->GetMaxEnergy(), (Float_t)MaxEnergyLimit) * 15.0);
                if (fVerbose > 3)
                    cout << "lineWidth: " << lineWidth << " for track " << track->GetPdg() << " - " << P->GetTitle()
                         << endl;

                if (lineWidth > 0)
                {
                    track->SetLineWidth(lineWidth);
                    if (P->GetMother(0) > -1)
                        track->SetLineStyle(2);
                }
                else
                {
                    track->SetLineStyle(1);
                    track->SetLineWidth(1);
                }
            }

            for (Int_t n = 0; n < Np; n++)
            {
                point = tr->GetPoint(n);
                track->SetPoint(n, point[0], point[1], point[2]);
                TEveVector pos = TEveVector(point[0], point[1], point[2]);
                TEvePathMark* path = new TEvePathMark();
                path->fV = pos;
                path->fTime = point[3];
                if (n == 0)
                {
                    TEveVector Mom = TEveVector(P->Px(), P->Py(), P->Pz());
                    path->fP = Mom;
                }
                if (fVerbose > 3)
                    cout << "Path marker added " << path << endl;

#if ROOT_VERSION_CODE <= ROOT_VERSION(5, 18, 0)
                track->AddPathMark(path);
#else
                track->AddPathMark(*path);
#endif
                if (fVerbose > 3)
                    cout << "Path marker added " << path << endl;
            }
            fTrList->AddElement(track);
            if (fVerbose > 3)
                cout << "track added " << track->GetName() << endl;
        }

        for (Int_t i = 0; i < fEveTrList->GetEntriesFast(); i++)
        {
            // TEveTrackList *TrListIn=( TEveTrackList *) fEveTrList->At(i);
            // TrListIn->FindMomentumLimits(TrListIn, kFALSE);
        }
        fEventManager->SetEvtMaxEnergy(MaxEnergyLimit);
        fEventManager->SetEvtMinEnergy(MinEnergyLimit);
        gEve->Redraw3D(kFALSE);
    }
}

R3BMCTracks::~R3BMCTracks() {}

void R3BMCTracks::SetParContainers() {}

void R3BMCTracks::Finish() {}

void R3BMCTracks::Reset()
{
    for (Int_t i = 0; i < fEveTrList->GetEntriesFast(); i++)
    {
        TEveTrackList* ele = static_cast<TEveTrackList*>(fEveTrList->At(i));
        gEve->RemoveElement(ele, fEventManager);
    }
    fEveTrList->Clear();
}

TEveTrackList* R3BMCTracks::GetTrGroup(TParticle* P)
{
    fTrList = 0;
    for (Int_t i = 0; i < fEveTrList->GetEntriesFast(); i++)
    {
        TEveTrackList* TrListIn = static_cast<TEveTrackList*>(fEveTrList->At(i));
        if (strcmp(TrListIn->GetName(), P->GetName()) == 0)
        {
            fTrList = TrListIn;
            break;
        }
    }
    if (fTrList == 0)
    {
        fTrPr = new TEveTrackPropagator();
        fTrList = new TEveTrackList(P->GetName(), fTrPr);
        fTrList->SetMainColor(fEventManager->Color(P->GetPdgCode()));
        fEveTrList->Add(fTrList);
        gEve->AddElement(fTrList, fEventManager);
        fTrList->SetRnrLine(kTRUE);
    }
    return fTrList;
}

ClassImp(R3BMCTracks);
