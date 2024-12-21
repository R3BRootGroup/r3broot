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

#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

#include "R3BCalifaCrystalCal2Cluster.h"
#include "R3BCalifaCrystalCalData.h"
#include "R3BCalifaGeometry.h"
#include "R3BCalifaMappingPar.h"
#include "R3BLogger.h"
#include "R3BTGeoPar.h"

#include <TClonesArray.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TRandom.h>
#include <TVector3.h>
#include <cmath>
#include <list>
#include <vector>

using namespace std;

struct califa_candidate
{
    Int_t motherId;
    vector<uint16_t> crystalList;
    Double_t energy;
    Double_t ns;
    Double_t nf;
    Double_t theta;
    Double_t phi;
    ULong64_t time;
};

void R3BCalifaCrystalCal2Cluster::SetRectangularWindow(Float_t thetaL, Float_t phiL)
{
    fWindowAlg = "Rectangular";

    fThetaLimit = thetaL;
    fPhiLimit = phiL;
};

void R3BCalifaCrystalCal2Cluster::SetRoundWindow(Double_t window)
{
    fWindowAlg = "Round";
    fRoundWindow = window;
};

bool compareByEnergy(R3BCalifaCrystalCalData* a, R3BCalifaCrystalCalData* b) { return a->GetEnergy() > b->GetEnergy(); }

bool isInside(vector<uint16_t>& vec, Int_t cryId) { return find(vec.begin(), vec.end(), cryId) != vec.end(); }

bool R3BCalifaCrystalCal2Cluster::InsideClusterWindow(TVector3 mother, TVector3 crystal)
{

    if (fWindowAlg == "Round")
        return mother.Angle(crystal) <= fRoundWindow;

    if (fWindowAlg == "Rectangular")
        return (abs(mother.Theta() - crystal.Theta()) < fThetaLimit) && (abs(mother.Phi() - crystal.Phi()) < fPhiLimit);

    else
        return 0;
}

void RemoveUsedCrystals(vector<uint16_t>& used,
                        vector<R3BCalifaCrystalCalData*>& all,
                        vector<R3BCalifaCrystalCalData*>& proton,
                        vector<R3BCalifaCrystalCalData*>& gamma,
                        vector<R3BCalifaCrystalCalData*>& saturated)
{

    for (int p = 0; p < used.size(); p++)
    {
        for (int s = 0; s < all.size(); s++)

            if (all.at(s)->GetCrystalId() == used.at(p))
                all.erase(all.begin() + s);
    }

    for (int p = 0; p < used.size(); p++)
    {
        for (int s = 0; s < gamma.size(); s++)

            if (gamma.at(s)->GetCrystalId() == used.at(p))
                gamma.erase(gamma.begin() + s);
    }

    for (int p = 0; p < used.size(); p++)
    {
        for (int s = 0; s < proton.size(); s++)

            if (proton.at(s)->GetCrystalId() == used.at(p))
                proton.erase(proton.begin() + s);
    }

    for (int p = 0; p < used.size(); p++)
    {
        for (int s = 0; s < saturated.size(); s++)

            if (saturated.at(s)->GetCrystalId() == used.at(p))
                saturated.erase(saturated.begin() + s);
    }
}

void addCrystal2Cluster(struct califa_candidate* cluster,
                        R3BCalifaCrystalCalData* crystalCal,
                        string range,
                        vector<uint16_t>* usedCrystals,
                        uint16_t totalCrystals)
{
    cluster->energy += crystalCal->GetEnergy();
    cluster->ns += crystalCal->GetNs();
    cluster->nf += crystalCal->GetNf();
    cluster->crystalList.push_back(crystalCal->GetCrystalId());

    usedCrystals->push_back(crystalCal->GetCrystalId());

    if (range == "gamma")
        usedCrystals->push_back(crystalCal->GetCrystalId() + totalCrystals);

    if (range == "proton")
        usedCrystals->push_back(crystalCal->GetCrystalId() - totalCrystals);
}

R3BCalifaCrystalCal2Cluster::R3BCalifaCrystalCal2Cluster()
    : FairTask("R3B CALIFA CrystalCal to Cluster Finder")
    , fCrystalCalData(NULL)
    , fCalifaClusterData(NULL)
    , fGeometryVersion(2024)
    , fCrystalThreshold(0)
    , fGammaClusterThreshold(0)
    , fRoundWindow(0.25)
    , fSimulation(kFALSE)
    , fTargetGeoPar(NULL)
    , fCalifaGeoPar(NULL)
    , fWindowAlg("Round")
    , fOnline(kFALSE)
    , fRand(0)
    , fTotalCrystals(2544)
    , fRandFile("")
{
}

R3BCalifaCrystalCal2Cluster::~R3BCalifaCrystalCal2Cluster()
{
    R3BLOG(debug1, "");
    if (fCalifaClusterData)
        delete fCalifaClusterData;
}

void R3BCalifaCrystalCal2Cluster::SetParContainers()
{
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(fatal, !rtdb, "FairRuntimeDb not found");

    fCalifaGeoPar = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer("CalifaGeoPar"));
    fTargetGeoPar = dynamic_cast<R3BTGeoPar*>(rtdb->getContainer("TargetGeoPar"));

    if (!fCalifaGeoPar)
    {
        R3BLOG_IF(warn, !fCalifaGeoPar, "Could not get access to CalifaGeoPar container. Setting nominal position.");
        fCalifaPos.SetXYZ(0.0, 0.0, 0.0);
    }
    else
    {
        R3BLOG(info, "Container CalifaGeoPar found.");
        fCalifaPos.SetXYZ(fCalifaGeoPar->GetPosX(), fCalifaGeoPar->GetPosY(), fCalifaGeoPar->GetPosZ());
    }

    if (!fCalifaGeoPar || !fTargetGeoPar)
    {
        R3BLOG_IF(warn, !fTargetGeoPar, "Could not get access to TargetGeoPar container. Setting nominal position. ");
        fTargetPos.SetXYZ(0.0, 0.0, 0.0);
    }

    else
    {
        R3BLOG(info, "Container TargetGeoPar found.");
        fTargetPos.SetXYZ(fTargetGeoPar->GetPosX(), fTargetGeoPar->GetPosY(), fTargetGeoPar->GetPosZ());
    }

    LOG(info) << "Califa position : X = " << fCalifaPos.X() << " Y = " << fCalifaPos.Y() << " Z = " << fCalifaPos.Z();
    LOG(info) << "Target position : X = " << fTargetPos.X() << " Y = " << fTargetPos.Y() << " Z = " << fTargetPos.Z();

    return;
}

InitStatus R3BCalifaCrystalCal2Cluster::Init()
{
    R3BLOG(info, "");
    assert(!fCalifaClusterData); // in case someone calls Init() twice.

    auto* rootman = FairRootManager::Instance();
    R3BLOG_IF(fatal, !rootman, "FairRootManager not found");

    fCrystalCalData = dynamic_cast<TClonesArray*>(rootman->GetObject("CalifaCrystalCalData"));

    fCalifaClusterData = new TClonesArray("R3BCalifaClusterData");
    rootman->Register("CalifaClusterData", "CALIFA Cluster", fCalifaClusterData, !fOnline);

    if (R3BCalifaGeometry::Instance()->IsSimulation() == true)
    {
        fTotalCrystals = R3BCalifaGeometry::Instance()->GetNbCrystals() / 2;
        R3BCalifaGeometry::Instance()->SetReferencePoint(fTargetPos);
        R3BLOG(info, "Working on a simulation framework");
        R3BLOG(info, "Position correction: ");
        fTargetPos.Print();
    }
    else
    {
        R3BCalifaGeometry::Instance()->Init(fGeometryVersion);
        fTotalCrystals = R3BCalifaGeometry::Instance()->GetNbCrystals() / 2;
        auto posCor = fTargetPos - fCalifaPos;
        R3BCalifaGeometry::Instance()->SetReferencePoint(posCor);
        R3BLOG(info, "Working on an analysis framework");
        R3BLOG(info, "Position correction: ");
        posCor.Print();
    }

    if (fRand)
    {
        R3BLOG_IF(fatal, !fHistoFile, "Randomization file not found");

        fAngularDistributions = new TH2F*[fTotalCrystals];

        char name[100];

        for (Int_t i = 0; i < fTotalCrystals; i++)
        {
            sprintf(name, "distributionCrystalID_%i", i + 1);
            fHistoFile->GetObject(name, fAngularDistributions[i]);
        }
    }

    return kSUCCESS;
}

InitStatus R3BCalifaCrystalCal2Cluster::ReInit()
{
    SetParContainers();
    return kSUCCESS;
}

void R3BCalifaCrystalCal2Cluster::Exec(Option_t* /*opt*/)
{
    Reset();

    const int numCrystalHits = fCrystalCalData->GetEntriesFast();

    if (numCrystalHits == 0)
        return;

    R3BLOG(debug1, "Crystal hits at start:" << numCrystalHits);

    vector<R3BCalifaCrystalCalData*> allCrystalVec;
    vector<R3BCalifaCrystalCalData*> protonCandidatesVec;
    vector<R3BCalifaCrystalCalData*> gammaCandidatesVec;
    vector<R3BCalifaCrystalCalData*> saturatedCandidatesVec;

    vector<uint16_t> usedCrystals;

    Double_t cryEnergy = 0.;
    Int_t cryId = 0;

    for (Int_t i = 0; i < numCrystalHits; i++)
    {
        cryId = dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalData->At(i))->GetCrystalId();
        cryEnergy = dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalData->At(i))->GetEnergy();

        if (cryEnergy >= fCrystalThreshold)
            allCrystalVec.push_back(dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalData->At(i)));

        /* -------- Real data ------- */
        if (!fSimulation)
        {
            if (cryEnergy >= fGammaClusterThreshold && !std::isnan(cryEnergy))
                gammaCandidatesVec.push_back(dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalData->At(i)));

            if (cryId > fTotalCrystals && cryEnergy >= fProtonClusterThreshold)
                protonCandidatesVec.push_back(dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalData->At(i)));

            if (cryId <= fTotalCrystals && std::isnan(cryEnergy))
                saturatedCandidatesVec.push_back(dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalData->At(i)));
        }

        /* -------- Simulation data ------- */
        else
        {
            if (cryEnergy >= fGammaClusterThreshold && cryEnergy < fProtonClusterThreshold)
                gammaCandidatesVec.push_back(dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalData->At(i)));

            if (cryEnergy >= fProtonClusterThreshold)
                protonCandidatesVec.push_back(dynamic_cast<R3BCalifaCrystalCalData*>(fCrystalCalData->At(i)));
        }
    }

    /* ----- Remove duplicate entries in gamma clusters ----- */
    for (int k = 0; k < gammaCandidatesVec.size(); k++)
    {
        for (int s = 0; s < gammaCandidatesVec.size(); s++)
        {
            if (gammaCandidatesVec.at(k)->GetCrystalId() == (gammaCandidatesVec.at(s)->GetCrystalId() - fTotalCrystals))
                gammaCandidatesVec.erase(gammaCandidatesVec.begin() + s);
        }
    }

    // Sort all vectors by energy
    std::sort(gammaCandidatesVec.begin(), gammaCandidatesVec.end(), compareByEnergy);
    std::sort(protonCandidatesVec.begin(), protonCandidatesVec.end(), compareByEnergy);
    std::sort(allCrystalVec.begin(), allCrystalVec.end(), compareByEnergy);

    TVector3 mother_angles, angles;
    Double_t fRandTheta = 0., fRandPhi = 0.;

    while (protonCandidatesVec.size())
    {
        Int_t motherId = protonCandidatesVec.at(0)->GetCrystalId();

        califa_candidate cluster = { motherId, vector<uint16_t>(), 0.0, 0.0, 0.0, 0.0, 0.0 };

        if (motherId <= fTotalCrystals)
        {
            mother_angles = R3BCalifaGeometry::Instance()->GetAngles(motherId);
        }
        else
        {
            mother_angles = R3BCalifaGeometry::Instance()->GetAngles(motherId - fTotalCrystals);
        }

        if (fRand)
        {

            if (motherId > fTotalCrystals)
                fAngularDistributions[protonCandidatesVec.at(0)->GetCrystalId() - 1 - fTotalCrystals]->GetRandom2(
                    fRandPhi, fRandTheta);

            else
                fAngularDistributions[protonCandidatesVec.at(0)->GetCrystalId() - 1]->GetRandom2(fRandPhi, fRandTheta);

            cluster.theta = TMath::DegToRad() * fRandTheta;
            cluster.phi = TMath::DegToRad() * fRandPhi;
        }

        else
        {
            cluster.theta = mother_angles.Theta();
            cluster.phi = mother_angles.Phi();
        }

        cluster.time = protonCandidatesVec.at(0)->GetTime();

        addCrystal2Cluster(&cluster, protonCandidatesVec.at(0), "proton", &usedCrystals, fTotalCrystals);

        for (Int_t j = 0; j < allCrystalVec.size(); j++)
        {
            Int_t thisCryId = allCrystalVec.at(j)->GetCrystalId();
            Float_t thisEnergy = allCrystalVec.at(j)->GetEnergy();

            if (thisCryId > fTotalCrystals && !isInside(usedCrystals, thisCryId) && !fSimulation)
            {

                angles = R3BCalifaGeometry::Instance()->GetAngles(thisCryId - fTotalCrystals);

                if (InsideClusterWindow(mother_angles, angles))
                    addCrystal2Cluster(&cluster, allCrystalVec.at(j), "proton", &usedCrystals, fTotalCrystals);
            }

            if (thisCryId <= fTotalCrystals && !isInside(usedCrystals, thisCryId) && !std::isnan(thisEnergy) &&
                !fSimulation)
            {

                angles = R3BCalifaGeometry::Instance()->GetAngles(thisCryId);

                if (InsideClusterWindow(mother_angles, angles))
                    addCrystal2Cluster(&cluster, allCrystalVec.at(j), "gamma", &usedCrystals, fTotalCrystals);
            }

            /* ------- Simulation ------- */
            if (!isInside(usedCrystals, thisCryId) && fSimulation)
            {

                angles = R3BCalifaGeometry::Instance()->GetAngles(thisCryId);

                if (InsideClusterWindow(mother_angles, angles))
                    addCrystal2Cluster(&cluster, allCrystalVec.at(j), "gamma", &usedCrystals, fTotalCrystals);
            }
        }

        AddCluster(
            cluster.crystalList, cluster.energy, cluster.nf, cluster.ns, cluster.theta, cluster.phi, cluster.time, 0);

        RemoveUsedCrystals(
            usedCrystals, allCrystalVec, protonCandidatesVec, gammaCandidatesVec, saturatedCandidatesVec);
    }

    /*------ Gamma Clusters ------- */
    while (gammaCandidatesVec.size())
    {
        Int_t motherId = gammaCandidatesVec.at(0)->GetCrystalId();

        califa_candidate cluster = { motherId, vector<uint16_t>(), 0.0, 0.0, 0.0, 0.0, 0.0 };

        if (motherId <= fTotalCrystals)
            mother_angles = R3BCalifaGeometry::Instance()->GetAngles(motherId);

        else
            mother_angles = R3BCalifaGeometry::Instance()->GetAngles(motherId - fTotalCrystals);

        if (fRand)
        {
            if (motherId > fTotalCrystals)
                fAngularDistributions[gammaCandidatesVec.at(0)->GetCrystalId() - 1 - fTotalCrystals]->GetRandom2(
                    fRandPhi, fRandTheta);

            else
                fAngularDistributions[gammaCandidatesVec.at(0)->GetCrystalId() - 1]->GetRandom2(fRandPhi, fRandTheta);

            cluster.theta = TMath::DegToRad() * fRandTheta;
            cluster.phi = TMath::DegToRad() * fRandPhi;
        }

        else
        {

            cluster.theta = mother_angles.Theta();
            cluster.phi = mother_angles.Phi();
        }

        cluster.time = gammaCandidatesVec.at(0)->GetTime();

        addCrystal2Cluster(&cluster, gammaCandidatesVec.at(0), "gamma", &usedCrystals, fTotalCrystals);

        for (Int_t j = 0; j < allCrystalVec.size(); j++)
        {
            Int_t thisCryId = allCrystalVec.at(j)->GetCrystalId();

            if (thisCryId <= fTotalCrystals && !isInside(usedCrystals, thisCryId) && !fSimulation)
            {
                angles = R3BCalifaGeometry::Instance()->GetAngles(thisCryId);

                if (InsideClusterWindow(mother_angles, angles))
                    addCrystal2Cluster(&cluster, allCrystalVec.at(j), "gamma", &usedCrystals, fTotalCrystals);
            }
        }

        for (Int_t j = 0; j < allCrystalVec.size(); j++)
        {
            Int_t thisCryId = allCrystalVec.at(j)->GetCrystalId();

            if (thisCryId > fTotalCrystals && !isInside(usedCrystals, thisCryId) && !fSimulation)
            {
                angles = R3BCalifaGeometry::Instance()->GetAngles(thisCryId - fTotalCrystals);

                if (InsideClusterWindow(mother_angles, angles))
                    addCrystal2Cluster(&cluster, allCrystalVec.at(j), "proton", &usedCrystals, fTotalCrystals);
            }
        }

        for (Int_t j = 0; j < allCrystalVec.size(); j++)
        {
            Int_t thisCryId = allCrystalVec.at(j)->GetCrystalId();

            if (!isInside(usedCrystals, thisCryId) && fSimulation)
            {

                angles = R3BCalifaGeometry::Instance()->GetAngles(thisCryId);

                if (InsideClusterWindow(mother_angles, angles))
                    addCrystal2Cluster(&cluster, allCrystalVec.at(j), "gamma", &usedCrystals, fTotalCrystals);
            }
        }

        AddCluster(
            cluster.crystalList, cluster.energy, cluster.nf, cluster.ns, cluster.theta, cluster.phi, cluster.time, 1);

        RemoveUsedCrystals(
            usedCrystals, allCrystalVec, protonCandidatesVec, gammaCandidatesVec, saturatedCandidatesVec);
    }

    /* ----------- Saturation Clusters ----------- */
    while (saturatedCandidatesVec.size())
    {
        Int_t motherId = saturatedCandidatesVec.at(0)->GetCrystalId();

        califa_candidate cluster = { motherId, vector<uint16_t>(), 0.0, 0.0, 0.0, 0.0, 0.0 };

        mother_angles = R3BCalifaGeometry::Instance()->GetAngles(motherId);

        if (fRand)
        {
            fAngularDistributions[saturatedCandidatesVec.at(0)->GetCrystalId() - 1]->GetRandom2(fRandPhi, fRandTheta);
            cluster.theta = TMath::DegToRad() * fRandTheta;
            cluster.phi = TMath::DegToRad() * fRandPhi;
        }
        else
        {
            cluster.theta = mother_angles.Theta();
            cluster.phi = mother_angles.Phi();
        }
        cluster.time = saturatedCandidatesVec.at(0)->GetTime();

        addCrystal2Cluster(&cluster, saturatedCandidatesVec.at(0), "gamma", &usedCrystals, fTotalCrystals);

        for (Int_t j = 0; j < allCrystalVec.size(); j++)
        {
            Int_t thisCryId = allCrystalVec.at(j)->GetCrystalId();

            if (thisCryId > fTotalCrystals && !isInside(usedCrystals, thisCryId))
            {
                angles = R3BCalifaGeometry::Instance()->GetAngles(thisCryId - fTotalCrystals);

                if (InsideClusterWindow(mother_angles, angles))
                    addCrystal2Cluster(&cluster, allCrystalVec.at(j), "proton", &usedCrystals, fTotalCrystals);
            }

            if (thisCryId <= fTotalCrystals && !isInside(usedCrystals, thisCryId))
            {
                angles = R3BCalifaGeometry::Instance()->GetAngles(thisCryId);

                if (InsideClusterWindow(mother_angles, angles))
                    addCrystal2Cluster(&cluster, allCrystalVec.at(j), "gamma", &usedCrystals, fTotalCrystals);
            }
        }

        AddCluster(
            cluster.crystalList, cluster.energy, cluster.nf, cluster.ns, cluster.theta, cluster.phi, cluster.time, 2);

        RemoveUsedCrystals(
            usedCrystals, allCrystalVec, protonCandidatesVec, gammaCandidatesVec, saturatedCandidatesVec);
    }

    return;
}

void R3BCalifaCrystalCal2Cluster::Reset()
{
    // Clear the CA structure
    R3BLOG(debug1, "Clearing CalifaClusterData Structure");
    if (fCalifaClusterData)
        fCalifaClusterData->Clear();
}

void R3BCalifaCrystalCal2Cluster::SelectGeometryVersion(Int_t version)
{
    fGeometryVersion = version;
    R3BLOG(info, "to " << fGeometryVersion);
}

R3BCalifaClusterData* R3BCalifaCrystalCal2Cluster::AddCluster(vector<uint16_t> crystalList,
                                                              Double_t ene,
                                                              Double_t Nf,
                                                              Double_t Ns,
                                                              Double_t pAngle,
                                                              Double_t aAngle,
                                                              ULong64_t time,
                                                              uint8_t clusterType)
{
    TClonesArray& clref = *fCalifaClusterData;
    Int_t size = clref.GetEntriesFast();
    return new (clref[size]) R3BCalifaClusterData(crystalList, ene, Nf, Ns, pAngle, aAngle, time, clusterType);
}

ClassImp(R3BCalifaCrystalCal2Cluster)
