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
// -----              R3BTttxCal2Hit source file                       -----
// -----               Created 18/01/24  by XXX                        -----
// -------------------------------------------------------------------------

// ROOT headers
#include <TClonesArray.h>
#include <TMath.h>
#include <TRandom.h>
#include <iomanip>

// Fair headers
#include <FairLogger.h>
#include <FairRootManager.h>
#include <FairRunAna.h>
#include <FairRuntimeDb.h>

// R3B headers
#include "R3BLogger.h"
#include "R3BTttxCal2Hit.h"
#include "R3BTttxCalData.h"
#include "R3BTttxHitData.h"
#include "R3BTttxHitPar.h"

// R3BTttxCal2Hit: Default Constructor --------------------------
R3BTttxCal2Hit::R3BTttxCal2Hit()
    : R3BTttxCal2Hit("R3BTttxCal2Hit", 1)
{
}

// R3BTttxCal2HitPar: Standard Constructor --------------------------
R3BTttxCal2Hit::R3BTttxCal2Hit(const TString& name, int iVerbose)
    : FairTask(name, iVerbose)
{
}

// Virtual R3BTttxCal2Hit: Destructor
R3BTttxCal2Hit::~R3BTttxCal2Hit()
{
    Reset();
    if (fTttxHitData)
        delete fTttxHitData;
}

// Constructor and destructor for HitDetector
R3BTttxCal2Hit::HitDetector::HitDetector(R3BTttxCalData* hit)
    : fX(StripID2X(hit->GetStripID()))
    , fE(static_cast<double>(hit->GetEnergy()))
    , fE2(static_cast<double>(hit->GetEnergy() * hit->GetEnergy()))
    , fT(static_cast<double>(hit->GetTime()))
    , fDetID(hit->GetDetID())
    , fTrefIncluded(false)
{
    R3BLOG(debug, "A hit constructed.");
}

R3BTttxCal2Hit::HitDetector::~HitDetector() { R3BLOG(debug, "Destructor is called."); }

void R3BTttxCal2Hit::SetParContainers()
{
    // Parameter Containers
    FairRuntimeDb* rtdb = FairRuntimeDb::instance();
    R3BLOG_IF(error, rtdb == nullptr, "FairRuntimeDb not found");

    fHit_Par = dynamic_cast<R3BTttxHitPar*>(rtdb->getContainer("tttxHitPar"));
    if (!fHit_Par)
    {
        R3BLOG(error, "Couldn't get handle on tttxHitPar container");
    }
    else
    {
        R3BLOG(info, "tttxHitPar container open");
    }
}

void R3BTttxCal2Hit::SetParameter()
{
    //--- Parameter Container ---
    NumDets = fHit_Par->GetNumDets();     // Number of Detectors
    NumPars = fHit_Par->GetNumParsZfit(); // Number of Parameters per each detector
    R3BLOG(info, "Nb detectors: " << NumDets << ", Nb of Parameters per det: " << NumPars);

    HitParams = new TArrayF();
    HitParams->Set(NumDets * NumPars);
    HitParams = const_cast<TArrayF*>(fHit_Par->GetZfitParams());

    vHitsDet.clear();
    vHitsDet.resize(NumDets);

    ResetVHitsDetector();
}

// -----   Public method Init   --------------------------------------------
InitStatus R3BTttxCal2Hit::Init()
{
    R3BLOG(info, "");

    FairRootManager* rootManager = FairRootManager::Instance();
    R3BLOG_IF(fatal, rootManager == nullptr, "FairRootManager not found");

    // INPUT DATA
    fTttxCalData = dynamic_cast<TClonesArray*>(rootManager->GetObject("tttxCalData"));
    R3BLOG_IF(fatal, fTttxCalData == nullptr, "tttxCalData not found");

    // OUTPUT DATA
    fTttxHitData = new TClonesArray("R3BTttxHitData");
    rootManager->Register("tttxHitData", "TTTX strip Hit", fTttxHitData, fOnline ? false : true);

    SetParameter();
    return kSUCCESS;
}

void R3BTttxCal2Hit::HitsDetector::Init()
{
    for (auto& fhit : fHits)
    {
        if (fhit)
            delete fhit;
    }
    fHits.clear();
    fDetID = 0;
}

void R3BTttxCal2Hit::InitVHitsDetector()
{
    for (auto& fhitdet : vHitsDet)
    {
        fhitdet->Init();
    }
}

// -----   Public method ReInit   ----------------------------------------------
InitStatus R3BTttxCal2Hit::ReInit()
{
    SetParContainers();
    SetParameter();
    return kSUCCESS;
}

// -----   Public method Reset   ------------------------------------------------
void R3BTttxCal2Hit::Reset()
{
    R3BLOG(debug, "Clearing tttxHitData structure");
    if (fTttxHitData)
        fTttxHitData->Clear();
    InitVHitsDetector();
    vUsed.clear();
}

void R3BTttxCal2Hit::ResetVHitsDetector()
{
    for (auto& fhits : vHitsDet)
    {
        fhits = new HitsDetector();
    }
    InitVHitsDetector();
}

// -----   Public method Execution   --------------------------------------------
void R3BTttxCal2Hit::Exec(Option_t* /*option*/)
{
    // Reset entries in output arrays, local arrays
    Reset();

    if (!fHit_Par)
    {
        LOG(error) << "NO Container Parameter";
    }

    // Reading the Input -- Cal Data --
    auto nHits = fTttxCalData->GetEntriesFast();
    R3BLOG(debug, "Number of Cal Data = " << nHits);
    if (!nHits)
        return;
    vUsed.resize(nHits, false);

    // Sort with energy in descending order
    fTttxCalData->Sort();

    // Clustering the hits
    for (int ihit = 0; ihit < nHits; ihit++)
    {
        auto* fhit = dynamic_cast<R3BTttxCalData*>(fTttxCalData->At(ihit));
        R3BLOG(debug, "Parent hit of the cluster: " << *fhit); // toString() is formatted for TttxCalData
        auto idet = fhit->GetDetID() - 1;
        vHitsDet.at(idet)->NewCluster(fhit);
        vUsed.at(ihit) = true;
        for (int jhit = ihit + 1; jhit < nHits; jhit++)
        {
            if (vUsed.at(jhit))
                continue;
            auto* ghit = dynamic_cast<R3BTttxCalData*>(fTttxCalData->At(jhit));
            if ((idet != ghit->GetDetID() - 1) || !IsInWindow(fhit, ghit))
                continue;
            R3BLOG(debug, "Child hit: " << *ghit);
            vHitsDet.at(idet)->AddHitToLast(ghit);
            vUsed.at(jhit) = true;
        }
    }

    // Merge hits of detectors
    // Calculate product of multiplicities
    int product_mult = 1;
    for (auto& fDet : vHitsDet)
    {
        product_mult *= fDet->GetMult();
    }

    // Following code assumes the number of detectors is up to two
    if (fStoreOnlyTwoDetHitEvent && product_mult == 0)
    {
        return;
    }

    if (product_mult == 0 || NumDets < 2)
    { // If at least one detector has 0 multiplicity, then use the available detector
        for (auto& fDet : vHitsDet)
        {
            if (fDet->GetMult() == 0)
                continue;
            for (int imult = 0; imult < fDet->GetMult(); imult++)
            {
                if (fRequireTref && !fDet->GetHit(imult)->IsTrefIn())
                {
                    continue;
                }
                AddHitData(fDet->GetHit(imult));
            }
        }
        //    }else if(product_mult == 1){
    }
    else
    {
        auto* fDet1 = vHitsDet.at(0);
        auto* fDet2 = vHitsDet.at(1);
        for (int imult1 = 0; imult1 < fDet1->GetMult(); imult1++)
        {
            if (fRequireTref && !fDet1->GetHit(imult1)->IsTrefIn())
            {
                continue;
            }
            for (int imult2 = 0; imult2 < fDet2->GetMult(); imult2++)
            {
                if (fRequireTref && !fDet2->GetHit(imult2)->IsTrefIn())
                {
                    continue;
                }
                if (TMath::Abs(fDet1->GetHit(imult1)->GetT() - fDet2->GetHit(imult2)->GetT()) > ftimewindow)
                {
                    continue;
                }
                // Maybe better to discard the used hits in the later analysis. Here, just filling HitData for any
                // combinations which satisfy the time window.
                AddHitData(fDet1->GetHit(imult1), fDet2->GetHit(imult2));
            }
        }
    }
    // For later update
    // RecalculateZ();
    return;
}

// -----   Miscellaneous public/private methods  --------------------------------------------
/*
void R3BTttxCal2Hit::RecalculateZ(double fbeta)
{ // Reconstruct Z based on the updated velocity.
    if (fbeta <= 0)
        fbeta = ffixedbeta;
    R3BLOG(warning, "Nothing will occur with the current code.");
    //
    // To be filled later
    //
    return;
}
*/
// Compare two hits whether they are within the windows
bool R3BTttxCal2Hit::IsInWindow(R3BTttxCalData* hit1, R3BTttxCalData* hit2)
{
    // Geometrical (Strip) Window
    // Include tref (E=0) anyway if it satisfys the criteria
    if ((hit2->GetStripID() < NStrips) && (TMath::Abs(hit1->GetStripID() - hit2->GetStripID()) > fneighbours))
        return false;
    // Time Window
    if (TMath::Abs(hit1->GetTime() - hit2->GetTime()) > ftimewindow)
        return false;
    return true;
}

double R3BTttxCal2Hit::CalculateZ(R3BTttxCal2Hit::HitDetector* hit, double fbeta)
{
    if (fbeta <= 0)
        fbeta = ffixedbeta;
    auto Zet = 0.;
    auto SqEBeta = sqrt(hit->GetE()) * fbeta;
    for (int ipar = 0; ipar < NumPars; ipar++)
    {
        Zet += HitParams->GetAt((hit->GetDetID() - 1) * NumPars + ipar) * TMath::Power(SqEBeta, ipar);
    }
    return Zet;
}

void R3BTttxCal2Hit::HitDetector::AddHit(R3BTttxCalData* hit)
{
    if (hit->GetStripID() >= NStrips)
    {
        fTrefIncluded = true;
        return;
    }
    auto E2this = hit->GetEnergy() * hit->GetEnergy();
    auto E2tmp = fE2 + E2this;
    // Calculate X with arithmetic average with w_i = E_i^2
    fX = (fE2 * fX + E2this * StripID2X(hit->GetStripID())) / E2tmp;
    fE += hit->GetEnergy();
    fE2 = E2tmp;
    // fT will not be updated. Take the time from the highest E hit.
    return;
}

void R3BTttxCal2Hit::HitsDetector::NewCluster(R3BTttxCalData* hit)
{
    fDetID = hit->GetDetID();
    fHits.push_back(new HitDetector(hit));
}

// Add hit to the existing cluster (last entry)
// Assume the CalData has been sorted in Energy in descending order, and window condition is satisfied.
void R3BTttxCal2Hit::HitsDetector::AddHitToLast(R3BTttxCalData* hit)
{
    R3BLOG_IF(error, fDetID == 0, "fDetID is not set in the cluster");
    fHits.back()->AddHit(hit);
}

// -----   Private method AddHitData  --------------------------------------------
R3BTttxHitData* R3BTttxCal2Hit::AddHitData(int8_t idet, // idet0: reconstructed with two dets, idet%: hit with only det%
                                           double xpos,
                                           double energy,
                                           double ang,
                                           double charge,
                                           double time)
{
    // It fills the R3BTttxHitData
    TClonesArray& clref = *fTttxHitData;
    int size = clref.GetEntriesFast();
    return new (clref[size]) R3BTttxHitData(idet, xpos, energy, ang, charge, time);
}

R3BTttxHitData* R3BTttxCal2Hit::AddHitData(HitDetector* hit1)
{
    if (!fStoreOnlyTwoDetHitEvent)
    {
        return AddHitData(hit1->GetDetID(), hit1->GetX(), hit1->GetE(), NAN, CalculateZ(hit1), hit1->GetT());
    }
    return nullptr;
}

R3BTttxHitData* R3BTttxCal2Hit::AddHitData(HitDetector* hit1, HitDetector* hit2)
{

    // For the later development, the angle should be updated with Z positions derived from GeoPar.
    auto angle = -hit1->GetX() + hit2->GetX();

    // Just taking average, but maybe it needs modification to weigh detectors if thicknesses are different.
    auto* hitdata0 = AddHitData(0,
                                0.5 * (hit1->GetX() + hit2->GetX()),
                                0.5 * (hit1->GetE() + hit2->GetE()),
                                angle,
                                0.5 * (CalculateZ(hit1) + CalculateZ(hit2)),
                                0.5 * (hit1->GetT() + hit2->GetT()));
    // First entry will be combined hit.

    // If necessary include them, then you can add data detector by detector. This will be useful for later analysis
    // with ReconstructZ();
    if (!fStoreOnlyTwoDetHitEvent)
    {
        AddHitData(hit1->GetDetID(), hit1->GetX(), hit1->GetE(), NAN, CalculateZ(hit1), hit1->GetT());
        AddHitData(hit2->GetDetID(), hit2->GetX(), hit2->GetE(), NAN, CalculateZ(hit2), hit2->GetT());
    }
    return hitdata0;
}

ClassImp(R3BTttxCal2Hit)
