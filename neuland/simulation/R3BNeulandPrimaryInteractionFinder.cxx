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

#include "R3BNeulandPrimaryInteractionFinder.h"
#include "FairLogger.h"
#include "TFile.h"
#include <FairRootManager.h>
#include <utility>
#include <vector>

bool IsPrimaryTrack(R3BMCTrack* track)
{
    // TODO: The test can be modified to rely on some other information,
    // e.g. if the neutrons were created in a reaction in the target simulated with Geant
    return track->GetMotherId() == -1 && track->GetPdgCode() == 2112;
}

// Map a hit to every point if possible (else -1)
std::map<R3BNeulandPoint*, R3BNeulandHit*> MapPointsToHits(const std::vector<R3BNeulandPoint*>& points,
                                                           const std::vector<R3BNeulandHit*>& hits)
{
    std::map<R3BNeulandPoint*, R3BNeulandHit*> m;
    for (const auto& point : points)
    {
        m[point] = nullptr;
        for (const auto hit : hits)
        {
            if (point->GetPaddle() == hit->GetPaddle())
            {
                m[point] = hit;
            }
        }
    }
    return m;
};

// Map a primary track to every point if possible (else -1)
std::map<R3BNeulandPoint*, R3BMCTrack*> MapPointsToPrimaryTracks(const std::vector<R3BNeulandPoint*>& points,
                                                                 const std::vector<R3BMCTrack*>& tracks)
{
    std::map<R3BNeulandPoint*, R3BMCTrack*> m;
    for (const auto& point : points)
    {
        m[point] = nullptr;
        auto iTrack = point->GetTrackID();
        while (iTrack > -1)
        {
            const auto track = tracks.at(iTrack);
            if (IsPrimaryTrack(track))
            {
                m[point] = track;
                break;
            }
            // Else, start tracing back:
            iTrack = track->GetMotherId();
        }
    }
    return m;
}

R3BNeulandPoint* FindFirstPoint(const R3BMCTrack* track, const std::map<R3BNeulandPoint*, R3BMCTrack*>& p2t)
{
    R3BNeulandPoint* minPoint = nullptr;
    //  Search for min ToF:
    double minToF = 1e99;

    for (const auto& pt : p2t)
    {
        // Only look at points traced back to this track:
        if (pt.second == track)
        {
            const auto ToF = pt.first->GetTime();
            if (ToF < minToF)
            {
                minToF = ToF;
                minPoint = pt.first;
            }
        }
    }
    return minPoint;
}

R3BNeulandHit* FindFirstHit(const R3BMCTrack* track,
                            const std::map<R3BNeulandPoint*, R3BMCTrack*>& p2t,
                            const std::map<R3BNeulandPoint*, R3BNeulandHit*>& p2h)
{
    R3BNeulandPoint* minHitPoint = nullptr;
    //  Search for min ToF:
    double minHitToF = 1e99;

    for (const auto& pt : p2t)
    {
        // Only look at points traced back to this track:
        // AND where a hit is registered
        if (pt.second == track && p2h.at(pt.first) != nullptr)
        {
            const auto ToF = pt.first->GetTime();
            if (ToF < minHitToF)
            {
                minHitToF = ToF;
                minHitPoint = pt.first;
            }
        }
    }
    return minHitPoint ? p2h.at(minHitPoint) : nullptr;
}

R3BNeulandPrimaryInteractionFinder::R3BNeulandPrimaryInteractionFinder(TString pointsIn,
                                                                       TString hitsIn,
                                                                       TString pointsOut,
                                                                       TString hitsOut,
                                                                       TString tracksOut)
    : FairTask("R3BNeulandPrimaryInteractionFinder")
    , fTracksIn("MCTrack")
    , fPointsIn(std::move(pointsIn))
    , fHitsIn(std::move(hitsIn))
    , fTracksOut(std::move(tracksOut))
    , fPointsOut(std::move(pointsOut))
    , fHitsOut(std::move(hitsOut))
    , fhDistance(new TH1D("fhDistance", "Distance firstPoint to firstHit", 10000, 0, 1000))
    , fhPointsVsHits(new TH2D("fhPointsVsHits", "Number of Primary Points vs Number of Primary Hits", 6, 0, 6, 6, 0, 6))
    , fhPointVsHitPaddle(
          new TH2D("fhPointVsHitPaddle", "First Point Paddle vs First Hit Paddle", 3001, -1, 3000, 3001, -1, 3000))
{
    fhDistance->GetXaxis()->SetTitle("Distance [cm]");
    fhPointsVsHits->GetXaxis()->SetTitle("# First Points");
    fhPointsVsHits->GetYaxis()->SetTitle("# First Hits");
    fhPointVsHitPaddle->GetXaxis()->SetTitle("First Point Paddle");
    fhPointVsHitPaddle->GetYaxis()->SetTitle("First Hit Paddle");
}

InitStatus R3BNeulandPrimaryInteractionFinder::Init()
{
    fTracksIn.Init();
    fPointsIn.Init();
    fHitsIn.Init();
    fTracksOut.Init();
    fPointsOut.Init();
    fHitsOut.Init();
    return kSUCCESS;
}

void R3BNeulandPrimaryInteractionFinder::Exec(Option_t*)
{
    const auto tracks = fTracksIn.Retrieve();
    const auto points = fPointsIn.Retrieve();
    const auto hits = fHitsIn.Retrieve();
    fTracksOut.Reset();
    fPointsOut.Reset();
    fHitsOut.Reset();

    const auto p2h = MapPointsToHits(points, hits);
    const auto p2t = MapPointsToPrimaryTracks(points, tracks);

    if (FairLogger::GetLogger()->IsLogNeeded(fair::Severity::debug))
    {
        LOG(debug) << "R3BNeulandPrimaryInteractionFinder: Points without Hit in: ";
        for (const auto& kv : p2h)
        {
            if (kv.second == nullptr)
            {
                const auto point = kv.first;
                LOG(debug) << point->GetDetectorID() << ":" << tracks.at(point->GetTrackID())->GetPdgCode() << ":"
                           << point->GetLightYield() << ":" << point->GetEnergyLoss() << "\t";
            }
        }
    }

    for (const auto& track : tracks)
    {
        if (IsPrimaryTrack(track))
        {
            fTracksOut.Insert(track);

            const auto firstPoint = FindFirstPoint(track, p2t);
            const auto firstHit = FindFirstHit(track, p2t, p2h);

            if (firstPoint)
            {
                fPointsOut.Insert(firstPoint);
            }

            if (firstHit)
            {
                // NOTE: This may add the same hit multiple times.
                fHitsOut.Insert(firstHit);
            }

            if (firstHit && firstPoint)
            {
                fhDistance->Fill((firstPoint->GetPosition() - firstHit->GetPosition()).Mag());
            }

            fhPointVsHitPaddle->Fill(firstPoint ? firstPoint->GetPaddle() : -1, firstHit ? firstHit->GetPaddle() : -1);
        }
    }

    fhPointsVsHits->Fill(fPointsOut.Size(), fHitsOut.Size());
}

void R3BNeulandPrimaryInteractionFinder::Finish()
{
    TDirectory* tmp = gDirectory;
    FairRootManager::Instance()->GetOutFile()->cd();
    gDirectory->mkdir("R3BNeulandPrimaryInteractionFinder");
    gDirectory->cd("R3BNeulandPrimaryInteractionFinder");

    fhDistance->Write();
    fhPointsVsHits->Write();
    fhPointVsHitPaddle->Write();

    gDirectory = tmp;
}
