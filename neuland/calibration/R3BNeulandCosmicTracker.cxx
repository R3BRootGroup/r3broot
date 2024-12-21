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

#include "R3BNeulandCosmicTracker.h"

#include "FairLogger.h"

#include "TCanvas.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "TGraphErrors.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <numeric>

namespace Neuland
{
    namespace Calibration
    {
        // Data to calculate the range and energyloss of cosmic muons
        // taken from "doi:10.1006/adnd.2001.0861"
        // { Kinetic Energies / MeV, Stopping Powers / MeV cm²/g, Ranges / cm}
        // const std::vector<Double_t> cosmicMuonData[3] = {
        //     { 10,   14,   20,   30,   40,    80,    100,   140,   200,   300,   400,    800,    1000,   1400,
        //       2000, 3000, 4000, 8000, 10000, 14000, 20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000 },
        //     { 7.917, 6.171, 4.816, 3.734, 3.187, 2.388, 2.237, 2.082, 1.992, 1.957, 1.962, 2.033, 2.066, 2.12,
        //       2.179, 2.246, 2.293, 2.4,   2.433, 2.48,  2.528, 2.58,  2.615, 2.697, 2.722, 2.76,  2.8,   2.845 },
        //     { 0.7062, 1.285, 2.398, 4.789, 7.707, 22.66, 31.34, 49.97, 79.55, 130.3, 181.4, 381.7, 479.3, 670.22,
        //       948.9,  1400,  1840,  3534,  4358,  4975,  8347,  12200, 15970, 30430, 37370, 50800, 69950, 99670 }
        // };

        using std::isfinite;
        using std::isnan;
        using DPair = std::array<Double_t, 2>;

        constexpr auto MinPoints = 3;
        constexpr auto TimeEps = 0.0001;
        constexpr auto MaxDistance = 2 * BarSize_XY;
        constexpr auto MaxSlope = 15.;

        constexpr bool WithinBounds(const double val, const double lower_bound, const double upper_bound)
        {
            return (val >= lower_bound && val <= upper_bound);
        }

        constexpr bool WithinBounds(const double val, const DPair& bounds)
        {
            return WithinBounds(val, bounds[0], bounds[1]);
        }

        bool WithinBounds(const TVector3& point, const DPair& xBounds, const DPair& yBounds, const DPair& zBounds)
        {
            return WithinBounds(point[0], xBounds) && WithinBounds(point[1], yBounds) &&
                   WithinBounds(point[2], zBounds);
        }

        CosmicTracker::CosmicTracker()
            : fDistances{ 0. }
            , fFit("CosmicTracker:fFit", "pol1")
        {
            fDistances.reserve(MaxNumberOfPlanes);
            while (fDistances.size() < MaxNumberOfPlanes)
            {
                fDistances.push_back(fDistances.back() + BarSize_Z);
            }
            fTrack.Interactions.reserve(256);
            fBarIDs.reserve(256);
        }

        void CosmicTracker::AddPoint(Int_t barID, const Double_t pos)
        {
            fBarIDs.push_back(barID);

            const auto plane = GetPlaneNumber(barID);
            const auto bar = barID % BarsPerPlane;

            const auto zPosition = fDistances[plane] + 0.5 * BarSize_Z;

            if (IsPlaneHorizontal(plane))
            {
                const auto barPoint = fYZ.GetN();
                // ig fYZ.SetPoint(barPoint, zPosition, (bar - 0.5 * BarsPerPlane + 0.66) * BarSize_XY);
                fYZ.SetPoint(barPoint, zPosition, (bar - 0.5 * BarsPerPlane + 0.5) * BarSize_XY);
                fYZ.SetPointError(barPoint, BarUncertainty_Z, BarUncertainty_XY);

                // ig if (!isnan(pos) && pos < 0.5 * BarLength + 10.) // + some margin
                /*ig if (!isnan(pos) && fabs(pos) < 0.5 * BarLength + 10.) // + some margin
                        {
                            const auto posPoint = fXZ.GetN();
                            fXZ.SetPoint(posPoint, zPosition, pos);
                            //ig fXZ.SetPointError(posPoint, BarUncertainty_Z, 3 * BarUncertainty_XY);
                            fXZ.SetPointError(posPoint, BarUncertainty_Z, BarUncertainty_XY);
                    }*/
            }
            else
            {
                const auto barPoint = fXZ.GetN();
                fXZ.SetPoint(barPoint, zPosition, (bar - 0.5 * BarsPerPlane + 0.5) * BarSize_XY);
                fXZ.SetPointError(barPoint, BarUncertainty_Z, BarUncertainty_XY);

                // ig if (!isnan(pos) && pos < 0.5 * BarLength + 10.) // + some margin
                /*ig if (!isnan(pos) && fabs(pos) < 0.5 * BarLength + 10.) // + some margin
                        {
                            const auto posPoint = fYZ.GetN();
                            fYZ.SetPoint(posPoint, zPosition, pos);
                            //ig fYZ.SetPointError(posPoint, BarUncertainty_Z, 3 * BarUncertainty_XY);
                            fYZ.SetPointError(posPoint, BarUncertainty_Z, BarUncertainty_XY);
                    }*/
            }
        }

        const R3BNeulandCosmicTrack& CosmicTracker::GetTrack()
        {
            LOG(debug) << "CosmicTracker::Fit : Number of Points: X-Z: " << fXZ.GetN() << "   Y-Z: " << fYZ.GetN();

            if (fYZ.GetN() < MinPoints || fXZ.GetN() < MinPoints)
            {
                LOG(debug) << "CosmicTracker::Fit : Not enough Points to make reasonable fit.";
                return fTrack;
            }

            filter(fYZ);
            if (fYZ.GetN() < MinPoints)
            {
                LOG(debug) << "CosmicTracker::Fit : Not enough Points to make reasonable fit after "
                              "horizontal filtering.";
                return fTrack;
            }

            filter(fXZ);
            if (fXZ.GetN() < MinPoints)
            {
                LOG(debug) << "CosmicTracker::Fit : Not enough Points to make reasonable fit after vertical "
                              "filtering.";
                return fTrack;
            }

            fXZ.Sort();
            fYZ.Sort();

            TVector3& direction = fTrack.Direction;
            TVector3& entryPoint = fTrack.EntryPoint;
            TVector3& invDirection = fTrack.InvDirection;

            const auto yFit = fit(fYZ);
            if (isnan(yFit[0]))
            {
                LOG(debug) << "CosmicTracker::Fit : Could not get a reasonable vertical fit.";
                return fTrack;
            }

            if (fabs(yFit[0]) < 0.1)
            {
                LOG(debug) << "CosmicTracker::Fit : Y-Slope too small. Better reject this Track.";
                return fTrack;
            }

            const auto xFit = fit(fXZ);
            if (isnan(xFit[0]))
            {
                LOG(debug) << "CosmicTracker::Fit : Could not get a reasonable horizontal fit.";
                return fTrack;
            }

            direction[0] = xFit[0];
            entryPoint[0] = xFit[1];

            direction[1] = yFit[0];
            entryPoint[1] = yFit[1];

            direction[2] = 1.0;
            entryPoint[2] = 0.;

            if (direction[1] > 0)
            {
                // Track is going up -> invert

                entryPoint[2] = fDistances.back() + BarSize_Z; // i.e. back of NeuLAND
                entryPoint[0] += direction[0] * entryPoint[2];
                entryPoint[1] += direction[1] * entryPoint[2];

                for (int i = 0; i < 3; ++i)
                    direction[i] *= -1.0;
            }

            // scale direction to c
            const auto scale = CLight / direction.Mag();
            for (int i = 0; i < 3; ++i)
            {
                direction[i] *= scale;
                invDirection[i] = 1.0 / direction[i];
            }

            // time until cosmic enters NeuLAND
            const auto flightTimeEntry = getCrossPointTime(entryPoint - direction * 2 * TimeEps,
                                                           direction,
                                                           invDirection,
                                                           { -0.5 * BarLength, 0.5 * BarLength },
                                                           { -0.5 * BarLength, 0.5 * BarLength },
                                                           { fDistances[0], fDistances.back() + BarSize_Z });

            if (isnan(flightTimeEntry))
                return fTrack;

            entryPoint += direction * (flightTimeEntry - TimeEps);

            LOG(debug) << "X: " << entryPoint[0] << " + t * " << direction[0];
            LOG(debug) << "Y: " << entryPoint[1] << " + t * " << direction[1];
            LOG(debug) << "Z: " << entryPoint[2] << " + t * " << direction[2];

            fillInteractions(fTrack);
            if (fTrack.Interactions.size() < 3) // ig || fTrack.Interactions.size() < fBarIDs.size() - 3)
            {
                fTrack.Interactions.clear();
                fTrack.TotalTrackLength = 0.;
            }

            LOG(debug) << "Interactions: " << fTrack.Interactions.size();

            return fTrack;
        }

        void CosmicTracker::filter(TGraphErrors& graph) const
        {
            std::array<Int_t, 256> remove;
            auto nRemove = 0;

            auto nPoints = graph.GetN();
            auto x = graph.GetX();
            auto y = graph.GetY();

            // first remove points which are far away from all others

            for (auto p = 0; p < nPoints; ++p)
            {
                Bool_t foundClose = false;
                for (auto op = 0; op < nPoints; ++op)
                {
                    if (p == op)
                        continue;
                    const auto dist2 = Sqr(x[p] - x[op]) + Sqr(y[p] - y[op]);
                    // ig if (dist2 < Sqr(MaxDistance))
                    if (dist2 < Sqr(2 * MaxDistance))
                    {
                        foundClose = true;
                        break;
                    }
                }
                if (!foundClose)
                {
                    remove[nRemove] = p;
                    if (++nRemove == remove.size())
                    {
                        graph.Set(0);
                        return;
                    }
                }
            }

            LOG(debug) << "   Removed : " << nRemove;

            if (nRemove > 0)
            {
                auto offset = 1;
                const auto totalPoints = nPoints - nRemove;
                for (auto p = remove[0]; p < totalPoints; ++p)
                {
                    if (offset < nRemove && p + offset == remove[offset])
                    {
                        ++offset;
                        --p;
                        continue;
                    }

                    x[p] = x[p + offset];
                    y[p] = y[p + offset];
                }

                graph.Set(totalPoints);
            }

            x = graph.GetX();
            y = graph.GetY();
            nPoints = graph.GetN();
            const auto linReg = linearRegression(x, y, nPoints);
            const auto factor = 1. / (Sqr(linReg[1]) + 1.);

            nRemove = 0;

            for (auto p = nPoints - 1; p >= 0; --p)
            {
                const auto dist2 = Sqr(linReg[1] * x[p] - y[p] + linReg[0]) * factor;
                // ig if (dist2 > Sqr(MaxDistance))
                if (dist2 > Sqr(MaxDistance))
                {
                    // seems like there is at least one cluster far away from the fit
                    // better reject this event
                    graph.Set(0);
                    return;
                }

                // ig if (dist2 > Sqr(1.5 * BarSize_XY))
                if (dist2 > Sqr(1.5 * BarSize_XY))
                {
                    // remove points which are a bit of
                    graph.RemovePoint(p);
                    if (++nRemove == 2)
                    {
                        // we removed to many points, better reject this event
                        graph.Set(0);
                        return;
                    }
                }
            }
        }

        void CosmicTracker::Reset()
        {
            fXZ.Set(0);
            fYZ.Set(0);

            fTrack.Interactions.clear();
            fTrack.TotalTrackLength = 0.;
            fTrack.Stopped = false;

            fBarIDs.clear();
        }

        std::array<Double_t, 2> CosmicTracker::fit(TGraphErrors& graph)
        {
            const auto linReg = linearRegression(graph.GetX(), graph.GetY(), graph.GetN());
            if (!isfinite(linReg[0]))
                return { NaN, NaN };

            fFit.SetParameter(0, linReg[0]);
            fFit.SetParameter(1, linReg[1]);

            TFitResultPtr resultptr = graph.Fit(&fFit, "NQS");

            const auto slope = fFit.GetParameter(1);
            const auto yIntercept = fFit.GetParameter(0);

            auto redchi2 = resultptr->Chi2() / resultptr->Ndf();

            // Throw away bad fits
            // ig if (redchi2 < 0.5 || redchi2 > 2. || fabs(slope) > MaxSlope)
            // ig   return { NaN, NaN };

            return { slope, yIntercept };
        }

        Double_t CosmicTracker::getCrossPointTime(const TVector3& point,
                                                  const TVector3& direction,
                                                  const TVector3& invDirection,
                                                  const DPair& xBounds,
                                                  const DPair& yBounds,
                                                  const DPair& zBounds) const
        {
            const auto invert = (WithinBounds(point, xBounds, yBounds, zBounds) ? -1. : 1.);

            TVector3 crossPoint;

            if (invert * invDirection[0] > 0)
            {
                crossPoint[0] = xBounds[0];
            }
            else
            {
                crossPoint[0] = xBounds[1];
            }

            if (invert * invDirection[1] > 0)
            {
                crossPoint[1] = yBounds[0];
            }
            else
            {
                crossPoint[1] = yBounds[1];
            }

            if (invert * invDirection[2] > 0)
            {
                crossPoint[2] = zBounds[0];
            }
            else
            {
                crossPoint[2] = zBounds[1];
            }

            auto wasHit = kFALSE;

            for (int i = 0; i < 3; ++i)
            {
                const auto flightTime = (crossPoint[i] - point[i]) * invDirection[i];
                if (flightTime > 0. &&
                    WithinBounds(point + direction * (flightTime + invert * TimeEps), xBounds, yBounds, zBounds))
                    return flightTime;
            }

            return NaN;
        }

        // Returns { YIntercept, Slope }
        DPair CosmicTracker::linearRegression(const Double_t* x, const Double_t* y, const Int_t points) const
        {
            const auto invPoints = 1. / points;
            const auto invRedPoints = 1. / (points - 1);
            const auto xMean = invPoints * std::accumulate(x, x + points, 0.);
            const auto yMean = std::accumulate(y, y + points, 0.) / points;
            const auto xVar = invRedPoints * std::accumulate(x,
                                                             x + points,
                                                             0.,
                                                             [xMean](const Double_t acc, const Double_t val)
                                                             { return acc + Sqr(val - xMean); });

            if (xVar == 0)
            {
                // we have a vertical line
                return { Inf, xMean };
            }
            const auto yVar = invRedPoints * std::accumulate(y,
                                                             y + points,
                                                             0.,
                                                             [yMean](const Double_t acc, const Double_t val)
                                                             { return acc + Sqr(val - yMean); });

            auto xyVar = 0.;
            for (auto p = 0; p < points; ++p)
                xyVar += (x[p] - xMean) * (y[p] - yMean);
            xyVar *= invRedPoints;

            const auto slope = (yVar - xVar + std::sqrt(Sqr(yVar - xVar) + 4 * Sqr(xyVar))) / (2 * xyVar);
            const auto yIntercept = yMean - slope * xMean;

            if (fabs(slope) > MaxSlope)
                return { Inf, NaN };

            return { yIntercept, slope };
        }

        void CosmicTracker::fillInteractions(R3BNeulandCosmicTrack& track) const
        {
            auto nMissedHits = 0;
            auto lastValidSize = 0;
            auto missedEnergy = 0.;

            std::array<Int_t, 3> barDirections;
            for (auto i = 0; i < 3; ++i)
            {
                barDirections[i] = (track.Direction[i] > 0. ? 1 : -1);
            }

            // first find the plane and bar where the cosmic entered NeuLAND
            Int_t currentPlane = 0, currentBar;
            track.EntryPoint += track.Direction * TimeEps;

            while (currentPlane < fDistances.size() - 1 &&
                   !WithinBounds(track.EntryPoint.Z(), fDistances[currentPlane], fDistances[currentPlane + 1]))
            {
                ++currentPlane;
            }

            // now we now that the cosmic is between currentPlane and currentPlane+1
            // there might be a gap so check if we are not in currentPlane
            if (!WithinBounds(track.EntryPoint.Z(), fDistances[currentPlane], fDistances[currentPlane] + BarSize_Z))
            {
                // Seems like we are in a gap
                if (track.Direction.Z() > 0.)
                {
                    // we will hit the next plane
                    ++currentPlane;
                }

                const auto crossTime =
                    getCrossPointTime(track.EntryPoint,
                                      track.Direction,
                                      track.InvDirection,
                                      { -0.5 * BarLength, 0.5 * BarLength },
                                      { -0.5 * BarLength, 0.5 * BarLength },
                                      { fDistances[currentPlane], fDistances[currentPlane] + BarSize_Z });

                if (isnan(crossTime))
                    return;

                track.EntryPoint += track.Direction * crossTime;
            }

            // we now know the plane

            auto planeChanged = kTRUE;
            DPair xBorder, yBorder;
            auto currentPoint = track.EntryPoint;
            auto tof = 0.;

            LOG(debug) << "Entry Plane: " << currentPlane << " at " << currentPoint.X() << " " << currentPoint.Y()
                       << " " << currentPoint.Z();

            while (kTRUE)
            {
                if (planeChanged)
                {
                    if (IsPlaneHorizontal(currentPlane))
                    {
                        currentBar = (currentPoint.Y() + BarsPerPlane * 0.5 * BarSize_XY) / BarSize_XY;
                        xBorder = { -0.5 * BarLength, 0.5 * BarLength };
                        yBorder = { -0.5 * BarLength + currentBar * BarSize_XY,
                                    -0.5 * BarLength + (currentBar + 1) * BarSize_XY };
                    }
                    else
                    {
                        currentBar = (currentPoint.X() + BarsPerPlane * 0.5 * BarSize_XY) / BarSize_XY;
                        xBorder = { -0.5 * BarLength + currentBar * BarSize_XY,
                                    -0.5 * BarLength + (currentBar + 1) * BarSize_XY };
                        yBorder = { -0.5 * BarLength, 0.5 * BarLength };
                    }
                }
                else
                {
                    if (IsPlaneHorizontal(currentPlane))
                    {
                        currentBar += barDirections[1];
                        yBorder[0] += barDirections[1] * BarSize_XY;
                        yBorder[1] += barDirections[1] * BarSize_XY;
                    }
                    else
                    {
                        currentBar += barDirections[0];
                        xBorder[0] += barDirections[0] * BarSize_XY;
                        xBorder[1] += barDirections[0] * BarSize_XY;
                    }
                }
                const auto crossTime =
                    getCrossPointTime(currentPoint,
                                      track.Direction,
                                      track.InvDirection,
                                      xBorder,
                                      yBorder,
                                      { fDistances[currentPlane], fDistances[currentPlane] + BarSize_Z });

                if (isnan(crossTime))
                    break;

                const auto barID = currentBar + currentPlane * BarsPerPlane;
                const auto trackLength = crossTime * CLight;
                if (IsPlaneHorizontal(currentPlane))
                {
                    track.Interactions.emplace_back(barID,
                                                    tof,
                                                    tof + crossTime,
                                                    currentPoint.X(),
                                                    currentPoint.X() + crossTime * track.Direction.X(),
                                                    trackLength * MIPStoppingPower,
                                                    trackLength);
                }
                else
                {
                    track.Interactions.emplace_back(barID,
                                                    tof,
                                                    tof + crossTime,
                                                    currentPoint.Y(),
                                                    currentPoint.Y() + crossTime * track.Direction.Y(),
                                                    trackLength * MIPStoppingPower,
                                                    trackLength);
                }

                if (std::find(fBarIDs.begin(), fBarIDs.end(), barID) != fBarIDs.end())
                {
                    missedEnergy = 0.;
                    nMissedHits = 0;
                    lastValidSize = track.Interactions.size();
                }
                else
                {
                    ++nMissedHits;
                    missedEnergy += trackLength * MIPStoppingPower;
                    if (missedEnergy > 10. || (missedEnergy > 7. && nMissedHits > 1))
                    {
                        if (lastValidSize == 0)
                        {
                            return;
                        }
                        track.Interactions.resize(lastValidSize);
                        track.TotalTrackLength = track.Interactions.back().ExitTime * CLight;
                        track.Stopped = kTRUE;
                        return;
                    }
                }

                tof += crossTime + TimeEps;
                currentPoint = track.EntryPoint + track.Direction * tof;

                // check if we are still in NeuLAND
                if (!WithinBounds(currentPoint,
                                  { -0.5 * BarLength, 0.5 * BarLength },
                                  { -0.5 * BarLength, 0.5 * BarLength },
                                  { fDistances[0], fDistances.back() + BarSize_Z }))
                {
                    LOG(debug) << "Left NeuLAND at Plane " << currentPlane << " with " << currentPoint.X() << " "
                               << currentPoint.Y() << " " << currentPoint.Z();
                    break;
                }

                // check if we left the plane
                if (!WithinBounds(currentPoint.Z(), fDistances[currentPlane], fDistances[currentPlane] + BarSize_Z))
                {
                    // we left the plane
                    planeChanged = true;
                    currentPlane += barDirections[2];

                    // check if we are already in the next plane
                    if (!WithinBounds(currentPoint.Z(), fDistances[currentPlane], fDistances[currentPlane] + BarSize_Z))
                    {
                        // propagate to next plane
                        const auto gapTime =
                            getCrossPointTime(currentPoint,
                                              track.Direction,
                                              track.InvDirection,
                                              { -0.5 * BarLength, 0.5 * BarLength },
                                              { -0.5 * BarLength, 0.5 * BarLength },
                                              { fDistances[currentPlane], fDistances[currentPlane] + BarSize_Z });

                        if (isnan(gapTime))
                            break;

                        tof += gapTime + TimeEps;
                        currentPoint = track.EntryPoint + track.Direction * tof;
                    }
                }
            }
            track.TotalTrackLength = tof * CLight;
            return;
        }
    } // namespace Calibration
} // namespace Neuland
