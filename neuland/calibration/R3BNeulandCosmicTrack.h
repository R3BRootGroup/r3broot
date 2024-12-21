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

#ifndef R3BNEULAND_COSMICTRACK_H
#define R3BNEULAND_COSMICTRACK_H

#include "TVector3.h"

#include "R3BNeulandCommon.h"

struct R3BNeulandCosmicTrack
{
    struct Interaction
    {
        Interaction(Int_t barID = -1,
                    Double_t entryTime = 0.,
                    Double_t exitTime = 0.,
                    Double_t entryPosition = 0.,
                    Double_t exitPosition = 0.,
                    Double_t energy = 0.,
                    Double_t trackLength = 0.)
            : BarID(barID)
            , EntryTime(entryTime)
            , ExitTime(exitTime)
            , EntryPosition(entryPosition)
            , ExitPosition(exitPosition)
            , Energy(energy)
            , TrackLength(trackLength)
        {
        }
        Int_t BarID;
        Double_t EntryTime;
        Double_t ExitTime;
        Double_t EntryPosition;
        Double_t ExitPosition;
        Double_t Energy;
        Double_t TrackLength;
    };

    std::vector<Interaction> Interactions;
    TVector3 EntryPoint;
    TVector3 Direction;
    TVector3 InvDirection;
    Double_t TotalTrackLength;
    Bool_t Stopped = kFALSE;
    // Double_t TimeInsideMaterial;
};
#endif
