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

#include "R3BPaddleCalData.h"

/* for the data analysis of the paddle-type detectors.
 * Introduced by Ralf, Apr 2016
 *
 *
 */

R3BPaddleCalData::R3BPaddleCalData()
    : fPlane(0)
    , fBar(0)
{
}

R3BPaddleCalData::R3BPaddleCalData(UInt_t plane, UInt_t bar)
    : fPlane(plane)
    , fBar(bar)
    , fTime1L_ns(0.0 / 0.0) // 0/0 produces a NAN. Macro NAN not available??
    , fTime1T_ns(0.0 / 0.0)
    , fTime2L_ns(0.0 / 0.0)
    , fTime2T_ns(0.0 / 0.0)
{
}

ClassImp(R3BPaddleCalData)
