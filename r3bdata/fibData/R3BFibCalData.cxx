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

#include "R3BFibCalData.h"

/* for the data analysis of the fiber-type detectors.
 * Introduced by M.Heil, Jan 2018
 *
 *
 */

R3BFibCalData::R3BFibCalData()
    : fPlane(0)
    , fFiber(0)
{
}

R3BFibCalData::R3BFibCalData(UInt_t plane, UInt_t fiber)
    : fPlane(plane)
    , fFiber(fiber)
    , fTime1L_ns(0.0 / 0.0) // 0/0 produces a NAN. Macro NAN not available??
    , fTime1T_ns(0.0 / 0.0)
    , fTime2L_ns(0.0 / 0.0)
    , fTime2T_ns(0.0 / 0.0)
{
}

ClassImp(R3BFibCalData)
