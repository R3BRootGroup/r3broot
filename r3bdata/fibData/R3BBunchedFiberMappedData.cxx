/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019 Members of R3B Collaboration                          *
 *                                                                            *
 *             This software is distributed under the terms of the            *
 *                 GNU General Public Licence (GPL) version 3,                *
 *                    copied verbatim in the file "LICENSE".                  *
 *                                                                            *
 * In applying this license GSI does not waive the privileges and immunities  *
 * granted to it by virtue of its status as an Intergovernmental Organization *
 * or submit itself to any jurisdiction.                                      *
 ******************************************************************************/

// ----------------------------------------------------------------
// -----              R3BBunchedFiberMappedData                -----
// -----             Created Jan 2018 by M.Heil        -----
// ----------------------------------------------------------------

#include "R3BBunchedFiberMappedData.h"

R3BBunchedFiberMappedData::R3BBunchedFiberMappedData()
  : fIsMAPMT()
  , fChannel(-1)
  , fIsLeading()
  , fCoarse(-1)
  , fFine(-1)
{
}

R3BBunchedFiberMappedData::R3BBunchedFiberMappedData(Bool_t a_is_mapmt, Int_t
    a_channel, Bool_t a_is_leading, Int_t a_coarse, Int_t a_fine)
  : fIsMAPMT(a_is_mapmt)
  , fChannel(a_channel)
  , fIsLeading(a_is_leading)
  , fCoarse(a_coarse)
  , fFine(a_fine)
{
}

R3BBunchedFiberMappedData::~R3BBunchedFiberMappedData()
{
}

Int_t R3BBunchedFiberMappedData::GetChannel() const
{
  return fChannel;
}

Int_t R3BBunchedFiberMappedData::GetCoarse() const
{
  return fCoarse;
}

Int_t R3BBunchedFiberMappedData::GetFine() const
{
  return fFine;
}

Bool_t R3BBunchedFiberMappedData::IsMAPMT() const
{
  return fIsMAPMT;
}

Bool_t R3BBunchedFiberMappedData::IsSPMT() const
{
  return !fIsMAPMT;
}

Bool_t R3BBunchedFiberMappedData::IsLeading() const
{
  return fIsLeading;
}

Bool_t R3BBunchedFiberMappedData::IsTrailing() const
{
  return !fIsLeading;
}

ClassImp(R3BBunchedFiberMappedData)
