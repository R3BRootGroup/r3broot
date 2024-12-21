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

#include <iomanip>

#include "TH1F.h"

#include "FairLogger.h"

#include "R3BTofUnpack.h"

R3BTofUnpack::R3BTofUnpack(Short_t type, Short_t subType, Short_t procId, Short_t subCrate, Short_t control)
    : FairUnpack(type, subType, procId, subCrate, control)
{
}

R3BTofUnpack::~R3BTofUnpack() {}

Bool_t R3BTofUnpack::Init() { return kTRUE; }

Bool_t R3BTofUnpack::DoUnpack(Int_t* data, Int_t size)
{
    LOG(info) << "R3BTofUnpack : Unpacking " << size << " double words of data ...";

    // Clear data map
    fmap_leading_t.clear();
    fmap_leading_c.clear();

    // Produce hex output
    LOG(info) << "!!!!!!! " << std::hex;
    for (Int_t i = 0; i < size; i++)
    {
        LOG(info) << data[i] << " ";
    }
    LOG(info) << std::dec;

    // Loop over double words of input data
    for (Int_t i = 1; i < size; i++)
    {
        UInt_t* p1 = (UInt_t*)(data + i);

        if ((p1[0] & 0xff000000) == 0x4c000000)
        {
            break;
        }

        if ((p1[0] & 0xff000000) != 0xab000000)
        {
            LOG(error) << "Wrong TFW sub-event header";
            return kFALSE;
        }

        UInt_t l_i = 0;

        UInt_t nr_cha = 0x1ff & (p1[l_i] >> 9);
        LOG(info) << "R3BTofUnpack : Number of raw channels to read: " << nr_cha;

        if (nr_cha)
        {
            l_i += 1;
            nr_cha -= 1;
        }

        if (nr_cha)
        {
            l_i += 1;
            nr_cha -= 1;
            while (l_i < size)
            {
                UInt_t raw_ch = (p1[l_i] & 0x7e000000) >> 25; // raw channel 0...31
                UInt_t phys_ch = raw_ch / 2;                  // physical channel 0...15
                UInt_t t = (p1[l_i] & 0x7ff);                 // tdc time
                UInt_t c = (p1[l_i] & 0xfff800) >> 11;        // clock

                // Store leading and trailing edge
                if (0 == (raw_ch % 2))
                {
                    if (0 == fmap_leading_t[phys_ch])
                    {
                        fmap_leading_t[phys_ch] = t;
                        fmap_leading_c[phys_ch] = c;
                    }
                }
                else
                {
                    if (0 == fmap_trailing_t[phys_ch])
                    {
                        fmap_trailing_t[phys_ch] = t;
                        fmap_trailing_c[phys_ch] = c;
                    }
                }
                if (0 == nr_cha)
                {
                    break;
                }
                l_i += 1;
                nr_cha -= 1;
            }
        }

        i += l_i;
    }

    // Set the data
    Int_t nhits = 0;
    Int_t t_leading, c_leading;
    Int_t t_trailing, c_trailing;
    for (Int_t i = 0; i < 16; i++)
    {
        if (0 != fmap_leading_t[i])
        {
            nhits += 1;
            t_leading = fmap_leading_t[i];
            c_leading = fmap_leading_c[i];
            t_trailing = fmap_trailing_t[i];
            c_trailing = fmap_trailing_c[i];
            LOG(info) << "######### ch:" << i << ", Leading: tdc=" << t_leading << ", clock=" << c_leading
                      << ";   Trailing: tdc=" << t_trailing << ",  clock=" << c_trailing;
        }
    }

    LOG(info) << "R3BTofUnpack : Number of hits in TFW: " << nhits;

    return kTRUE;
}

void R3BTofUnpack::Reset() {}

void R3BTofUnpack::Register() {}

ClassImp(R3BTofUnpack)
