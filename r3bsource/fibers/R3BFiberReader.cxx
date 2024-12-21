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

#include "R3BFiberReader.h"
#include "R3BFiberMappedData.h"

#include <FairLogger.h>
#include <FairRootManager.h>
#include <TClonesArray.h>

/**
 ** ext_h101_fib30.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:FIBEO,id=h101_FIBEO,NOTRIGEVENTNO,ext_h101_fib30.h
 **
 ** ext_h101_fib31.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:FIBEI,id=h101_FIBEI,NOTRIGEVENTNO,ext_h101_fib31.h
 **
 ** ext_h101_fib32.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:FIBEZ,id=h101_FIBEZ,NOTRIGEVENTNO,ext_h101_fib32.h
 **
 ** ext_h101_fib33.h was created by running
 ** $unpacker --ntuple=STRUCT_HH,RAW:FIBEE,id=h101_FIBEE,NOTRIGEVENTNO,ext_h101_fib33.h
 **/

extern "C"
{
#include "ext_data_client.h"
#include "ext_h101_fib.h"
#include "ext_h101_fib23a.h"
#include "ext_h101_fib23b.h"
#include "ext_h101_fib30.h"
#include "ext_h101_fib31.h"
#include "ext_h101_fib32.h"
#include "ext_h101_fib33.h"
#include "ext_h101_fibeight.h"
#include "ext_h101_fibeleven.h"
#include "ext_h101_fibseven.h"
#include "ext_h101_fibten.h"
#include "ext_h101_fibthirteen.h"
#include "ext_h101_fibthreea.h"
#include "ext_h101_fibthreeb.h"
#include "ext_h101_fibtwelve.h"
}

R3BFiberReader::R3BFiberReader(char const* a_name, UInt_t fiber_num, EXT_STR_h101_FIBZEA_onion* data, size_t offset)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData23a(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(fiber_num)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
}

R3BFiberReader::R3BFiberReader(char const* a_name, UInt_t fiber_num, EXT_STR_h101_FIBZEB_onion* data, size_t offset)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData23b(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(fiber_num)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
}

R3BFiberReader::R3BFiberReader(char const* a_name,
                               EXT_STR_h101_FIBTHREEA_onion* data,
                               size_t offset,
                               UInt_t a_sub_num,
                               UInt_t a_mapmt_channel_num,
                               UInt_t a_spmt_channel_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData3a(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(0)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
    fChannelNum[0] = a_sub_num * a_mapmt_channel_num;
    fChannelNum[1] = a_sub_num * a_spmt_channel_num;
}

R3BFiberReader::R3BFiberReader(char const* a_name,
                               EXT_STR_h101_FIBTHREEB_onion* data,
                               size_t offset,
                               UInt_t a_sub_num,
                               UInt_t a_mapmt_channel_num,
                               UInt_t a_spmt_channel_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData3b(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(0)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
    fChannelNum[0] = a_sub_num * a_mapmt_channel_num;
    fChannelNum[1] = a_sub_num * a_spmt_channel_num;
}

R3BFiberReader::R3BFiberReader(char const* a_name,
                               EXT_STR_h101_FIBSEVEN_onion* data,
                               size_t offset,
                               UInt_t a_sub_num,
                               UInt_t a_mapmt_channel_num,
                               UInt_t a_spmt_channel_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData7(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(0)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
    fChannelNum[0] = a_sub_num * a_mapmt_channel_num;
    fChannelNum[1] = a_sub_num * a_spmt_channel_num;
}

R3BFiberReader::R3BFiberReader(char const* a_name,
                               EXT_STR_h101_FIBEIGHT_onion* data,
                               size_t offset,
                               UInt_t a_sub_num,
                               UInt_t a_mapmt_channel_num,
                               UInt_t a_spmt_channel_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData8(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(0)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
    fChannelNum[0] = a_sub_num * a_mapmt_channel_num;
    fChannelNum[1] = a_sub_num * a_spmt_channel_num;
}

R3BFiberReader::R3BFiberReader(char const* a_name,
                               EXT_STR_h101_FIBTEN_onion* data,
                               size_t offset,
                               UInt_t a_sub_num,
                               UInt_t a_mapmt_channel_num,
                               UInt_t a_spmt_channel_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData10(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(0)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
    fChannelNum[0] = a_sub_num * a_mapmt_channel_num;
    fChannelNum[1] = a_sub_num * a_spmt_channel_num;
}

R3BFiberReader::R3BFiberReader(char const* a_name,
                               EXT_STR_h101_FIBELEVEN_onion* data,
                               size_t offset,
                               UInt_t a_sub_num,
                               UInt_t a_mapmt_channel_num,
                               UInt_t a_spmt_channel_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData11(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(0)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
    fChannelNum[0] = a_sub_num * a_mapmt_channel_num;
    fChannelNum[1] = a_sub_num * a_spmt_channel_num;
}

R3BFiberReader::R3BFiberReader(char const* a_name,
                               EXT_STR_h101_FIBTWELVE_onion* data,
                               size_t offset,
                               UInt_t a_sub_num,
                               UInt_t a_mapmt_channel_num,
                               UInt_t a_spmt_channel_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData12(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(0)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
    fChannelNum[0] = a_sub_num * a_mapmt_channel_num;
    fChannelNum[1] = a_sub_num * a_spmt_channel_num;
}

R3BFiberReader::R3BFiberReader(char const* a_name,
                               EXT_STR_h101_FIBTHIRTEEN_onion* data,
                               size_t offset,
                               UInt_t a_sub_num,
                               UInt_t a_mapmt_channel_num,
                               UInt_t a_spmt_channel_num)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData13(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(0)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
    fChannelNum[0] = a_sub_num * a_mapmt_channel_num;
    fChannelNum[1] = a_sub_num * a_spmt_channel_num;
}

R3BFiberReader::R3BFiberReader(char const* a_name, UInt_t fiber_num, EXT_STR_h101_FIBEO_onion* data, size_t offset)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData30(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(fiber_num)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
}

R3BFiberReader::R3BFiberReader(char const* a_name, UInt_t fiber_num, EXT_STR_h101_FIBEI_onion* data, size_t offset)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData31(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(fiber_num)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
}

R3BFiberReader::R3BFiberReader(char const* a_name, UInt_t fiber_num, EXT_STR_h101_FIBEZ_onion* data, size_t offset)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData32(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(fiber_num)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
}

R3BFiberReader::R3BFiberReader(char const* a_name, UInt_t fiber_num, EXT_STR_h101_FIBEE_onion* data, size_t offset)
    : R3BReader(TString("R3B") + a_name + "Reader")
    , fData33(data)
    , fOffset(offset)
    , fShortName(a_name)
    , fFiberNum(fiber_num)
    , fDataSPMTTrig(NULL)
    , fMappedArray(new TClonesArray("R3BFiberMappedData"))
    , fOnline(kFALSE)
{
}

void R3BFiberReader::SetUcesbSPMTTrigStructure(EXT_STR_h101_FIB_onion* data, size_t offset)
{
    fDataSPMTTrig = data;
    fOffsetspmttrig = offset;
    return;
}

R3BFiberReader::~R3BFiberReader()
{
    R3BLOG(debug1, "");
    if (fMappedArray)
        delete fMappedArray;
}

Bool_t R3BFiberReader::Init(ext_data_struct_info* a_struct_info)
{
    R3BLOG(info, "For fiber " << fShortName);
    if (fShortName == "Fi23a")
    {
        R3B_MAPMT_INIT_TRIG(FIBZEA, fData23a);
    }
    else if (fShortName == "Fi23b")
    {
        R3B_MAPMT_INIT_TRIG(FIBZEB, fData23b);
    }
    else if (fShortName == "Fi3a")
    {
        R3BBUNCHED_FIBER_INIT_MAPMT_TRIG(FIBTHREEA, fData3a);
    }
    else if (fShortName == "Fi3b")
    {
        R3BBUNCHED_FIBER_INIT_MAPMT_TRIG(FIBTHREEB, fData3b);
    }
    else if (fShortName == "Fi7")
    {
        R3BBUNCHED_FIBER_INIT(FIBSEVEN, fData7);
    }
    else if (fShortName == "Fi8")
    {
        R3BBUNCHED_FIBER_INIT(FIBEIGHT, fData8);
    }
    else if (fShortName == "Fi10")
    {
        R3BBUNCHED_FIBER_INIT_MAPMT_TRIG(FIBTEN, fData10);
    }
    else if (fShortName == "Fi11")
    {
        R3BBUNCHED_FIBER_INIT_MAPMT_TRIG(FIBELEVEN, fData11);
    }
    else if (fShortName == "Fi12")
    {
        R3BBUNCHED_FIBER_INIT_MAPMT_TRIG(FIBTWELVE, fData12);
    }
    else if (fShortName == "Fi13")
    {
        R3BBUNCHED_FIBER_INIT_MAPMT_TRIG(FIBTHIRTEEN, fData13);
    }
    else if (fShortName == "Fi30")
    {
        R3B_MAPMT_INIT_TRIG(FIBEO, fData30);
    }
    else if (fShortName == "Fi31")
    {
        R3B_MAPMT_INIT_TRIG(FIBEI, fData31);
    }
    else if (fShortName == "Fi32")
    {
        R3B_MAPMT_INIT_TRIG(FIBEZ, fData32);
    }
    else if (fShortName == "Fi33")
    {
        R3B_MAPMT_INIT_TRIG(FIBEE, fData33);
    }

    if (fShortName == "Fi3a" || fShortName == "Fi3b" || fShortName == "Fi7" || fShortName == "Fi8" ||
        fShortName == "Fi10" || fShortName == "Fi11" || fShortName == "Fi12" || fShortName == "Fi13")
        for (size_t side_i = 0; side_i < 2; ++side_i)
        {
            for (size_t edge_i = 0; edge_i < 2; ++edge_i)
            {
                for (size_t prec_i = 0; prec_i < 2; ++prec_i)
                {
                    auto const& ch = fMHL[side_i][edge_i][prec_i];
                    if (fChannelNum[side_i] != ch._MI_len)
                    {
                        LOG(fatal) << "Multi-hit array sizes mismatch (fChannelNum[" << side_i
                                   << "]=" << fChannelNum[side_i] << " != MI-len=" << ch._MI_len << ").";
                        return kFALSE;
                    }
                }
            }
            if (0 == fChannelNum[1])
            {
                break;
            }
        }

    if (fShortName == "Fi23a" || fShortName == "Fi23b" || fShortName == "Fi30" || fShortName == "Fi31" ||
        fShortName == "Fi32" || fShortName == "Fi33")
        for (size_t side_i = 0; side_i < 2; ++side_i)
        {
            for (size_t edge_i = 0; edge_i < 2; ++edge_i)
            {
                for (size_t prec_i = 0; prec_i < 2; ++prec_i)
                {
                    auto const& ch = fMHL[side_i][edge_i][prec_i];
                    if (fFiberNum != ch._MI_len)
                    {
                        R3BLOG(fatal,
                               "Multi-hit array sizes mismatch (fFiberNum=" << fFiberNum << " != MI-len=" << ch._MI_len
                                                                            << ").");
                        return kFALSE;
                    }
                }
            }
        }

    if ((fShortName == "Fi10" || fShortName == "Fi11" || fShortName == "Fi12" || fShortName == "Fi13") && fDataSPMTTrig)
    {
        Int_t okt;
        EXT_STR_h101_FIB_ITEMS_INFO(okt, *a_struct_info, fOffsetspmttrig, EXT_STR_h101_FIB, 0);
        if (!okt)
        {
            R3BLOG(fatal, "Failed to setup structure information of EXT_STR_h101_FIB_ITEMS_INFO.");
            return kFALSE;
        }
        else
            R3BLOG(info, "SPMT trigger times for " << fShortName);
    }

    // Register of fiber mapped data
    FairRootManager::Instance()->Register(fShortName + "Mapped", fShortName + " mapped data", fMappedArray, !fOnline);

    return kTRUE;
}

Bool_t R3BFiberReader::R3BRead()
{
    R3BLOG(debug1, "Read BEGIN for fib " << fShortName);
    for (size_t side_i = 0; side_i < 2; ++side_i)
    {
        for (size_t edge_i = 0; edge_i < 2; ++edge_i)
        {
            auto const& e = fMHL[side_i][edge_i];

            // Check that coarse and fine list are the same size.
            uint32_t c_M = *e[0]._M;
            uint32_t f_M = *e[1]._M;
            uint32_t c_ = *e[0]._;
            uint32_t f_ = *e[1]._;

            if (c_M != f_M || c_ != f_)
            {
                R3BLOG(warn,
                       "Coarse and fine multi-hit list counts mismatch "
                       "(edge="
                           << edge_i << ";M{c=" << c_M << ",f=" << f_M << "};_{c=" << c_ << ",f=" << f_ << "}).");
                return kFALSE;
            }
            if (c_M > e[0]._MI_len || c_M > e[0]._ME_len || c_ > e[0]._v_len)
            {
                R3BLOG(warn,
                       "Multi-hit indexing out of range "
                       "(edge="
                           << edge_i << ";(M=" << c_M << ")>=(MI=" << e[0]._MI_len << ",ME=" << e[0]._ME_len
                           << ");(_=" << c_ << ")>=(v=" << e[0]._v_len << ")).");
                return kFALSE;
            }

            // Simply dump edges, can only be sorted perfectly after calibration,
            // and hits should be somewhat sorted already, that helps many sorting
            // algos.
            uint32_t cur_entry = 0;
            for (uint32_t i = 0; i < c_M; i++)
            {
                uint32_t c_MI = e[0]._MI[i];
                uint32_t f_MI = e[1]._MI[i];
                uint32_t c_ME = e[0]._ME[i];
                uint32_t f_ME = e[1]._ME[i];

                if (c_MI != f_MI || c_ME != f_ME)
                {
                    R3BLOG(warn,
                           "Coarse and fine multi-hit data mismatch "
                           "(edge="
                               << edge_i << ";MI{c=" << c_MI << ",f=" << f_MI << "};ME{c=" << c_ME << ",f=" << f_ME
                               << "}).");
                    return kFALSE;
                }
                for (; cur_entry < c_ME; cur_entry++)
                {
                    new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
                        R3BFiberMappedData(side_i + 1, c_MI, 0 == edge_i, e[0]._v[cur_entry], e[1]._v[cur_entry]);
                }
            }
        }
        if (0 == fChannelNum[1] &&
            (fShortName == "Fi3a" || fShortName == "Fi3b" || fShortName == "Fi7" || fShortName == "Fi8" ||
             fShortName == "Fi10" || fShortName == "Fi11" || fShortName == "Fi12" || fShortName == "Fi13"))
        {
            break;
        }
    }

    // MAPMT trigger times.
    {
        auto const& e = fMHL[2][0];

        // Check that coarse and fine list are the same size.
        uint32_t c_M = *e[0]._M;
        uint32_t f_M = *e[1]._M;
        uint32_t c_ = *e[0]._;
        uint32_t f_ = *e[1]._;

        if (c_ != f_)
        {
            R3BLOG(warn, "Coarse and fine single-hit list counts mismatch (_{c=" << c_ << ",f=" << f_ << "}).");
            return kFALSE;
        }

        // Simply dump edges, can only be sorted perfectly after calibration,
        // and hits should be somewhat sorted already, that helps many sorting
        // algos.
        uint32_t cur_entry = 0;
        for (uint32_t i = 0; i < c_M; i++)
        {
            uint32_t c_MI = e[0]._MI[i];
            uint32_t f_MI = e[1]._MI[i];
            uint32_t c_ME = e[0]._ME[i];
            uint32_t f_ME = e[1]._ME[i];
            // Only take the first hit per trigger channel, seems like it can
            // have multi-hits.
            if (cur_entry < c_ME)
            {
                new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
                    R3BFiberMappedData(3, c_MI, true, e[0]._v[cur_entry], e[1]._v[cur_entry]);
            }
            cur_entry = c_ME;
        }
    }

    if ((fShortName == "Fi10" || fShortName == "Fi11" || fShortName == "Fi12" || fShortName == "Fi13") && fDataSPMTTrig)
    {
        // Leading TAMEX trigger times.
        auto numChannels = fDataSPMTTrig->FIB_TRIGSLF;
        for (uint32_t i = 0; i < numChannels; i++)
        {
            uint32_t channel = fDataSPMTTrig->FIB_TRIGSLFI[i];
            new ((*fMappedArray)[fMappedArray->GetEntriesFast()])
                R3BFiberMappedData(4, channel, true, fDataSPMTTrig->FIB_TRIGSLCv[i], fDataSPMTTrig->FIB_TRIGSLFv[i]);
        }
    }

    return kTRUE;
}

void R3BFiberReader::Reset() { fMappedArray->Clear(); }

ClassImp(R3BFiberReader)
