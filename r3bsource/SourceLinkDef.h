// clang-format off

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

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class R3BUcesbSource+;
#pragma link C++ class R3BReader+;
#pragma link C++ class R3BUnpackReader+;
#pragma link C++ class R3BWhiterabbitNeulandReader+;
#pragma link C++ class R3BWhiterabbitMasterReader+;
#pragma link C++ class R3BWhiterabbitCalifaReader+;
#pragma link C++ class R3BWhiterabbitAmsReader+;
#pragma link C++ class R3BWhiterabbitFootReader+;
#pragma link C++ class R3BWhiterabbitS2Reader+;
#pragma link C++ class R3BWhiterabbitS8Reader+;
#pragma link C++ class R3BWhiterabbitPspReader+;
#pragma link C++ class R3BWhiterabbitLosReader+;
#pragma link C++ class R3BWhiterabbitPspReader+;
#pragma link C++ class R3BWhiterabbitMusicReader+;
#pragma link C++ class R3BTrloiiTpatReader+;
#pragma link C++ class R3BTrloiiSampReader+;
#pragma link C++ class R3BTrloiiScalerReader+;
#pragma link C++ class R3BTrloiiScalerReader_s467+;
#pragma link C++ class R3BTimestampMasterReader+;
#pragma link C++ class R3BBeamMonitorReader+;
#pragma link C++ class R3BS515PspxReader+;
#pragma link C++ class R3BBunchedFiberReader+;
#pragma link C++ class R3BBunchedFiberSPMTTrigReader+;
#pragma link C++ class R3BSfibReader+;
#pragma link C++ class R3BFi0Reader+;
#pragma link C++ class R3BFi1aReader+;
#pragma link C++ class R3BFi1bReader+;
#pragma link C++ class R3BFi3aReader+;
#pragma link C++ class R3BFi3bReader+;
#pragma link C++ class R3BFi4Reader+;
#pragma link C++ class R3BFi5Reader+;
#pragma link C++ class R3BFi6Reader+;
#pragma link C++ class R3BFi7Reader+;
#pragma link C++ class R3BFi8Reader+;
#pragma link C++ class R3BFiberReader+;
#pragma link C++ class R3BLosReader+;
#pragma link C++ class R3BRoluReader+;
#pragma link C++ class R3BRoluReaderOnline+;
#pragma link C++ class R3BFrsSciReader+;
#pragma link C++ class R3BSci2Reader+;
#pragma link C++ class R3BSci8Reader+;
#pragma link C++ class R3BTofdReader+;
#pragma link C++ class R3BTofiReader+;
#pragma link C++ class R3BPdcReader+;
#pragma link C++ class R3BCalifaFebexReader+;
#pragma link C++ class R3BPtofReader+;
#pragma link C++ class R3BNeulandTamexReader+;
#pragma link C++ class R3BNeulandTacquilaReader+;
#pragma link C++ class R3BStrawtubesReader+;
#pragma link C++ class R3BAmsReader+;
#pragma link C++ class R3BMusicReader+;
#pragma link C++ class R3BFootReader+;
#pragma link C++ class R3BFootSiReader+;
#pragma link C++ class R3BAlpideReader+;
#pragma link C++ class R3BCalifaJulichReader+;
#pragma link C++ class R3BRpcReader+;
#pragma link C++ class R3BTwimReader+;
#pragma link C++ class R3BTttxReader+;
#pragma link C++ class R3BMusliReader+;
#pragma link C++ class R3BMwpcReader+;
#pragma link C++ class R3BSyncCheckReader+;
#pragma link C++ class R3BMosaicReader+;
#pragma link C++ class R3BActafReader+;


#pragma link C++ class EXT_STR_h101_unpack_t;
#pragma link C++ class EXT_STR_h101_whiterabbit_onion_t;
#pragma link C++ class EXT_STR_h101_TPAT_t;
#pragma link C++ class EXT_STR_h101_timestamp_master_t;
#pragma link C++ class EXT_STR_h101_BMON_onion_t;
#pragma link C++ class EXT_STR_h101_FIB_onion_t;
#pragma link C++ class EXT_STR_h101_FIBZERO_onion_t;
#pragma link C++ class EXT_STR_h101_FIBONEA_onion_t;
#pragma link C++ class EXT_STR_h101_FIBONEB_onion_t;
#pragma link C++ class EXT_STR_h101_FIBTHREEA_onion_t;
#pragma link C++ class EXT_STR_h101_FIBTHREEB_onion_t;
#pragma link C++ class EXT_STR_h101_FIBFOUR_onion_t;
#pragma link C++ class EXT_STR_h101_FIBFIVE_onion_t;
#pragma link C++ class EXT_STR_h101_FIBSIX_onion_t;
#pragma link C++ class EXT_STR_h101_FIBSEVEN_onion_t;
#pragma link C++ class EXT_STR_h101_FIBEIGHT_onion_t;
#pragma link C++ class EXT_STR_h101_FIBTEN_onion_t;
#pragma link C++ class EXT_STR_h101_FIBELEVEN_onion_t;
#pragma link C++ class EXT_STR_h101_FIBTWELVE_onion_t;
#pragma link C++ class EXT_STR_h101_FIBTHIRTEEN_onion_t;
#pragma link C++ class EXT_STR_h101_FIBZEA_onion_t;
#pragma link C++ class EXT_STR_h101_FIBZEB_onion_t;
#pragma link C++ class EXT_STR_h101_FIBEO_onion_t;
#pragma link C++ class EXT_STR_h101_FIBEI_onion_t;
#pragma link C++ class EXT_STR_h101_FIBEZ_onion_t;
#pragma link C++ class EXT_STR_h101_FIBEE_onion_t;
#pragma link C++ class EXT_STR_h101_SCI2_onion_t;
#pragma link C++ class EXT_STR_h101_SCI8_onion_t;
#pragma link C++ class EXT_STR_h101_ROLU_onion_t;
#pragma link C++ class EXT_STR_h101_ROLU_NOTRIG_onion_t;
#pragma link C++ class EXT_STR_h101_TOFD_onion_t;
#pragma link C++ class EXT_STR_h101_TOFI_onion_t;
#pragma link C++ class EXT_STR_h101_PTOF_onion_t;
#pragma link C++ class EXT_STR_h101_PSP_onion_t;
#pragma link C++ class EXT_STR_h101_CALIFA_onion_t;
#pragma link C++ class EXT_STR_h101_WRNEULAND_onion_t;
#pragma link C++ class EXT_STR_h101_raw_nnp_tacquila_onion_t;
#pragma link C++ class EXT_STR_h101_raw_nnp_tamex_onion_t;
#pragma link C++ class EXT_STR_h101_AMS_onion_t;
#pragma link C++ class EXT_STR_h101_FOOT_onion_t;
#pragma link C++ class EXT_STR_h101_ALPIDE_onion_t;
#pragma link C++ class EXT_STR_h101_WRMASTER_onion_t;
#pragma link C++ class EXT_STR_h101_WRCALIFA_onion_t;
#pragma link C++ class EXT_STR_h101_WRAMS_onion_t;
#pragma link C++ class EXT_STR_h101_WRFOOT_onion_t;
#pragma link C++ class EXT_STR_h101_WRS2_onion_t;
#pragma link C++ class EXT_STR_h101_WRS8_onion_t;
#pragma link C++ class EXT_STR_h101_WRLOS_onion_t;
#pragma link C++ class EXT_STR_h101_WRMUSIC_onion_t;
#pragma link C++ class EXT_STR_h101_TRLO_onion_t;
#pragma link C++ class EXT_STR_h101_TIMESTAMP_PSPX_onion_t;
#pragma link C++ class EXT_STR_h101_MUSIC_onion_t;
#pragma link C++ class EXT_STR_h101_SOFTWIM_onion_t;
#pragma link C++ class EXT_STR_h101_TTTX_onion_t;
#pragma link C++ class EXT_STR_h101_SOFMWPC_onion_t;
#pragma link C++ class EXT_STR_h101_MUSLI_onion_t;
#pragma link C++ class EXT_STR_h101_SAMP_onion_t;
#pragma link C++ class EXT_STR_h101_SAMPLOS_onion_t;
#pragma link C++ class EXT_STR_h101_PDC_onion_t;
#pragma link C++ class EXT_STR_h101_LOS_t;
#pragma link C++ class EXT_STR_h101_RPC_t;
#pragma link C++ class EXT_STR_h101_SYNC_CHECK_t;
#pragma link C++ class EXT_STR_h101_MOSAIC_onion_t;
#pragma link C++ class EXT_STR_h101_ACTAF_onion_t;

#endif
