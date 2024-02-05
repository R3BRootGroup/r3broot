// clang-format off

/******************************************************************************
 *   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    *
 *   Copyright (C) 2019-2024 Members of R3B Collaboration                     *
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

#pragma link C++ class R3BStack+;
#pragma link C++ class R3BMCTrack+;
#pragma link C++ class R3BHit+;
#pragma link C++ class R3BDouble+;
#pragma link C++ class R3BTrack+;

#pragma link C++ class R3BXBallPoint+;
#pragma link C++ class R3BXBallCrystalHit+;
#pragma link C++ class R3BXBallCrystalHitSim+;
#pragma link C++ class R3BCalifaMappedData+;
#pragma link C++ class R3BCalifaCrystalCalData+;
#pragma link C++ class R3BCalifaClusterData+;
#pragma link C++ class R3BCalifaPoint+;
#pragma link C++ class R3BDchPoint+;
#pragma link C++ class R3BDchFullPoint+;
#pragma link C++ class R3BDchDigi+;
#pragma link C++ class R3BDch2pDigi+;
#pragma link C++ class R3BGfiPoint+;
#pragma link C++ class R3BGfiDigi+;
#pragma link C++ class R3BmTofPoint+;
#pragma link C++ class R3BmTofDigi+;
#pragma link C++ class R3BTofdPoint+;
#pragma link C++ class R3BTofdDigi+;
#pragma link C++ class R3BTofDigi+;
#pragma link C++ class R3BTofPoint+;
#pragma link C++ class R3BTof2pDigi+;
#pragma link C++ class R3BTofiPoint+;
#pragma link C++ class R3BTofiDigi+;
#pragma link C++ class R3BTofiMappedData+;
#pragma link C++ class R3BTofiCalData+;
#pragma link C++ class R3BTofiHitData+;
#pragma link C++ class R3BVetoPoint+;
#pragma link C++ class R3BAmsMappedData+;
#pragma link C++ class R3BAmsStripCalData+;
#pragma link C++ class R3BAmsHitData+;
#pragma link C++ class R3BFootMappedData+;
#pragma link C++ class R3BFootCalData+;
#pragma link C++ class R3BFootHitData+;
#pragma link C++ class R3BAlpidePoint+;
#pragma link C++ class R3BAlpideMappedData+;
#pragma link C++ class R3BAlpideCalData+;
#pragma link C++ class R3BAlpideHitData+;
#pragma link C++ class R3BRpcMappedData+;
#pragma link C++ class R3BRpcPreCalData+;
#pragma link C++ class R3BRpcCalData+;
#pragma link C++ class R3BRpcHitData+;
#pragma link C++ class R3BRpcPoint+;

#pragma link C++ class R3BWRData+;
#pragma link C++ class R3BTrloiiData+;

#pragma link C++ class R3BMusicPoint+;
#pragma link C++ class R3BMusicMappedData+;
#pragma link C++ class R3BMusicCalData+;
#pragma link C++ class R3BMusicHitData+;

#pragma link C++ class R3BMwpcMappedData+;
#pragma link C++ class R3BMwpcCalData+;
#pragma link C++ class R3BMwpcHitData+;
#pragma link C++ class R3BMwpcPoint+;

#pragma link C++ class R3BTwimMappedData+;
#pragma link C++ class R3BTwimCalData+;
#pragma link C++ class R3BTwimHitData+;
#pragma link C++ class R3BTwimPoint+;

#pragma link C++ class R3BTttxMappedData+;
#pragma link C++ class R3BTttxCalData+;
#pragma link C++ class R3BTttxHitData+;

#pragma link C++ class R3BMusliMappedData+;
#pragma link C++ class R3BMusliCalData+;
#pragma link C++ class R3BMusliHitData+;

#pragma link C++ class R3BTraPoint+;
#pragma link C++ class R3BTraHit+;
#pragma link C++ class R3BTrackerHit+;
#pragma link C++ class R3BTraDigi+;
#pragma link C++ class R3BTra2pDigi+;
#pragma link C++ class R3BTraFraDigi+;
#pragma link C++ class R3BTargetDigi+;
#pragma link C++ class R3BTarget2pDigi+;
#pragma link C++ class R3BSitrackerHit+;
#pragma link C++ class R3BLandPoint+;
#pragma link C++ class R3BNeutronTrack+;
#pragma link C++ class R3BLandDigi+;
#pragma link C++ class R3BLandFirstHits+;
#pragma link C++ class R3BLandHit+;
#pragma link C++ class R3BPaddleCrystalHit+;
#pragma link C++ class ELILuMonPoint+;
#pragma link C++ class R3BMfiPoint+;
#pragma link C++ class R3BMfiDigi+;
#pragma link C++ class R3BPspPoint+;
#pragma link C++ class R3BPspDigi+;
#pragma link C++ class R3BNeulandTacquilaMappedData+;
#pragma link C++ class R3BNeulandPoint+;
#pragma link C++ class std::vector<R3BNeulandPoint>+;
#pragma link C++ class R3BNeulandHit+;
#pragma link C++ class R3BNeulandMultiplicity+;
#pragma link C++ class R3BNeulandCluster+;
#pragma link C++ class R3BNeulandNeutron+;
#pragma link C++ class R3BPspxMappedData+;
#pragma link C++ class R3BPspxPrecalData+;
#pragma link C++ class R3BPspxCalData+;
#pragma link C++ class R3BPspxHitData+;
#pragma link C++ class R3BLosMappedData+;
#pragma link C++ class R3BLosCalData+;
#pragma link C++ class R3BLosTCalData+;
#pragma link C++ class R3BLosHitData+;
#pragma link C++ class R3BRoluMappedData+;
#pragma link C++ class R3BRoluCalData+;
#pragma link C++ class R3BRoluHitData+;
#pragma link C++ class R3BFrsSciMappedData+;
#pragma link C++ class R3BFrsSciTcalData+;
#pragma link C++ class R3BFrsSciPosCalData+;
#pragma link C++ class R3BFrsSciTofCalData+;
#pragma link C++ class R3BSci2MappedData+;
#pragma link C++ class R3BSci2TcalData+;
#pragma link C++ class R3BSci2CalData+;
#pragma link C++ class R3BSci2HitData+;
#pragma link C++ class R3BSci8MappedData+;
#pragma link C++ class R3BSci8CalData+;
#pragma link C++ class R3BSci8HitData+;
#pragma link C++ class R3BPaddleTamexMappedData+;
#pragma link C++ class R3BPaddleCalData+;
#pragma link C++ class R3BTofdMappedData+;
#pragma link C++ class R3BTofdCalData+;
#pragma link C++ class R3BTofdHitData+;
#pragma link C++ class R3BPdcMappedData+;
#pragma link C++ class R3BPdcCalData+;
#pragma link C++ class R3BPdcHitData+;
#pragma link C++ class R3BNeulandCalData+;
#pragma link C++ class R3BActarPoint+;
#pragma link C++ class R3BStrawtubesMappedData+;
#pragma link C++ class R3BStrawtubesCalData+;
#pragma link C++ class R3BStrawtubesHitData+;
#pragma link C++ class R3BPtofHit+;
#pragma link C++ class R3BBunchedFiberCalData+;
#pragma link C++ class R3BBunchedFiberHitData+;
#pragma link C++ class R3BBunchedFiberMappedData+;
#pragma link C++ class R3BFiberMappedData+;
#pragma link C++ class R3BFiberMAPMTCalData+;
#pragma link C++ class R3BFiberMAPMTHitData+;
#pragma link C++ class R3BFiberMAPMTMappedData+;
#pragma link C++ class R3BFi4Point+;
#pragma link C++ class R3BFi4CalItem+;
#pragma link C++ class R3BFi4HitItem+;
#pragma link C++ class R3BFibMappedData+;
#pragma link C++ class R3BFibCalData+;
#pragma link C++ class R3BFibPoint+;
#pragma link C++ class R3BSci8MappedData+;
#pragma link C++ class R3BSci8CalData+;
#pragma link C++ class R3BSci8HitData+;
#pragma link C++ class R3BBeamMonitorMappedData+;
#pragma link C++ class R3BSamplerMappedData+;
#pragma link C++ class R3BSampLosMappedData+;
#pragma link C++ class R3BSfibMappedData+;
#pragma link C++ class R3BSfibCalData+;
#pragma link C++ class R3BSfibHitData+;
#pragma link C++ class R3BFrsData+;
#pragma link C++ class R3BSyncCheckData+;

#endif
