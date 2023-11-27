##############################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2019-2024 Members of R3B Collaboration                     #
#                                                                            #
#             This software is distributed under the terms of the            #
#                 GNU General Public Licence (GPL) version 3,                #
#                    copied verbatim in the file "LICENSE".                  #
#                                                                            #
# In applying this license GSI does not waive the privileges and immunities  #
# granted to it by virtue of its status as an Intergovernmental Organization #
# or submit itself to any jurisdiction.                                      #
##############################################################################

set(SRCS
    R3BNeulandMappedHist.cxx
    R3BNeulandMapped2Cal.cxx
    R3BNeulandMapped2CalPar.cxx
    R3BNeulandTacquilaMapped2Cal.cxx
    R3BNeulandTacquilaMapped2CalPar.cxx
    R3BNeulandTacquilaMapped2QCalPar.cxx
    R3BNeulandHitHist.cxx
    # R3BNeulandCalTest.cxx
    R3BNeulandHitPar.cxx
    R3BNeulandCosmicTracker.cxx
    R3BNeulandHitCalibrationEngine.cxx
    R3BNeulandHitCalibrationBar.cxx
    R3BNeulandTSyncer.cxx
    LSQR.cxx
    R3BNeulandCal2HitPar.cxx
    R3BNeulandParFact.cxx
    R3BNeulandCal2Hit.cxx
    R3BNeulandHitModulePar.cxx
    R3BNeulandQCalPar.cxx
    R3BNeulandQCalFiller.cxx
    R3BNeulandProvideTStart.cxx
    R3BNeulandMappingPar.cxx)

set(HEADERS
    R3BNeulandMappedHist.h
    R3BNeulandMapped2Cal.h
    R3BNeulandMapped2CalPar.h
    R3BNeulandTacquilaMapped2Cal.h
    R3BNeulandTacquilaMapped2CalPar.h
    R3BNeulandTacquilaMapped2QCalPar.h
    R3BNeulandHitHist.h
    # R3BNeulandCalTest.h
    R3BNeulandHitPar.h
    R3BNeulandCosmicTracker.h
    R3BNeulandHitCalibrationEngine.h
    R3BNeulandHitCalibrationBar.h
    R3BNeulandTSyncer.h
    LSQR.h
    R3BNeulandCal2HitPar.h
    R3BNeulandParFact.h
    R3BNeulandCal2Hit.h
    R3BNeulandHitModulePar.h
    R3BNeulandQCalPar.h
    R3BNeulandQCalFiller.h
    R3BNeulandProvideTStart.h
    R3BNeulandMappingPar.h)

add_library(R3BNeulandCalibration SHARED ${SRCS})
target_include_directories(R3BNeulandCalibration PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(R3BNeulandCalibration PUBLIC R3BTCal R3BNeulandShared)
root_generate_dictionary(R3BNeulandCalibration_dict ${HEADERS} MODULE R3BNeulandCalibration LINKDEF
                         NeulandCalibrationLinkDef.h)
