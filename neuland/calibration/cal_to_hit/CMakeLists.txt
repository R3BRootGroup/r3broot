# cmake-format: off
##############################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2019-2023 Members of R3B Collaboration                     #
#                                                                            #
#             This software is distributed under the terms of the            #
#                 GNU General Public Licence (GPL) version 3,                #
#                    copied verbatim in the file "LICENSE".                  #
#                                                                            #
# In applying this license GSI does not waive the privileges and immunities  #
# granted to it by virtue of its status as an Intergovernmental Organization #
# or submit itself to any jurisdiction.                                      #
##############################################################################
# cmake-format: on

set(LIBRARY_NAME R3BNeulandCalToHit)
set(LINKDEF NeulandCalToHitLinkDef.h)

set(SRCS
    # cmake-format: sortable
    engine/BarPositionFilter.cxx
    engine/R3BNeulandLSQREngineAdaptor.cxx
    engine/R3BNeulandMille.cxx
    engine/R3BNeulandMillepede.cxx
    engine/R3BNeulandPede.cxx
    engine/R3BNeulandPositionCal.cxx
    R3BNeulandCalToHitPar.cxx
    R3BNeulandCalToHitParTask.cxx
    R3BNeulandCalToHitTask.cxx)

set(HEADERS
    # cmake-format: sortable
    engine/BarPositionFilter.h
    engine/R3BNeulandCosmicEngine.h
    engine/R3BNeulandLSQREngineAdaptor.h
    engine/R3BNeulandMille.h
    engine/R3BNeulandMillepede.h
    engine/R3BNeulandPede.h
    engine/R3BNeulandPositionCal.h
    engine/Utilities.h
    R3BNeulandCalToHitPar.h
    R3BNeulandCalToHitParTask.h
    R3BNeulandCalToHitTask.h)

add_library(R3BNeulandCalToHit SHARED ${SRCS})
target_include_directories(R3BNeulandCalToHit PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} engine)
target_link_libraries(R3BNeulandCalToHit PUBLIC R3BNeulandMapToCal R3BNeulandCalibrationLegacy
                                                R3BMille)
root_generate_dictionary(
    R3BNeulandCalToHit_dict
    ${HEADERS}
    MODULE
    R3BNeulandCalToHit
    LINKDEF
    NeulandCalToHitLinkDef.h)
