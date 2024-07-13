# cmake-format: off
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
# cmake-format: on

set(SRCS
    ElasticScattering.cxx
    IsElastic.cxx
    R3BNeulandGeoPar.cxx
    R3BNeulandContFact.cxx
    R3BProgramOptions.cxx
    R3BNeulandSignalMatcher.cxx
    R3BNeulandVisualizer.cxx)

set(HEADERS
    # cmake-format: sortable
    ClusteringEngine.h
    ElasticScattering.h
    Filterable.h
    IsElastic.h
    R3BNeulandCommon.h
    R3BNeulandContFact.h
    R3BNeulandGeoPar.h
    R3BNeulandSignalMatcher.h
    R3BNeulandVisualizer.h
    R3BProgramOptions.h
    TCAConnector.h
    Validated.h)

add_library(R3BNeulandShared SHARED ${SRCS})
target_include_directories(R3BNeulandShared PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(R3BNeulandShared PUBLIC R3BData Boost::program_options)

root_generate_dictionary(
    R3BNeulandShared_dict
    ${HEADERS}
    MODULE
    R3BNeulandShared
    LINKDEF
    NeulandSharedLinkDef.h)
