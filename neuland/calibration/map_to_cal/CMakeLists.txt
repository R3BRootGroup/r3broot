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

set(SRCS
    # cmake-format: sortable
    R3BFTCalEngine.cxx
    R3BNeulandMapToCalPar.cxx
    R3BNeulandMapToCalParTask.cxx
    R3BNeulandMapToCalTask.cxx
    R3BTrigIDMappingFinder.cxx)

set(HEADERS
    # cmake-format: sortable
    R3BFTCalEngine.h
    R3BNeulandCalData2.h
    R3BNeulandMapToCalPar.h
    R3BNeulandMapToCalParTask.h
    R3BNeulandMapToCalTask.h
    R3BTrigIDMappingFinder.h)

add_library(R3BNeulandMapToCal SHARED ${SRCS})
target_include_directories(R3BNeulandMapToCal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
    R3BNeulandMapToCal
    PUBLIC R3BNeulandCalibrationShare range-v3::range-v3
    PRIVATE nlohmann_json::nlohmann_json)

root_generate_dictionary(
    R3BNeulandMapToCal_dict
    ${HEADERS}
    MODULE
    R3BNeulandMapToCal
    LINKDEF
    NeulandMapToCalLinkDef.h)
