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

set(SRCS R3BNeulandOnlineReconstruction.cxx R3BNeulandOnlineSpectra.cxx)
set(HEADERS R3BNeulandOnlineReconstruction.h R3BNeulandOnlineSpectra.h)

if(CONAN_ENABLED)
    list(
        APPEND
        SRCS
        # cmake-format: sortable
        canvases/R3BNeulandCalCanvas.cxx
        canvases/R3BNeulandCountRateCanvas.cxx
        canvases/R3BNeulandEventHeaderCanvas.cxx
        canvases/R3BNeulandHitCanvas.cxx
        canvases/R3BNeulandHitCosmicCanvas.cxx
        canvases/R3BNeulandHitXYCanvas.cxx
        canvases/R3BNeulandMappedCanvas.cxx
        canvases/R3BNeulandOnlineCanvas.cxx
        canvases/R3BNeulandTimingCanvas.cxx
        canvases/R3BNeulandTJumpCanvas.cxx
        R3BNeulandOnlineReconstruction.cxx
        R3BNeulandOnlineSpectra.cxx
        R3BNeulandOnlineSpectra2.cxx)
    list(
        APPEND
        HEADERS
        # cmake-format: sortable
        canvases/R3BNeulandCalCanvas.h
        canvases/R3BNeulandCountRateCanvas.h
        canvases/R3BNeulandEventHeaderCanvas.h
        canvases/R3BNeulandHitCanvas.h
        canvases/R3BNeulandHitCosmicCanvas.h
        canvases/R3BNeulandHitXYCanvas.h
        canvases/R3BNeulandMappedCanvas.h
        canvases/R3BNeulandOnlineCanvas.h
        canvases/R3BNeulandTimingCanvas.h
        canvases/R3BNeulandTJumpCanvas.h
        R3BNeulandOnlineReconstruction.h
        R3BNeulandOnlineSpectra.h
        R3BNeulandOnlineSpectra2.h)

endif()

add_library(R3BNeulandOnline SHARED ${SRCS})
target_include_directories(R3BNeulandOnline PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} canvases)
target_link_libraries(R3BNeulandOnline PUBLIC R3BNeulandCalibration)
if(CONAN_ENABLED)
    target_link_libraries(R3BNeulandOnline PRIVATE range-v3::range-v3)
endif()
root_generate_dictionary(
    R3BNeulandOnline_dict
    ${HEADERS}
    MODULE
    R3BNeulandOnline
    LINKDEF
    NeulandOnlineLinkDef.h)
