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
    R3BDigitizingChannel.cxx
    R3BDigitizingPaddle.cxx
    R3BDigitizingTacQuila.cxx
    R3BDigitizingTamex.cxx
    R3BDigitizingPaddleNeuland.cxx
    R3BNeulandHitMon.cxx
    R3BNeulandDigitizer.cxx
    R3BNeulandHitMon.cxx
    NeulandPointFilter.cxx)

set(HEADERS
    R3BDigitizingChannel.h
    R3BDigitizingChannelMock.h
    R3BDigitizingEngine.h
    R3BDigitizingPaddle.h
    R3BDigitizingPaddleMock.h
    R3BDigitizingPaddleNeuland.h
    R3BDigitizingTacQuila.h
    R3BDigitizingTamex.h
    R3BNeulandDigitizer.h
    R3BNeulandHitMon.h
    NeulandPointFilter.h)

add_library_with_dictionary(
    LIBNAME
    R3BNeulandDigitizing
    LINKDEF
    NeulandDigitizingLinkDef.h
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    HEADERS
    ${HEADERS}
    DEPENDENCIES
    R3BNeulandCalibration)
