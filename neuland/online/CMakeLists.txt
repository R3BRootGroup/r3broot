##############################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2019-2025 Members of R3B Collaboration                     #
#                                                                            #
#             This software is distributed under the terms of the            #
#                 GNU General Public Licence (GPL) version 3,                #
#                    copied verbatim in the file "LICENSE".                  #
#                                                                            #
# In applying this license GSI does not waive the privileges and immunities  #
# granted to it by virtue of its status as an Intergovernmental Organization #
# or submit itself to any jurisdiction.                                      #
##############################################################################
set(SRCS R3BNeulandOnlineReconstruction.cxx R3BNeulandOnlineSpectra.cxx)
set(HEADERS R3BNeulandOnlineReconstruction.h R3BNeulandOnlineSpectra.h)

add_library(R3BNeulandOnline SHARED ${SRCS})
target_include_directories(R3BNeulandOnline PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(R3BNeulandOnline PUBLIC R3BTCal R3BNeulandShared)
root_generate_dictionary(R3BNeulandOnline_dict ${HEADERS} MODULE R3BNeulandOnline LINKDEF
                         NeulandOnlineLinkDef.h)
