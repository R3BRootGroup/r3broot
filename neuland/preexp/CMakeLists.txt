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

set(SRCS R3BNeulandCheckMapping.cxx R3BNeulandGainMatching.cxx R3BNeulandTimeRes.cxx)
set(HEADERS R3BNeulandCheckMapping.h R3BNeulandGainMatching.h R3BNeulandTimeRes.h)

add_library(R3BNeulandPreexp SHARED ${SRCS})
target_include_directories(R3BNeulandPreexp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(R3BNeulandPreexp PUBLIC R3BData R3BChannelAccess)
root_generate_dictionary(R3BNeulandPreexp_dict ${HEADERS} MODULE R3BNeulandPreexp LINKDEF
                         R3BNeulandPreexpLinkDef.h)
