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

set(SRCS R3BNeulandTcalFill.cxx R3BNeulandTcal.cxx R3BNeulandPmt.cxx)
set(SRCS R3BNeulandTcalFill.h R3BNeulandTcal.h R3BNeulandPmt.h)

add_library(R3BNeulandUnpack SHARED ${SRCS})
target_include_directories(R3BNeulandUnpack PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(R3BNeulandUnpack PUBLIC R3BNeulandShared)
root_generate_dictionary(R3BNeulandUnpack_dict ${HEADERS} MODULE R3BNeulandUnpack LINKDEF
                         NeulandUnpackLinkDef.h)
