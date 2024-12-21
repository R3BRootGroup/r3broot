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
set(SRCS R3BNeulandClusterMon.cxx R3BNeulandClusterFinder.cxx)
set(HEADERS R3BNeulandClusterMon.h R3BNeulandClusterFinder.h)

add_library(R3BNeulandClustering SHARED ${SRCS})
target_include_directories(R3BNeulandClustering PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(R3BNeulandClustering PUBLIC R3BNeulandShared)
root_generate_dictionary(R3BNeulandClustering_dict ${HEADERS} MODULE R3BNeulandClustering LINKDEF
                         NeulandClusteringLinkDef.h)
