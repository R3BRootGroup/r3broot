# ##################################################################################################
# Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    # Copyright (C) 2019-2024
# Members of R3B Collaboration                     # # This software is distributed under the terms
# of the            # GNU General Public Licence (GPL) version 3,                # copied verbatim
# in the file "LICENSE".                  # # In applying this license GSI does not waive the
# privileges and immunities  # granted to it by virtue of its status as an Intergovernmental
# Organization # or submit itself to any jurisdiction.                                      #
# ##################################################################################################

set(SRCS
    R3BNeuland.cxx
    R3BNeulandMCMon.cxx
    R3BNeulandPrimaryClusterFinder.cxx
    R3BNeulandPrimaryInteractionFinder.cxx)

set(HEADERS
    R3BNeuland.h
    R3BNeulandMCMon.h
    R3BNeulandPrimaryClusterFinder.h
    R3BNeulandPrimaryInteractionFinder.h)

add_library(R3BNeulandSimulation SHARED ${SRCS})
target_include_directories(R3BNeulandSimulation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(R3BNeulandSimulation PUBLIC R3BNeulandShared R3BPassive)
root_generate_dictionary(
    R3BNeulandSimulation_dict
    ${HEADERS}
    MODULE
    R3BNeulandSimulation
    LINKDEF
    NeulandSimulationLinkDef.h)
