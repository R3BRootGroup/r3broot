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

message(STATUS "Detected c3w configuration. Setting default target to c3w tar file.")
message(STATUS "Run make all to build locally only.")
execute_process(COMMAND ln -vs ${CMAKE_SOURCE_DIR}/util/c3w/c3w.make ${CMAKE_CURRENT_BINARY_DIR}/../../GNUmakefile)
execute_process(COMMAND ln -vs ${CMAKE_SOURCE_DIR}/util/c3w/local_env.sh ${CMAKE_CURRENT_BINARY_DIR}/../../local_c3w.sh)
