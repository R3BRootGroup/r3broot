################################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH      #
#   Copyright (C) 2019-2025 Members of R3B Collaboration                       #
#                                                                              #
#             This software is distributed under the terms of the              #
#                 GNU General Public Licence (GPL) version 3,                  #
#                    copied verbatim in the file "LICENSE".                    #
#                                                                              #
# In applying this license GSI does not waive the privileges and immunities    #
# granted to it by virtue of its status as an Intergovernmental Organization   #
# or submit itself to any jurisdiction.                                        #
################################################################################

# FindAtima.cmake
#
# Finds the Atima library
#
# This will define the following variables
#
# Atima_FOUND Atima_LIBRARY
#
# and the following imported targets
#
# Atima::Atima
#

message(STATUS "Looking for Atima...")

find_library(Atima_LIBRARY
             NAMES Atima
             HINTS $ENV{ATIMAPATH})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Atima
                                  FOUND_VAR
                                  Atima_FOUND
                                  REQUIRED_VARS
                                  Atima_LIBRARY)

if(Atima_FOUND AND NOT TARGET Atima::Atima)
    message(STATUS "Found Atima: ${Atima_LIBRARY}")
    add_library(Atima::Atima UNKNOWN IMPORTED)
    set_target_properties(Atima::Atima PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${Atima_LIBRARY}"
    )
endif()
