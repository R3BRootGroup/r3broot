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

set(SRCS R3BGladMagnet.cxx)
set(HEADERS R3BGladMagnet.h)

add_library_with_dictionary(
    LIBNAME
    R3BGlad
    LINKDEF
    GladLinkDef.h
    SRCS
    R3BGladMagnet.cxx
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    HEADERS
    R3BGladMagnet.h
    DEPENDENCIES
    R3BTracking)

if(BUILD_GEOMETRY)
    add_subdirectory(geobase)
endif()
