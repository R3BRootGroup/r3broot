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

file(
    GLOB
    SRCS
    ./sim/*.cxx
    ./digi/*.cxx
    ./pars/*.cxx
    ./calibration/*.cxx
    ./online/*.cxx)
change_file_extension(*.cxx *.h HEADERS "${SRCS}")

add_library_with_dictionary(
    LIBNAME
    R3BMusic
    LINKDEF
    MusicLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    sim
    digi
    pars
    calibration
    online
    DEPENDENCIES
    R3BTracking)

if(BUILD_GEOMETRY)
    add_subdirectory(geobase)
endif()
add_subdirectory(test)
