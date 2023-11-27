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

file(
    GLOB
    SRCS
    ./sim/*.cxx
    ./digi/*.cxx
    ./pars/*.cxx
    ./online/*.cxx
    ./calibration/*.cxx
    ./ana/*.cxx)

change_file_extension(*.cxx *.h HEADERS "${SRCS}")

add_library_with_dictionary(
    LIBNAME
    R3BTwim
    LINKDEF
    TwimLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    sim
    digi
    pars
    online
    calibration
    ana
    DEPENDENCIES
    R3BTracking)

if(WITH_SOFIA)
    target_link_libraries(R3BTwim PRIVATE R3BSofData)
endif(WITH_SOFIA)

if(BUILD_GEOMETRY)
    add_subdirectory(geobase)
endif()

add_subdirectory(test)
