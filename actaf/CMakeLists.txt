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

# Create a library called "libR3BActaf" which includes the source files given 
# in the array. The extension is already found. Any number of sources could be
# listed here.

file(
    GLOB
    SRCS
    ./sim/*.cxx
    ./pars/*.cxx
    ./online/*.cxx
#    ./calibration/*.cxx
#    ./ana/*.cxx
    )

change_file_extension(*.cxx *.h HEADERS "${SRCS}")

add_library_with_dictionary(
    LIBNAME
    R3BActaf
    LINKDEF
    ActafLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    sim
    pars
    online
    calibration
    ana
    DEPENDENCIES
    R3BBase 
    R3BTracking
    R3BData)

