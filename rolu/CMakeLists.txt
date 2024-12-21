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

# Create a library called "libR3BRolu" which includes the source files given in
# the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS R3BRoluMapped2Cal.cxx R3BRoluMapped2CalPar.cxx R3BRoluCal2Hit.cxx R3BRoluOnlineSpectra.cxx)

set(HEADERS R3BRoluMapped2Cal.h R3BRoluMapped2CalPar.h R3BRoluCal2Hit.h R3BRoluOnlineSpectra.h)

add_library_with_dictionary(
    LIBNAME
    R3BRolu
    LINKDEF
    RoluLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData
    R3BTCal)
