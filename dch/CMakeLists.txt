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

# Create a library called "libR3BDch" which includes the source files given in
# the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    R3BDch.cxx
    R3BGeoDch.cxx
    R3BGeoDchPar.cxx
    R3BDchContFact.cxx
    R3BDchDigiPar.cxx
    R3BDchDigitizer.cxx
    R3BDch2pDigiPar.cxx
    R3BDch2pDigitizer.cxx)

set(HEADERS
    R3BDch.h
    R3BGeoDch.h
    R3BGeoDchPar.h
    R3BDchContFact.h
    R3BDchDigiPar.h
    R3BDchDigitizer.h
    R3BDch2pDigiPar.h
    R3BDch2pDigitizer.h)

add_library_with_dictionary(
    LIBNAME
    R3BDch
    LINKDEF
    DchLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData)
