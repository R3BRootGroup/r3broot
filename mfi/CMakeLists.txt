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

# Create a library called "libR3BMfi" which includes the source files given in
# the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    R3BMfi.cxx
    R3BGeoMfi.cxx
    #R3BGeoMfiPar.cxx
    #R3BMfiContFact.cxx
    R3BMfiDigiPar.cxx
    R3BMfiDigitizer.cxx
    R3BKeepEventProbabilityMfi.cxx
    #unpack/R3BMfiRawHit.cxx
    #unpack/R3BMfiUnpack.cxx
)

set(HEADERS
    R3BMfi.h
    R3BGeoMfi.h
    #R3BGeoMfiPar.h
    #R3BMfiContFact.h
    R3BMfiDigiPar.h
    R3BMfiDigitizer.h
    R3BKeepEventProbabilityMfi.h
    #unpack/R3BMfiRawHit.h
    #unpack/R3BMfiUnpack.h
)

add_library_with_dictionary(
    LIBNAME
    R3BMfi
    LINKDEF
    MfiLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData)
