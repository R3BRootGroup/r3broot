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

# Create a library called "libR3BSfib" which includes the source files given
# in the array. The extension is already found. Any number of sources could
# be listed here.

set(SRCS
    sim/R3Bsfi.cxx sim/R3BsfiDigitizer.cxx R3BSfibMapped2Cal.cxx R3BSfibMapped2CalPar.cxx
    R3BSfibCal2Hit.cxx
    #R3BSfibCal2HitPar.cxx
)

set(HEADERS
    sim/R3Bsfi.h sim/R3BsfiDigitizer.h R3BSfibMapped2Cal.h R3BSfibMapped2CalPar.h R3BSfibCal2Hit.h
    #R3BSfibCal2HitPar.h
)
add_library_with_dictionary(
    LIBNAME
    R3BSfib
    LINKDEF
    SfibLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    sim
    DEPENDENCIES
    R3BTCal
    R3BData)
