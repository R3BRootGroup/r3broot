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

# Create a library called "libR3BSci2" which includes the source files given in
# the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    calib/R3BSci2Mapped2Tcal.cxx
    calib/R3BSci2Mapped2Cal.cxx
    params/R3BSci2Mapped2CalPar.cxx
    calib/R3BSci2Tcal2Hit.cxx
    params/R3BSci2HitPar.cxx
    params/R3BSci2ContFact.cxx
    online/R3BOnlineSpectraSci2.cxx)

set(HEADERS
    calib/R3BSci2Mapped2Tcal.h
    calib/R3BSci2Mapped2Cal.h
    params/R3BSci2Mapped2CalPar.h
    calib/R3BSci2Tcal2Hit.h
    params/R3BSci2HitPar.h
    params/R3BSci2ContFact.h
    online/R3BOnlineSpectraSci2.h)

add_library_with_dictionary(
    LIBNAME
    R3BSci2
    LINKDEF
    Sci2LinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    calib
    online
    params
    DEPENDENCIES
    R3BData
    R3BTCal)
