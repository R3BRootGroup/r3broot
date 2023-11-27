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

# Create a library called "libR3BLos" which includes the source files given in
# the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    calib/R3BLosMapped2Cal.cxx
    calib/R3BLosMapped2TCal.cxx
    params/R3BLosMapped2CalPar.cxx
    calib/R3BLosCal2Hit.cxx
    calib/R3BLosTCal2Hit.cxx
    R3BLosProvideTStart.cxx
    params/R3BLosHitPar.cxx
    params/R3BLosContFact.cxx
    online/R3BLosOnlineSpectra.cxx)

set(HEADERS
    calib/R3BLosMapped2Cal.h
    calib/R3BLosMapped2TCal.h
    params/R3BLosMapped2CalPar.h
    calib/R3BLosCal2Hit.h
    calib/R3BLosTCal2Hit.h
    R3BLosProvideTStart.h
    params/R3BLosHitPar.h
    params/R3BLosContFact.h
    online/R3BLosOnlineSpectra.h)

set(LINKDEF LosLinkDef.h)
set(LIBRARY_NAME R3BLos)
set(DEPENDENCIES R3BBase R3BData R3BTracking R3BTCal)

add_library_with_dictionary(
    LIBNAME
    R3BLos
    LINKDEF
    LosLinkDef.h
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
    R3BTracking
    R3BTCal
    # TODO: only due to TCAConnector.
    R3BNeulandShared)
