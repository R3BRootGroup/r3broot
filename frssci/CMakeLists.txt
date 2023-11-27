#############################################################################
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

# Create a library called "libR3BFrsSci" which includes the source files given in
# the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    parameters/R3BFrsSciContFact.cxx
    parameters/R3BFrsSciTcalPar.cxx
    parameters/R3BFrsSciCalPar.cxx
    calibration/R3BFrsSciMapped2Tcal.cxx
    calibration/R3BFrsSciMapped2TcalPar.cxx
    calibration/R3BFrsSciTcal2Cal.cxx
    online/R3BOnlineSpectraFrsSci.cxx
    online/R3BOnlineSpectraFrsSciVsMusli.cxx
    online/R3BOnlineSpectraFrsSciVsLos.cxx)

set(HEADERS
    parameters/R3BFrsSciContFact.h
    parameters/R3BFrsSciTcalPar.h
    parameters/R3BFrsSciCalPar.h
    calibration/R3BFrsSciMapped2Tcal.h
    calibration/R3BFrsSciMapped2TcalPar.h
    calibration/R3BFrsSciTcal2Cal.h
    online/R3BOnlineSpectraFrsSci.h
    online/R3BOnlineSpectraFrsSciVsMusli.h
    online/R3BOnlineSpectraFrsSciVsLos.h)

add_library_with_dictionary(
    LIBNAME
    R3BFrsSci
    LINKDEF
    FrsSciLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    calibration
    online
    parameters
    DEPENDENCIES
    R3BData
    R3BTCal)
