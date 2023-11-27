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

# Create a library called "libR3BLand" which includes the source files given in
# the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    R3BLand.cxx
    R3BGeoLand.cxx
    R3BGeoLandPar.cxx
    R3BLandContFact.cxx
    R3BLandDigitizer.cxx
    R3BLandDigitizerQA.cxx
    R3BNeutronTracker.cxx
    R3BLandDigiPar.cxx
    R3BNeuLandCluster.cxx
    R3BNeuLandClusterFinder.cxx
    R3BNeutronCalibr2D.cxx
    R3BPrimPart.cxx
    R3BNeutronTracker2D.cxx
    R3BNeutHit.cxx)

set(HEADERS
    R3BLand.h
    R3BGeoLand.h
    R3BGeoLandPar.h
    R3BLandContFact.h
    R3BLandDigitizer.h
    R3BLandDigitizerQA.h
    R3BNeutronTracker.h
    R3BLandDigiPar.h
    R3BNeuLandCluster.h
    R3BNeuLandClusterFinder.h
    R3BNeutronCalibr2D.h
    R3BPrimPart.h
    R3BNeutronTracker2D.h
    R3BNeutHit.h)

add_library_with_dictionary(
    LIBNAME
    R3BLand
    LINKDEF
    LandLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData
    R3BTCal)
