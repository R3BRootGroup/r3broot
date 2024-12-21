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

# Create a library called "libR3BPassive" which includes the source files given
# in the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    R3BCave.cxx
    R3BGeoCave.cxx
    R3BGeoPassivePar.cxx
    R3BGeoTarget.cxx
    R3BAladinMagnet.cxx
    R3BNeutronWindowAndSomeAir.cxx
    R3BPassiveContFact.cxx
    R3BPipe.cxx
    R3BVacVesselCool.cxx
    R3BTarget.cxx
    R3BCollimator.cxx)

set(HEADERS
    R3BCave.h
    R3BGeoCave.h
    R3BGeoPassivePar.h
    R3BGeoTarget.h
    R3BAladinMagnet.h
    R3BNeutronWindowAndSomeAir.h
    R3BPassiveContFact.h
    R3BPipe.h
    R3BVacVesselCool.h
    R3BTarget.h
    R3BCollimator.h)

add_library_with_dictionary(
    LIBNAME
    R3BPassive
    LINKDEF
    PassiveLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BTracking)
