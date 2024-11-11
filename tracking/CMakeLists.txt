##############################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH    #
#   Copyright (C) 2019 Members of R3B Collaboration                          #
#                                                                            #
#             This software is distributed under the terms of the            #
#                 GNU General Public Licence (GPL) version 3,                #
#                    copied verbatim in the file "LICENSE".                  #
#                                                                            #
# In applying this license GSI does not waive the privileges and immunities  #
# granted to it by virtue of its status as an Intergovernmental Organization #
# or submit itself to any jurisdiction.                                      #
##############################################################################

# Create a library called "libR3BTracking" which includes the source files given in
# the array .
# The extension is already found.  Any number of sources could be listed here.

Set(SYSTEM_INCLUDE_DIRECTORIES
${SYSTEM_INCLUDE_DIRECTORIES}
${BASE_INCLUDE_DIRECTORIES}
)

set(INCLUDE_DIRECTORIES
#put here all directories where header files are located
${R3BROOT_SOURCE_DIR}/atima
${R3BROOT_SOURCE_DIR}/tracking
${R3BROOT_SOURCE_DIR}/field
${R3BROOT_SOURCE_DIR}/tof
${R3BROOT_SOURCE_DIR}/fiber
${R3BROOT_SOURCE_DIR}/passive
${R3BROOT_SOURCE_DIR}/r3bdata
${R3BROOT_SOURCE_DIR}/r3bbase
${R3BROOT_SOURCE_DIR}/r3bdata/tofData
${R3BROOT_SOURCE_DIR}/r3bdata/pspData
${R3BROOT_SOURCE_DIR}/r3bdata/fibData
${R3BROOT_SOURCE_DIR}/r3bdata/califaData
${Geant4_INCLUDE_DIRS}
)

include_directories( ${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES
${ROOT_LIBRARY_DIR}
${FAIRROOT_LIBRARY_DIR}
)

link_directories( ${LINK_DIRECTORIES})

set(SRCS
R3BTPropagator.cxx
R3BTGeoPar.cxx
R3BFragmentFitterGeneric.cxx
R3BFragmentTrackerS494.cxx
R3BFragmentTrajectoryS494.cxx
#R3BFragmentTrackerPAS.cxx
R3BFragmentFitterChi2S494.cxx
#R3BFragmentFitterChi2PAS.cxx
R3BTrackingDetector.cxx
R3BTrackingParticle.cxx
R3BTrackingSetup.cxx
R3BOptimizeGeometryS494.cxx
#R3BOptimizeMagfieldS494.cxx
)

# fill list of header files from list of source files
# by exchanging the file extension
CHANGE_FILE_EXTENSION(*.cxx *.h HEADERS "${SRCS}")

Set(LINKDEF TrackingLinkDef.h)
Set(LIBRARY_NAME R3BTracking)
Set(DEPENDENCIES
    Base ParBase TrkBase R3Bbase Minuit R3BData)

GENERATE_LIBRARY()

