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

# Create a library called "libR3BAnalysis" which includes the source files given
# in the array. The extension is already found. Any number of sources could be
# listed here.

Set(SYSTEM_INCLUDE_DIRECTORIES
${SYSTEM_INCLUDE_DIRECTORIES}
${BASE_INCLUDE_DIRECTORIES}
)

set(INCLUDE_DIRECTORIES
#put here all directories where header files are located
${R3BROOT_SOURCE_DIR}/r3bbase
${R3BROOT_SOURCE_DIR}/passive
${R3BROOT_SOURCE_DIR}/tracking
${R3BROOT_SOURCE_DIR}/los
${R3BROOT_SOURCE_DIR}/rolu
${R3BROOT_SOURCE_DIR}/sci2
${R3BROOT_SOURCE_DIR}/sci8
${R3BROOT_SOURCE_DIR}/music/pars
${R3BROOT_SOURCE_DIR}/tof
${R3BROOT_SOURCE_DIR}/tofi
${R3BROOT_SOURCE_DIR}/tofd/online
${R3BROOT_SOURCE_DIR}/tcal
${R3BROOT_SOURCE_DIR}/fiber
${R3BROOT_SOURCE_DIR}/fiber/pars
${R3BROOT_SOURCE_DIR}/fibers
${R3BROOT_SOURCE_DIR}/fibers/fi30
${R3BROOT_SOURCE_DIR}/fibers/fi31
${R3BROOT_SOURCE_DIR}/fibers/fi32
${R3BROOT_SOURCE_DIR}/fibers/fi33
${R3BROOT_SOURCE_DIR}/fibers/fi23a
${R3BROOT_SOURCE_DIR}/fibers/fi23b
${R3BROOT_SOURCE_DIR}/ssd/online
${R3BROOT_SOURCE_DIR}/los/online
${R3BROOT_SOURCE_DIR}/califa/online
${R3BROOT_SOURCE_DIR}/music/online
${R3BROOT_SOURCE_DIR}/twim/online
${R3BROOT_SOURCE_DIR}/mwpc/online
${R3BROOT_SOURCE_DIR}/r3bdata
${R3BROOT_SOURCE_DIR}/r3bdata/landData
${R3BROOT_SOURCE_DIR}/r3bdata/beammonitorData
${R3BROOT_SOURCE_DIR}/r3bdata/neulandData
${R3BROOT_SOURCE_DIR}/r3bdata/califaData
${R3BROOT_SOURCE_DIR}/r3bdata/alpideData
${R3BROOT_SOURCE_DIR}/r3bdata/tttxData
${R3BROOT_SOURCE_DIR}/r3bdata/startrackData
${R3BROOT_SOURCE_DIR}/r3bdata/losData
${R3BROOT_SOURCE_DIR}/r3bdata/roluData
${R3BROOT_SOURCE_DIR}/r3bdata/sci2Data
${R3BROOT_SOURCE_DIR}/r3bdata/sci8Data
${R3BROOT_SOURCE_DIR}/r3bdata/tofData
${R3BROOT_SOURCE_DIR}/r3bdata/tofiData
${R3BROOT_SOURCE_DIR}/r3bdata/pdcData
${R3BROOT_SOURCE_DIR}/r3bdata/pspData
${R3BROOT_SOURCE_DIR}/r3bdata/fibData
${R3BROOT_SOURCE_DIR}/r3bdata/frssciData
${R3BROOT_SOURCE_DIR}/r3bdata/ptofData
${R3BROOT_SOURCE_DIR}/r3bdata/sampData
${R3BROOT_SOURCE_DIR}/r3bdata/sfibData
${R3BROOT_SOURCE_DIR}/r3bdata/musicData
${R3BROOT_SOURCE_DIR}/r3bdata/musliData
${R3BROOT_SOURCE_DIR}/r3bdata/twimData
${R3BROOT_SOURCE_DIR}/r3bdata/mwpcData
${R3BROOT_SOURCE_DIR}/r3bdata/frsData
${R3BROOT_SOURCE_DIR}/r3bdata/footData
${R3BROOT_SOURCE_DIR}/r3bdata/rpcData
${R3BROOT_SOURCE_DIR}/r3bdata/wrData
${R3BROOT_SOURCE_DIR}/r3bdata/synccheckData
${R3BROOT_SOURCE_DIR}/r3bdata/trloiiData
${R3BROOT_SOURCE_DIR}/tracker_rene
${R3BROOT_SOURCE_DIR}/analysis
${R3BROOT_SOURCE_DIR}/analysis/pars
${R3BROOT_SOURCE_DIR}/analysis/offline
${R3BROOT_SOURCE_DIR}/analysis/online
${R3BROOT_SOURCE_DIR}/neuland/shared
)

include_directories( ${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

link_directories(${FAIRROOT_LIBRARY_DIR} ${ROOT_LIBRARY_DIR})

file(GLOB SRCS ./online/*.cxx ./offline/*.cxx ./pars/*.cxx)

# fill list of header files from list of source files
# by exchanging the file extension
CHANGE_FILE_EXTENSION(*.cxx *.h HEADERS "${SRCS}")

Set(LINKDEF AnaLinkDef.h)

Set(DEPENDENCIES
GeoBase ParBase MbsAPI Base FairTools R3BPassive R3BData Core Geom GenVector Physics Matrix MathCore R3BTraRene R3BTracking
R3BSsd R3BLos R3BCalifa R3BMwpc R3BTofD R3BTwim R3BMusic R3BBunchedFiber
)

Set(LIBRARY_NAME R3BAnalysis)

GENERATE_LIBRARY()
