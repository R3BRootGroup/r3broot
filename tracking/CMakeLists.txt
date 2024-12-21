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

# Create a library called "libR3BTracking" which includes the source files given
# in the array. The extension is already found. Any number of sources could be
# listed here.

set(SRCS
    R3BTPropagator.cxx
    R3BTGeoPar.cxx
    R3BFragmentTracker.cxx
    R3BFragmentFitterGeneric.cxx
    R3BFragmentFitterChi2.cxx
    R3BTrackingDetector.cxx
    R3BTrackingParticle.cxx
    R3BTrackingSetup.cxx
    R3BMDFWrapper.cxx
    R3BTrackingS515.cxx)

set(HEADERS
    R3BTPropagator.h
    R3BTGeoPar.h
    R3BFragmentTracker.h
    R3BFragmentFitterGeneric.h
    R3BFragmentFitterChi2.h
    R3BTrackingDetector.h
    R3BTrackingParticle.h
    R3BTrackingSetup.h
    R3BMDFWrapper.h
    R3BTrackingS515.h)

add_library_with_dictionary(
    LIBNAME
    R3BTracking
    LINKDEF
    TrackingLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BData
    Field
    FairRoot::TrackBase)
