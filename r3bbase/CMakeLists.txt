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

set(SRCS
    # cmake-format: sortable
    data_monitor/R3BDataMonitor.cxx
    data_monitor/R3BDataMonitorCanvas.cxx
    R3BCoarseTimeStitch.cxx
    R3BDataPropagator.cxx
    R3BDetector.cxx
    R3BEventHeader.cxx
    R3BEventHeaderPropagator.cxx
    R3BException.cxx
    R3BFileSource.cxx
    R3BFileSource2.cxx
    R3BLogger.cxx
    R3BModule.cxx
    R3BTcutPar.cxx
    R3BTsplinePar.cxx
    R3BWhiterabbitPropagator.cxx)

set(HEADERS
    # cmake-format: sortable
    data_monitor/R3BDataMonitor.h
    data_monitor/R3BDataMonitorCanvas.h
    R3BCoarseTimeStitch.h
    R3BDataPropagator.h
    R3BDetector.h
    R3BEventHeader.h
    R3BEventHeaderPropagator.h
    R3BException.h
    R3BFileSource.h
    R3BFileSource2.h
    R3BIOConnector.h
    R3BLogger.h
    R3BModule.h
    R3BShared.h
    R3BTcutPar.h
    R3BTsplinePar.h
    R3BWhiterabbitPropagator.h)

set(DEPENDENCIES FairRoot::Generators FairRoot::Online)

add_library_with_dictionary(
    LIBNAME
    R3BBase
    LINKDEF
    BaseLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    data_monitor
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    ${DEPENDENCIES})

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    message(STATUS "Linking -lstdc++fs due to the old version of gcc")
    target_link_libraries(R3BBase PUBLIC -lstdc++fs)
endif()

if("${FairRoot_VERSION}" VERSION_LESS 19.0.0)
    target_link_directories(R3BBase PUBLIC ${ROOT_LIBRARY_DIR})
endif()
