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
    R3BEventManager.cxx
    R3BEventManagerEditor.cxx
    R3BCalifaEventDisplay.cxx
    R3BCalifaClusterEventDisplay.cxx
    R3BMCTracks.cxx
    R3BIonName.cxx)

set(HEADERS
    R3BEventManager.h
    R3BEventManagerEditor.h
    R3BCalifaEventDisplay.h
    R3BCalifaClusterEventDisplay.h
    R3BMCTracks.h
    R3BIonName.h)

add_library_with_dictionary(
    LIBNAME
    R3BEvtVis
    LINKDEF
    R3BEventDisplayLinkDef.h
    HEADERS
    ${HEADERS}
    SRCS
    ${SRCS}
    INCLUDEDIRS
    ${CMAKE_CURRENT_SOURCE_DIR}
    DEPENDENCIES
    R3BCalifa)

# set(DEPENDENCIES
#     Field
#     R3BData
#     R3BCalifa
#     Geane
#     EventDisplay
#     Eve
#     Gui
#     Ged
#     RGL
#     Graf3d
#     Geom
#     GeomPainter)
