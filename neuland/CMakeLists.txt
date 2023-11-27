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

add_subdirectory(calibration)
add_subdirectory(clustering)
add_subdirectory(digitizing)
add_subdirectory(online)
if(WITH_EPICS)
    add_subdirectory(preexp)
endif(WITH_EPICS)
add_subdirectory(reconstruction)
add_subdirectory(shared)
add_subdirectory(simulation)
add_subdirectory(test)
add_subdirectory(unpack)
add_subdirectory(executables)
if(BUILD_GEOMETRY)
    add_subdirectory(geobase)
endif()
