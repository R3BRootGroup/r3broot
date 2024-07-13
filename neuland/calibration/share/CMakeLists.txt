# cmake-format: off
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
# cmake-format: on

set(SRCS # cmake-format: sortable
         R3BNeulandBasePar.cxx R3BNeulandCalibrationTask.cxx R3BNeulandParSet.cxx)

set(HEADERS
    # cmake-format: sortable
    R3BNeulandBasePar.h
    R3BNeulandCalibrationTask.h
    R3BNeulandParDirCreator.h
    R3BNeulandParSet.h
    R3BNeulandTriggerTypes.h)

add_library(R3BNeulandCalibrationShare SHARED ${SRCS})
target_include_directories(R3BNeulandCalibrationShare PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(
    R3BNeulandCalibrationShare
    PUBLIC R3BSource R3BNeulandShared
    PRIVATE Microsoft.GSL::GSL range-v3::range-v3)
root_generate_dictionary(
    R3BNeulandCalibrationShare_dict
    ${HEADERS}
    MODULE
    R3BNeulandCalibrationShare
    LINKDEF
    NeulandCalShareLinkDef.h)
