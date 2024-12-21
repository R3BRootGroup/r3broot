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

#cmake_minimum_required(VERSION 3.0)

#enable_testing()
#set(PROJECT_TEST_NAME AtimaUnitTests)
#set(GTEST_ROOT ${SIMPATH})
#find_package(GTest)

#if(GTEST_FOUND)
#file(GLOB TEST_SRC_FILES ${PROJECT_SOURCE_DIR}/atima/test/*.cxx)

#include_directories(${GTEST_INCLUDE_DIRS}
#                    ${SYSTEM_INCLUDE_DIRECTORIES}
#                    ${BASE_INCLUDE_DIRECTORIES}
#                    ${R3BROOT_SOURCE_DIR}/atima)

#link_directories(${GTEST_LIBS_DIR}
#                 ${ROOT_LIBRARY_DIR}
#                 ${FAIRROOT_LIBRARY_DIR}
#                 ${Boost_LIBRARY_DIRS})

#set(TEST_DEPENDENCIES
#    ${GTEST_BOTH_LIBRARIES}
#    ${ROOT_LIBRARIES}
#    R3BAtima)

#add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
#target_link_libraries(${PROJECT_TEST_NAME} ${TEST_DEPENDENCIES})
#add_test(${PROJECT_TEST_NAME} ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_TEST_NAME})
#endif(GTEST_FOUND)
