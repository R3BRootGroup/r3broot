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

include(CMakePackageConfigHelpers)

# Configure/Install CMake package
function(install_cmake_package)
  # Correct spelling, in local scope
  set(PROJECT_NAME R3BRoot)

  set(PACKAGE_INSTALL_DESTINATION
    ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}-${PROJECT_VERSION})
  if(PROJECT_EXPORT_SET)
    install(EXPORT ${PROJECT_EXPORT_SET}
      NAMESPACE ${PROJECT_NAME}::
      DESTINATION ${PACKAGE_INSTALL_DESTINATION}
      EXPORT_LINK_INTERFACE_LIBRARIES)
  endif()
  configure_package_config_file(
    "cmake/private/${PROJECT_NAME}Config.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION "${PACKAGE_INSTALL_DESTINATION}"
    PATH_VARS
      CMAKE_INSTALL_PREFIX
      CMAKE_INSTALL_FULL_BINDIR
      PACKAGE_INSTALL_DESTINATION)
  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    VERSION "${PROJECT_VERSION}"
    COMPATIBILITY AnyNewerVersion)
  install(
    FILES
      "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
    DESTINATION "${PACKAGE_INSTALL_DESTINATION}")
endfunction()
