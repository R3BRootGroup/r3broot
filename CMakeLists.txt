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

# The name of our project is “R3BROOT”.  CMakeLists files in this project can refer to the root
# source directory of the project as ${R3BROOT_SOURCE_DIR} or as ${CMAKE_SOURCE_DIR} and to the root
# binary directory of the project as ${R3BROOT_BINARY_DIR} or ${CMAKE_BINARY_DIR}. This difference
# is is important for the base classes which are in R3BROOT and PANDAROOT.

if(CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "In-source builds are not supported, please, create a folder build.")
endif()

# Check if cmake has the required version
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# set cmake policies
cmake_policy(VERSION 3.15...3.31)

# Set name of our project to "R3BROOT". Has to be done after check of cmake version
project(R3BROOT LANGUAGES C CXX Fortran)

# use default gnu install convention:
include(GNUInstallDirs)

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/macros")
    message(STATUS "  ")
    message(
        ERROR
        " macros repository is missing in the R3BRoot source directory\ncd ${CMAKE_SOURCE_DIR}\ngit clone -b dev https://github.com/R3BRootGroup/macros.git"
    )
    message(STATUS "  ")
endif(NOT EXISTS "${PROJECT_SOURCE_DIR}/macros")

# Check for needed environment variables
if(NOT DEFINED ENV{SIMPATH})
    message(
        FATAL_ERROR
            "You did not define the environment variable SIMPATH which is nedded to find the external packages. Please set this variable and execute cmake again."
    )
endif(NOT DEFINED ENV{SIMPATH})

if(NOT DEFINED ENV{FAIRROOTPATH})
    message(
        FATAL_ERROR
            "You did not define the environment variable FAIRROOTPATH which is needed to find FairRoot. Please set this variable and execute cmake again."
    )
endif()

set(SIMPATH $ENV{SIMPATH})
set(FAIRROOTPATH $ENV{FAIRROOTPATH})

include(cmake/utils/create_library_with_dictionary.cmake)
include(cmake/scripts/set_ld_library_dir.cmake)
set(GRAPHVIZ_EXTERNAL_LIBS FALSE)

list(APPEND CMAKE_PREFIX_PATH ${SIMPATH})
list(APPEND CMAKE_PREFIX_PATH ${FAIRROOTPATH})

list(APPEND CMAKE_MODULE_PATH "${FAIRROOTPATH}/share/fairbase/cmake/modules")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/scripts")

# For CMakeLists that still use the old way
set(ROOT_FOUND_VERSION 600000)
set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")
set(FairRoot_ROOT "${FAIRROOTPATH}")

set(CheckSrcDir "${FAIRROOTPATH}/share/fairbase/cmake/checks")

# Set option(s)
option(BUILD_GEOMETRY "Build all ROOT geometries" ON)
option(USE_PATH_INFO "Information from PATH and LD_LIBRARY_PATH are used." OFF)
option(WITH_UCESB "Build UCESB" OFF)
option(WITH_EPICS "Build EPICS" OFF)
option(WITH_ACTAF "Build ACTAF" OFF)
option(WITH_SOFIA "Build SOFIA" OFF)
option(WITH_GTPC "Build GTPC" OFF)
option(WITH_FRS "Build FRS" OFF)
option(WITH_ASYEOS "Build ASYEOS" OFF)
option(BUILD_C3W "Build C3W" OFF)
option(CONAN_ENABLED "conan package is enabled. (Use CMakePreset to set it)" OFF)
option(FETCH_GTEST "Fetch GTest from CMake" OFF)

find_package(FairCMakeModules CONFIG REQUIRED)
if(NOT FairCMakeModules_FOUND)
    message(FATAL_ERROR "Could not find FairCMakeModules! ")
endif()

include(FairFindPackage2)
include(FairFormattedOutput)
include(FetchAllContents)

find_package2(PUBLIC FairRoot CONFIG REQUIRED)

# FairRoot_VERSION is eventually set by find_package, for now let's set it explicitly
execute_process(
    COMMAND ${FAIRROOTPATH}/bin/fairroot-config --version
    COMMAND sed "s/^v//"
    OUTPUT_VARIABLE FairRoot_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
if("${FairRoot_VERSION}" VERSION_LESS 18.2.1)
    message(
        FATAL_ERROR
            "Incompatible FairRoot version detected: ${FairRoot_VERSION}. R3BRoot requires FairRoot v18.2.1 or later"
    )
endif()

# added a variable to point to r3broot path in the system
set(R3BROOTPATH ${CMAKE_INSTALL_PREFIX})

# Load some basic macros which are needed later on
include(FairMacros)

set(required_root_components
    RIO
    Geom
    Physics
    Matrix
    Hist
    Graf
    Gpad
    Net
    Imt
    Core
    MathCore
    RHTTP
    Spectrum)
list(APPEND optional_root_components opengl)
list(APPEND optional_root_components gdml vmc)

find_package2(
    PUBLIC
    ROOT
    CONFIG
    REQUIRED
    COMPONENTS
    ${required_root_components}
    OPTIONAL_COMPONENTS
    ${optional_root_components})

if("${ROOT_vmc_FOUND}" MATCHES "no")
    find_package2(PUBLIC VMC REQUIRED)
    set(VMCLIB VMCLibrary)
else()
    set(VMCLIB VMC)
endif()

include(R3BMacros)
include(WriteConfigFile_r3b)
include(WriteEnvCache)
include(CheckCompiler)

find_package2(
    PUBLIC
    FairLogger
    VERSION
    1.6.0
    REQUIRED)

foreach(dep IN LISTS FairLogger_PACKAGE_DEPENDENCIES)
    if(NOT dep STREQUAL "Boost")
        find_package2(
            PUBLIC
            ${dep}
            REQUIRED
            VERSION
            ${FairLogger_${dep}_VERSION})
    endif()
endforeach()

if("${FairRoot_VERSION}" VERSION_GREATER 18.6.4)

    find_package2(
        PUBLIC
        FairMQ
        VERSION
        1.4.0
        REQUIRED)
    find_package2(PUBLIC DDS)

    if(FAIRSOFT_EXTERN)
        set(Boost_ROOT ${SIMPATH})
        set(GSL_DIR ${SIMPATH})
    else(FAIRSOFT_EXTERN)
        set(Boost_ROOT ${SIMPATH}/basics/boost)
        set(GSL_DIR ${SIMPATH}/basics/gsl)
    endif(FAIRSOFT_EXTERN)

    if(NOT DEFINED Boost_NO_SYSTEM_PATHS)
        set(Boost_NO_SYSTEM_PATHS TRUE)
    endif()
    if(NOT DEFINED Boost_NO_BOOST_CMAKE AND CMAKE_VERSION VERSION_LESS 3.15)
        set(Boost_NO_BOOST_CMAKE TRUE)
    endif()
    if(Boost_NO_BOOST_CMAKE)
        # If an older version of boost is found both of the variables below are
        # cached and in a second cmake run, a good boost version is found even
        # if the version is to old.
        # To overcome this problem both variables are cleared before checking
        # for boost.
        unset(Boost_INCLUDE_DIR CACHE)
        unset(Boost_LIBRARY_DIRS CACHE)
    endif()

    list(
        APPEND
        boost_dependencies
        thread
        system
        timer
        program_options
        random
        filesystem
        chrono
        exception
        iostreams
        regex
        serialization
        log
        atomic)

    find_package2(
        PUBLIC
        Boost
        VERSION
        1.67
        COMPONENTS
        ${boost_dependencies}
        ADD_REQUIREMENTS_OF
        FairMQ)
    #  If (FairMQ_FOUND AND Boost_FOUND)
    #    add_subdirectory(base/MQ)
    #    add_subdirectory(parmq)
    #  Else ()
    #    Message(STATUS "base/MQ and parmq will not be built, because FairMQ and/or Boost was not found.")
    #  EndIf ()
else()
    set(Boost_NO_SYSTEM_PATHS TRUE)
    set(Boost_NO_BOOST_CMAKE TRUE)

    if(${ROOT_LIBRARY_DIR} MATCHES /lib/root)
        set(Boost_ROOT ${SIMPATH})
        set(GSL_DIR ${SIMPATH})
    else(${ROOT_LIBRARY_DIR} MATCHES /lib/root)
        set(Boost_ROOT ${SIMPATH}/basics/boost)
        set(GSL_DIR ${SIMPATH}/basics/gsl)
    endif(${ROOT_LIBRARY_DIR} MATCHES /lib/root)

    unset(Boost_INCLUDE_DIR CACHE)
    unset(Boost_LIBRARY_DIRS CACHE)
    find_package2(
        PUBLIC
        Boost
        VERSION
        1.67
        ${FairMQ_Boost_VERSION}
        COMPONENTS
        thread
        system
        timer
        program_options
        random
        filesystem
        chrono
        exception
        iostreams
        regex
        serialization
        log
        log_setup
        atomic
        date_time
        ${FairMQ_Boost_COMPONENTS}
        REQUIRED)
endif()

if(CONAN_ENABLED)
    find_package(range-v3 REQUIRED)
    find_package(Microsoft.GSL REQUIRED)
    find_package(nlohmann_json REQUIRED)
    set(JSON_ImplicitConversions
        OFF
        CACHE BOOL "disable implicit conversion of JSON" FORCE)
endif()

find_package2(PRIVATE GTest)
find_package2(PUBLIC Pythia6)
find_package2(PUBLIC Pythia8)
find_package2(PUBLIC Protobuf)
find_package2(PUBLIC msgpack)
find_package2(PUBLIC Flatbuffers)

find_package2(PUBLIC Geant3)
find_package2(PUBLIC Geant4)

find_package2(PUBLIC VGM)

find_package2(PUBLIC Geant4VMC)

find_package2(PUBLIC IWYU)

find_package2(PUBLIC GSL)

find_package2(PUBLIC Git)
set(Git_VERSION ${GIT_VERSION_STRING})
string(REPLACE "git" " " Git_ROOT ${GIT_EXECUTABLE})

# Check the compiler and set the compile and link flags
check_compiler()

if(DEFINED ENV{WERROR} AND "$ENV{WERROR}")
    message(STATUS "Will compile with -Werror. ")
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -Wall -Wformat-security -Werror -Wno-sign-compare -Wno-reorder -Wno-unused-variable -Wno-unused-but-set-variable -Wfatal-errors"
    )
else()
    if(NOT APPLE)
        message(
            STATUS
                "Set env WERROR to 1 to enable -Werror. If origin/dev compiles on your platform with that option, it is definitly a good idea to do that."
        )
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -Wall -Wformat-security -Wno-sign-compare -Wno-reorder -Wno-unused-variable -Wno-unused-but-set-variable"
        )
    endif()
endif()
set(CMAKE_CXX_FLAGS "-D_GLIBCXX_ASSERTIONS ${CMAKE_CXX_FLAGS}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # used by clang-tidy

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(VMCWORKDIR ${CMAKE_SOURCE_DIR})

# Display SUPPORTED_CXX_STANDARDS
set(SUPPORTED_CXX_STANDARDS "")
foreach(version 14 17 20 23)
    string(CONCAT flag "cxx_std_" ${version})
    if(flag IN_LIST CMAKE_CXX_COMPILE_FEATURES)
        list(APPEND SUPPORTED_CXX_STANDARDS ${version})
    endif()
endforeach()

if(SUPPORTED_CXX_STANDARDS)
    message(STATUS "${Cyan}Supported C++ standards:${CR} ${SUPPORTED_CXX_STANDARDS}")
    message(STATUS "${Cyan}CMake is using C++ standard:${CR} ${BGreen}${CMAKE_CXX_STANDARD}${CR}")
else()
    message(WARNING "No supported C++ standards found!")
endif()

# Display configurations
message(STATUS "CMake C++ Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Library Output Directory: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS "Runtime Output Directory: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

if(USE_PATH_INFO)
    if(APPLE)
        set(LD_LIBRARY_PATH $ENV{DYLD_LIBRARY_PATH})
    else(APPLE)
        set(LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})
    endif(APPLE)
endif()

# Check if the user wants to build the project in the source directory
check_out_of_source_build()

# Check if we are on an UNIX system. If not stop with an error message
if(NOT UNIX)
    message(
        FATAL_ERROR
            "You're not on an UNIX system. The project was up to now only tested on UNIX systems, so we break here. If you want to go on please edit the CMakeLists.txt in the source directory."
    )
endif(NOT UNIX)

# Check if the external packages are installed into a separate install directory
if(COMMAND check_external_package_install_dir)
    check_external_package_install_dir()
endif()

if(Boost_FOUND)
    set(Boost_Avail 1)
else(Boost_FOUND)
    set(Boost_Avail 0)
endif(Boost_FOUND)

find_package2(PUBLIC yaml-cpp)
# Workaround missing exported include directories
# Upstream has fixed this in https://github.com/jbeder/yaml-cpp/commit/ab5f9259a4e67d3fe0e4832bd407a9e596e2d884 (since 0.6.3)
if(yaml-cpp_FOUND)
    get_filename_component(
        YAML_CPP_INCLUDE_DIR
        "${YAML_CPP_INCLUDE_DIR}"
        REALPATH
        BASE_DIR
        "/")
endif()
if(yaml-cpp_FOUND
   AND TARGET yaml-cpp
   AND EXISTS YAML_CPP_INCLUDE_DIR
   AND yaml-cpp_VERSION VERSION_LESS 0.6.3)
    set_target_properties(yaml-cpp PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                              "${YAML_CPP_INCLUDE_DIR}")
endif()

# include(${CMAKE_SOURCE_DIR}/cmake/scripts/fetchGTest.cmake)

setbasicvariables()

# Set the library version in the main CMakeLists.txt
set(R3BROOT_MAJOR_VERSION 0)
set(R3BROOT_MINOR_VERSION 0)
set(R3BROOT_PATCH_VERSION 0)
set(R3BROOT_VERSION "${FairRoot_VERSION}")

if(NOT DEFINED ROOT_FOUND_VERSION)
    set(ROOT_FOUND_VERSION "${found_vers}")
endif()
message("ROOT_FOUND_VERSION = ${ROOT_FOUND_VERSION}")

set(FAIRROOT_LIBRARY_PROPERTIES
    ${FAIRROOT_LIBRARY_PROPERTIES}
    VERSION
    "${R3BROOT_VERSION}"
    SOVERSION
    "${R3BROOT_VERSION}")

string(COMPARE EQUAL "${CMAKE_GENERATOR}" "Xcode" IS_XCODE)
set(R3BLIBDIR ${CMAKE_BINARY_DIR}/lib)
if(${IS_XCODE})
    set(LD_LIBRARY_PATH ${R3BLIBDIR}/RELWITHDEBINFO ${R3BLIBDIR} ${LD_LIBRARY_PATH})
else(${IS_XCODE})
    set(LD_LIBRARY_PATH ${R3BLIBDIR} ${LD_LIBRARY_PATH})
endif(${IS_XCODE})

find_package(ucesb)
if(ucesb_FOUND)
    set(WITH_UCESB ON)
    set(UCESBPATH $ENV{UCESB_DIR})
    set(ucesb_summary "${BGreen}YES${CR}       ${BMagenta}${UCESBPATH}${CR}")
    set(LD_LIBRARY_PATH ${ucesb_LIBRARY_DIR} ${LD_LIBRARY_PATH})
else(ucesb_FOUND)
    message(
        WARNING
            "${BMagenta}No UCESB configured. You will not be able to unpack lmd files. Set UCESB_DIR"
            " to a ucesb installation to fix this.${CR}")
endif(ucesb_FOUND)

find_package(EPICS)
if(EPICS_FOUND)
    set(WITH_EPICS ON)
    set(epics_summary "${BGreen}YES${CR}       ${BMagenta}${EPICS_LIBRARY_DIR}${CR}")
    set(LD_LIBRARY_PATH ${EPICS_LIBRARY_DIR} ${LD_LIBRARY_PATH})
else(EPICS_FOUND)
    set(epics_summary "${BYellow}NO${CR}")
endif(EPICS_FOUND)

find_package(Atima)
if(Atima_FOUND)
    set(atima_summary "${BGreen}YES${CR}       ${BMagenta}${R3BROOT_SOURCE_DIR}/atima/${CR}")
else(Atima_FOUND)
    set(atima_summary "${BYellow}NO${CR}")
endif(Atima_FOUND)

find_package(Garfield)

if(GARFIELD_FOUND)
    message("Building ACTAF framework...")
    set(WITH_ACTAF ON)
    set(actaf_summary "${BGreen}YES${CR}       ${BMagenta}${R3BROOT_SOURCE_DIR}/actaf/${CR}")
    set(LD_LIBRARY_PATH ${GARFIELD_LIBRARY_DIR} ${LD_LIBRARY_PATH})
else(GARFIELD_FOUND)
    set(actaf_summary "${BYellow}NO${CR}")
endif(GARFIELD_FOUND)

if(EXISTS "${R3BROOT_SOURCE_DIR}/sofia/")
    message("Building SOFIA framework...")
    set(WITH_SOFIA ON)
    add_definitions(-DSOFIA)
    set(sofia_summary "${BGreen}YES${CR}       ${BMagenta}${R3BROOT_SOURCE_DIR}/sofia/${CR}")
else(EXISTS "${R3BROOT_SOURCE_DIR}/sofia/")
    message("${BYellow}SOFIA framework not found${CR}")
    set(sofia_summary "${BYellow}NO${CR}")
endif(EXISTS "${R3BROOT_SOURCE_DIR}/sofia/")

if(EXISTS "${R3BROOT_SOURCE_DIR}/NeuLAND_DNN/")
    message("Building NeuLAND DNN framework...")
    set(WITH_NeuLAND_DNN ON)
    set(neuland_dnn_summary
        "${BGreen}YES${CR}       ${BMagenta}${R3BROOT_SOURCE_DIR}/NeuLAND_DNN/${CR}")
    add_definitions(-DNEULAND_DNN)
else(EXISTS "${R3BROOT_SOURCE_DIR}/NeuLAND_DNN/")
    message("${BYellow}NeuLAND_DNN framework not found${CR}")
    set(neuland_dnn_summary "${BYellow}NO${CR}")
endif(EXISTS "${R3BROOT_SOURCE_DIR}/NeuLAND_DNN/")

if(EXISTS "${R3BROOT_SOURCE_DIR}/glad-tpc/")
    message("Building GLAD-TPC detector...")
    set(WITH_GTPC ON)
    set(gladtpc_summary "${BGreen}YES${CR}       ${BMagenta}${R3BROOT_SOURCE_DIR}/glad-tpc/${CR}")
    add_definitions(-DGTPC)
else(EXISTS "${R3BROOT_SOURCE_DIR}/glad-tpc/")
    message("${BYellow}GLAD-TPC framework not found${CR}")
    set(gladtpc_summary "${BYellow}NO${CR}")
endif(EXISTS "${R3BROOT_SOURCE_DIR}/glad-tpc/")

if(EXISTS "${R3BROOT_SOURCE_DIR}/frs/")
    message("Building FRS framework...")
    set(WITH_FRS ON)
    set(frs_summary "${BGreen}YES${CR}       ${BMagenta}${R3BROOT_SOURCE_DIR}/frs/${CR}")
    add_definitions(-DFRS)
else(EXISTS "${R3BROOT_SOURCE_DIR}/frs/")
    message("${BYellow}FRS framework not found${CR}")
    set(frs_summary "${BYellow}NO${CR}")
endif(EXISTS "${R3BROOT_SOURCE_DIR}/frs/")

if(EXISTS "${R3BROOT_SOURCE_DIR}/asyeos/")
    message("Building ASYEOS framework...")
    set(WITH_ASYEOS ON)
    set(asyeos_summary "${BGreen}YES${CR}       ${BMagenta}${R3BROOT_SOURCE_DIR}/asyeos/${CR}")
    add_definitions(-DASYEOS)
else(EXISTS "${R3BROOT_SOURCE_DIR}/asyeos/")
    message("${BYellow}ASYEOS framework not found${CR}")
    set(asyeos_summary "${BYellow}NO${CR}")
endif(EXISTS "${R3BROOT_SOURCE_DIR}/asyeos/")

configure_file(${CMAKE_SOURCE_DIR}/cmake/scripts/R3BRootConfigVersion.cmake.in
               "${CMAKE_BINARY_DIR}/R3BRootConfigVersion.cmake" @ONLY)

r3b_generate_version_info()

# Recurse into the given subdirectories.  This does not actually cause another cmake executable to
# run.  The same process will walk through the project's entire directory structure.

# FIXME in FairSoft / FairRoot
find_package(VMC QUIET)
set(BASE_INCLUDE_DIRECTORIES ${BASE_INCLUDE_DIRECTORIES} ${VMC_INCLUDE_DIRS})
set(ROOT_CINT_EXECUTABLE ${ROOT_rootcling_CMD})

include(CTest)
include(RootCTest)

if(MODULE)
    add_subdirectory(${MODULE})
endif(MODULE)

if(NOT MODULE)
    if(BUILD_GEOMETRY)
        set(ROOTSYS ${SIMPATH})
        set(unitgeo_summary
            "${BGreen} YES${CR}       (default, disable with ${BMagenta}-DBUILD_GEOMETRY=OFF${CR})")
    else()
        set(unitgeo_summary
            "${BRed} NO${CR}        (enable with ${BMagenta}-DBUILD_GEOMETRY=ON${CR})")
    endif()
    message(STATUS "  ${BWhite}ROOT geometries${CR}          ${unitgeo_summary}")
    add_subdirectory(r3bbase)
    add_subdirectory(r3bdata)
    add_subdirectory(field)
    add_subdirectory(tracking)
    add_subdirectory(passive)
    add_subdirectory(tcal)
    if(WITH_UCESB)
        add_subdirectory(r3bsource)
    endif()
    add_subdirectory(neuland)
    add_subdirectory(glad)
    add_subdirectory(dch)
    add_subdirectory(fiber)
    add_subdirectory(gfi)
    add_subdirectory(land)
    add_subdirectory(los)
    add_subdirectory(rolu)
    add_subdirectory(frssci)
    add_subdirectory(sci2)
    add_subdirectory(sci8)
    add_subdirectory(sfib)
    add_subdirectory(mtof)
    add_subdirectory(tof)
    add_subdirectory(tofi)
    add_subdirectory(tofd)
    add_subdirectory(ssd)
    add_subdirectory(tttx)
    add_subdirectory(strawtubes)
    add_subdirectory(music)
    add_subdirectory(r3bgen)
    add_subdirectory(califa)
    add_subdirectory(evtvis)
    add_subdirectory(psp)
    add_subdirectory(compilehelper)
    add_subdirectory(tracker_rene)
    add_subdirectory(pdc)
    add_subdirectory(rpc)
    add_subdirectory(alpide)
    add_subdirectory(mwpc)
    add_subdirectory(twim)
    add_subdirectory(musli)
    add_subdirectory(analysis)
    if(WITH_SOFIA)
        add_subdirectory(sofia)
    endif(WITH_SOFIA)
    if(WITH_GTPC)
        add_subdirectory(glad-tpc)
    endif(WITH_GTPC)
    if(WITH_EPICS)
        add_subdirectory(epics)
    endif(WITH_EPICS)
    if(WITH_FRS)
        add_subdirectory(frs)
    endif(WITH_FRS)
    if(WITH_ASYEOS)
        add_subdirectory(asyeos)
    endif(WITH_ASYEOS)

    #================================================================================
    # TODO: Since other dependencies still use the old cmake configurations, include(ROOTMacros) has to be used. Please change them using CMake targets ASAP.
    include(ROOTMacros)
    if(Atima_FOUND)
        add_subdirectory(atima)
    endif(Atima_FOUND)
    if(WITH_ACTAF)
        add_subdirectory(actaf)
    endif(WITH_ACTAF)
    if(EXISTS "${PROJECT_SOURCE_DIR}/R3BFileSource")
        add_subdirectory(R3BFileSource)
    endif(EXISTS "${PROJECT_SOURCE_DIR}/R3BFileSource")
    if(EXISTS "${PROJECT_SOURCE_DIR}/macros")
        add_subdirectory(macros)
    endif(EXISTS "${PROJECT_SOURCE_DIR}/macros")
    #================================================================================

endif(NOT MODULE)

# add_subdirectory(gconfig)
add_subdirectory(cmake)

if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    execute_process(
        COMMAND ${GIT_EXECUTABLE} describe --tags
        OUTPUT_VARIABLE R3BROOT_GIT_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
else()
    set(R3BROOT_GIT_VERSION
        v${R3BROOT_MAJOR_VERSION}.${R3BROOT_MINOR_VERSION}.${R3BROOT_PATCH_VERSION})
endif()

write_config_file(config.sh)
write_config_file(config.csh)
write_config_file(config.sh_install)
write_config_file(config.csh_install)

write_env_cache(env_cache.sh)

# Get rid of leftover c3w configuration
execute_process(COMMAND rm -f ${CMAKE_CURRENT_BINARY_DIR}/GNUmakefile
                        ${CMAKE_CURRENT_BINARY_DIR}/c3w.sh)

if(EXISTS ${CMAKE_SOURCE_DIR}/c3w.conf)
    set(BUILD_C3W ON)
    set(c3w_summary "${BGreen}YES${CR}       ${BMagenta}${R3BROOT_SOURCE_DIR}/util/c3w/${CR}")
    add_subdirectory(util/c3w)
else()
    set(c3w_summary "${BYellow}NO${CR}")
endif()

# Summary ######################################################################
if(CMAKE_CXX_FLAGS)
    message(STATUS "  ")
    message(STATUS "  ${Cyan}GLOBAL CXX FLAGS${CR}  ${BGreen}${CMAKE_CXX_FLAGS}${CR}")
endif()
if(CMAKE_CONFIGURATION_TYPES)
    message(STATUS "  ")
    message(STATUS "  ${Cyan}BUILD TYPE         CXX FLAGS${CR}")
    string(TOUPPER "${CMAKE_BUILD_TYPE}" selected_type)
    foreach(type IN LISTS CMAKE_CONFIGURATION_TYPES)
        string(TOUPPER "${type}" type_upper)
        if(type_upper STREQUAL selected_type)
            fair_pad("${type}" 19 " " type_padded)
            message(STATUS "${BGreen}* ${type_padded}${CMAKE_CXX_FLAGS_${type_upper}}${CR}")
        else()
            fair_pad("${type}" 19 " " type_padded)
            message(STATUS "  ${BWhite}${type_padded}${CR}${CMAKE_CXX_FLAGS_${type_upper}}")
        endif()
        unset(type_padded)
        unset(type_upper)
    endforeach()
    message(STATUS "  ")
    message(STATUS "  (Change the build type with ${BMagenta}-DCMAKE_BUILD_TYPE=...${CR})")
endif()

if(PROJECT_PACKAGE_DEPENDENCIES)
    message(STATUS "  ")
    message(STATUS "  ${Cyan}DEPENDENCY FOUND     VERSION                   PREFIX${CR}")
    foreach(dep IN LISTS PROJECT_PACKAGE_DEPENDENCIES)
        if(${dep}_VERSION)
            if(${dep} STREQUAL DDS)
                set(version_str "${BGreen}${${dep}_MAJOR_VERSION}.${${dep}_MINOR_VERSION}${CR}")
            elseif(${dep} STREQUAL Boost)
                if(Boost_VERSION_MAJOR)
                    set(version_str "${BGreen}${${dep}_VERSION_MAJOR}.${${dep}_VERSION_MINOR}${CR}")
                else()
                    set(version_str "${BGreen}${${dep}_MAJOR_VERSION}.${${dep}_MINOR_VERSION}${CR}")
                endif()
            else()
                set(version_str "${BGreen}${${dep}_VERSION}${CR}")
            endif()
        else()
            set(version_str "${BYellow}unknown${CR}")
        endif()
        if(PROJECT_${dep}_VERSION)
            set(version_req_str " (>= ${PROJECT_${dep}_VERSION})")
        endif()
        fair_pad(${dep} 21 " " dep_padded)
        if(DISABLE_COLOR)
            fair_pad("${version_str}${version_req_str}" 26 " " version_padded)
        else()
            fair_pad(
                "${version_str}${version_req_str}"
                26
                " "
                version_padded
                COLOR
                1)
        endif()
        set(prefix ${${dep}_PREFIX})
        if(${dep} STREQUAL Boost)
            if(TARGET Boost::headers)
                get_target_property(boost_include Boost::headers INTERFACE_INCLUDE_DIRECTORIES)
            else()
                get_target_property(boost_include Boost::boost INTERFACE_INCLUDE_DIRECTORIES)
            endif()
            get_filename_component(prefix ${boost_include}/.. ABSOLUTE)
        elseif(${dep} STREQUAL Protobuf)
            get_filename_component(prefix ${Protobuf_INCLUDE_DIRS}/.. ABSOLUTE)
        elseif(${dep} STREQUAL msgpack)
            get_target_property(msgpack_include msgpackc INTERFACE_INCLUDE_DIRECTORIES)
            get_filename_component(prefix ${msgpack_include}/.. ABSOLUTE)
        elseif(${dep} STREQUAL Pythia6)
            get_filename_component(prefix ${PYTHIA6_LIBRARY_DIR}/.. ABSOLUTE)
        elseif(${dep} STREQUAL Pythia8)
            get_filename_component(prefix ${PYTHIA8_LIB_DIR}/.. ABSOLUTE)
        elseif(${dep} STREQUAL FairLogger)
            if(FairLogger_PREFIX)
                set(prefix ${FairLogger_PREFIX})
            else()
                set(prefix ${FairLogger_ROOT})
            endif()
        elseif(${dep} STREQUAL FairMQ)
            if(FairMQ_PREFIX)
                set(prefix ${FairMQ_PREFIX})
            else()
                set(prefix ${FairMQ_ROOT})
            endif()
        elseif(${dep} STREQUAL yaml-cpp)
            get_filename_component(prefix ${YAML_CPP_INCLUDE_DIR}/.. ABSOLUTE)
        elseif(${dep} STREQUAL Geant4VMC)
            string(REPLACE ":" ";" geant4vmc_include ${Geant4VMC_INCLUDE_DIRS})
            list(GET geant4vmc_include 0 geant4vmc_include)
            get_filename_component(prefix ${geant4vmc_include}/../.. ABSOLUTE)
        elseif(${dep} STREQUAL Geant3)
            string(REPLACE ":" ";" geant3_include ${Geant3_INCLUDE_DIRS})
            list(GET geant3_include 0 geant3_include)
            get_filename_component(prefix ${geant3_include}/../.. ABSOLUTE)
        elseif(${dep} STREQUAL Geant4)
            list(GET Geant4_INCLUDE_DIRS 0 geant4_include)
            get_filename_component(prefix ${geant4_include}/../.. ABSOLUTE)
        elseif(${dep} STREQUAL VGM)
            string(REPLACE ":" ";" vgm_include ${VGM_INCLUDE_DIRS})
            list(GET vgm_include 0 vgm_include)
            get_filename_component(prefix ${vgm_include}/.. ABSOLUTE)
        elseif(${dep} STREQUAL ROOT)
            set(prefix ${ROOT_INSTALL_DIR})
        elseif(${dep} STREQUAL IWYU)
            get_filename_component(prefix ${IWYU_BINARY}/.. ABSOLUTE)
        elseif(${dep} STREQUAL yaml-cpp)
            get_filename_component(prefix ${YAML_CPP_INCLUDE_DIR}/.. ABSOLUTE)
        endif()

        message(STATUS "  ${BWhite}${dep_padded}${CR}${version_padded}${prefix}")

        #    if(${dep} STREQUAL Geant3)
        #      message(STATUS "                                          G3SYS: ${Geant3_SYSTEM_DIR}")
        #    endif()
        if(${dep} STREQUAL Geant4)
            foreach(dataset IN LISTS Geant4_DATASETS)
                fair_pad(
                    ${Geant4_DATASET_${dataset}_ENVVAR}
                    21
                    " "
                    envvar_padded
                    LEFT)
                message(
                    STATUS
                        "                          ${envvar_padded}: ${Geant4_DATASET_${dataset}_PATH}"
                )
            endforeach()
        endif()

        unset(version_str)
        unset(version_padded)
        unset(version_req_str)
    endforeach()
endif()

# Summary frameworks ###########################################################
message(STATUS "  ")
message(STATUS "  ${Cyan}Building framework   BUILT?    BASEPATH${CR}")
message(STATUS "  ${BWhite}UCESB${CR}                ${ucesb_summary}")
message(STATUS "  ${BWhite}ATIMA${CR}                ${atima_summary}")
message(STATUS "  ${BWhite}EPICS${CR}                ${epics_summary}")
message(STATUS "  ${BWhite}ACTAF${CR}                ${actaf_summary}")
message(STATUS "  ${BWhite}SOFIA${CR}                ${sofia_summary}")
message(STATUS "  ${BWhite}FRS${CR}                  ${frs_summary}")
message(STATUS "  ${BWhite}NeuLAND_DNN${CR}          ${neuland_dnn_summary}")
message(STATUS "  ${BWhite}GLAD-TPC${CR}             ${gladtpc_summary}")
message(STATUS "  ${BWhite}ASYEOS${CR}               ${asyeos_summary}")
message(STATUS "  ${BWhite}C3W${CR}                  ${c3w_summary}")

message(STATUS "  ")
message(STATUS "  ${Cyan}COMPONENT            BUILT?    INFO${CR}")
message(STATUS "  ${BWhite}ROOT geometries${CR}     ${unitgeo_summary}")
message(STATUS "  ")
################################################################################
