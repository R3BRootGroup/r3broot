################################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH      #
#   Copyright (C) 2019-2025 Members of R3B Collaboration                       #
#                                                                              #
#             This software is distributed under the terms of the              #
#                 GNU General Public Licence (GPL) version 3,                  #
#                    copied verbatim in the file "LICENSE".                    #
#                                                                              #
# In applying this license GSI does not waive the privileges and immunities    #
# granted to it by virtue of its status as an Intergovernmental Organization   #
# or submit itself to any jurisdiction.                                        #
################################################################################
#
# Plox to install EPICS base and set EPICS_BASE and EPICS_HOST_ARCH.
#

message(STATUS "Looking for EPICS (base)...")

find_path(
    EPICS_INCLUDE_DIR
    NAMES cadef.h
    PATHS $ENV{EPICS_BASE}/include
    NO_DEFAULT_PATH
    DOC "Path to EPICS include header files.")

if(EPICS_INCLUDE_DIR)
    set(EPICS_HOST_ARCH $ENV{EPICS_HOST_ARCH})
    if(NOT EPICS_HOST_ARCH)
        set(EPICS_HOST_ARCH linux-x86_64)
    endif(NOT EPICS_HOST_ARCH)

    # Extract ARCH_CLASS from EPICS config/* and look for OS-specific headers.
    file(STRINGS $ENV{EPICS_BASE}/config/CONFIG.Host.${EPICS_HOST_ARCH} EPICS_ARCH_CLASS
         REGEX "ARCH_CLASS[ \t]*=.*")
    string(REGEX REPLACE "[^=]*=[ \t]*(.*)" "\\1" EPICS_ARCH_CLASS ${EPICS_ARCH_CLASS})
    find_path(
        EPICS_OS_INCLUDE_DIR
        NAMES osdEvent.h
        PATHS $ENV{EPICS_BASE}/include/os/${EPICS_ARCH_CLASS}
        NO_DEFAULT_PATH
        DOC "Path to EPICS OS include header files.")

    find_path(
        EPICS_LIBRARY_DIR
        NAMES "libca.so"
        PATHS $ENV{EPICS_BASE}/lib/${EPICS_HOST_ARCH}
        NO_DEFAULT_PATH
        DOC "Path to EPICS library files.")
endif(EPICS_INCLUDE_DIR)

if(EPICS_INCLUDE_DIR
   AND EPICS_OS_INCLUDE_DIR
   AND EPICS_LIBRARY_DIR)
    set(EPICS_FOUND true)
else(
    EPICS_INCLUDE_DIR
    AND EPICS_OS_INCLUDE_DIR
    AND EPICS_LIBRARY_DIR)
    set(EPICS_FOUND false)
endif(
    EPICS_INCLUDE_DIR
    AND EPICS_OS_INCLUDE_DIR
    AND EPICS_LIBRARY_DIR)

message("    - ${Cyan}EPICS_BASE${CR} = ${BGreen}$ENV{EPICS_BASE}${CR}")
message("    - ${Cyan}EPICS_INCLUDE_DIR${CR}   = ${BGreen}${EPICS_INCLUDE_DIR}${CR}")
message("    - ${Cyan}EPICS_OSINCLUDE_DIR${CR} = ${BGreen}${EPICS_OS_INCLUDE_DIR}${CR}")
message("    - ${Cyan}EPICS_LIBRARY_DIR${CR}   = ${BGreen}${EPICS_LIBRARY_DIR}${CR}")

if(EPICS_FOUND)
    add_definitions(-DWITH_EPICS)
    add_library(epicslib SHARED IMPORTED)
    add_library(epics::epics ALIAS epicslib)
    set_property(TARGET epicslib PROPERTY IMPORTED_LOCATION "${EPICS_LIBRARY_DIR}/libca.so")
    target_include_directories(epicslib INTERFACE ${EPICS_INCLUDE_DIR} ${EPICS_OS_INCLUDE_DIR})

    message(STATUS "${BGreen}EPICS was FOUND${CR}")
else(EPICS_FOUND)
    message(STATUS "${BYellow}Could not find package EPICS${CR}")
endif(EPICS_FOUND)
