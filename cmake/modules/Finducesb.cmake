################################################################################
#   Copyright (C) 2019 GSI Helmholtzzentrum für Schwerionenforschung GmbH      #
#   Copyright (C) 2019-2024 Members of R3B Collaboration                       #
#                                                                              #
#             This software is distributed under the terms of the              #
#                 GNU General Public Licence (GPL) version 3,                  #
#                    copied verbatim in the file "LICENSE".                    #
#                                                                              #
# In applying this license GSI does not waive the privileges and immunities    #
# granted to it by virtue of its status as an Intergovernmental Organization   #
# or submit itself to any jurisdiction.                                        #
################################################################################

################################################################################
# For this script to succeed, you need to install ucesb:                       #
# Check out the git repository and compile:                                    #
#  git clone http://fy.chalmers.se/~f96hajo/ucesb/ucesb.git                    #
#  cd ucesb                                                                    #
#  make -j empty/empty								  #
# Make sure to use the same ROOT version as for R3BROOT			  #
# After installation export an environment variable UCESB_DIR                  #
# which points to the 'unpacker' directory you have checked out.               #
# More information: http://fy.chalmers.se/~f96hajo/ucesb/                      #
################################################################################

message(STATUS "Looking for ucesb...")

set(ucesb_LIBRARY_SHARED ext_data_clnt)

set(ucesb_LIBRARY_STATIC ext_data_clnt.o)

find_path(
    ucesb_INCLUDE_DIR
    NAMES ext_data_clnt.hh ext_h101.h
    PATHS $ENV{UCESB_DIR} $ENV{UCESB_DIR}/hbook
    NO_DEFAULT_PATH
    DOC "Path to ucesb include header files.")

find_path(
    ucesb_LIBRARY_DIR
    NAMES "lib${ucesb_LIBRARY_SHARED}.so"
    PATHS $ENV{UCESB_DIR} $ENV{UCESB_DIR}/hbook
    NO_DEFAULT_PATH
    DOC "Path to ucesb library files.")
#find_library(ucesb_LIBRARY_SHARED NAMES ${ucesb_LIBRARY_SHARED}
#  PATHS $ENV{UCESB_DIR}/hbook
#  NO_DEFAULT_PATH
#  DOC   "Path to ${LIBucesb_SHARED}."
#)

find_library(
    ucesb_LIBRARY_STATIC
    NAMES ${ucesb_LIBRARY_STATIC}
    PATHS $ENV{UCESB_DIR}/hbook
    NO_DEFAULT_PATH
    DOC "Path to ${LIBucesb_STATIC}.")

add_library(Ucesb SHARED IMPORTED GLOBAL)
set_property(TARGET Ucesb PROPERTY IMPORTED_LOCATION
                                   "${ucesb_LIBRARY_DIR}/lib${ucesb_LIBRARY_SHARED}.so")
target_include_directories(Ucesb INTERFACE $ENV{UCESB_DIR}/hbook)

if(ucesb_INCLUDE_DIR
   AND ucesb_LIBRARY_SHARED
   AND ucesb_LIBRARY_STATIC
   AND ucesb_LIBRARY_DIR)
    set(ucesb_FOUND true)
else(
    ucesb_INCLUDE_DIR
    AND ucesb_LIBRARY_SHARED
    AND ucesb_LIBRARY_STATIC
    AND ucesb_LIBRARY_DIR)
    set(ucesb_FOUND false)
endif(
    ucesb_INCLUDE_DIR
    AND ucesb_LIBRARY_SHARED
    AND ucesb_LIBRARY_STATIC
    AND ucesb_LIBRARY_DIR)

message("    - ${Cyan}UCESB_DIR${CR} = ${BGreen}$ENV{UCESB_DIR}${CR}")
message("    - ${Cyan}ucesb_INCLUDE_DIR${CR} = ${BGreen}${ucesb_INCLUDE_DIR}${CR}")
message("    - ${Cyan}ucesb_LIBRARY_SHARED${CR} = ${BGreen}${ucesb_LIBRARY_SHARED}${CR}")
message("    - ${Cyan}ucesb_LIBRARY_STATIC${CR} = ${BGreen}${ucesb_LIBRARY_STATIC}${CR}")
message("    - ${Cyan}ucesb_LIBRARY_DIR${CR} = ${BGreen}${ucesb_LIBRARY_DIR}${CR}")

if(ucesb_FOUND)
    add_definitions(-DWITH_UCESB)
    set(ucesb_LIBRARIES "${ucesb_LIBRARY_STATIC};${ucesb_LIBRARY_SHARED}")
    if(NOT ucesb_FIND_QUIETLY)
        set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${ucesb_LIBRARY_DIR})
    endif(NOT ucesb_FIND_QUIETLY)
    message(STATUS "${BBlue}UCESB was FOUND${CR}")
else(ucesb_FOUND)
    if(NOT ucesb_FIND_QUIETLY)
        if(ucesb_FIND_REQUIRED)
            message(FATAL_ERROR "${BRed}Could not find package UCESB${CR}")
        else(ucesb_FIND_REQUIRED)
            message(STATUS "${BYellow}Could not find package UCESB${CR}")
        endif(ucesb_FIND_REQUIRED)
    endif(NOT ucesb_FIND_QUIETLY)
endif(ucesb_FOUND)

mark_as_advanced(ucesb_INCLUDE_DIR ucesb_LIBRARIES ucesb_LIBRARY_SHARED ucesb_LIBRARY_STATIC)
