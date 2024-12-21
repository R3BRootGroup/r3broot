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

# Remove trailing slash from a path passed as argument
Macro(Remove_Trailing_Slash _variable)
  String(FIND ${_variable} "/" _pos_last_slash REVERSE)
  STRING(LENGTH ${_variable} _length)
  Math(EXPR _last_pos ${_pos_last_slash}+1)
  If(${_last_pos} EQUAL ${_length})  
    String(SUBSTRING ${_variable} 0 ${_pos_last_slash} _ret_val)
  Else()
    Set(_ret_val ${_variable})
  EndIf()
EndMacro()


MACRO (ROOT_GENERATE_DICTIONARY_OLD_EXTRA INFILES LINKDEF_FILE OUTFILE INCLUDE_DIRS_IN)

  set(INCLUDE_DIRS)

  foreach (_current_FILE ${INCLUDE_DIRS_IN})
    set(INCLUDE_DIRS ${INCLUDE_DIRS} -I${_current_FILE})
  endforeach (_current_FILE ${INCLUDE_DIRS_IN})

#  Message("Definitions: ${DEFINITIONS}")
#  MESSAGE("INFILES: ${INFILES}")
#  MESSAGE("OutFILE: ${OUTFILE}")
#  MESSAGE("LINKDEF_FILE: ${LINKDEF_FILE}")
#  MESSAGE("INCLUDE_DIRS: ${INCLUDE_DIRS}")

  STRING(REGEX REPLACE "^(.*)\\.(.*)$" "\\1.h" bla "${OUTFILE}")
#  MESSAGE("BLA: ${bla}")
  SET (OUTFILES ${OUTFILE} ${bla})

  if (CMAKE_SYSTEM_NAME MATCHES Linux)
    ADD_CUSTOM_COMMAND(OUTPUT ${OUTFILES}
       COMMAND LD_LIBRARY_PATH=${ROOT_LIBRARY_DIR}:${_intel_lib_dirs} ROOTSYS=${ROOTSYS} ${ROOT_CINT_EXECUTABLE}
       ARGS -f ${OUTFILE} -c -p -DHAVE_CONFIG_H ${INCLUDE_DIRS} ${INFILES} ${LINKDEF_FILE} DEPENDS ${INFILES} ${LINKDEF_FILE})
  else (CMAKE_SYSTEM_NAME MATCHES Linux)
    if (CMAKE_SYSTEM_NAME MATCHES Darwin)
      ADD_CUSTOM_COMMAND(OUTPUT ${OUTFILES}
       COMMAND DYLD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS} ${ROOT_CINT_EXECUTABLE}
       ARGS -f ${OUTFILE} -c -p -DHAVE_CONFIG_H ${INCLUDE_DIRS} ${INFILES} ${LINKDEF_FILE} DEPENDS ${INFILES} ${LINKDEF_FILE})
    endif (CMAKE_SYSTEM_NAME MATCHES Darwin)
  endif (CMAKE_SYSTEM_NAME MATCHES Linux)

ENDMACRO (ROOT_GENERATE_DICTIONARY_OLD_EXTRA)


Macro (R3B_Generate_Version_Info)
  Add_Custom_Target(svnr3bheader ALL)

  Add_Custom_Command(TARGET svnr3bheader
                     COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${CMAKE_SOURCE_DIR}
                     -DBINARY_DIR=${CMAKE_BINARY_DIR}
                     -DINCLUDE_OUTPUT_DIRECTORY=${INCLUDE_OUTPUT_DIRECTORY}
                     -P ${CMAKE_SOURCE_DIR}/cmake/modules/R3BGenerateVersionInfo.cmake
                    )
EndMacro (R3B_Generate_Version_Info)
