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

macro(WRITE_CONFIG_FILE filename)

    string(REGEX REPLACE "^.*(install).*$" "\\1" INSTALL_VERSION ${filename})
    string(COMPARE EQUAL "install" "${INSTALL_VERSION}" INSTALL_TRUE)

    list(REMOVE_ITEM LD_LIBRARY_PATH ${CMAKE_BINARY_DIR}/lib)

    if(INSTALL_TRUE)
        set(_INSTALLDIR ${CMAKE_INSTALL_PREFIX})
        set(_BINDIR ${CMAKE_INSTALL_PREFIX}/bin)
        #    SET(VMCWORKDIR ${CMAKE_INSTALL_PREFIX}/share/fairbase)
        set(FAIRLIBDIR ${CMAKE_INSTALL_PREFIX}/lib)
        set(_LD_LIBRARY_PATH ${FAIRLIBDIR} ${LD_LIBRARY_PATH})
    else(INSTALL_TRUE)
        set(_INSTALLDIR ${CMAKE_BINARY_DIR})
        set(_BINDIR ${CMAKE_BINARY_DIR})
        #    SET(VMCWORKDIR ${CMAKE_SOURCE_DIR})
        set(FAIRLIBDIR ${CMAKE_BINARY_DIR}/lib)
        set(_LD_LIBRARY_PATH ${FAIRLIBDIR} ${LD_LIBRARY_PATH})
    endif(INSTALL_TRUE)

    if(NOT DEFINED FULL_CONFIG_FILE)
        set(FULL_CONFIG_FILE "true")
    endif(NOT DEFINED FULL_CONFIG_FILE)

    if(CMAKE_SYSTEM_NAME MATCHES Linux)
        if(FAIRROOTPATH)
            configure_file(${FAIRROOTPATH}/share/fairbase/cmake/scripts/check_system.sh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/check_system.sh)
            configure_file(${FAIRROOTPATH}/share/fairbase/cmake/scripts/check_system.csh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/check_system.csh)
        else(FAIRROOTPATH)
            configure_file(${PROJECT_SOURCE_DIR}/cmake/scripts/check_system.sh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/check_system.sh)
            configure_file(${PROJECT_SOURCE_DIR}/cmake/scripts/check_system.csh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/check_system.csh)
        endif(FAIRROOTPATH)

        execute_process(
            COMMAND lsb_release -sd
            OUTPUT_VARIABLE _linux_flavour
            OUTPUT_STRIP_TRAILING_WHITESPACE)

        if(_linux_flavour)
            string(REGEX REPLACE "^\"" "" _linux_flavour ${_linux_flavour})
            string(REGEX REPLACE "\"$" "" _linux_flavour ${_linux_flavour})
        endif(_linux_flavour)

        execute_process(
            COMMAND uname -m
            OUTPUT_VARIABLE _system
            OUTPUT_STRIP_TRAILING_WHITESPACE)

    elseif(CMAKE_SYSTEM_NAME MATCHES Darwin)
        if(FAIRROOTPATH)
            configure_file(${FAIRROOTPATH}/share/fairbase/cmake/scripts/check_system_mac.sh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/check_system.sh)
            configure_file(${FAIRROOTPATH}/share/fairbase/cmake/scripts/check_system_mac.csh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/check_system.csh)
        else(FAIRROOTPATH)
            configure_file(${PROJECT_SOURCE_DIR}/cmake/scripts/check_system_mac.sh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/check_system.sh)
            configure_file(${PROJECT_SOURCE_DIR}/cmake/scripts/check_system_mac.csh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/check_system.csh)
        endif(FAIRROOTPATH)
        execute_process(
            COMMAND uname -sr
            OUTPUT_VARIABLE _linux_flavour
            OUTPUT_STRIP_TRAILING_WHITESPACE)
        execute_process(
            COMMAND uname -m
            OUTPUT_VARIABLE _system
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    endif(CMAKE_SYSTEM_NAME MATCHES Linux)

    convert_list_to_string(${Geant4_INCLUDE_DIRS})
    set(Geant4_INCLUDE_DIRS ${output})

    convert_list_to_string(${Geant4VMC_INCLUDE_DIRS})
    set(Geant4VMC_INCLUDE_DIRS ${output})

    convert_list_to_string(${Geant4VMC_CMAKE_INSTALL_LIBDIR})
    set(Geant4VMC_LIBRARY_DIR ${output})

    convert_list_to_string(${Geant4VMC_MACRO_DIR})
    set(Geant4VMC_MACRO_DIR ${output})

    write_geant4data_variables_sh()
    write_geant4data_variables_csh()

    convert_list_to_string(${PLUTO_LIBRARY_DIR})
    set(PLUTO_LIBRARY_DIR ${output})

    convert_list_to_string(${PLUTO_INCLUDE_DIR})
    set(PLUTO_INCLUDE_DIR ${output})

    convert_list_to_string(${PYTHIA6_LIBRARY_DIR})
    set(PYTHIA6_LIBRARY_DIR ${output})

    convert_list_to_string(${Geant3_SYSTEM_DIR})
    set(G3SYS ${output})

    convert_list_to_string(${Geant3_INCLUDE_DIRS})
    set(Geant3_INCLUDE_DIRS ${output})

    convert_list_to_string(${Geant3_LIBRARY_DIR})
    set(Geant3_LIBRARY_DIR ${output})

    convert_list_to_string(${Geant3_LIBRARIES})
    set(Geant3_LIBRARIES ${output})

    convert_list_to_string(${ROOT_LIBRARY_DIR})
    set(ROOT_LIBRARY_DIR ${output})

    convert_list_to_string(${ROOT_LIBRARIES})
    set(ROOT_LIBRARIES ${output})

    convert_list_to_string(${ROOT_INCLUDE_DIR})
    set(ROOT_INCLUDE_DIR ${output})

    #  Set(VMCWORKDIR ${C})

    set(FAIRLIBDIR ${FAIRLIBDIR})

    list(REMOVE_DUPLICATES _LD_LIBRARY_PATH)
    convert_list_to_string(${_LD_LIBRARY_PATH})

    if(CMAKE_SYSTEM_NAME MATCHES Linux)
        set(MY_LD_LIBRARY_PATH ${output})
    else(CMAKE_SYSTEM_NAME MATCHES Linux)
        if(CMAKE_SYSTEM_NAME MATCHES Darwin)
            set(MY_DYLD_LIBRARY_PATH ${output})
        endif(CMAKE_SYSTEM_NAME MATCHES Darwin)
    endif(CMAKE_SYSTEM_NAME MATCHES Linux)

    set(USE_VGM 1)

    set(PYTHONPATH
        ${CMAKE_SOURCE_DIR}/python
        ${SIMPATH}/lib
        ${SIMPATH}/lib/root
        ${SIMPATH}/lib/Geant4
        ${SIMPATH}/lib/g4py
        ${PYTHONPATH})
    unique(PYTHONPATH "${PYTHONPATH}")
    convert_list_to_string(${PYTHONPATH})
    set(MY_PYTHONPATH ${output})

    # TODO: this should change to CMAKE_INSTALL_PREFIX/bin in the future
    list(APPEND PATH ${ROOT_BINDIR})
    list(APPEND PATH ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    unique(PATH "${PATH}")
    convert_list_to_string(${PATH})
    set(MY_PATH ${output})

    if(FAIRSOFT_EXTERN)
        set(PYTHIA8DATA "${SIMPATH}/share/pythia8/xmldoc")
    else()
        set(PYTHIA8DATA "${SIMPATH}/generators/pythia8/xmldoc")
    endif()

    convert_list_to_string($ENV{NEW_CLASSPATH})
    set(MY_CLASSPATH ${output})

    convert_list_to_string(${ROOT_INCLUDE_PATH})
    set(ROOT_INCLUDE_PATH ${output})

    if(${filename} MATCHES "[.]csh.*$")
        if(R3BROOTPATH)
            configure_file(${PROJECT_SOURCE_DIR}/cmake/scripts/config.csh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/${filename})
        elseif(R3BROOTPATH)
            configure_file(${FAIRROOTPATH}/share/fairbase/cmake/scripts/config.csh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/${filename})
        endif(R3BROOTPATH)

    else()
        if(R3BROOTPATH)
            configure_file(${PROJECT_SOURCE_DIR}/cmake/scripts/config.sh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/${filename})
        elseif(R3BROOTPATH)
            configure_file(${FAIRROOTPATH}/share/fairbase/cmake/scripts/config.sh.in
                           ${CMAKE_CURRENT_BINARY_DIR}/${filename})
        endif(R3BROOTPATH)

    endif()

endmacro(WRITE_CONFIG_FILE)

macro(CONVERT_LIST_TO_STRING)
    set(tmp "")
    foreach(_current ${ARGN})

        set(tmp1 ${tmp})
        set(tmp "")
        set(tmp ${tmp1}:${_current})

    endforeach(_current ${ARGN})
    if(tmp)
        string(REGEX REPLACE "^:(.*)" "\\1" output ${tmp})
    else(tmp)
        set(output "")
    endif(tmp)
endmacro(CONVERT_LIST_TO_STRING LIST)
