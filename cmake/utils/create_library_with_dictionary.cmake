# generate root dictioanry from cmake target
function(add_library_with_dictionary)
    set(options)
    set(oneValueArgs LIBNAME LINKDEF)
    set(multiValueArgs SRCS INCLUDEDIRS HEADERS PRIVATE_INCLUDEDIRS DEPENDENCIES
                       PRIVATE_DEPENDENCIES)
    cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    add_library(${ARGS_LIBNAME} SHARED ${ARGS_SRCS})
    target_include_directories(
        ${ARGS_LIBNAME}
        PUBLIC ${ARGS_INCLUDEDIRS}
        PRIVATE ${ARGS_PRIVATE_INCLUDEDIRS})
    target_link_libraries(
        ${ARGS_LIBNAME}
        PUBLIC ${ARGS_DEPENDENCIES}
        PRIVATE ${ARGS_PRIVATE_DEPENDENCIES})
    root_generate_dictionary(${ARGS_LIBNAME}_dict ${ARGS_HEADERS} MODULE ${ARGS_LIBNAME} LINKDEF
                             ${ARGS_LINKDEF})
endfunction()

# generate macros using availble variables. Please use create_library_with_dictionary instead.
macro(ADD_LIB_WITH_DIC)
    add_library_with_dictionary(
        LIBNAME
        ${LIBRARY_NAME}
        LINKDEF
        ${LINKDEF}
        HEADERS
        ${HEADERS}
        SRCS
        ${SRCS}
        INCLUDEDIRS
        ${INCLUDE_DIRS}
        DEPENDENCIES
        ${DEPENDENCIES})
endmacro()
