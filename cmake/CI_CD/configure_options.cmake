find_program(CCACHE "ccache")
if(CCACHE)
    message("enable ccache for cmake build")
    set(CMAKE_C_COMPILER_LAUNCHER ${CCACHE} CACHE PATH "set c launcher to ccache")
    set(CMAKE_CXX_COMPILER_LAUNCHER ${CCACHE} CACHE PATH "set c++ launcher to ccache")
endif()
