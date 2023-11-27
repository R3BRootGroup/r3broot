# Note the Fortran project mention in the main CMakeLists.txt.
include(CheckFortran)
add_library(R3BTraRene SHARED tracker_routines.f95)
target_include_directories(R3BTraRene PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
