find_package(Eigen3)

add_library(
    R3BMille SHARED
    Mille.cxx
    ParResultReader.cxx
    PedeLauncher.cxx
    SteerWriter.cxx)
target_include_directories(R3BMille PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${Boost_INCLUDE_DIR})
target_link_libraries(
    R3BMille
    PUBLIC ${Boost_LIBRARIES}
    PRIVATE fmt::fmt range-v3::range-v3)
if(DEFINED ENV{PEDEDIR})
    message(STATUS "Millepede directory is set to be $ENV{PEDEDIR}")
    target_compile_definitions(R3BMille PUBLIC PEDEDIR="$ENV{PEDEDIR}")
endif()
