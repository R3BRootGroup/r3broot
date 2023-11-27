# Input parameters: -DBUILD_J=1 -DTEST_NAME="R3BRoot CI/CD" -DTEST_MODEL="Experimental" -DBUILD_TYPE= -DEXTRA_FLAGS= -DSITE_NAME=

# mandatory cmake variable definition
set(CTEST_SOURCE_DIRECTORY "${CTEST_SCRIPT_DIRECTORY}/../..")
set(CTEST_BINARY_DIRECTORY "${CTEST_SOURCE_DIRECTORY}/build")
set(CTEST_CONFIGURE_COMMAND
    "${CMAKE_COMMAND} -S ${CTEST_SOURCE_DIRECTORY} -B ${CTEST_BINARY_DIRECTORY}")
set(CTEST_BUILD_COMMAND "${CMAKE_COMMAND} --build ${CTEST_BINARY_DIRECTORY}")

#include macros and function definition
include(${CTEST_SOURCE_DIRECTORY}/cmake/CI_CD/Macros.cmake)

set_input_options()

check_test_model(${ctest_model})

enable_launcher(1)

cdash_start()

cdash_configure()

cdash_build()

cdash_test()

cdash_submit()

print_error()
