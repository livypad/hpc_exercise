cmake_minimum_required(VERSION 3.14)

enable_testing()

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)
include(FetchContent)
cmake_policy(SET CMP0135 NEW)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

include(GoogleTest)

file(GLOB_RECURSE HPC_TEST_SOURCES "${PROJECT_SOURCE_DIR}/test/*test.cpp")
message(STATUS "HPC_TEST_SOURCES: ${HPC_TEST_SOURCES}")
#add_custom_target(build-tests COMMAND ${CMAKE_CTEST_COMMAND} --show-only)
#add_custom_target(check-tests COMMAND ${CMAKE_CTEST_COMMAND} --verbose)
foreach (hpc_test_source ${HPC_TEST_SOURCES})
    # Create a human readable name.
    get_filename_component(hpc_test_filename ${hpc_test_source} NAME)
    string(REPLACE ".cpp" "" hpc_test_name ${hpc_test_filename})

    # Add the test target separately and as part of "make check-tests".
    add_executable(${hpc_test_name} EXCLUDE_FROM_ALL ${hpc_test_source})
    #    add_dependencies(build-tests ${hpc_test_name})
    #    add_dependencies(check-tests ${hpc_test_name})

    gtest_discover_tests(${hpc_test_name}
            EXTRA_ARGS
            --gtest_color=auto
            --gtest_catch_exceptions=0
            DISCOVERY_TIMEOUT 120
            PROPERTIES
            TIMEOUT 120
            )

    #    target_include_directories(${hpc_test_name} src/include)
    target_link_libraries(${hpc_test_name} gtest gmock_main hpc)

    # Set test target properties and dependencies.
    set_target_properties(${hpc_test_name}
            PROPERTIES
            RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/test"
            COMMAND ${hpc_test_name}
            )
endforeach ()