cmake_minimum_required(VERSION 3.14)
project(HPC_Exercises)

# GoogleTest requires at least C++14
set(CMAKE_CXX_STANDARD 14)

include_directories(${PROJECT_SOURCE_DIR}/src/include)
include_directories(BEFORE src) # This is needed for gtest.
add_subdirectory(src)

enable_testing()
add_subdirectory(test)
