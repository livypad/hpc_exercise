file(GLOB_RECURSE HPC_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_library(hpc
        STATIC
        ${HPC_SOURCES})

target_link_libraries(hpc)


find_package(MPI)
if (MPI_CXX_FOUND)
    message(STATUS "MPI version: ${MPI_CXX_VERSION}")
    target_link_libraries(hpc PUBLIC MPI::MPI_CXX)
endif ()


find_package(Eigen3)
target_link_libraries(hpc PUBLIC Eigen3::Eigen)

if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
    message(STATUS "current platform MacOS")

    find_package(OpenMP REQUIRED)
    if (OpenMP_CXX_FOUND)
        message(STATUS "OpenMP_CXX_VERSION: ${OpenMP_CXX_VERSION}")
        set(OpenMP_CXX_LIBRARY "/opt/homebrew/opt/libomp/lib")
        set(OpenMP_CXX_INCLUDE_DIRS "/opt/homebrew/opt/libomp/include")
        target_link_libraries(hpc PUBLIC "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    endif ()

elseif (CMAKE_SYSTEM_NAME MATCHES "Linux")
    message(STATUS "current platform Linux")

    find_package(OpenMP REQUIRED)
    if (OpenMP_CXX_FOUND)
        message(STATUS "OpenMP_CXX_VERSION: ${OpenMP_CXX_VERSION}")
        target_link_libraries(hpc PUBLIC OpenMP::OpenMP_CXX)
    endif ()
endif ()

target_include_directories(
        hpc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
