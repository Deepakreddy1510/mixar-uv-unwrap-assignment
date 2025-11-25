cmake_minimum_required(VERSION 3.15)
project(uvunwrap)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Eigen library - try to find system Eigen first
find_package(Eigen3 3.3 QUIET)
if(EIGEN3_FOUND)
    include_directories(${EIGEN3_INCLUDE_DIR})
    message(STATUS "Using system Eigen: ${EIGEN3_INCLUDE_DIR}")
else()
    # Fall back to bundled version
    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)
        include_directories(${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen)
        message(STATUS "Using bundled Eigen")
    else()
        message(FATAL_ERROR "Eigen not found. Please install Eigen3 (apt-get install libeigen3-dev)")
    endif()
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Source files
set(SOURCES
    src/mesh_io.cpp
    src/math_utils.cpp
    src/topology.cpp
    src/seam_detection.cpp
    src/lscm.cpp
    src/packing.cpp
    src/unwrap.cpp
)

# Main library
add_library(uvunwrap SHARED ${SOURCES})

# Test executable
add_executable(test_unwrap tests/test_unwrap.cpp)
target_link_libraries(test_unwrap uvunwrap)

# Enable warnings
if(MSVC)
    target_compile_options(uvunwrap PRIVATE /W4)
else()
    target_compile_options(uvunwrap PRIVATE -Wall -Wextra -Wpedantic)
endif()

# Install targets
install(TARGETS uvunwrap
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        RUNTIME DESTINATION bin)

install(DIRECTORY include/
        DESTINATION include
        FILES_MATCHING PATTERN "*.h")
