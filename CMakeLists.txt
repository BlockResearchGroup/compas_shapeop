cmake_minimum_required(VERSION 3.15...3.26)

project(compas_shapeop LANGUAGES CXX)

# Force CMake to use a specific build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Fast compilation option (ON = faster compilation, OFF = faster execution)
option(FAST_COMPILE "Optimize for faster compilation" ON)
if(FAST_COMPILE)
  add_compile_options(-O0)
else()
  add_compile_options(-O3)
endif()

# Ninja is faster for incremental builds
if (CMAKE_GENERATOR MATCHES "Ninja")
  message(STATUS "Using Ninja generator")
endif()

# Create external directory if it doesn't exist
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external)

# External dependencies setup (Eigen and ShapeOp)
set(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)
set(SHAPEOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/shapeop)
set(SHAPEOP_SRC_DIR ${SHAPEOP_SOURCE_DIR}/libShapeOp/src)
set(SHAPEOP_API_DIR ${SHAPEOP_SOURCE_DIR}/libShapeOp/api)
set(SHAPEOP_INCLUDE_DIR ${SHAPEOP_SOURCE_DIR}/libShapeOp)

# Download Eigen
if(NOT EXISTS ${EIGEN_INCLUDE_DIR})
  message(STATUS "Downloading Eigen...")
  file(DOWNLOAD
    https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen.zip
    SHOW_PROGRESS
  )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/eigen.zip
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external
  )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E rename 
      ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen-3.4.0 
      ${EIGEN_INCLUDE_DIR}
  )
  file(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/eigen.zip)
endif()

# Download ShapeOp
if(NOT EXISTS ${SHAPEOP_SOURCE_DIR})
  message(STATUS "Downloading ShapeOp...")
  file(DOWNLOAD
    https://github.com/EPFL-LGG/ShapeOp/archive/refs/heads/master.zip
    ${CMAKE_CURRENT_SOURCE_DIR}/shapeop.zip
    SHOW_PROGRESS
  )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E tar xf ${CMAKE_CURRENT_SOURCE_DIR}/shapeop.zip
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external
  )
  execute_process(
    COMMAND ${CMAKE_COMMAND} -E rename 
      ${CMAKE_CURRENT_SOURCE_DIR}/external/ShapeOp-master 
      ${SHAPEOP_SOURCE_DIR}
  )
  file(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/shapeop.zip)
endif()

# Verify that required files exist
if(NOT EXISTS ${SHAPEOP_SRC_DIR}/Constraint.cpp)
  message(FATAL_ERROR "ShapeOp source files missing. Expected: ${SHAPEOP_SRC_DIR}/Constraint.cpp")
endif()

# Build ShapeOp as a library
add_library(shapeop STATIC
  ${SHAPEOP_SRC_DIR}/Constraint.cpp
  ${SHAPEOP_SRC_DIR}/Force.cpp
  ${SHAPEOP_SRC_DIR}/LSSolver.cpp
  ${SHAPEOP_SRC_DIR}/Solver.cpp
)

# Ensure the library is compiled with -fPIC
# To build static library to shared library.
target_compile_options(shapeop PRIVATE -fPIC)

target_include_directories(shapeop PRIVATE
  ${EIGEN_INCLUDE_DIR}
  ${SHAPEOP_INCLUDE_DIR}
  ${SHAPEOP_SRC_DIR}
  ${SHAPEOP_API_DIR}
)

# Enable precompiled headers for shapeop
target_precompile_headers(shapeop PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/compas.h)


target_compile_options(shapeop PRIVATE -O3)

# Build the Python extension module
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(
  _shapeop
  STABLE_ABI
  NB_STATIC
  src/shapeop.cpp

)

target_include_directories(_shapeop PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${EIGEN_INCLUDE_DIR}
  ${SHAPEOP_INCLUDE_DIR}
  ${SHAPEOP_SRC_DIR}
  ${SHAPEOP_API_DIR}
)

target_link_libraries(_shapeop PRIVATE shapeop)

# Install directive for scikit-build-core
install(TARGETS _shapeop LIBRARY DESTINATION compas_shapeop)

# Build configuration info
message(STATUS "============= Build Configuration =============")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Optimization: ${FAST_COMPILE} (O0 if ON, O3 if OFF)")
message(STATUS "=============================================")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)