cmake_minimum_required(VERSION 3.15...3.26)
project(compas_shapeop LANGUAGES CXX)

# ==============================================================================
# Build configuration
# ==============================================================================
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build options
option(ENABLE_PRECOMPILED_HEADERS "Enable precompiled headers" ON)
option(FAST_COMPILE "Optimize for faster compilation (-O0) vs execution (-O3)" OFF)
option(USE_OPENMP "Enable OpenMP support for parallel processing" ON)

# Apply optimization flags
if(FAST_COMPILE)
  add_compile_options(-O0)
else()
  add_compile_options(-O3)
endif()

# ==============================================================================
# Dependencies
# ==============================================================================

# Setup Eigen
set(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
  # Download Eigen if not found
  if(NOT EXISTS ${EIGEN_INCLUDE_DIR})
    message(STATUS "Downloading Eigen...")
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/external)
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
  set(EIGEN3_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})
endif()

# Find Python and nanobind
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)
find_package(Threads REQUIRED)

# Setup OpenMP
if(USE_OPENMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - enabling parallel processing")
    add_definitions(-DSHAPEOP_OPENMP)
  else()
    message(STATUS "OpenMP not found - parallel processing disabled")
  endif()
endif()

# ==============================================================================
# ShapeOp library
# ==============================================================================

# Set directories
set(SHAPEOP_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/shapeop)
set(SHAPEOP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# ShapeOp library with core and custom components
add_library(shapeop STATIC
  # Core ShapeOp files
  ${SHAPEOP_SRC_DIR}/Constraint.cpp
  ${SHAPEOP_SRC_DIR}/Force.cpp
  ${SHAPEOP_SRC_DIR}/LSSolver.cpp
  ${SHAPEOP_SRC_DIR}/Solver.cpp
  # Custom constraints/forces
  ${SHAPEOP_SRC_DIR}/custom_constraints/normalforce.cpp
)

target_include_directories(shapeop PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${EIGEN3_INCLUDE_DIRS}
)

target_compile_options(shapeop PRIVATE -fPIC)
if(OpenMP_CXX_FOUND)
  target_compile_options(shapeop PUBLIC ${OpenMP_CXX_FLAGS})
  target_link_libraries(shapeop PUBLIC ${OpenMP_CXX_LIBRARIES})
endif()

# ==============================================================================
# Precompiled headers
# ==============================================================================

if(ENABLE_PRECOMPILED_HEADERS)
  add_library(compas_pch INTERFACE)
  target_precompile_headers(compas_pch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/compas.h)
  target_include_directories(compas_pch INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIRS}
    ${nanobind_INCLUDE_DIRS}
  )
endif()

# ==============================================================================
# Python module
# ==============================================================================

# Add the ShapeOpSolver Python module
nanobind_add_module(
  _shapeop
  STABLE_ABI
  NB_STATIC
  src/shapeop.cpp
)

target_include_directories(_shapeop PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${EIGEN3_INCLUDE_DIRS}
  ${nanobind_INCLUDE_DIRS}
)

target_link_libraries(_shapeop PRIVATE shapeop)
if(ENABLE_PRECOMPILED_HEADERS)
  target_link_libraries(_shapeop PRIVATE compas_pch)
endif()

install(TARGETS _shapeop LIBRARY DESTINATION compas_shapeop)

# ==============================================================================
# Summary
# ==============================================================================

message(STATUS "============= Build Configuration =============")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Optimization: ${FAST_COMPILE} (O0 if ON, O3 if OFF)")
message(STATUS "Precompiled Headers: ${ENABLE_PRECOMPILED_HEADERS}")
message(STATUS "=============================================")