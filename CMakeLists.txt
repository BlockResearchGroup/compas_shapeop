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
option(USE_OPENMP "Enable OpenMP support for parallel processing" OFF)

# Apply optimization flags in a cross-platform way
if(FAST_COMPILE)
  if(MSVC)
    add_compile_options(/Od)  # Disable optimization (MSVC)
  else()
    add_compile_options(-O0)  # Disable optimization (GCC/Clang)
  endif()
else()
  if(MSVC)
    add_compile_options(/O2)  # Optimize for speed (MSVC)
  else()
    add_compile_options(-O3)  # Optimize for speed (GCC/Clang)
  endif()
endif()

# ==============================================================================
# Dependencies
# ==============================================================================
include(ExternalProject)

# Define source directories for external dependencies
set(EXTERNAL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(EIGEN_INCLUDE_DIR "${EXTERNAL_DIR}/eigen")

# Create a custom target for all external dependencies
add_custom_target(external_downloads ALL)

# ========================================================================
# Setup Eigen
# ========================================================================
if(NOT EXISTS "${EIGEN_INCLUDE_DIR}")
  message(STATUS "Downloading Eigen...")
  ExternalProject_Add(
      eigen_download
      URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.zip
      SOURCE_DIR "${EIGEN_INCLUDE_DIR}"
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      LOG_DOWNLOAD ON
      UPDATE_COMMAND ""
      PATCH_COMMAND ""
      TLS_VERIFY ON
  )
  add_dependencies(external_downloads eigen_download)
endif()
set(EIGEN3_INCLUDE_DIRS "${EIGEN_INCLUDE_DIR}")

# Print include directories for verification
message(STATUS "============= Dependencies =============")
message(STATUS "EIGEN INCLUDE DIRECTORY: ${EIGEN_INCLUDE_DIR}")
message(STATUS "=========================================")

# Find Python and nanobind
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)
find_package(Threads REQUIRED)

# Setup OpenMP - Properly handle macOS and make it optional
if(USE_OPENMP)
  # Check what platform we're on
  if(APPLE)
    # On macOS, OpenMP might not be available with the default Clang
    # Check if we have Homebrew's LLVM/Clang with OpenMP support
    find_program(BREW NAMES brew)
    if(BREW)
      execute_process(
        COMMAND ${BREW} --prefix llvm
        OUTPUT_VARIABLE LLVM_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
      )
      if(LLVM_DIR)
        # Use Homebrew's LLVM for OpenMP support
        set(OpenMP_CXX_FLAGS "-fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LLVM_DIR}/lib/libomp.dylib")
        message(STATUS "Using Homebrew LLVM for OpenMP support: ${LLVM_DIR}")
        add_definitions(-DSHAPEOP_OPENMP)
      else()
        message(STATUS "Homebrew LLVM not found. OpenMP support disabled on macOS.")
        set(USE_OPENMP OFF)
      endif()
    else()
      message(STATUS "Homebrew not found. OpenMP support disabled on macOS.")
      set(USE_OPENMP OFF)
    endif()
  else()
    # For non-macOS platforms, use standard FindOpenMP
    find_package(OpenMP QUIET)
    if(OpenMP_CXX_FOUND)
      message(STATUS "OpenMP found - enabling parallel processing")
      add_definitions(-DSHAPEOP_OPENMP)
    else()
      message(STATUS "OpenMP not found - parallel processing disabled")
      set(USE_OPENMP OFF)
    endif()
  endif()
endif()

# Add include directories
include_directories(
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${EIGEN3_INCLUDE_DIRS}
)

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
  ${SHAPEOP_SRC_DIR}/custom_constraints/NormalForce.cpp
)

# Make shapeop depend on the external downloads
add_dependencies(shapeop external_downloads)

target_compile_options(shapeop PRIVATE -fPIC)
target_compile_definitions(shapeop PUBLIC SHAPEOP_EXPORT)

# Apply OpenMP flags only if enabled and available
if(USE_OPENMP)
  if(DEFINED OpenMP_CXX_FLAGS)
    target_compile_options(shapeop PUBLIC ${OpenMP_CXX_FLAGS})
    if(DEFINED OpenMP_CXX_LIBRARIES)
      target_link_libraries(shapeop PUBLIC ${OpenMP_CXX_LIBRARIES})
    endif()
    if(APPLE AND DEFINED OpenMP_omp_LIBRARY)
      target_link_libraries(shapeop PUBLIC ${OpenMP_omp_LIBRARY})
    endif()
  endif()
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

# Function to add a nanobind module with include directories
function(add_nanobind_module module_name source_file)
  nanobind_add_module(
    ${module_name}
    STABLE_ABI
    NB_STATIC
    ${source_file}
  )

  target_include_directories(${module_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIRS}
    ${nanobind_INCLUDE_DIRS}
  )

  # Make module depend on the external downloads
  add_dependencies(${module_name} external_downloads)

  # Apply OpenMP flags only if enabled and available
  if(USE_OPENMP)
    if(DEFINED OpenMP_CXX_FLAGS)
      target_compile_options(${module_name} PRIVATE ${OpenMP_CXX_FLAGS})
      if(DEFINED OpenMP_CXX_LIBRARIES)
        target_link_libraries(${module_name} PRIVATE ${OpenMP_CXX_LIBRARIES})
      endif()
      if(APPLE AND DEFINED OpenMP_omp_LIBRARY)
        target_link_libraries(${module_name} PRIVATE ${OpenMP_omp_LIBRARY})
      endif()
    endif()
  endif()

  if(ENABLE_PRECOMPILED_HEADERS)
    target_link_libraries(${module_name} PRIVATE compas_pch)
  endif()

  install(TARGETS ${module_name} LIBRARY DESTINATION compas_shapeop)
endfunction()

# Add ShapeOp module
add_nanobind_module(_shapeop src/shapeop.cpp)
target_link_libraries(_shapeop PRIVATE shapeop)

# ==============================================================================
# Summary
# ==============================================================================

message(STATUS "============= Build Configuration =============")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Optimization: ${FAST_COMPILE} (O0 if ON, O3 if OFF)")
message(STATUS "Precompiled Headers: ${ENABLE_PRECOMPILED_HEADERS}")
message(STATUS "OpenMP Support: ${USE_OPENMP}")
message(STATUS "Eigen Include Dir: ${EIGEN3_INCLUDE_DIRS}")
message(STATUS "External Dir: ${EXTERNAL_DIR}")
message(STATUS "===============================================")