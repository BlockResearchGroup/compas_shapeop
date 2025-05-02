cmake_minimum_required(VERSION 3.15...3.26)

project(compas_shapeop LANGUAGES CXX)

# Force CMake to use a specific build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set C++ standard to C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option for precompiled headers
option(ENABLE_PRECOMPILED_HEADERS "Enable precompiled headers for the build" ON)

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
# First set up Eigen directory to maintain compatibility with existing code
set(EIGEN_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen)
find_package(Eigen3 QUIET)
if(NOT Eigen3_FOUND)
  # Use our downloaded Eigen if system Eigen not found
  set(EIGEN3_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})
endif()

# ShapeOp is now integrated directly into our source tree
set(SHAPEOP_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(SHAPEOP_API_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src/api)
set(SHAPEOP_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Add debug message for SHAPEOP_API_DIR
message(STATUS "SHAPEOP_API_DIR: ${SHAPEOP_API_DIR}")

# OpenMP support for parallel processing (especially for large meshes)
option(USE_OPENMP "Enable OpenMP support for parallel processing" ON)
if(USE_OPENMP)
  find_package(OpenMP QUIET)
  if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found - enabling parallel processing")
    # Define OpenMP flag for ShapeOp 
    add_definitions(-DSHAPEOP_OPENMP)
    # We'll apply OpenMP flags to the libraries that need them below
  else()
    message(STATUS "OpenMP not found - parallel processing disabled")
  endif()
endif()

# Platform specific settings
if(APPLE)
  # macOS has historically had OpenMP issues with default compilers
  message(STATUS "On macOS: If you have OpenMP issues, you can disable with -DUSE_OPENMP=OFF")
endif()

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
  
  # When using downloaded Eigen, ensure the include dirs point to it
  set(EIGEN3_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR})
endif()

# Find Python and nanobind
find_package(Python 3.8 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(nanobind CONFIG REQUIRED)
find_package(Threads REQUIRED)

# Debugging: Print nanobind include directories
message(STATUS "nanobind_INCLUDE_DIRS: ${nanobind_INCLUDE_DIRS}")

# Build ShapeOp as a library
add_library(shapeop STATIC
  ${SHAPEOP_SRC_DIR}/Constraint.cpp
  ${SHAPEOP_SRC_DIR}/Force.cpp
  ${SHAPEOP_SRC_DIR}/LSSolver.cpp
  ${SHAPEOP_SRC_DIR}/Solver.cpp
)

# Create a separate library for our custom ShapeOp extensions
add_library(shapeop_extensions STATIC 
  ${CMAKE_CURRENT_SOURCE_DIR}/src/api/API.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/custom_constraints/normalforce.cpp
)

target_include_directories(shapeop_extensions PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/custom_constraints
  ${EIGEN3_INCLUDE_DIRS}
)

target_link_libraries(shapeop_extensions PUBLIC shapeop)

if(OpenMP_CXX_FOUND)
  target_compile_options(shapeop_extensions PRIVATE ${OpenMP_CXX_FLAGS})
  target_link_libraries(shapeop_extensions PRIVATE ${OpenMP_CXX_LIBRARIES})
endif()

# Add include directories for ShapeOp
target_include_directories(shapeop PRIVATE
  ${EIGEN3_INCLUDE_DIRS}
  ${SHAPEOP_INCLUDE_DIR}
  ${SHAPEOP_SRC_DIR}
  ${SHAPEOP_API_DIR}
  ${nanobind_INCLUDE_DIRS}
)

# Make sure position independent code is used
target_compile_options(shapeop PRIVATE -fPIC)
target_compile_options(shapeop_extensions PRIVATE -fPIC)

# Apply OpenMP flags if available
if(OpenMP_CXX_FOUND)
  target_compile_options(shapeop PUBLIC ${OpenMP_CXX_FLAGS})
  target_link_libraries(shapeop PUBLIC ${OpenMP_CXX_LIBRARIES})
endif()

# Force rebuild every time by updating timestamp
add_definitions(-DBUILD_TIMESTAMP="${CMAKE_CURRENT_TIMESTAMP}")
add_definitions(-DFORCE_REBUILD=1)

# Create a precompiled header library (properly implemented)
if(ENABLE_PRECOMPILED_HEADERS)
  add_library(compas_pch INTERFACE)
  target_precompile_headers(compas_pch INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/compas.h)
  target_include_directories(compas_pch INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIRS}
    ${SHAPEOP_INCLUDE_DIR}
    ${SHAPEOP_SRC_DIR}
    ${SHAPEOP_API_DIR}
    ${nanobind_INCLUDE_DIRS}
  )
endif()

# Function to add our modules with common settings
function(add_shapeop_module module_name source_file)
  nanobind_add_module(
    ${module_name}
    STABLE_ABI
    NB_STATIC
    ${source_file}
  )
  
  # Add include directories for the module
  target_include_directories(${module_name} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIRS}
    ${SHAPEOP_INCLUDE_DIR}
    ${SHAPEOP_SRC_DIR}
    ${SHAPEOP_API_DIR}
    ${nanobind_INCLUDE_DIRS}
  )
  
  # Link shapeop and potentially PCH
  target_link_libraries(${module_name} PRIVATE shapeop shapeop_extensions)
  
  if(ENABLE_PRECOMPILED_HEADERS)
    target_link_libraries(${module_name} PRIVATE compas_pch)
  endif()
  
  # Install directive
  install(TARGETS ${module_name} LIBRARY DESTINATION compas_shapeop)
endfunction()

# Add our module using the function
add_shapeop_module(_shapeop src/shapeop.cpp)

# Add our optimized zero-copy module for high-performance applications
add_shapeop_module(_shapeoplibrary src/shapeoplibrary.cpp)

# Build configuration info
message(STATUS "============= Build Configuration =============")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Standard: C++${CMAKE_CXX_STANDARD}")
message(STATUS "Optimization: ${FAST_COMPILE} (O0 if ON, O3 if OFF)")
message(STATUS "Precompiled Headers: ${ENABLE_PRECOMPILED_HEADERS}")
message(STATUS "=============================================")

# Export compile commands for debugging
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)