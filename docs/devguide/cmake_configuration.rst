********************************************************************************
CMake Configuration
********************************************************************************

The COMPAS ShapeOp extension is built with CMake, which configures and generates the C++ build system. The main ``CMakeLists.txt`` file contains all the necessary configuration to build the ShapeOp library and its Python bindings.

Project Structure
-----------------

The CMake configuration is organized as follows:

1. Basic project configuration and compiler options
2. Dependency handling (Eigen, nanobind, OpenMP)
3. ShapeOp library compilation
4. Precompiled headers setup
5. Python module definition and linking

Build Options
-------------

The CMakeLists.txt file provides several options to customize the build:

.. code-block:: cmake

    option(ENABLE_PRECOMPILED_HEADERS "Enable precompiled headers" ON)
    option(FAST_COMPILE "Optimize for faster compilation (-O0) vs execution (-O3)" OFF)
    option(USE_OPENMP "Enable OpenMP support for parallel processing" ON)

- ``ENABLE_PRECOMPILED_HEADERS``: Enable precompiled headers for faster compilation (default: ON)
- ``FAST_COMPILE``: Optimize for faster compilation time with -O0 instead of faster execution with -O3 (default: OFF)
- ``USE_OPENMP``: Enable OpenMP for parallel constraint solving (default: ON)

Dependencies
------------

The project has several external dependencies:

- **Eigen**: A header-only C++ library for linear algebra that is automatically downloaded if not already installed
- **nanobind**: A lightweight library for creating Python bindings, used to expose the C++ ShapeOp implementation to Python
- **OpenMP**: An optional dependency for parallel processing, which significantly improves performance for large meshes

ShapeOp Library
---------------

The core ShapeOp library is built as a static library:

.. code-block:: cmake

    add_library(shapeop STATIC
      # Core ShapeOp files
      ${SHAPEOP_SRC_DIR}/Constraint.cpp
      ${SHAPEOP_SRC_DIR}/Force.cpp
      ${SHAPEOP_SRC_DIR}/LSSolver.cpp
      ${SHAPEOP_SRC_DIR}/Solver.cpp
      # Custom constraints/forces
      ${SHAPEOP_SRC_DIR}/custom_constraints/normalforce.cpp
    )

This library contains both the core ShapeOp implementation and custom constraints/forces specific to COMPAS ShapeOp.

Python Module
-------------

The Python module is built using nanobind:

.. code-block:: cmake

    nanobind_add_module(
      _shapeop
      STABLE_ABI
      NB_STATIC
      src/shapeop.cpp
    )

The module is named ``_shapeop`` and is linked against the ShapeOp static library. It uses the stable ABI to ensure compatibility across Python versions.

Installation
------------

The installation target installs the module to the ``compas_shapeop`` package directory:

.. code-block:: cmake

    install(TARGETS _shapeop LIBRARY DESTINATION compas_shapeop)

This allows the Python frontend to import the C++ module as ``from compas_shapeop import _shapeop``.

Adding New Constraints or Forces
--------------------------------

To add a new constraint or force:

1. Add the C++ implementation to ``src/shapeop/`` or ``src/shapeop/custom_constraints/``
2. Add the file to the ``shapeop`` library in ``CMakeLists.txt``
3. Expose the new constraint/force in ``src/shapeop.cpp``
4. Create a Python wrapper in ``src/compas_shapeop/shapeop.py``

Build Process
-------------

The build process is handled by scikit-build-core, which manages the CMake configuration and build process during Python package installation. When you run ``pip install -e .``, scikit-build-core automatically:

1. Runs CMake to configure the build
2. Builds the C++ extension
3. Installs it alongside the Python files
