********************************************************************************
Overview
********************************************************************************

`ShapeOp <https://shapeop.org/>`_ integration uses `nanobind <https://github.com/wjakob/nanobind>`_ for Python bindings in COMPAS ShapeOp. This guide explains how to contribute to the project, implement new constraints or forces, and understand the overall architecture of the codebase.

File and Folder Structure
-------------------------

The compas_shapeop package is organized with a clear separation between C++ backend and Python frontend:

Source Code
^^^^^^^^^^^
* ``src/shapeop/`` - Core ShapeOp C++ library files (Solver, Constraints, Forces)
* ``src/shapeop.cpp`` - Main C++ binding file that connects ShapeOp to Python via nanobind
* ``src/compas_shapeop/`` - Python frontend code that wraps the C++ bindings
* ``src/compas.h`` - Precompiled header with common includes for faster compilation

Build & Dependencies
^^^^^^^^^^^^^^^^^^^^
* ``build/`` - Build artifacts from CMake
* ``CMakeLists.txt`` - C++ project configuration for ShapeOp and nanobind integration
* ``external/`` - External C++ dependencies like Eigen (automatically downloaded if needed)
* ``pyproject.toml`` - Python project configuration (pip install -e .)
* ``requirements.txt`` - Runtime requirements
* ``requirements-dev.txt`` - Development requirements for building docs, testing, etc.
* ``tasks.py`` - Development tasks (invoke test, invoke docs, invoke format, invoke lint)

Examples and Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^
* ``examples/`` - Example scripts demonstrating various constraints and applications
* ``data/`` - Sample mesh files for examples
* ``docs/`` - Documentation source files

Code Structure
--------------

The compas_shapeop package has three main components:

1. **ShapeOp Core Library**: The C++ implementation of the ShapeOp solver, constraints, and forces (in ``src/shapeop/``)
2. **C++ Binding Layer**: The nanobind wrapper that exposes ShapeOp to Python (in ``src/shapeop.cpp``)
3. **Python Interface**: The Python classes that provide a user-friendly API (in ``src/compas_shapeop/``)

The C++ binding layer is implemented as the ``SolverWrapper`` class in ``shapeop.cpp``, which wraps the ``ShapeOp::Solver`` class and provides methods for adding constraints, applying forces, and solving the system.

Zero-Copy Integration
---------------------

One of the key features of compas_shapeop is the zero-copy integration between Python and C++. This means that when the solver modifies point positions, the NumPy array is automatically updated without any data copying or conversion.

This is achieved through nanobind's Eigen integration, which allows the C++ Eigen matrices to be directly exposed as NumPy arrays with shared memory. The Python ``Solver`` class maintains a direct view of the C++ solver's internal points matrix, enabling high-performance simulations.

Constraints and Forces
----------------------

ShapeOp provides several constraints and forces for geometry processing:

**Constraints**:

* Closeness - Keeps vertices close to their original or target positions
* EdgeStrain - Maintains edge lengths within a specified range
* Circle - Makes vertices lie on a circle
* Plane - Makes vertices lie on a plane
* Similarity - Transforms vertices to match a target shape
* Regular Polygon - Makes vertices form a regular polygon

**Forces**:

* VertexForce - Applies a force vector to specific vertices
* NormalForce - Applies forces along face normals for inflation/deflation
