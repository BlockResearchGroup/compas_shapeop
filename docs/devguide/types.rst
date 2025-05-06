********************************************************************************
Types
********************************************************************************

Binding Approach by Exposing C++ Class
======================================

The compas_shapeop extension uses nanobind to expose C++ classes and methods to Python. The binding approach in compas_shapeop has two main components:

1. The C++ ``SolverWrapper`` class which wraps the ShapeOp C++ library and provides methods for Python
2. The Python ``Solver`` class which imports and provides a user-friendly interface to the C++ functionality

Memory Management with std::unique_ptr
--------------------------------------

The ``SolverWrapper`` class uses ``std::unique_ptr`` to manage the lifetime of the ``ShapeOp::Solver`` instance:

.. code-block:: cpp

    class SolverWrapper {
    private:
        std::unique_ptr<ShapeOp::Solver> solver;
        
        // Helper method to check if solver is valid
        bool is_valid() const {
            return solver != nullptr;
        }
        
    public:
        SolverWrapper() : solver(std::make_unique<ShapeOp::Solver>()) {
            if (!solver) {
                throw std::runtime_error("Failed to create ShapeOp solver");
            }
        }
        
        ~SolverWrapper() {
            // The unique_ptr will automatically release the solver
        }
        // ...
    };

There are several important reasons for using ``std::unique_ptr`` here:

1. **Automatic Memory Management**: The ``unique_ptr`` automatically deallocates the solver when the ``SolverWrapper`` instance is destroyed, preventing memory leaks
2. **Ownership Semantics**: The ``unique_ptr`` indicates that the ``SolverWrapper`` has exclusive ownership of the ``ShapeOp::Solver`` instance
3. **Exception Safety**: If an exception occurs during construction or operation, the solver will still be properly cleaned up
4. **Explicit Lifetime Management**: The solver's lifetime is explicitly tied to the ``SolverWrapper`` instance
5. **Python Binding Compatibility**: This approach works well with nanobind's memory management system, ensuring proper cleanup when Python objects are garbage collected

Using a raw pointer or a direct class member would require manual memory management or additional copying operations that could affect performance. With ``unique_ptr``, we get safe, efficient memory management with clear ownership semantics.

Class Binding
-------------

The C++ ``SolverWrapper`` class is bound to Python in ``src/shapeop.cpp`` using nanobind:

.. code-block:: cpp

    NB_MODULE(_shapeop, m) {
        // Give a clear docstring about this module
        m.doc() = "ShapeOp dynamic solver binding";
        
        // Define the solver class with a unique name
        nb::class_<SolverWrapper>(m, "SolverWrapper")
            .def(nb::init<>())
            .def("set_points", &SolverWrapper::set_points)
            .def("get_points_ref", &SolverWrapper::get_points_ref)
            .def("add_closeness_constraint", &SolverWrapper::add_closeness_constraint)
            // ... other methods
            .def("initialize", &SolverWrapper::initialize)
            .def("solve", &SolverWrapper::solve);
    }

The ``nb::class_`` template creates a Python binding for the C++ ``SolverWrapper`` class, and the ``.def()`` calls expose individual methods.

Method Binding
--------------

Methods are bound with appropriate type signatures to handle conversions between C++ and Python types:

.. code-block:: cpp

    // Example method binding with Python list input
    int add_closeness_constraint(nb::list indices, float weight = 1.0) {
        if (!is_valid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        // Convert Python list to std::vector<int>
        std::vector<int> ids;
        for (size_t i = 0; i < len(indices); i++) {
            int idx = nb::cast<int>(indices[i]);
            ids.push_back(idx);
        }
        
        // Call the ShapeOp method with the converted vector
        auto constraint = ShapeOp::Constraint::shapeConstraintFactory(
            "Closeness", ids, weight
        );
        return solver->addConstraint(constraint);
    }

Each method follows this pattern:
1. Accept Python-friendly types as input
2. Convert inputs to C++ types using nanobind's casting functions
3. Call the underlying ShapeOp C++ methods
4. Return results as Python-friendly types

Zero-Copy Memory Sharing
------------------------

One of the most important aspects of the binding is the zero-copy memory sharing between Eigen matrices and NumPy arrays:

.. code-block:: cpp

    // Direct access to ShapeOp's internal points matrix with zero-copy
    Eigen::Ref<Eigen::MatrixXd> get_points_ref() {
        if (!is_valid()) {
            throw std::runtime_error("Invalid solver");
        }
        // Direct access to the solver's points matrix
        return solver->points_;
    }

This method returns an ``Eigen::Ref`` to the internal points matrix, which nanobind automatically converts to a NumPy array view without copying the data. This allows Python code to directly read and write to the C++ memory. In python we access the coordinates from the solver.init() method. The return value is a numpy array that is a view of the C++ memory. This value is updated every time you call solver.solve(n) where n is the number of iterations.

Type Conversion
===============

Matching C++/Python types often takes the most of the time and requires careful attention. When implementing C++/Python bindings, follow these key patterns from the existing files or implement your own. If there are specific types you want to implement, review the `nanobind tests <https://github.com/wjakob/nanobind/tree/master/tests>`_ . Ask questions in discussion section for nanobind typing or follow previous issues. Current implementation provides examples for the following types:


* C++:
    * Use ``Eigen::Ref`` for matrix parameters, e.g. to transfer mesh vertex coordinates.
    * Return complex data as ``std::tuple<type, ...>`` types.
    * Use ``std::vector<type>`` for list copies otherwise use ``const std::vector<type> &``.
    * Use Eigen Matrix types in vectors ``const std::vector<Eigen::Matrix<type, ...>> &`` instead of reference type ``const std::vector<Eigen::Ref<...>> &``.
    * On Windows, ensure ``NOMINMAX`` is defined before including any Windows headers to prevent max/min macro conflicts.

* Python:
    * Use ``float64`` for vertices and ``int32`` for faces in numpy arrays
    * Enforce row-major (C-contiguous) order for matrices
    * Use libigl's matrix types (e.g., ``Eigen::MatrixXd``, ``Eigen::MatrixXi``)


Type Conversion Patterns
========================

When implementing C++/Python bindings, follow these established patterns:

Matrix Operations
-----------------

Use ``Eigen::Ref`` for efficient matrix passing:

.. code-block:: cpp

    void my_function(const Eigen::Ref<const Eigen::MatrixXd>& vertices,
                    const Eigen::Ref<const Eigen::MatrixXi>& faces);

Return complex mesh data as tuples:

.. code-block:: cpp

    return std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> my_function();

Enforce proper numpy array types using float64 and int32 in C-contiguous order:

.. code-block:: python

    import numpy as np
    from compas_libigl._nanobind import my_submodule

    # Convert mesh vertices and faces to proper numpy arrays
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    faces = np.asarray(mesh.faces, dtype=np.int32)

    # Pass to C++ function
    V, F = my_submodule.my_function(vertices, faces)


Vector Types
------------

For list data, choose between ``std::vector`` for value copies, ``const std::vector&`` for references, and ``std::vector<Eigen::Matrix<type, ...>>`` for matrix vectors.

Bind vector types explicitly:

.. code-block:: cpp

    // In module initialization
    nb::bind_vector<std::vector<double>>(m, "VectorDouble");

Access in Python:

.. code-block:: python

    # Get vector result
    vector_result = my_function()
    # Access elements by index
    x, y, z = vector_result[0], vector_result[1], vector_result[2]


Type Conversion Best Practices
==============================

When implementing new functionality:

* Matrix Operations:

  .. code-block:: cpp

      // GOOD: Use Eigen::Ref for matrix parameters
      void my_function(Eigen::Ref<const Eigen::MatrixXd> vertices);

      // BAD: Don't use raw matrices
      void my_function(Eigen::MatrixXd vertices);

* Return Types:

  .. code-block:: cpp

      // GOOD: Return complex data as tuples
      std::tuple<Eigen::MatrixXd, Eigen::MatrixXi> my_mesh_operation();

      // BAD: Don't use output parameters
      void my_mesh_operation(Eigen::MatrixXd& out_vertices);

* Vector Handling:

  .. code-block:: cpp

      // GOOD: Use const references for input vectors
      void my_function(const std::vector<double>& input);

      // GOOD: Return vectors by value
      std::vector<double> MyOperation();

      // BAD: Don't use non-const references
      void my_function(std::vector<double>& input);

* Matrix Vectors:

  .. code-block:: cpp

      // GOOD: Use Matrix types in vectors
      std::vector<Eigen::Matrix<double, 3, 1>> points;

      // BAD: Don't use Ref types in vectors
      std::vector<Eigen::Ref<Eigen::Vector3d>> points;

* Python Integration:

  .. code-block:: python

      # GOOD: Enforce proper types
      vertices = np.array(points, dtype=np.float64)
      faces = np.array(indices, dtype=np.int32)

      # BAD: Don't rely on automatic conversion
      vertices = points  # type not enforced
      faces = indices   # type not enforced

* Windows-Specific:

  .. code-block:: cpp

      // GOOD: Define NOMINMAX before Windows headers
      #ifdef _WIN32
      #define NOMINMAX
      #endif
      #include <windows.h>

      // BAD: Don't use Windows headers without NOMINMAX
      #include <windows.h>  # May cause conflicts with std::min/max