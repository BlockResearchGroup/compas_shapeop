********************************************************************************
Contribute
********************************************************************************

This guide explains how to contribute to the COMPAS ShapeOp project, including adding new constraints, forces, or improving the existing codebase.

Setting Up for Development
--------------------------

1. Fork the repository on GitHub
2. Clone your fork locally
3. Set up the development environment as described in :doc:`conda_environment`
4. Create a new branch for your feature or bug fix

Binding Existing ShapeOp Constraints
------------------------------------

Not all constraints from the ShapeOp C++ library are currently bound to Python in COMPAS ShapeOp. Here's how to expose an existing ShapeOp constraint:

1. **Identify the Constraint**:
   
   * Locate the constraint class in the ShapeOp C++ library (typically in ``ext/ShapeOp/src/Constraints.h``)
   * Understand the parameters required by its constructor and any additional methods it may have

2. **Add a Method to SolverWrapper**:
   
   * In ``src/shapeop.cpp``, add a new method to the ``SolverWrapper`` class:

   .. code-block:: cpp
   
       int add_your_constraint(nb::list indices, float weight = 1.0) {
           if (!is_valid()) {
               throw std::runtime_error("Invalid solver");
           }
           
           // Convert Python list to std::vector<int>
           std::vector<int> ids;
           for (size_t i = 0; i < len(indices); i++) {
               ids.push_back(nb::cast<int>(indices[i]));
           }
           
           // Create the constraint using ShapeOp's factory
           auto constraint = ShapeOp::Constraint::shapeConstraintFactory(
               "YourConstraintType", ids, weight
           );
           
           // Additional configuration if needed
           // constraint->setSpecificParameter(...);
           
           // Add to solver and return ID
           return solver->addConstraint(constraint);
       }

3. **Expose the Method to Python**:
   
   * In the ``NB_MODULE`` section at the bottom of ``src/shapeop.cpp``, add your method:
   
   .. code-block:: cpp
   
       nb::class_<SolverWrapper>(m, "SolverWrapper")
           // ... existing methods
           .def("add_your_constraint", &SolverWrapper::add_your_constraint, 
                nb::arg("indices"), nb::arg("weight") = 1.0)
           // ... more methods

4. **Add Python Wrapper**:
   
   * In ``src/compas_shapeop/shapeop.py``, add a method to the ``Solver`` class:
   
   .. code-block:: python
   
       def add_your_constraint(self, indices, weight=1.0):
           """Add a your_constraint constraint to the solver.
           
           Parameters
           ----------
           indices : list[int]
               List of point indices to constrain.
           weight : float, optional
               Weight of the constraint, by default 1.0
               
           Returns
           -------
           int
               Constraint ID
           """
           return self._solver.add_your_constraint(indices, weight)

5. **Document in API**:
   
   * Ensure your method has a proper docstring as shown above
   * Update API documentation in ``docs/api/compas_shapeop.rst`` if needed

6. **Build and Test**:
   
   * Rebuild the extension: ``pip install -e .``
   * Create a simple test example to verify the constraint works

Common ShapeOp constraints you might want to bind include ``TriangleStrainConstraint``, ``TetrahedronStrainConstraint``, ``AreaConstraint``, or ``VolumeConstraint``.

Adding New Constraints
----------------------

To add a new constraint to COMPAS ShapeOp:

1. **C++ Implementation**: 
   
   * Add the constraint class to ``src/shapeop/custom_constraints/`` 
   * The class should inherit from ``ShapeOp::Constraint``
   * Implement required methods like ``project()`` and ``addConstraint()``

2. **Update C++ Binding**:
   
   * Add a method to the ``SolverWrapper`` class in ``src/shapeop.cpp``
   * Expose the new method in the Python module definition

3. **Python Wrapper**:

   * Add a corresponding method to the ``Solver`` class in ``src/compas_shapeop/shapeop.py``
   * Document the method with proper docstrings

Example: Adding a New Constraint
--------------------------------

Here's a simplified example of how to add a new constraint:

1. Create the C++ constraint implementation in ``src/shapeop/custom_constraints/myconstraint.cpp``:

.. code-block:: cpp

    #include "ShapeOp/Constraint.h"
    
    namespace ShapeOp {
    
    class MyConstraint : public Constraint {
    public:
        MyConstraint(const std::vector<int> &idI, Scalar weight) 
            : Constraint(idI, weight) {}
        
        void project(Matrix3X &positions, Matrix3X &projections) override {
            // Implement your constraint logic here
        }
        
        static std::shared_ptr<Constraint> create(const std::vector<int> &idI, Scalar weight) {
            return std::make_shared<MyConstraint>(idI, weight);
        }
    };
    
    } // namespace ShapeOp

2. Add the binding method to ``SolverWrapper`` in ``src/shapeop.cpp``:

.. code-block:: cpp

    // Add a new custom constraint
    int add_my_constraint(nb::list indices, float weight = 1.0) {
        if (!is_valid()) {
            throw std::runtime_error("Invalid solver");
        }
        
        // Convert Python list to std::vector<int>
        std::vector<int> ids;
        for (size_t i = 0; i < len(indices); i++) {
            int idx = nb::cast<int>(indices[i]);
            ids.push_back(idx);
        }
        
        // Create the constraint
        auto constraint = ShapeOp::MyConstraint::create(ids, weight);
        return solver->addConstraint(constraint);
    }

3. Expose the method in the Python module definition:

.. code-block:: cpp

    nb::class_<SolverWrapper>(m, "SolverWrapper")
        // ... existing methods
        .def("add_my_constraint", &SolverWrapper::add_my_constraint);

4. Add the Python wrapper in ``src/compas_shapeop/shapeop.py``:

.. code-block:: python

    def add_my_constraint(self, indices, weight=1.0):
        """Add a custom constraint to the solver.
        
        Parameters
        ----------
        indices : list
            List of vertex indices to constrain.
        weight : float, optional
            Weight of the constraint, by default 1.0
            
        Returns
        -------
        int
            The ID of the newly added constraint.
        """
        return self._solver.add_my_constraint(indices, weight)

Adding New Forces
-----------------

The process for adding new forces is similar to adding constraints:

1. Create a C++ force implementation inheriting from ``ShapeOp::Force``
2. Add binding methods to ``SolverWrapper``
3. Add Python wrapper methods to the ``Solver`` class

Testing Your Contributions
--------------------------

When adding new features, always include tests and examples:

1. Create an example script in the ``examples/`` directory
2. Add documentation for the new feature
3. Run linting checks with ``invoke lint``
4. Ensure all existing tests pass with ``invoke test``

Pull Request Process
--------------------

1. Push your changes to your fork
2. Create a pull request to the main repository
3. Describe your changes clearly in the PR description
4. Make sure your code follows the project's style guidelines
5. Address any feedback from code reviewers

Style Guidelines
----------------

- Follow the existing code style
- Use clear, descriptive variable and function names
- Document all public methods with docstrings
- Keep PR scope focused on a single feature or bug fix
- Write clean, maintainable code

Documentation
-------------

When adding new features, always update the documentation:

1. Add docstrings to all new methods
2. Update relevant tutorial sections or create new ones
3. Include examples that demonstrate the new functionality
4. Build and check the documentation with ``invoke docs``

For additional guidance, please reach out to the project maintainers or open an issue on GitHub.
