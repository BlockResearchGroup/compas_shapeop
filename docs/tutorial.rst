********************************************************************************
Tutorial
********************************************************************************

Introduction to COMPAS ShapeOp
==============================

COMPAS ShapeOp provides Python bindings for the ShapeOp C++ physics solver through a zero-copy integration using nanobind. This tutorial will guide you through the basics of using COMPAS ShapeOp for geometry processing tasks like mesh planarization, 
regularization, and dynamic simulations.

Getting Started
===============

The core of COMPAS ShapeOp is the ``Solver`` class, which provides an interface to the C++ ShapeOp solver with efficient memory sharing between Python and C++.

Basic Workflow
--------------

The basic workflow for using COMPAS ShapeOp involves:

1. Creating or loading a mesh
2. Creating a solver with the mesh
3. Adding constraints and forces
4. Running the solver to optimize the mesh

The solver automatically handles initialization and mesh vertex updates.

Here's a simple example:

.. code-block:: python

    from compas.datastructures import Mesh
    from compas_shapeop import MeshSolver

    # 1. Create a mesh
    mesh = Mesh.from_obj("data/m0.obj")

    # 2. Initialize the solver
    solver = MeshSolver(mesh)

    # 3. Add constraints and forces
    solver.constrain_edge_lengths(weight=1.0)
    solver.constrain_face_planarity(weight=10.0)

    # 4. Run the solver
    solver.solve(iterations=10)

    # The mesh vertices are automatically updated
    for i, vertex in enumerate(mesh.vertices()):
        mesh.vertex_attributes(vertex, "xyz", points_ref[i])


Understanding Constraints
=========================

Constraints are the core of ShapeOp's functionality. They define goals for how the geometry should behave during the optimization process.

Geometric Constraints
---------------------

* **Closeness**: Keeps vertices close to their original positions
* **Edge Strain**: Maintains edge lengths within a specified range
* **Plane**: Makes vertices of a face lie on a plane
* **Circle**: Makes vertices lie on a circle
* **Similarity**: Transforms vertices to match a target shape
* **Regular Polygon**: Makes vertices form a regular polygon

Each constraint has a weight that determines its influence during the solving process. Higher weights make the constraint stronger.

Example: Adding Planarization Constraints
-----------------------------------------

.. code-block:: python

    # Add plane constraints to all faces
    for fkey in mesh.faces():
        face_vertices = mesh.face_vertices(fkey)
        if len(face_vertices) > 3:  # Triangles are already planar
            solver.add_plane_constraint(face_vertices, weight=10.0)

Working with Forces
===================

Forces provide external influences on the geometry during solving:

* **Vertex Force**: Apply a force vector to a specific vertex
* **Normal Force**: Apply forces along face normals (e.g., for inflation)

Example: Adding Gravity
-----------------------

.. code-block:: python

    # Add downward force (gravity) to all vertices
    solver.add_gravity(fz=-0.001)

Mesh Integration
================

COMPAS ShapeOp provides convenience methods for working with COMPAS meshes:

* ``MeshSolver(mesh)``: Create a solver from a COMPAS mesh
* ``MeshSolver.from_obj(path)``: Create a solver from an OBJ file
* ``MeshSolver.from_grid(dx, nx, dy, ny)``: Create a solver from a grid mesh
* ``constrain_edge_lengths()``: Add edge length constraints
* ``constrain_face_planarity()``: Add face planarity constraints
* ``constrain_face_regularization()``: Add face regularization constraints
* ``constrain_triface_bending()``: Add bending constraints between triangular faces
* ``fix_vertices()``: Fix vertices in place or to target positions
* ``add_gravity()``: Add gravity force
* ``inflate()``: Add inflation force

These methods simplify the process of setting up constraints for complex meshes.

Zero-Copy Integration
=====================

One of the key features of COMPAS ShapeOp is its zero-copy integration between Python and C++. When setting and getting point data, the library provides efficient memory handling:

.. code-block:: python

    # The solver's points property provides direct access to the C++ solver's memory
    solver = MeshSolver(mesh)
    points = solver.points  # This is a zero-copy view into the C++ solver's memory
    
    # Points are directly modified in the solver's memory
    # No need to call set_points() again!
    points[0, 2] += 1.0  # Modify Z-coordinate of first point
    
    # Solve to apply constraints with the modified points
    solver.solve(10)

Interactive Visualization
=========================

COMPAS ShapeOp works seamlessly with COMPAS Viewer for interactive visualization:

.. code-block:: python

    from compas_viewer import Viewer
    
    viewer = Viewer()
    mesh_obj = viewer.scene.add(mesh)
    
    @viewer.on(interval=1)
    def update(frame):
        # Run solver iteration
        solver.solve(10)
        
        # Update mesh
        for i, vertex in enumerate(mesh.vertices()):
            mesh.vertex_attributes(vertex, "xyz", points_ref[i])
        
        # Update viewer
        mesh_obj.update(update_data=True)
    
    viewer.show()

Next Steps
==========

Check out the :doc:`examples` section for more advanced usage scenarios, and the :doc:`api` reference for detailed documentation of all available classes and methods.
